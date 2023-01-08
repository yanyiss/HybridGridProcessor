#include "OccReader.h"
#include "BRepBndLib.hxx"
#include "Bnd_Box.hxx"

namespace CADMesher
{
	OccReader::OccReader(QString &file)
	{
		fileName = file;
		SetShape();
		SetCADWireFrame();
	}

	void OccReader::SetShape()
	{
		std::string filetype;
		if (fileName.endsWith(".stp") || fileName.endsWith(".step") || fileName.endsWith(".STP") || fileName.endsWith(".STEP")) {
			reader = new STEPControl_Reader();
			dprint("CAD model from STEP file");
			filetype = "STEP";
		}
		else if (fileName.endsWith(".igs") || fileName.endsWith(".iges") || fileName.endsWith(".IGS") || fileName.endsWith(".IGES")) {
			reader = new IGESControl_Reader();
			dprint("CAD model from IGES file");
			filetype = "IGES";
		}
		else
		{
			dprinterror("Is anything wrong? It can't be other file types");
			exit(0);
		}

		dprint(filetype + "file read beginning!\n");
		reader->ReadFile(fileName.toLatin1().data());
		Standard_Integer NbTrans = reader->TransferRoots();
		globalmodel.aShape = reader->OneShape();
		dprint(filetype + "file read finished\n");
	}

	void OccReader::SetCADWireFrame()
	{
		ComputeFaceAndEdge();
		Discrete_Edge();
	}

	void OccReader::Set_TriMesh()
	{
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;
		Surface_TriMeshes.resize(faceshape.size());
		for (int i = 0; i < faceshape.size(); i++)
		{
			//if (i != 8)
				//continue;
			TriMesh &aMesh = Surface_TriMeshes[i];
			auto &wires = faceshape[i].wires;
			if (wires.empty() || !faceshape[i].if_exisited)
			{
				continue;
			}

			auto &aface = faceshape[i].face;
			TopLoc_Location loc;
			Handle(Geom_Surface) asurface = BRep_Tool::Surface(aface, loc);
			Standard_Real x, y, z, w;
			asurface->Bounds(x, y, z, w);

			int pointsnumber = 0;
			vector<Matrix2Xi> bnd;
			double x_step = 0;
			double y_step = 0;
			for (int m = 0; m < wires.size(); m++)
			{
				auto &edges = wires[m];
				int boundsum = 0;
				auto &end_paras = edgeshape[*edges.rbegin()].parameters;
				auto end_para = end_paras.col(end_paras.cols() - 1);
				for (int j = 0; j < edges.size(); j++)
				{
					auto &aedge = edgeshape[edges[j]];
					auto &boundpos = aedge.parameters;
					int cols = boundpos.cols() - 1;
					pointsnumber += cols;
					boundsum += cols;

					double temp;
					for (int r = 0; r < cols; r++)
					{
						temp = fabs(boundpos(0, r + 1) - boundpos(0, r));
						x_step = std::max(x_step, temp);
						temp = fabs(boundpos(1, r + 1) - boundpos(1, r));
						y_step = std::max(y_step, temp);
					}
				}
				Matrix2Xi bound(2, boundsum);
				for (int j = 0; j < boundsum - 1; j++)
				{
					bound.col(j) << j + pointsnumber - boundsum, j + pointsnumber - boundsum + 1;
				}
				bound.col(boundsum - 1) << pointsnumber - 1, pointsnumber - boundsum;
				bnd.push_back(bound);
			}

			Matrix2Xd all_pnts(2, pointsnumber);
			double if_reverse = aface.Orientation() ? -1.0 : 1.0;
			double ra = (y_step < epsilonerror ? 10 : x_step / y_step) * if_reverse;

			int s = 0;
			pointsnumber = 0;
			for (int j = 0; j < wires.size(); j++)
			{
				auto &edges = wires[j];
				for (int k = 0; k < edges.size(); k++)
				{
					auto &boundpos = edgeshape[edges[k]].parameters;
					int cols = boundpos.cols() - 1;
					all_pnts.block(0, s, 2, cols) = boundpos.block(0, 0, 2, cols);
#if 0
					dprint();
					dprint(k);
					for (int pp = 0; pp < boundpos.cols(); ++pp)
					{
						dprintwithprecision(15, boundpos(0, pp), boundpos(1, pp));
					}
#endif
					s += cols;
				}
				pointsnumber = s;
			}
			//处理内外边界相切的情况
			if (ProcessTangentialBoundary(i, outerFlag(all_pnts, bnd)))
			{
				s = 0;
				pointsnumber = 0;
				for (int j = 0; j < wires.size(); j++)
				{
					auto &edges = wires[j];
					for (int k = 0; k < edges.size(); k++)
					{
						auto &boundpos = edgeshape[edges[k]].parameters;
						int cols = boundpos.cols() - 1;
						all_pnts.block(0, s, 2, cols) = boundpos.block(0, 0, 2, cols);
						s += cols;
					}
					pointsnumber = s;
				}
			}

			Matrix2Xd boundary = all_pnts.block(0, 0, 2, pointsnumber);
			all_pnts.block(1, 0, 1, all_pnts.cols()) *= ra;
			triangulate(all_pnts, bnd, sqrt(3) * x_step * x_step *0.25 * mu, aMesh);
			double ra_inv = 1.0 / ra;
			for (auto tv : aMesh.vertices())
			{
				if (tv.idx() < pointsnumber)
				{
					aMesh.set_point(tv, TriMesh::Point(boundary(0, tv.idx()), boundary(1, tv.idx()), 0));
					continue;
				}
				auto pos = aMesh.point(tv);
				aMesh.set_point(tv, TriMesh::Point(pos[0], pos[1] * ra_inv, 0));
			}

			/*ClearBoundary(aMesh);
			continue;*/
			auto &Surface = faceshape[i].Surface;

			//Remesh in domain		
			if (1)
			{
				Riemannremesh Remesh(Surface, &aMesh);
				Remesh.remesh();
				//dprint("domain remesh done!");
			}

			ClearBoundary(aMesh);
			if (1)
			{
				Set_Curvature(Surface, aMesh);
				//dprint("GaussCurvature compute done!");
			}
			else
			{
				for (auto tv : aMesh.vertices())
				{
					aMesh.data(tv).GaussCurvature = 0;
				}
			}


			for (auto tv : aMesh.vertices())
			{
				auto pos = aMesh.point(tv);
				auto v = asurface->Value(pos[0], pos[1]);
				aMesh.set_point(tv, TriMesh::Point(v.X(), v.Y(), v.Z()));
			}

			//set flag1 and flag2
			int id = 0, begin, endid;
			for (int j = 0; j < wires.size(); j++)
			{
				begin = id;
				auto &edges = wires[j];
				for (int k = 0; k < edges.size(); k++)
				{
					auto &aedge = edgeshape[edges[k]];
					int cols = aedge.parameters.cols() - 1;
					if (aedge.if_C0)
					{
						if (k == edges.size() - 1) endid = begin;
						else endid = id + cols;
						for (int m = id; m < id + cols - 1; m++)
						{
							auto e = aMesh.edge_handle(aMesh.find_halfedge(aMesh.vertex_handle(m), aMesh.vertex_handle(m + 1)));
							aMesh.data(e).flag1 = true;
						}
						auto e = aMesh.edge_handle(aMesh.find_halfedge(aMesh.vertex_handle(id + cols - 1), aMesh.vertex_handle(endid)));
						aMesh.data(e).flag1 = true;
					}
					id += cols;
				}
			}
			id = 0;
			for (int j = 0; j < wires.size(); j++)
			{
				begin = id;
				auto &edges = wires[j];
				for (int k = 0; k < edges.size(); k++)
				{
					auto &aedge = edgeshape[edges[k]];
					int cols = aedge.parameters.cols() - 1;
					if (aedge.if_curvature)
					{
						if (k == edges.size() - 1) endid = begin;
						else endid = id + cols;
						for (int m = id; m < id + cols - 1; m++)
						{
							auto e = aMesh.edge_handle(aMesh.find_halfedge(aMesh.vertex_handle(m), aMesh.vertex_handle(m + 1)));
							aMesh.data(e).flag2 = true;
						}
						auto e = aMesh.edge_handle(aMesh.find_halfedge(aMesh.vertex_handle(id + cols - 1), aMesh.vertex_handle(endid)));
						aMesh.data(e).flag2 = true;
					}
					id += cols;
				}
			}
		}
		dprint("Piecewise TriMesh Done!");
	}

	void OccReader::Set_PolyMesh()
	{
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;

		Set_Offset_Grid();

		Surface_PolyMeshes.resize(faceshape.size());
		for (int i = 0; i < faceshape.size(); i++)
		{
			//if (i != 68) continue;
			Mesh &newmesh = Surface_PolyMeshes[i];
			TriMesh aMesh;
			auto &wires = faceshape[i].wires;
			if (wires.empty() || !faceshape[i].if_exisited)
			{
				continue;
			}

			auto &aface = faceshape[i].face;
			TopLoc_Location loc;
			Handle(Geom_Surface) asurface = BRep_Tool::Surface(aface, loc);
			Standard_Real x, y, z, w;
			asurface->Bounds(x, y, z, w);

			int pointsnumber = 0;
			vector<Matrix2Xi> bnd;
			double x_step = 0;
			double y_step = 0;
			for (int m = 0; m < wires.size(); m++)
			{
				auto &edges = wires[m];
				int boundsum = 0;
				auto &end_paras = edgeshape[*edges.rbegin()].parameters;
				auto end_para = end_paras.col(end_paras.cols() - 1);
				for (int j = 0; j < edges.size(); j++)
				{
					auto &aedge = edgeshape[edges[j]];
					auto &boundpos = aedge.parameters;
					int cols = boundpos.cols() - 1;
					pointsnumber += cols;
					boundsum += cols;

					double temp;
					for (int r = 0; r < cols; r++)
					{
						temp = fabs(boundpos(0, r + 1) - boundpos(0, r));
						x_step = std::max(x_step, temp);
						temp = fabs(boundpos(1, r + 1) - boundpos(1, r));
						y_step = std::max(y_step, temp);
					}
				}
				Matrix2Xi bound(2, boundsum);
				for (int j = 0; j < boundsum - 1; j++)
				{
					bound.col(j) << j + pointsnumber - boundsum, j + pointsnumber - boundsum + 1;
				}
				bound.col(boundsum - 1) << pointsnumber - 1, pointsnumber - boundsum;
				bnd.push_back(bound);
			}

			Matrix2Xd all_pnts(2, pointsnumber);
			pointsnumber = 0;
			int s = 0;
			double if_reverse = aface.Orientation() ? -1.0 : 1.0;


			double ra = (y_step < epsilonerror ? 10 : x_step / y_step) * if_reverse;

			for (int j = 0; j < wires.size(); j++)
			{
				auto &edges = wires[j];
				for (int k = 0; k < edges.size(); k++)
				{
					auto &boundpos = edgeshape[edges[k]].parameters;
					int cols = boundpos.cols() - 1;
					all_pnts.block(0, s, 2, cols) = boundpos.block(0, 0, 2, cols);
					s += cols;
				}
				pointsnumber = s;
			}
			Matrix2Xd boundary = all_pnts.block(0, 0, 2, pointsnumber);
			auto &Surface = faceshape[i].Surface;
			double ra_inv = 1.0 / ra;

			//该面是否需生成四边形
			if (!faceshape[i].if_quad)
			{
				all_pnts.block(1, 0, 1, all_pnts.cols()) *= ra;
				triangulate(all_pnts, bnd, sqrt(3) * x_step * x_step *0.25 * mu, aMesh);
				for (auto tv : aMesh.vertices())
				{
					if (tv.idx() < pointsnumber)
					{
						aMesh.set_point(tv, TriMesh::Point(boundary(0, tv.idx()), boundary(1, tv.idx()), 0));
						continue;
					}
					auto pos = aMesh.point(tv);
					aMesh.set_point(tv, TriMesh::Point(pos[0], pos[1] * ra_inv, 0));
				}

				//if (!OpenMesh::IO::write_mesh(aMesh, "1.obj"))
				//{
				//	std::cerr << "fail";
				//}
				//Remesh in domain		
				Riemannremesh Remesh(Surface, &aMesh);
				Remesh.remesh();
				tri2poly(aMesh, newmesh);
				//if (!OpenMesh::IO::write_mesh(newmesh, "2.obj"))
				//{
				//	std::cerr << "fail";
				//}
			}
			else
			{
				int quad_num = faceshape[i].quad_num;
				all_pnts.block(1, 0, 1, all_pnts.cols()) *= if_reverse;
				int beginid = 0, endid, offset_pnt_num;
				auto &edges = wires.front();
				for (int j = 0; j < edges.size(); j++)
				{
					if (!edgeshape[edges[j]].if_curvature)
					{
						beginid += (edgeshape[edges[j]].parameters).cols() - 1;
						continue;
					}
					offset_pnt_num = (edgeshape[edges[j]].parameters).cols();
					break;
				}
				endid = (beginid + offset_pnt_num - 1 >= pointsnumber) ? 0 : beginid + offset_pnt_num - 1;
				if (endid + quad_num > pointsnumber || (beginid && beginid < quad_num))
				{
					dprint("There exist a too short edge that cannot construct offset!");
					system("pause");
				}
				vector<Matrix2Xd> offset_pnts;
				Offset_lines(all_pnts, offset_pnts, beginid, offset_pnt_num, quad_num);
				Matrix2Xd sub_all_pnts(2, pointsnumber - 2 * quad_num);
				if (!beginid || !endid)
				{
					sub_all_pnts.block(0, 0, 2, pointsnumber - 2 * quad_num - offset_pnt_num + 2) = all_pnts.block(0, endid + quad_num, 2, pointsnumber - 2 * quad_num - offset_pnt_num + 2);
				}
				else
				{
					sub_all_pnts.block(0, 0, 2, pointsnumber - endid - quad_num) = all_pnts.block(0, endid + quad_num, 2, pointsnumber - endid - quad_num);
					sub_all_pnts.block(0, pointsnumber - endid - quad_num, 2, beginid - quad_num + 1) = all_pnts.block(0, 0, 2, beginid - quad_num + 1);
				}
				sub_all_pnts.block(0, pointsnumber - 2 * quad_num - offset_pnt_num + 2, 2, offset_pnt_num - 2) = offset_pnts.back();
				vector<Matrix2Xi> sub_bnd;
				Matrix2Xi sub_bnd_(2, pointsnumber - 2 * quad_num);
				sub_bnd_.block(0, 0, 2, pointsnumber - 2 * quad_num - 1) = (bnd.front()).block(0, 0, 2, pointsnumber - 2 * quad_num - 1);
				sub_bnd_(0, pointsnumber - 2 * quad_num - 1) = pointsnumber - 2 * quad_num - 1;
				sub_bnd_(1, pointsnumber - 2 * quad_num - 1) = 0;
				sub_bnd.push_back(sub_bnd_);
				sub_all_pnts.block(1, 0, 1, sub_all_pnts.cols()) *= if_reverse;
				Matrix2Xd aboundary(2, sub_all_pnts.cols());
				aboundary = sub_all_pnts;
				sub_all_pnts.block(1, 0, 1, sub_all_pnts.cols()) *= ra;
				triangulate(sub_all_pnts, sub_bnd, sqrt(3) * x_step * x_step *0.25 * mu, aMesh);
				for (auto tv : aMesh.vertices())
				{
					if (tv.idx() < sub_all_pnts.cols())
					{
						aMesh.set_point(tv, TriMesh::Point(aboundary(0, tv.idx()), aboundary(1, tv.idx()), 0));
						continue;
					}
					auto pos = aMesh.point(tv);
					aMesh.set_point(tv, TriMesh::Point(pos[0], pos[1] * ra_inv, 0));
				}

				//Remesh in domain		
				Riemannremesh Remesh(Surface, &aMesh);
				Remesh.remesh();

				//******construct hybridmesh*****//
				//add boundary points
				for (int j = 0; j < pointsnumber; j++)
				{
					auto newp = Mesh::Point(boundary(0, j), boundary(1, j), 0);
					newmesh.add_vertex(newp);
				}

				//if (!OpenMesh::IO::write_mesh(newmesh, "1.obj"))
				//{
				//	std::cerr << "fail";
				//}

				//add offset points and internal points
				for (int j = 0; j < quad_num; j++)
				{
					auto &offsetline = offset_pnts[j];
					for (int k = 0; k < offset_pnt_num - 2; k++)
					{
						auto newp = Mesh::Point(offsetline(0, k), offsetline(1, k)* if_reverse, 0);
						newmesh.add_vertex(newp);
					}
				}
				for (int j = pointsnumber - 2 * quad_num; j < aMesh.n_vertices(); j++)
				{
					auto newp = aMesh.point(aMesh.vertex_handle(j));
					newmesh.add_vertex(Mesh::Point(newp[0], newp[1], 0));
				}

				//add offset faces
				vector<Mesh::VertexHandle> facevhandle;
				Mesh::VertexHandle v1, v2, v3, v4;
				v1 = newmesh.vertex_handle(beginid);
				v2 = newmesh.vertex_handle(beginid + 1);
				int newid = beginid ? beginid : pointsnumber;
				for (int k = 0; k < quad_num; k++)
				{
					v3 = newmesh.vertex_handle(pointsnumber + k * (offset_pnt_num - 2));
					v4 = newmesh.vertex_handle(newid - k - 1);
					facevhandle = { v1, v2, v3, v4 };
					newmesh.add_face(facevhandle);
					v1 = v4;
					v2 = v3;
				}
				for (int j = 1; j < offset_pnt_num - 2; j++)
				{
					v1 = newmesh.vertex_handle(j + beginid);
					v2 = newmesh.vertex_handle(j + 1 + beginid);
					for (int k = 0; k < quad_num; k++)
					{
						v3 = newmesh.vertex_handle(pointsnumber + k * (offset_pnt_num - 2) + j);
						v4 = newmesh.vertex_handle(pointsnumber + k * (offset_pnt_num - 2) + j - 1);
						facevhandle = { v1, v2, v3, v4 };
						newmesh.add_face(facevhandle);
						v1 = v4;
						v2 = v3;
					}
				}
				v1 = newmesh.vertex_handle(beginid + offset_pnt_num - 2);
				v2 = newmesh.vertex_handle(endid);
				for (int k = 0; k < quad_num; k++)
				{
					v3 = newmesh.vertex_handle(endid + k + 1);
					v4 = newmesh.vertex_handle(pointsnumber + (k + 1) * (offset_pnt_num - 2) - 1);
					facevhandle = { v1, v2, v3, v4 };
					newmesh.add_face(facevhandle);
					v1 = v4;
					v2 = v3;
				}

				//if (!OpenMesh::IO::write_mesh(newmesh, "2.obj"))
				//{
				//	std::cerr << "fail";
				//}

				//add internal faces
				for (auto f : aMesh.faces())
				{
					facevhandle = {};
					for (auto fv : aMesh.fv_range(f))
					{
						if (fv.idx() < pointsnumber - 2 * quad_num - offset_pnt_num + 2)
							facevhandle.push_back(newmesh.vertex_handle((endid + quad_num + fv.idx()) % pointsnumber));
						else
							facevhandle.push_back(newmesh.vertex_handle(fv.idx() + quad_num * offset_pnt_num));
					}
					newmesh.add_face(facevhandle);
				}
				//if (!OpenMesh::IO::write_mesh(newmesh, "3.obj"))
				//{
				//	std::cerr << "fail";
				//}
			}

			if (0)
			{
				Set_Curvature(Surface, newmesh);
				//dprint("GaussCurvature compute done!");
			}
			else
			{
				for (auto tv : newmesh.vertices())
				{
					newmesh.data(tv).GaussCurvature = 0;
				}
			}

			//Set_Curvature(Surface, newmesh);

			for (auto tv : newmesh.vertices())
			{
				auto pos = newmesh.point(tv);
				auto v = asurface->Value(pos[0], pos[1]);
				newmesh.set_point(tv, Mesh::Point(v.X(), v.Y(), v.Z()));
			}

			//if (!OpenMesh::IO::write_mesh(newmesh, "4.obj"))
			//{
			//	std::cerr << "fail";
			//}

			//set flag1
			int id = 0, begin, endid;
			for (int j = 0; j < wires.size(); j++)
			{
				begin = id;
				auto &edges = wires[j];
				for (int k = 0; k < edges.size(); k++)
				{
					auto &aedge = edgeshape[edges[k]];
					int cols = aedge.parameters.cols() - 1;
					if (aedge.if_C0)
					{
						if (k == edges.size() - 1) endid = begin;
						else endid = id + cols;
						for (int m = id; m < id + cols - 1; m++)
						{
							auto e = newmesh.edge_handle(newmesh.find_halfedge(newmesh.vertex_handle(m), newmesh.vertex_handle(m + 1)));
							newmesh.data(e).flag1 = true;
						}
						auto e = newmesh.edge_handle(newmesh.find_halfedge(newmesh.vertex_handle(id + cols - 1), newmesh.vertex_handle(endid)));
						newmesh.data(e).flag1 = true;
					}
					id += cols;
				}
			}
		}
		dprint("Piecewise PolyMesh Done!");
	}

	void OccReader::ComputeFaceAndEdge()
	{
		TopoDS_Shape &aShape = globalmodel.aShape;

		bbmin[0] = DBL_MAX; bbmin[1] = DBL_MAX; bbmin[2] = DBL_MAX;
		bbmax[0] = -DBL_MAX; bbmax[1] = -DBL_MAX; bbmax[2] = -DBL_MAX;
		for (TopExp_Explorer vertexExp(aShape, TopAbs_VERTEX); vertexExp.More(); vertexExp.Next())
		{
			gp_Pnt v = BRep_Tool::Pnt(TopoDS::Vertex(vertexExp.Current()));
			bbmin[0] = std::min(bbmin[0], v.X()); bbmax[0] = std::max(bbmax[0], v.X());
			bbmin[1] = std::min(bbmin[1], v.Y()); bbmax[1] = std::max(bbmax[1], v.Y());
			bbmin[2] = std::min(bbmin[2], v.Z()); bbmax[2] = std::max(bbmax[2], v.Z());
		}
		//initialRate = 0.05;
		epsratio *= (bbmin - bbmax).norm();

		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		faceshape.clear();
		int faceSize = 0;
		for (TopExp_Explorer faceExp(aShape, TopAbs_FACE); faceExp.More(); faceExp.Next())
		{
			faceSize++;
		}
		faceshape.reserve(faceSize);

		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;
		edgeshape.clear();
		int edgeSize = 0;
		for (TopExp_Explorer edgeExp(aShape, TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
		{
			edgeSize++;
		}
		edgeshape.reserve(edgeSize);

		//a new method to order the faces and edges
		faceSize = 0;
		edgeSize = 0;
		//double vertexThreshold = expected_edge_length * degeneratedRate;
		for (TopExp_Explorer faceExp(aShape, TopAbs_FACE); faceExp.More(); faceExp.Next())
		{
			auto &aface = TopoDS::Face(faceExp.Current());
			faceshape.emplace_back(faceSize, aface);
			faceshape[faceSize].wires.reserve(aface.NbChildren());
			for (TopExp_Explorer wireExp(aface, TopAbs_WIRE); wireExp.More(); wireExp.Next())
			{
				TopoDS_Wire awire = TopoDS::Wire(wireExp.Current());
				vector<int> edges;
				edges.reserve(awire.NbChildren());

				for (TopExp_Explorer edgeExp(awire, TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
				{
					auto &aedge = TopoDS::Edge(edgeExp.Current());
					if (BRep_Tool::Degenerated(aedge))
					{
						continue;
					}
					/*if (!BRep_Tool::IsClosed(aedge))
					{
						gp_Pnt p0 = BRep_Tool::Pnt(TopExp::FirstVertex(aedge));
						gp_Pnt p1 = BRep_Tool::Pnt(TopExp::LastVertex(aedge));
						if (p0.IsEqual(p1, vertexThreshold))
						{
							continue;
						}
					}*/

					edgeshape.emplace_back(edgeSize, aedge);
					edgeshape[edgeSize].main_face = faceSize;
					edges.push_back(edgeSize);
					edgeSize++;
				}
				if (awire.Orientation() == TopAbs_REVERSED)
				{
					edges = vector<int>(edges.rbegin(), edges.rend());
				}
				if (!edges.empty())
				{
					//在大量测试模型时，发现有少量模型出现边并非完全逆时针排序，而是有“插队”现象，特此根据边的首尾节点重新排序
					auto getVertex = [&](TopoDS_Edge &aedge, bool isLast)
					{
						if ((aedge.Orientation() == TopAbs_FORWARD) == isLast)
							return TopExp::LastVertex(aedge);
						else
							return TopExp::FirstVertex(aedge);
					};
					for (int i = 0; i < edges.size() - 1; ++i)
					{
						if (getVertex(edgeshape[edges[i]].edge, true).IsSame(getVertex(edgeshape[edges[i + 1]].edge, false)))
							continue;
						for (int j = i + 2; j < edges.size(); ++j)
						{
							if (getVertex(edgeshape[edges[i]].edge, true).IsSame(getVertex(edgeshape[edges[j]].edge, false)))
							{
								std::swap(edges[i + 1], edges[j]);
								continue;
							}
						}
					}
					faceshape[faceSize].wires.push_back(edges);
				}
			}
			faceSize++;
		}

		//测试 AR15-milSpec-lower-receiver-scale.stp 时产生未合并的边，推测是模型文件问题
		/*auto e0 = edgeshape[575].edge;
		auto e1 = edgeshape[1600].edge;
		auto pnt0 = BRep_Tool::Pnt(TopExp::FirstVertex(e0));
		auto pnt1 = BRep_Tool::Pnt(TopExp::LastVertex(e0));
		dprint(pnt0.X(), pnt0.Y(), pnt0.Z(), pnt1.X(), pnt1.Y(), pnt1.Z());
		pnt0 = BRep_Tool::Pnt(TopExp::FirstVertex(e1));
		pnt1 = BRep_Tool::Pnt(TopExp::LastVertex(e1));
		dprint(pnt0.X(), pnt0.Y(), pnt0.Z(), pnt1.X(), pnt1.Y(), pnt1.Z());
		dprint(TopExp::FirstVertex(e0).IsPartner(TopExp::LastVertex(e1)));
		dprint(TopExp::FirstVertex(e0).IsSame(TopExp::LastVertex(e1)));
		dprint(TopExp::FirstVertex(e0).IsEqual(TopExp::LastVertex(e1)));
		dprint(TopExp::FirstVertex(e1).IsPartner(TopExp::LastVertex(e0)));
		dprint(TopExp::FirstVertex(e1).IsSame(TopExp::LastVertex(e0)));
		dprint(TopExp::FirstVertex(e1).IsEqual(TopExp::LastVertex(e0)));
		dprint(e0.IsEqual(e1));
		dprint(e0.IsSame(e1));
		dprint(e0.IsPartner(e1));
		auto t0 = e0.TShape();
		auto t1 = e1.TShape();*/


		auto edgeEnd = edgeshape.end();
		for (auto aedge = edgeshape.begin(); aedge != edgeEnd; ++aedge)
		{
			if (aedge->secondary_face == -1)
			{
				for (auto redge = aedge + 1; redge != edgeEnd; ++redge)
				{
					if (aedge->edge.IsSame(redge->edge))
					{
						aedge->reversed_edge = redge->id;
						redge->reversed_edge = aedge->id;
						aedge->secondary_face = redge->main_face;
						redge->secondary_face = aedge->main_face;
						break;
					}
				}
			}
		}

		dprint("CAD Info:\nModel Size: [", bbmin[0], bbmax[0], "],[", bbmin[1], bbmax[1], "],[", bbmin[2], bbmax[2], "]"
			"\nSurface Number: ", faceshape.size(),
			"\nEdge Number: ", edgeshape.size(),
			"\n\nCompute Topology Done!");
	}

	void OccReader::Discrete_Edge()
	{
		expected_edge_length = initialRate * (bbmin - bbmax).norm();
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;

		const double GL_c[8] = { 0.1012285363,0.2223810345,0.3137066459,0.3626837834,0.3626837834 ,0.3137066459,0.2223810345,0.1012285363 };
		const double GL_x[8] = { -0.9602898565,-0.7966664774,-0.5255354099,-0.1834346425,0.1834346425,0.5255354099,0.7966664774,0.9602898565 };

		for (auto edge = edgeshape.begin(); edge != edgeshape.end(); ++edge)
		{
			/*static int ti = 0;
			dprint(ti++);*/
			if (!edge->if_exisited) continue;
			auto &aedge = edge->edge;
			if (!edge->if_exisited) continue;
			auto &mainface = faceshape[edge->main_face].face;
			TopLoc_Location loc;
			Handle(Geom_Surface) asurface = BRep_Tool::Surface(mainface, loc);
			Standard_Real first = 0;
			Standard_Real last = 0;
			Handle(Geom_Curve) acurve = BRep_Tool::Curve(aedge, first, last);
			double curve_length = 0;
			gp_Pnt P;
			gp_Vec d1;
			//Gauss-Lengred积分公式
			double a = (last - first) / 2;
			double b = (last + first) / 2;
			for (int i = 0; i < 8; i++)
			{
				acurve->D1(a*GL_x[i] + b, P, d1);
				curve_length += a * GL_c[i] * d1.Magnitude();
			}
			edge->length = curve_length;
			int segment_number = std::max(/*BRep_Tool::IsClosed(aedge) ? 4 : 4*/2, int(curve_length / expected_edge_length));
			edge->parameters.resize(2, segment_number + 1);
			Handle_Geom2d_Curve thePCurve = BRep_Tool::CurveOnSurface(aedge, mainface, first, last);
			double step = (last - first) / segment_number;
			gp_Pnt2d uv;
			if (aedge.Orientation() == TopAbs_FORWARD)
			{
				for (int i = 0; i < segment_number; i++)
				{
					uv = thePCurve->Value(first + i * step);
					edge->parameters.col(i) << uv.X(), uv.Y();
				}
				uv = thePCurve->Value(last);
				edge->parameters.col(segment_number) << uv.X(), uv.Y();
			}
			else
			{
				for (int i = 0; i < segment_number; i++)
				{
					uv = thePCurve->Value(last - i * step);
					edge->parameters.col(i) << uv.X(), uv.Y();
				}
				uv = thePCurve->Value(first);
				edge->parameters.col(segment_number) << uv.X(), uv.Y();
			}
		}
		dprint("Discrete Edges Done!");
	}

	bool OccReader::ProcessTangentialBoundary(int fid, int bid)
	{
		//return false;
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;
		auto &wires = faceshape[fid].wires;
		bool flag = false;
		if (wires.size() > 1)
		{
			std::vector<TopoDS_Vertex> boundVertex;
			boundVertex.reserve(wires[bid].size());
			for (int j = 0; j < wires[bid].size(); ++j)
			{
				auto &ae = edgeshape[wires[bid][j]].edge;
				if (ae.Orientation() == TopAbs_FORWARD)
					boundVertex.push_back(TopExp::LastVertex(ae));
				else
					boundVertex.push_back(TopExp::FirstVertex(ae));
			}
			for (int m = 0; m < wires.size(); m++)
			{
				if (m == bid)
					continue;
				auto &edges = wires[m];
				for (int j = 0; j < edges.size(); j++)
				{
					auto &ae = edgeshape[edges[j]].edge;
					TopoDS_Vertex tv = ae.Orientation() == TopAbs_FORWARD ? TopExp::LastVertex(ae) : TopExp::FirstVertex(ae);
					for (int k = 0; k < boundVertex.size(); ++k)
					{
						if (tv.IsSame(boundVertex[k]))
						{
							auto &e0 = edgeshape[edges[j]].parameters;
							auto &e1 = edgeshape[edges[(j + 1) % edges.size()]].parameters;
							Vector2d v0 = e0.col(e0.cols() - 2) - e0.col(e0.cols() - 1);
							Vector2d v1 = e1.col(1) - e1.col(0);
							double alpha = acos(v0.dot(v1))*(v0(0)*v1(1) < v0(1)*v1(0) ? 1 : -1);
							alpha = alpha > 0 ? alpha : alpha + 2 * PI;
							Matrix2d M; M << cos(alpha), -sin(alpha), cos(alpha), sin(alpha);
							Vector2d v = (M * v1).normalized()*v1.norm()*0.1*(faceshape[fid].face.Orientation() ? -1 : 1) + e1.col(0);
							e0.col(e0.cols() - 1) = v;
							e1.col(0) = v;
							flag = true;
						}
					}
				}
			}
		}
		return flag;
	}

	void OccReader::ClearBoundary(TriMesh &tm)
	{
		for (auto &te : tm.edges())
		{
			if (te.v0().is_boundary() && te.v1().is_boundary() && !te.is_boundary())
			{
				auto newvert = tm.add_vertex(tm.calc_edge_midpoint(te));
				tm.split_edge(te, newvert);
			}
		}
		tm.garbage_collection();
	}

	void OccReader::Face_type()
	{
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;

		for (int i = 0; i < faceshape.size(); i++)
		{
			//dprint("face", i);
			if (!faceshape[i].if_exisited) continue;
			if (faceshape[i].wires.empty()) continue;
			//if (i >= 6 && i <= 28) continue;
			TopLoc_Location loca;
			opencascade::handle<Geom_Surface> geom_surface = BRep_Tool::Surface(faceshape[i].face, loca);
			opencascade::handle<Standard_Type> type = geom_surface->DynamicType();
			if (type == STANDARD_TYPE(Geom_Plane))
			{
				//dprint(i, "plane");
				opencascade::handle<Geom_Plane> geom_plane = Handle(Geom_Plane)::DownCast(geom_surface);
				auto local_coordinate = &geom_plane->Position();
				gp_Pnt oringe = local_coordinate->Location();
				gp_Dir xdir = local_coordinate->XDirection();
				gp_Dir ydir = local_coordinate->YDirection();
				faceshape[i].Surface = new PlaneType(Point(oringe.X(), oringe.Y(), oringe.Z()), Point(xdir.X(), xdir.Y(), xdir.Z()), Point(ydir.X(), ydir.Y(), ydir.Z()));
			}
			else if (type == STANDARD_TYPE(Geom_CylindricalSurface))
			{
				//dprint(i, "CylindricalSurface");
				opencascade::handle<Geom_CylindricalSurface> geom_cylindricalsurface = opencascade::handle<Geom_CylindricalSurface>::DownCast(geom_surface);;
				auto local_coordinate = &geom_cylindricalsurface->Position();
				gp_Pnt oringe = local_coordinate->Location();
				gp_Dir xdir = local_coordinate->XDirection();
				gp_Dir ydir = local_coordinate->YDirection();
				gp_Dir zdir = (geom_cylindricalsurface->Axis()).Direction();
				double r = geom_cylindricalsurface->Radius();
				faceshape[i].Surface = new CylindricalType(Point(oringe.X(), oringe.Y(), oringe.Z()), Point(xdir.X(), xdir.Y(), xdir.Z()), Point(ydir.X(), ydir.Y(), ydir.Z()), Point(zdir.X(), zdir.Y(), zdir.Z()), r);
			}
			else if (type == STANDARD_TYPE(Geom_ConicalSurface))
			{
				//dprint(i, "ConicalSurface");
				opencascade::handle<Geom_ConicalSurface> geom_ConicalSurface = opencascade::handle<Geom_ConicalSurface>::DownCast(geom_surface);
				auto local_coordinate = &(geom_ConicalSurface->Position());
				gp_Pnt oringe = local_coordinate->Location();
				gp_Dir xdir = local_coordinate->XDirection();
				gp_Dir ydir = local_coordinate->YDirection();
				gp_Dir zdir = (geom_ConicalSurface->Axis()).Direction();
				double r = geom_ConicalSurface->RefRadius();
				double ang = geom_ConicalSurface->SemiAngle();
				faceshape[i].Surface = new ConicalType(Point(oringe.X(), oringe.Y(), oringe.Z()), Point(xdir.X(), xdir.Y(), xdir.Z()), Point(ydir.X(), ydir.Y(), ydir.Z()), Point(zdir.X(), zdir.Y(), zdir.Z()), r, ang);
			}
			else if (type == STANDARD_TYPE(Geom_SphericalSurface))
			{
				//dprint(i, "SphericalSurface");
				opencascade::handle<Geom_SphericalSurface> geom_SphericalSurface = opencascade::handle<Geom_SphericalSurface>::DownCast(geom_surface);
				auto local_coordinate = &(geom_SphericalSurface->Position());
				gp_Pnt oringe = local_coordinate->Location();
				gp_Dir xdir = local_coordinate->XDirection();
				gp_Dir ydir = local_coordinate->YDirection();
				gp_Dir zdir = (geom_SphericalSurface->Axis()).Direction();
				double r = geom_SphericalSurface->Radius();
				faceshape[i].Surface = new SphericalType(Point(oringe.X(), oringe.Y(), oringe.Z()), Point(xdir.X(), xdir.Y(), xdir.Z()), Point(ydir.X(), ydir.Y(), ydir.Z()), Point(zdir.X(), zdir.Y(), zdir.Z()), r);
			}
			else if (type == STANDARD_TYPE(Geom_ToroidalSurface))
			{
				//dprint(i, "ToroidalSurface");
				opencascade::handle<Geom_ToroidalSurface> geom_ToroidalSurface = opencascade::handle<Geom_ToroidalSurface>::DownCast(geom_surface);
				auto local_coordinate = &(geom_ToroidalSurface->Position());
				gp_Pnt oringe = local_coordinate->Location();
				gp_Dir xdir = local_coordinate->XDirection();
				gp_Dir ydir = local_coordinate->YDirection();
				gp_Dir zdir = (geom_ToroidalSurface->Axis()).Direction();
				double r1 = geom_ToroidalSurface->MajorRadius();
				double r2 = geom_ToroidalSurface->MinorRadius();
				faceshape[i].Surface = new ToroidalType(Point(oringe.X(), oringe.Y(), oringe.Z()), Point(xdir.X(), xdir.Y(), xdir.Z()), Point(ydir.X(), ydir.Y(), ydir.Z()), Point(zdir.X(), zdir.Y(), zdir.Z()), r1, r2);
			}
			else if (type == STANDARD_TYPE(Geom_BezierSurface))
			{
				//dprint(i, "BezierSurface");
				opencascade::handle<Geom_BezierSurface> geom_beziersurface = opencascade::handle<Geom_BezierSurface>::DownCast(geom_surface);

				Standard_Real U1, U2, V1, V2;
				geom_beziersurface->Bounds(U1, U2, V1, V2);
				TColgp_Array2OfPnt controlpoints = geom_beziersurface->Poles();
				vector<vector<Point>> cp(controlpoints.NbRows());
				for (int r = 1; r <= controlpoints.NbRows(); r++)
					for (int c = 1; c <= controlpoints.NbColumns(); c++) {
						gp_Pnt pos = controlpoints.Value(r, c);
						cp[r - 1].emplace_back(pos.X(), pos.Y(), pos.Z());
					}

				const TColStd_Array2OfReal* weights = geom_beziersurface->Weights();
				if (weights) {
					vector<vector<double>> w(weights->NbRows());
					for (int r = 1; r <= weights->NbRows(); r++)
					{
						w[r - 1].reserve(weights->NbColumns());
						for (int c = 1; c <= weights->NbColumns(); c++)
							w[r - 1].push_back(weights->Value(r, c));
					}
					faceshape[i].Surface = new BezierSurface(geom_beziersurface->UDegree(), geom_beziersurface->VDegree(), U1, U2, V1, V2, w, cp);
				}
				else
					faceshape[i].Surface = new BezierSurface(geom_beziersurface->UDegree(), geom_beziersurface->VDegree(), U1, U2, V1, V2, cp);
			}
			else if (type == STANDARD_TYPE(Geom_BSplineSurface))
			{
				//dprint(i, "BSplineSurface");
				opencascade::handle<Geom_BSplineSurface> geom_bsplinesurface = Handle(Geom_BSplineSurface)::DownCast(geom_surface);
				TColStd_Array1OfReal uknotsequence = geom_bsplinesurface->UKnotSequence();
				TColStd_Array1OfReal vknotsequence = geom_bsplinesurface->VKnotSequence();
				vector<double> u;
				vector<double> v;
				for (auto itr = uknotsequence.begin(); itr != uknotsequence.end(); itr++)
					u.push_back(*itr);
				for (auto itr = vknotsequence.begin(); itr != vknotsequence.end(); itr++)
					v.push_back(*itr);

				auto isuclosed = geom_bsplinesurface->IsUClosed(), isvclosed = geom_bsplinesurface->IsVClosed();
				TColgp_Array2OfPnt controlpoints = geom_bsplinesurface->Poles();
				vector<vector<Point>> cp(controlpoints.NbRows());
				for (int r = 1; r <= controlpoints.NbRows(); r++)
				{
					cp[r - 1].reserve(controlpoints.NbColumns());
					for (int c = 1; c <= controlpoints.NbColumns(); c++) {
						gp_Pnt pos = controlpoints.Value(r, c);
						cp[r - 1].emplace_back(pos.X(), pos.Y(), pos.Z());
					}
					if (isvclosed) cp[r - 1].push_back(cp[r - 1][0]);
				}
				if (isuclosed) cp.push_back(cp.front());
				const TColStd_Array2OfReal* weights = geom_bsplinesurface->Weights();
				if (weights) {
					vector<vector<double>> w(weights->NbRows());
					for (int r = 1; r <= weights->NbRows(); r++)
					{
						w[r - 1].reserve(weights->NbColumns());
						for (int c = 1; c <= weights->NbColumns(); c++)
							w[r - 1].push_back(weights->Value(r, c));
						if (isvclosed) w[r - 1].push_back(w[r - 1][0]);
					}
					if (isuclosed) w.push_back(w.front());
					faceshape[i].Surface = new BSplineSurface(geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), u, v, w, cp);
				}
				else
					faceshape[i].Surface = new BSplineSurface(geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), u, v, cp);
				//dprint(i, cp.size(), cp[0].size(), geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree());
				//gp_Vec u2, v2, uu, uv, vv;
				//gp_Pnt oringe;
				//geom_bsplinesurface->D2(ut, vt, oringe, u2, v2, uu, vv, uv);
				////double k1, k2;
				////faceshape[i].Surface->PrincipalCurvature(ut, vt, k1, k2);
				//auto v1 = faceshape[i].Surface->PartialDerivativeV(ut, vt);
				//auto u1 = faceshape[i].Surface->PartialDerivativeU(ut, vt);
				//auto uu1 = faceshape[i].Surface->PartialDerivativeUU(ut, vt);
				//auto vv1 = faceshape[i].Surface->PartialDerivativeVV(ut, vt);
				//auto uv1 = faceshape[i].Surface->PartialDerivativeUV(ut, vt);
				//dprint("D1U:", u2.X(),u2.Y(),u2.Z(), u1(0), u1(1), u1(2));
				//dprint("D1V:", v2.X(), v2.Y(), v2.Z(), v1(0), v1(1), v1(2));
				//dprint("D2U:", uu.X(), uu.Y(), uu.Z(), uu1(0), uu1(1), uu1(2));
				//dprint("D2V:", vv.X(), vv.Y(), vv.Z(), vv1(0), vv1(1), vv1(2));
				//dprint("D2UV:", uv.X(), uv.Y(), uv.Z(), uv1(0), uv1(1), uv1(2));
				//dprint(i, geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), geom_bsplinesurface->Weights());
			}
			else if (type == STANDARD_TYPE(Geom_SurfaceOfRevolution))
			{
				//dprint(i, "SurfaceOfRevolution");
				opencascade::handle<Geom_SurfaceOfRevolution> geom_RevolutionSurface = opencascade::handle<Geom_SurfaceOfRevolution>::DownCast(geom_surface);
				auto basecurve = geom_RevolutionSurface->BasisCurve();
				auto typecurve = basecurve->DynamicType();
				if (typecurve == STANDARD_TYPE(Geom_BSplineCurve))
				{
					//获取B样条曲线的信息
					BSplineCurve *B;
					opencascade::handle<Geom_BSplineCurve> geom_bsplinecurve = Handle(Geom_BSplineCurve)::DownCast(basecurve);
					TColStd_Array1OfReal knot = geom_bsplinecurve->KnotSequence();
					TColgp_Array1OfPnt ctrpnts = geom_bsplinecurve->Poles();
					std::vector<double> u;
					std::vector<Point> ctr;
					ctr.reserve(ctrpnts.Size());
					for (auto itr = knot.begin(); itr != knot.end(); itr++)
					{
						u.push_back(*itr);
					}
					for (int r = 1; r <= ctrpnts.Size(); r++)
					{
						gp_Pnt pos = ctrpnts.Value(r);
						ctr.emplace_back(pos.X(), pos.Y(), pos.Z());
					}
					const TColStd_Array1OfReal* weights = geom_bsplinecurve->Weights();
					if (weights)
					{
						vector<double> w(weights->Size());
						for (int r = 1; r <= weights->Size(); r++) w.push_back(weights->Value(r));
						B = new BSplineCurve(geom_bsplinecurve->Degree(), u, w, ctr);
					}
					else
					{
						B = new BSplineCurve(geom_bsplinecurve->Degree(), u, ctr);
					}
					auto axis = geom_RevolutionSurface->Axis();
					auto p = axis.Location();
					auto dir = axis.Direction();
					faceshape[i].Surface = new SurfaceRevolutionType(Point(p.X(), p.Y(), p.Z()), Point(dir.X(), dir.Y(), dir.Z()), B);
					//double ut = 2.5, vt = 5 ,k1, k2;
					//faceshape[i].Surface->PrincipalCurvature(ut, vt, k1, k2);
					//gp_Pnt origin;
					//gp_Vec u2, v2, uu, uv, vv;
					//geom_RevolutionSurface->D2(ut, vt, origin, u2, v2, uu, vv, uv);
					//auto v1 = faceshape[i].Surface->PartialDerivativeV(ut, vt);
					//auto u1 = faceshape[i].Surface->PartialDerivativeU(ut, vt);
					//auto uu1 = faceshape[i].Surface->PartialDerivativeUU(ut, vt);
					//auto vv1 = faceshape[i].Surface->PartialDerivativeVV(ut, vt);
					//auto uv1 = faceshape[i].Surface->PartialDerivativeUV(ut, vt);
					//dprint("D1U:", u2.X(),u2.Y(),u2.Z(), u1(0), u1(1), u1(2));
					//dprint("D1V:", v2.X(), v2.Y(), v2.Z(), v1(0), v1(1), v1(2));
					//dprint("D2U:", uu.X(), uu.Y(), uu.Z(), uu1(0), uu1(1), uu1(2));
					//dprint("D2V:", vv.X(), vv.Y(), vv.Z(), vv1(0), vv1(1), vv1(2));
					//dprint("D2UV:", uv.X(), uv.Y(), uv.Z(), uv1(0), uv1(1), uv1(2));

					//auto p1 = geom_RevolutionSurface->Value(0, 0.1);
					//auto p2 = geom_RevolutionSurface->Value(0, 0.2);
					//auto p3 = geom_RevolutionSurface->Value(0, 0.4);
					//auto p4 = geom_RevolutionSurface->Value(0, 0.6);
					//int y = 0;
				}
				else
				{
					dprint(typecurve);
					dprint("Cannot recognize the type of the curve of this revolution surface!");
					system("pause");
				}
			}
			else
			{
				dprint(type);
				dprint("Cannot recognize the type of this surface!");
				system("pause");
			}
		}
		for (int i = 0; i < edgeshape.size(); i++)
		{
			auto &edge = edgeshape[i];
			//if (edge.main_face >= 6 && edge.main_face <= 28) continue;
			if (!edge.if_exisited) continue;
			auto &facer = faceshape[edge.main_face].Surface;
			auto &pnts = edge.parameters;
			Point4 UV = facer->Getbounds();
			GeneralMathMethod::DataSet(UV, 10e-8, pnts);
		}
		dprint("face type done!");
	}

	void OccReader::C0_Feature()//通过计算相交线的两个面的法向内积来判断是否光滑
	{
		auto& edge = globalmodel.edgeshape;
		auto& face = globalmodel.faceshape;

		int discrete_num = 20;
		for (int i = 0; i < edge.size(); i++)
		{
			auto& aedge = edge[i];
			if (aedge.if_C0 || aedge.if_visited) continue;

			//将边界曲线设置为C0
			if (aedge.reversed_edge == -1)
			{
				aedge.if_C0 = true;
				continue;
			}

			double C0 = 0.94;     //两个法向量内积的阈值
			int single_flag = 0;
			auto& face1 = face[aedge.main_face].Surface;
			auto& face2 = face[aedge.secondary_face].Surface;

			//得到边的等距离散点
			Standard_Real first = 0, last = 0;
			Matrix2Xd pnts1(2, discrete_num), pnts2(2, discrete_num);
			double step;
			gp_Pnt2d uv;
			Handle_Geom2d_Curve thePCurve1 = BRep_Tool::CurveOnSurface(aedge.edge, face[aedge.main_face].face, first, last);
			step = (last - first) / (discrete_num - 1);
			if (aedge.edge.Orientation() == TopAbs_FORWARD)
			{
				for (int i = 0; i < discrete_num - 1; i++)
				{
					uv = thePCurve1->Value(first + i * step);
					pnts1.col(i) << uv.X(), uv.Y();
				}
				uv = thePCurve1->Value(last);
				pnts1.col(discrete_num - 1) << uv.X(), uv.Y();
			}
			else
			{
				for (int i = 0; i < discrete_num - 1; i++)
				{
					uv = thePCurve1->Value(last - i * step);
					pnts1.col(i) << uv.X(), uv.Y();
				}
				uv = thePCurve1->Value(first);
				pnts1.col(discrete_num - 1) << uv.X(), uv.Y();
			}
			Point4 UV = face1->Getbounds();
			GeneralMathMethod::DataSet(UV, 10e-8, pnts1);
			Handle_Geom2d_Curve thePCurve2 = BRep_Tool::CurveOnSurface(edge[aedge.reversed_edge].edge, face[aedge.secondary_face].face, first, last);
			step = (last - first) / (discrete_num - 1);
			if (edge[aedge.reversed_edge].edge.Orientation() == TopAbs_FORWARD)
			{
				for (int i = 0; i < discrete_num - 1; i++)
				{
					uv = thePCurve2->Value(first + i * step);
					pnts2.col(i) << uv.X(), uv.Y();
				}
				uv = thePCurve2->Value(last);
				pnts2.col(discrete_num - 1) << uv.X(), uv.Y();
			}
			else
			{
				for (int i = 0; i < discrete_num - 1; i++)
				{
					uv = thePCurve2->Value(last - i * step);
					pnts2.col(i) << uv.X(), uv.Y();
				}
				uv = thePCurve2->Value(first);
				pnts2.col(discrete_num - 1) << uv.X(), uv.Y();
			}
			UV = face2->Getbounds();
			GeneralMathMethod::DataSet(UV, 10e-8, pnts2);

			//计算采样点处的法向量内积
			Vector3d n1, n2;
			for (int j = 0; j < pnts1.cols(); j++)
			{
				if (j - single_flag > pnts1.cols() * 0.4) break;
				double u = pnts1(0, j), v = pnts1(1, j);
				if (!face1->NormalValue(u, v, n1)) continue;
				if ((face[aedge.main_face].face).Orientation() == TopAbs_REVERSED) n1 = -n1;
				u = pnts2(0, pnts1.cols() - 1 - j), v = pnts2(1, pnts1.cols() - 1 - j);
				if (!face2->NormalValue(u, v, n2)) continue;
				if ((face[aedge.secondary_face].face).Orientation() == TopAbs_REVERSED) n2 = -n2;
				if (n1.dot(n2) < C0) single_flag++;
			}

			//若有超过60%采样点是single_flag，则将该曲线标记为flag
			if (single_flag > pnts1.cols() * 0.6)
			{
				aedge.if_C0 = true;
				edge[aedge.reversed_edge].if_C0 = true;
			}
		}
		Reset_feature();
		dprint("C0 feature done!");
	}

	void OccReader::Reset_feature()
	{
		auto& edge = globalmodel.edgeshape;
		auto& face = globalmodel.faceshape;

		//计算每个面的bounding box长度
		double bbox = 0;
		for (int i = 0; i < face.size(); i++)
		{
			Bnd_Box aBound;
			BRepBndLib::Add(face[i].face, aBound);
			Standard_Real xmin, ymin, zmin, xmax, ymax, zmax;
			aBound.Get(xmin, ymin, zmin, xmax, ymax, zmax);
			double boundingbox = sqrt((xmin - xmax) * (xmin - xmax) + (ymin - ymax) * (ymin - ymax) +
				(zmin - zmax) * (zmin - zmax));
			face[i].bbox = boundingbox;
			bbox += boundingbox;
		}
		bbox /= face.size();

		//set pre edge
		for (int i = 0; i < face.size(); i++)
		{
			auto &wires = face[i].wires;
			if (wires.empty() || !face[i].if_exisited)
			{
				continue;
			}
			for (int m = 0; m < wires.size(); m++)
			{
				auto &edges = wires[m];
				int pre = edges.front();
				for (int j = 1; j < edges.size(); j++)
				{
					edge[edges[j]].prev_edge = pre;
					edge[pre].next_edge = edges[j];
					pre = edges[j];
				}
				edge[edges.front()].prev_edge = pre;
				edge[pre].next_edge = edges.front();
			}
		}

		//处理孤立C0特征边
		for (int i = 0; i < edge.size(); i++)
		{
			auto& aedge = edge[i];
			if (!aedge.if_C0 || aedge.reversed_edge < -0.5) continue;
			if (aedge.length > 0.2 * bbox) continue;
			if (edge[aedge.prev_edge].if_C0 || edge[aedge.next_edge].if_C0) continue;
			if (edge[edge[edge[aedge.next_edge].reversed_edge].next_edge].if_C0) continue;
			auto &reaedge = edge[aedge.reversed_edge];
			if (edge[reaedge.prev_edge].if_C0 || edge[reaedge.next_edge].if_C0) continue;
			if (edge[edge[edge[reaedge.next_edge].reversed_edge].next_edge].if_C0) continue;
			aedge.if_C0 = false;
			reaedge.if_C0 = false;
		}
	}

	//找到曲率变化较大的曲线并标记（类似机翼前端曲线，后期需要加密或做各向异性、四边形）
	void OccReader::Curvature_Feature()
	{
		auto& edge = globalmodel.edgeshape;
		auto& face = globalmodel.faceshape;

		//边离散点数目，和沿(stepu, stepv)方向的踩点数
		int discrete_num = 20, pnt_num = 15;

		for (int i = 0; i < edge.size(); i++)
		{
			//if (i != 255) continue;
			auto& aedge = edge[i];
			if (!aedge.if_curvature) continue;

			//C0特征边不再考虑是否为曲率特征边
			if (aedge.if_C0)
			{
				aedge.if_curvature = false;
				continue;
			}

			//得到边的等距离散点
			Matrix2Xd pnts(2, discrete_num);
			Standard_Real first = 0, last = 0;
			double step0;
			gp_Pnt2d uv;
			Handle_Geom2d_Curve thePCurve = BRep_Tool::CurveOnSurface(aedge.edge, face[aedge.main_face].face, first, last);
			step0 = (last - first) / (discrete_num - 1);
			if (aedge.edge.Orientation() == TopAbs_FORWARD)
			{
				for (int i = 0; i < discrete_num - 1; i++)
				{
					uv = thePCurve->Value(first + i * step0);
					pnts.col(i) << uv.X(), uv.Y();
				}
				uv = thePCurve->Value(last);
				pnts.col(discrete_num - 1) << uv.X(), uv.Y();
			}
			else
			{
				for (int i = 0; i < discrete_num - 1; i++)
				{
					uv = thePCurve->Value(last - i * step0);
					pnts.col(i) << uv.X(), uv.Y();
				}
				uv = thePCurve->Value(first);
				pnts.col(discrete_num - 1) << uv.X(), uv.Y();
			}
			auto & facer = face[aedge.main_face].Surface;
			Point4 UV = facer->Getbounds();
			GeneralMathMethod::DataSet(UV, 10e-8, pnts);

			//获取面的参数与边界和初始步长step
			double u1 = pnts(0, 0), v1 = pnts(1, 0), u2 = pnts(0, discrete_num - 1), v2 = pnts(1, discrete_num - 1);
			Eigen::Vector2d step(v1 - v2, u2 - u1);
			step /= (discrete_num - 1) * 4;
			if ((face[aedge.main_face].face).Orientation() == TopAbs_REVERSED) step = -step;
			double base_stepu = step(0), base_stepv = step(1);
			double umin = UV(0), umax = UV(1), vmin = UV(2), vmax = UV(3);

			//判断单个点是否满足曲率特征
			int single_flag = 0;
			for (int j = 1; j < pnts.cols() - 1; j++)
			{
				if (j - 1 - single_flag > (pnts.cols() - 2) * 0.4) break;

				//更新步长，使得所有采样点均落在参数域内
				double stepu = base_stepu, stepv = base_stepv;
				int count = 0;
				while (true)
				{
					u1 = pnts(0, j);
					v1 = pnts(1, j);
					for (int m = 0; m < pnt_num; m++)
					{
						u1 += stepu;
						v1 += stepv;
						if (m == pnt_num - 1 || u1 < umin || u1 > umax || v1 < vmin || v1 > vmax)
						{
							count = m;
							break;
						}
					}
					if (count) break;
					else
					{
						stepu *= 0.5;
						stepv *= 0.5;
					}
				}
				stepu = (stepu * count) / pnt_num;
				stepv = (stepv * count) / pnt_num;

				//沿(stepu, stepv)方向按固定步长采样并计算相应的法曲率
				int maxtime = 3;
				std::vector<double> curvature;
				double K, decrease_ratio = 0.8;
				Vector3d start_normal, end_normal;
				bool is_normal_ok = facer->NormalValue(pnts(0, j), v1 = pnts(1, j), start_normal);
				while (maxtime)
				{
					u1 = pnts(0, j);
					v1 = pnts(1, j);
					curvature.clear();
					for (int m = 0; m < pnt_num; m++)
					{
						if (!facer->NormalCurvature(u1, v1, stepu, stepv, K)) K = -1;
						curvature.push_back(K);
						if (is_normal_ok && m == pnt_num - 1)
							is_normal_ok = facer->NormalValue(u1, v1, end_normal);
						u1 += stepu;
						v1 += stepv;
					}
					if (is_normal_ok)
					{
						if (start_normal.dot(end_normal) < 0.85 && If_decline(curvature))
						{
							single_flag++;
							break;
						}
					}
					else if (If_decline(curvature))
					{
						single_flag++;
						break;
					}
					stepu *= decrease_ratio;
					stepv *= decrease_ratio;
					maxtime--;
				}
			}

			//如果有超过60%的采样点是曲率特征点，则将该边标记为曲率特征边
			if (single_flag <= (pnts.cols() - 2) * 0.6)
			{
				aedge.if_curvature = false;
				edge[aedge.reversed_edge].if_curvature = false;
			}
		}

		for (int i = 0; i < edge.size(); i++)
		{
			auto& aedge = edge[i];
			if (!aedge.if_curvature) continue;
			//dprint("edge", i);

			//目前只考虑在无洞面片内做offset
			if (face[aedge.main_face].wires.size() > 1) continue;

			//曲率特征边有一定长度才做offset
			if ((aedge.length) > 0.1 * face[aedge.main_face].bbox)
				face[aedge.main_face].if_quad = true;
		}

		dprint("Curvature feature done!");
	}

	bool OccReader::If_decline(vector<double>& curvature)
	{
		int n = (int)(curvature.size() * 0.5);

		//前半段的曲率平均值比后半段的平均值要大一定倍数
		double front_curv = 0, back_curv = 0;
		for (int i = 0; i < n; i++)
		{
			front_curv += curvature[i];
		}
		for (int i = curvature.size() - n; i < curvature.size(); i++)
		{
			back_curv += curvature[i];
		}
		double ration;
		if (front_curv > 1000 * n) ration = 1.2;
		else if (front_curv > n) ration = 1.5;
		else ration = 2;
		if (front_curv < ration * back_curv) return false;

		//连续下降一定次数或通过Cox stuart趋势性检验
		front_curv = curvature.front();
		int T = 0;
		for (int i = 1; i < curvature.size(); i++)
		{
			if (curvature[i] >= front_curv) break;
			front_curv = curvature[i];
			T++;
		}
		if (T < (curvature.size() - 1) * 0.25)
		{
			//Cox stuart趋势性检验
			double p0 = 0.2;
			T = 0;
			if (curvature.size() % 2)
			{
				for (int i = 0; i < n; i++)
				{
					if (curvature[i] > curvature[n + i + 1]) T++;
				}
			}
			else
			{
				for (int i = 0; i < n; i++)
				{
					if (curvature[i] > curvature[n + i]) T++;
				}
			}
			if (2 * T <= n) return false;
			if (GeneralMathMethod::Binomial(n, 0.5, T) > p0) return false;   //拒绝原假设，认为不存在下降趋势
		}

		//存在一定数量的有效数据
		std::sort(curvature.begin(), curvature.end());
		if (curvature[(int)(curvature.size() * 0.8)] < 10e-4) return false;
		else return true;
	}

	void OccReader::Set_Offset_Grid()
	{
		if (offset_quad_num <= 0)
		{
			dprint("the offset quad number must be positive!");
			system("pause");
		}
		if (offset_increase_ratio < 1)
		{
			dprint("the offset increase ratio must more than 1!");
			system("pause");
		}

		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;
		int discrete_num = 10;

		//determine each edge's discrete points 
		vector<vector<int>> edge_prev_info(edgeshape.size());
		vector<vector<int>> edge_next_info(edgeshape.size());
		for (int i = 0; i < faceshape.size(); i++)
		{
			if (!faceshape[i].if_quad) continue;
			auto &edges = faceshape[i].wires.front();

			//find curvature edge
			int now, prev, next;
			for (int j = 0; j < edges.size(); j++)
			{
				if (edgeshape[edges[j]].if_curvature)
				{
					now = j;
					break;
				}
			}
			if (now) prev = edges[now - 1];
			else prev = edges.back();
			if (now < edges.size() - 1) next = edges[now + 1];
			else next = edges.front();
			auto &curv = edgeshape[edges[now]];

			//determain the offset height
			int prev_id = -1, next_id = -1;
			Vector2d p0, dir1, dir2;
			double ave_len = 0;
			p0 = curv.parameters.col(0);
			for (int j = 1; j < curv.parameters.cols(); j++)
			{
				ave_len += (curv.parameters.col(j) - p0).norm();
				p0 = curv.parameters.col(j);
			}
			ave_len /= curv.parameters.cols() - 1;
			double offset_height;
			//int offset_quad_num = curv.quad_num;
			//double offset_initial_ratio = curv.initial_ratio;
			//double offset_increase_ratio = curv.increase_ratio;
			//注意下行代码更改了邻边的increase_ratio，在一个面做多个四边形会有问题
			//edgeshape[prev].increase_ratio = edgeshape[next].increase_ratio = offset_increase_ratio;

			if (abs(offset_increase_ratio - 1) < DBL_EPSILON)
				offset_height = ave_len * offset_initial_ratio * offset_quad_num;
			else
				offset_height = ave_len * offset_initial_ratio *(1 - pow(offset_increase_ratio, offset_quad_num)) / (1 - offset_increase_ratio);

			//whether the two nearby edges can construct the expect offet
			bool flag = false;
			auto &prev_para = edgeshape[prev].parameters;
			auto &next_para = edgeshape[next].parameters;
			GeneralMathMethod::Find_Span(curv.parameters, prev_para, false, prev_id, discrete_num, offset_height);
			if (prev_id >= 0)
			{
				GeneralMathMethod::Find_Span(curv.parameters, next_para, true, next_id, discrete_num, offset_height);
				if (next_id >= 0) flag = true;
			}

			//compute offset grid to each nearby edge
			int prev_quad_num, next_quad_num;
			bool prev_end_flag, next_end_flag;
			if (flag)
			{
				prev_end_flag = !prev_id;
				next_end_flag = next_id == (next_para.cols() - 1) * discrete_num;
				if (prev_end_flag && next_end_flag)
				{
					prev_id++;
					next_id--;
				}
				prev_quad_num = next_quad_num = std::max(4, offset_quad_num);
			}
			else
			{
				p0 = curv.parameters.col(0);
				dir1 = (p0 - curv.parameters.col(1)).normalized();
				dir2 = prev_para.col(0) - p0;
				double prev_height = abs(dir1[0] * dir2[1] - dir1[1] * dir2[0]);
				p0 = curv.parameters.col(curv.parameters.cols() - 1);
				dir1 = (p0 - curv.parameters.col(curv.parameters.cols() - 2)).normalized();
				dir2 = next_para.col(next_para.cols() - 1) - p0;
				double next_height = abs(dir1[0] * dir2[1] - dir1[1] * dir2[0]);
				if (prev_height > next_height)
				{
					if (abs(offset_increase_ratio - 1) < DBL_EPSILON)
						prev_quad_num = next_quad_num = std::max(4, (int)(next_height / (ave_len * offset_initial_ratio) + 0.5));
					else
						prev_quad_num = next_quad_num = std::max(4, (int)(log(1 - (next_height * (1 - offset_increase_ratio) / (ave_len * offset_initial_ratio))) / log(offset_increase_ratio) + 0.5));
					GeneralMathMethod::Find_Span(curv.parameters, prev_para, false, prev_id, discrete_num, next_height);
					prev_end_flag = !prev_id;
					if (!prev_end_flag) next_id = (next_para.cols() - 1) * discrete_num;
					else
					{
						prev_id = 1;
						next_id = (next_para.cols() - 1) * discrete_num - 1;
					}
				}
				else
				{
					if (abs(offset_increase_ratio - 1) < DBL_EPSILON)
						prev_quad_num = next_quad_num = std::max(4, (int)(prev_height / (ave_len * offset_initial_ratio) + 0.5));
					else
						prev_quad_num = next_quad_num = std::max(4, (int)(log(1 - (prev_height * (1 - offset_increase_ratio) / (ave_len * offset_initial_ratio))) / log(offset_increase_ratio) + 0.5));
					GeneralMathMethod::Find_Span(curv.parameters, next_para, true, next_id, discrete_num, prev_height);
					next_end_flag = next_id == (next_para.cols() - 1) * discrete_num;
					if (!next_end_flag) prev_id = 0;
					else
					{
						next_id = (next_para.cols() - 1) * discrete_num - 1;
						prev_id = 1;
					}
				}
			}

			//save the offset grid information 
			auto &prev_info = edge_prev_info[prev];
			if (prev_info.empty())
			{
				prev_info = { prev_id, prev_quad_num };
				if (edgeshape[prev].reversed_edge != -1)
					edge_next_info[edgeshape[prev].reversed_edge] = { (int)((prev_para.cols() - 1) * discrete_num - prev_id), prev_quad_num };
			}
			else
			{
				prev_id = (prev_id + prev_info.front()) / 2;
				prev_quad_num = std::min(prev_quad_num, prev_info.back());
				prev_info = { prev_id, prev_quad_num };
				if (edgeshape[prev].reversed_edge != -1)
					edge_next_info[edgeshape[prev].reversed_edge] = { (int)((prev_para.cols() - 1) * discrete_num - prev_id), prev_quad_num };
			}
			auto &next_info = edge_next_info[next];
			if (next_info.empty())
			{
				next_info = { next_id, next_quad_num };
				if (edgeshape[next].reversed_edge != -1)
					edge_prev_info[edgeshape[next].reversed_edge] = { (int)((next_para.cols() - 1) * discrete_num - next_id), next_quad_num };
			}
			else
			{
				next_id = (next_id + next_info.front()) / 2;
				next_quad_num = std::min(next_quad_num, next_info.back());
				next_info = { next_id,  next_quad_num };
				if (edgeshape[next].reversed_edge != -1)
					edge_prev_info[edgeshape[next].reversed_edge] = { (int)((next_para.cols() - 1) * discrete_num - next_id), next_quad_num };
			}
		}

		//re-discrete point in edge which is nearby curvature edge
		for (int i = 0; i < edgeshape.size(); i++)
		{
			auto &prev_info = edge_prev_info[i];
			auto &next_info = edge_next_info[i];
			if (prev_info.empty() && next_info.empty()) continue;
			if (!prev_info.empty() && !next_info.empty())
			{
				int t1 = prev_info.front(), t2 = next_info.front();
				if (t1 >= t2)
				{
					Re_discrete(edgeshape[i], t1, discrete_num, prev_info.back(), false);
					Re_discrete(edgeshape[i], t2, discrete_num, next_info.back(), true);
				}
				else
				{
					Re_discrete(edgeshape[i], (t1 + t2) / 2, discrete_num, prev_info.back(), false);
					Re_discrete(edgeshape[i], (t1 + t2) / 2, discrete_num, next_info.back(), true);
				}
				continue;
			}
			if (!prev_info.empty())
			{
				int t1 = prev_info.front();
				Re_discrete(edgeshape[i], t1, discrete_num, prev_info.back(), false);
			}
			else
			{
				int t2 = next_info.front();
				Re_discrete(edgeshape[i], t2, discrete_num, next_info.back(), true);
			}
		}

		//set offset quad num for each surface which has curvture feature edge
		for (int i = 0; i < faceshape.size(); i++)
		{
			if (!faceshape[i].if_quad) continue;
			auto &edges = faceshape[i].wires.front();

			//find curvature edge
			int now, prev, next;
			for (int j = 0; j < edges.size(); j++)
			{
				if (edgeshape[edges[j]].if_curvature)
				{
					now = j;
					break;
				}
			}
			if (now) prev = edges[now - 1];
			else prev = edges.back();
			if (now < edges.size() - 1) next = edges[now + 1];
			else next = edges.front();
			faceshape[i].quad_num = std::min(edge_prev_info[prev].back(), edge_next_info[next].back());
		}
		for (int i = 0; i < faceshape.size(); i++)
		{
			if (!faceshape[i].if_quad) continue;
			auto &edges = faceshape[i].wires.front();

			//find curvature edge
			int now, prev, next;
			for (int j = 0; j < edges.size(); j++)
			{
				if (edgeshape[edges[j]].if_curvature)
				{
					now = j;
					break;
				}
			}
			if (now) prev = edges[now - 1];
			if (edgeshape[edges[now]].reversed_edge < -0.5) continue;
			faceshape[i].quad_num = std::min(faceshape[i].quad_num, faceshape[edgeshape[edges[now]].secondary_face].quad_num);
		}

		for (int i = 0; i < edgeshape.size(); i++)
		{
			if (edgeshape[i].reversed_edge < 0) continue;
			int a = edgeshape[i].parameters.cols();
			int b = edgeshape[edgeshape[i].reversed_edge].parameters.cols();
			if (a != b)
			{
				dprint("big bug!!!!!", i, edgeshape[i].reversed_edge);
				system("pause");
			}
		}
	}

	void OccReader::Re_discrete(ShapeEdge &edge, int id, int discrete_num, int quad_num, bool direction)
	{
		if (!quad_num) return;
		int num = edge.parameters.cols(), id1 = id / discrete_num, id2 = id % discrete_num, left_segment;
		Eigen::Matrix2Xd para;
		if (direction) left_segment = num - id1 - 1;
		else
		{
			if (id2) left_segment = id1 + 1;
			else left_segment = id1;
		}

		para.resize(2, quad_num + 1 + left_segment);
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		gp_Pnt2d uv;
		Standard_Real first = 0, last = 0;
		Handle_Geom2d_Curve thePCurve = BRep_Tool::CurveOnSurface(edge.edge, faceshape[edge.main_face].face, first, last);
		double t, t1, initial_step;
		//double offset_increase_ratio = edge.increase_ratio;

		if (edge.edge.Orientation() == TopAbs_FORWARD)
		{
			t = (last - first) * id / ((num - 1.0) * discrete_num) + first;
			if (direction)
			{
				if (abs(offset_increase_ratio - 1) < DBL_EPSILON)
					initial_step = (t - first) / quad_num;
				else
					initial_step = (t - first) * (1 - offset_increase_ratio) / (1 - std::pow(offset_increase_ratio, quad_num));
				para.col(0) = edge.parameters.col(0);
				t1 = initial_step;
				for (int i = 1; i <= quad_num; i++)
				{
					uv = thePCurve->Value(t1);
					para.col(i) << uv.X(), uv.Y();
					t1 += initial_step * std::pow(offset_increase_ratio, i);
				}
				for (int i = 1; i <= left_segment; i++)
				{
					para.col(i + quad_num) = edge.parameters.col(id1 + i);
				}
			}
			else
			{
				if (abs(offset_increase_ratio - 1) < DBL_EPSILON)
					initial_step = (last - t) / quad_num;
				else
					initial_step = (last - t) * (1 - offset_increase_ratio) / (1 - std::pow(offset_increase_ratio, quad_num));
				t1 = last - initial_step;
				para.col(para.cols() - 1) = edge.parameters.col(num - 1);
				for (int i = 1; i <= quad_num; i++)
				{
					uv = thePCurve->Value(t1);
					para.col(para.cols() - 1 - i) << uv.X(), uv.Y();
					t1 -= initial_step * std::pow(offset_increase_ratio, i);
				}
				for (int i = 0; i < left_segment; i++)
				{
					para.col(i) = edge.parameters.col(i);
				}
			}
		}
		else
		{
			t = last - (last - first) * id / ((num - 1.0) * discrete_num);
			if (direction)
			{
				if (abs(offset_increase_ratio - 1) < DBL_EPSILON)
					initial_step = (last - t) / quad_num;
				else
					initial_step = (last - t) * (1 - offset_increase_ratio) / (1 - std::pow(offset_increase_ratio, quad_num));
				para.col(0) = edge.parameters.col(0);
				t1 = last - initial_step;
				for (int i = 1; i <= quad_num; i++)
				{
					uv = thePCurve->Value(t1);
					para.col(i) << uv.X(), uv.Y();
					t1 -= initial_step * std::pow(offset_increase_ratio, i);
				}
				for (int i = 1; i <= left_segment; i++)
				{
					para.col(i + quad_num) = edge.parameters.col(id1 + i);
				}
			}
			else
			{
				if (abs(offset_increase_ratio - 1) < DBL_EPSILON)
					initial_step = (t - first) / quad_num;
				else
					initial_step = (t - first) * (1 - offset_increase_ratio) / (1 - std::pow(offset_increase_ratio, quad_num));
				t1 = initial_step;
				para.col(para.cols() - 1) = edge.parameters.col(num - 1);
				for (int i = 1; i <= quad_num; i++)
				{
					uv = thePCurve->Value(t1);
					para.col(para.cols() - 1 - i) << uv.X(), uv.Y();
					t1 += initial_step * std::pow(offset_increase_ratio, i);
				}
				for (int i = 0; i < left_segment; i++)
				{
					para.col(i) = edge.parameters.col(i);
				}
			}
		}
		edge.parameters.setZero();
		edge.parameters.resize(2, para.cols());
		edge.parameters = para;
	}

	void OccReader::Offset_lines(Matrix2Xd &parameters, vector<Matrix2Xd> &offset_pnts, int beginid, int pntnum, int quadnum)
	{
		offset_pnts.reserve(quadnum);
		int totalpnts = parameters.cols();
		int endid = (beginid + pntnum - 1) >= totalpnts ? 0 : (beginid + pntnum - 1);
		if (!beginid) beginid = totalpnts;
		for (int i = 1; i <= quadnum; i++)
		{
			Matrix2Xd singleline(2, pntnum - 2);
			auto start_pnt = parameters.col(beginid - i), end_pnt = parameters.col(endid + i);
			for (int j = 1; j <= pntnum - 2; j++)
			{
				singleline.col(j - 1) = (j * end_pnt + (pntnum - 1 - j)*start_pnt) / (pntnum - 1);
			}
			offset_pnts.push_back(singleline);
		}
	}

	int OccReader::EndId(vector<ShapeEdge> &edgeshape, int edge_id)
	{
		int id = edgeshape[edge_id].next_reversed_edge;
		if (id == -1)
		{
			return edgeshape[edge_id].end_id;
		}
		else
		{
			return EndId(edgeshape, id);
		}
	}

	void OccReader::SetTriFeature()
	{
		TriMesh& model_mesh = globalmodel.initial_trimesh;
		for (auto te : model_mesh.edges())
		{
			if (!model_mesh.data(te).flag1 && !model_mesh.data(te).flag2 && !model_mesh.is_boundary(te)) continue;
			if (model_mesh.is_boundary(te))
				model_mesh.data(te).flag1 = true;
			model_mesh.data(te.v0()).set_vertflag(true);
			model_mesh.data(te.v1()).set_vertflag(true);
		}

		dprint("Reset Feature Done!");
	}

	void OccReader::SetPolyFeature()
	{
		Mesh& model_mesh = globalmodel.initial_polymesh;
		for (auto te : model_mesh.edges())
		{
			if (!model_mesh.data(te).flag1 && !model_mesh.is_boundary(te)) continue;
			if (model_mesh.is_boundary(te))
				model_mesh.data(te).flag1 = true;
			model_mesh.data(te.v0()).set_vertflag(true);
			model_mesh.data(te.v1()).set_vertflag(true);
		}
		//���ı��ε�4�����Ϊ�����
		for (auto f : model_mesh.faces())
		{
			if (f.valence() < 4) continue;
			for (auto fv : model_mesh.fv_range(f))
			{
				model_mesh.data(fv).set_vertflag(true);
			}
			for (auto fe : model_mesh.fe_range(f))
			{
				model_mesh.data(fe).flag2 = true;
			}
		}

		dprint("Reset Feature Done!");
	}

}


#include "OccReader.h"

namespace CADMesher
{
	void OccReader::Set_TriMesh()
	{
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;
		Surface_TriMeshes.resize(faceshape.size());
		for (int i = 0; i < faceshape.size(); i++)
		{
			dprint(i,"*******");
			//if (i != 24) continue;
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
			auto &Surface = faceshape[i].Surface;

			//calculate the target edge length with respect to the error between mesh and surface		
			//double target_h = 0.0;
			//double pntgap = (double)aMesh.n_vertices() / std::min(200.0, (double)aMesh.n_vertices()), pntid=0;
			//std::vector<double> errorvector;
			//int testnum = 0, sumnum = 0;
			//while (pntid <= (int)aMesh.n_vertices()-1)
			//{
			//	auto p = aMesh.point(aMesh.vertex_handle((int)pntid));
			//	errorvector.push_back((Surface->PartialDerivativeUU(p[0], p[1])).norm() + 2 * (Surface->PartialDerivativeUV(p[0], p[1])).norm() + (Surface->PartialDerivativeVV(p[0], p[1])).norm());
			//	pntid += pntgap;
			//	testnum++;
			//	double k1, k2;
			//}
			//std::sort(errorvector.begin(), errorvector.end());
			//for (int j = (int)(0.2*testnum); j < (int)(0.9*testnum); j++)
			//{
			//	sumnum++;
			//	target_h += std::sqrt(2*epsratio / errorvector[j]);
			//}
			//target_h /= sumnum;
			//dprint(aMesh.n_vertices(), target_h, x_step);
			//if (target_h < 2.5 * x_step && target_h > 0.4 * x_step)
			//{
			//	triangulate(all_pnts, bnd, sqrt(3) * target_h * target_h *0.25, aMesh);
			//	for (auto tv : aMesh.vertices())
			//	{
			//		if (tv.idx() < pointsnumber)
			//		{
			//			aMesh.set_point(tv, TriMesh::Point(boundary(0, tv.idx()), boundary(1, tv.idx()), 0));
			//			continue;
			//		}
			//		auto pos = aMesh.point(tv);
			//		aMesh.set_point(tv, TriMesh::Point(pos[0], pos[1] * ra_inv, 0));
			//	}
			//}

			dprint("initial vertices:", aMesh.n_vertices());
			//Remesh in domain		
			//Riemannremesh Remesh(Surface, &aMesh);
			//Remesh.remesh();
			dprint("domain remesh done!");

			//double k1, k2;
			//for (auto v : aMesh.vertices())
			//{
			//	auto p = aMesh.point(v);
			//	Surface->PrincipalCurvature(p[0], p[1], k1, k2);
			//	//dprint("v.idx()", k1, k2);
			//	if (isnan(k1*k2))
			//	{
			//		dprint("bug", v.idx());
			//	}
			//	aMesh.data(v).GaussCurvature = std::max(std::fabs(k1), std::fabs(k2));
			//}
			dprint("GaussCurvature compute done!");
			if (!OpenMesh::IO::write_mesh(aMesh, "one.obj"))
			{
				std::cerr << "fail";
			}
			for (auto tv : aMesh.vertices())
			{
				auto pos = aMesh.point(tv);
				auto v = asurface->Value(pos[0], pos[1]);
				aMesh.set_point(tv, TriMesh::Point(v.X(), v.Y(), v.Z()));
			}
			if (!OpenMesh::IO::write_mesh(aMesh, "two.obj"))
			{
				std::cerr << "fail";
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
						for (int m = id; m < id + cols-1; m++)
						{
							auto e = aMesh.edge_handle(aMesh.find_halfedge(aMesh.vertex_handle(m), aMesh.vertex_handle(m+1)));
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

		Surface_PolyMeshes.resize(faceshape.size());
		for (int i = 0; i < faceshape.size(); i++)
		{
			dprint(i);
			//if (i != 69) continue;
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
			all_pnts.block(1, 0, 1, all_pnts.cols()) *= if_reverse;

			int startid = 0, endid = -1, beginid;
			Matrix2Xd offsetline;
			int wireoffset;
			for (int j = 0; j < wires.size(); j++)
			{
				beginid = startid;
				auto &edges = wires[j];
				for (int k = 0; k < edges.size(); k++)
				{
					auto &aedge = edgeshape[edges[k]];
					if (aedge.if_curvature)
					{
						int id1, id2;
						if (k == edges.size() - 1) endid = beginid;
						else endid = startid + (aedge.parameters).cols()-1;
						wireoffset = bnd[j].cols();
						if (startid == beginid)
						{
							if (startid) id1 = beginid + wireoffset - 1;
							else id1 = pointsnumber - 1;
						}
						else id1 = startid - 1;
						id2 = endid + 1;
						offsetline = Subdomain(aedge.parameters, all_pnts.col(id1), all_pnts.col(id2));
						Matrix2Xi bnd_offset(2, wireoffset - 2);
						bnd_offset.block(0, 0, 2, wireoffset - 3) = bnd[j].block(0, 0, 2, wireoffset - 3);
						bnd_offset.col(wireoffset - 3) << beginid + wireoffset - 3, beginid;
						bnd[j] = bnd_offset;
						for (int m = j + 1; m < wires.size(); m++)
						{
							bnd[m] -= 2 * Matrix2i::Ones(2, bnd[m].cols());
						}
						break;
					}
					startid += (aedge.parameters).cols() - 1;
				}
				if (endid != -1) break;
			}
			if (endid == -1)
			{
				all_pnts.block(1, 0, 1, all_pnts.cols()) *= ra * if_reverse;
				triangulate(all_pnts, bnd, sqrt(3) * x_step * x_step / 4 * mu, aMesh);

				double ra_inv = 1.0 / ra;
				for (auto tv : aMesh.vertices())
				{
					if (tv.idx() < pointsnumber)
					{
						aMesh.set_point(tv, TriMesh::Point(boundary(0, tv.idx()), boundary(1, tv.idx()), 0));
						continue;
					}
					auto p = aMesh.point(tv);
					aMesh.set_point(tv, Mesh::Point{ p[0],p[1] * ra_inv,0 });
				}

				//calculate the target edge length with respect to the error between mesh and surface		
				//double errorbounds = 0.0, testnum = std::min(200, (int)aMesh.n_vertices());
				//for (int j = 0; j < testnum; j++)
				//{
				//	auto p = aMesh.point(aMesh.vertex_handle(j));
				//	errorbounds += (Surface->PartialDerivativeUU(p[0], p[1])).norm() + 2 * (Surface->PartialDerivativeUV(p[0], p[1])).norm() + (Surface->PartialDerivativeVV(p[0], p[1])).norm();
				//}
				//errorbounds /= 2 * testnum;
				//errorbounds = std::max(errorbounds, 0.00001);
				//double target_h = std::sqrt(epsratio / errorbounds);
				//dprint(aMesh.n_vertices(), target_h, x_step);
				//if (target_h < x_step)
				//{
				//	triangulate(all_pnts, bnd, sqrt(3) * target_h * target_h *0.25, aMesh);
				//	for (auto tv : aMesh.vertices())
				//	{
				//		if (tv.idx() < pointsnumber)
				//		{
				//			aMesh.set_point(tv, TriMesh::Point(boundary(0, tv.idx()), boundary(1, tv.idx()), 0));
				//			continue;
				//		}
				//		auto pos = aMesh.point(tv);
				//		aMesh.set_point(tv, TriMesh::Point(pos[0], pos[1] * ra_inv, 0));
				//	}
				//}

				//Remesh in domain
				Riemannremesh Remesh(Surface, &aMesh);
				Remesh.remesh();

				newmesh = Mesh(aMesh);
			}
			else
			{
				Matrix2Xd newall_pnts(2, pointsnumber-2);
				int offsetpntnum = offsetline.cols();
				if (endid < startid)
				{
					newall_pnts.block(0, 0, 2, beginid) = all_pnts.block(0, 1, 2, beginid);
					newall_pnts.block(0, beginid, 2, startid - beginid - 1) = all_pnts.block(0, beginid+1, 2, startid - beginid - 1);
					newall_pnts.block(0, startid - 1, 2, offsetpntnum) = offsetline;
					newall_pnts.block(0, startid + offsetpntnum - 1, 2, pointsnumber - startid - offsetpntnum - 1) = all_pnts.block(0, startid + offsetpntnum + 1, 2, pointsnumber - startid - offsetpntnum - 1);
				}
				else
				{
					newall_pnts.block(0, 0, 2, startid) = all_pnts.block(0, 0, 2, startid);
					newall_pnts.block(0, startid, 2, offsetpntnum) = offsetline;
					newall_pnts.block(0, endid - 1, 2, pointsnumber - endid - 1) = all_pnts.block(0, endid + 1, 2, pointsnumber - endid - 1);
				}
				newall_pnts.block(1, 0, 1, newall_pnts.cols()) *= ra * if_reverse;
				triangulate(newall_pnts, bnd, sqrt(3) * x_step * x_step *0.25 * mu, aMesh);
				double ra_inv = 1.0 / ra;
				for (auto v = aMesh.vertices_begin(); v != aMesh.vertices_end(); v++)
				{
					auto p = aMesh.point(*v);
					aMesh.set_point(*v, Mesh::Point{ p[0],p[1] * ra_inv,0 });
				}
				auto &Surface = faceshape[i].Surface;

				//calculate the target edge length with respect to the error between mesh and surface		
				//double errorbounds = 0.0, testnum = std::min(200, (int)aMesh.n_vertices());
				//for (int j = 0; j < testnum; j++)
				//{
				//	auto p = aMesh.point(aMesh.vertex_handle(j));
				//	errorbounds += (Surface->PartialDerivativeUU(p[0], p[1])).norm() + 2 * (Surface->PartialDerivativeUV(p[0], p[1])).norm() + (Surface->PartialDerivativeVV(p[0], p[1])).norm();
				//}
				//errorbounds /= 2 * testnum;
				//errorbounds = std::max(errorbounds, 0.00001);
				//double target_h = std::sqrt(epsratio / errorbounds);
				//dprint(aMesh.n_vertices(), target_h, x_step);
				//if (target_h < x_step)
				//{
				//	triangulate(newall_pnts, bnd, sqrt(3) * target_h * target_h *0.25, aMesh);
				//	for (auto tv : aMesh.vertices())
				//	{
				//		auto pos = aMesh.point(tv);
				//		aMesh.set_point(tv, TriMesh::Point(pos[0], pos[1] * ra_inv, 0));					
				//	}
				//}

				//Remesh in domain
				Riemannremesh Remesh(Surface, &aMesh);
				Remesh.remesh();

				Mesh::Point newp;
				Mesh::VertexHandle vh1, vh2, vh3, vh4;
				std::vector<Mesh::VertexHandle> facevhandle;

				//add boundary points
				for (int j = 0; j < pointsnumber; j++)
				{
					newp = Mesh::Point(all_pnts(0, j), all_pnts(1, j)*if_reverse, 0);
					newmesh.add_vertex(newp);
				}

				//add offset points and internal points
				for (int j = 0; j < offsetline.cols(); j++)
				{
					newp = Mesh::Point(offsetline(0, j), offsetline(1, j), 0);
					newmesh.add_vertex(newp);
				}
				for (int j = pointsnumber - 2; j < aMesh.n_vertices(); j++)
				{
					auto p = aMesh.point(aMesh.vertex_handle(j));
					newp = Mesh::Point(p[0], p[1], p[2]);
					newmesh.add_vertex(newp);
				}

				//add offset faces
				if (startid == beginid)
				{
					if(startid) vh1 = newmesh.vertex_handle(beginid + wireoffset - 1);
					else vh1 = newmesh.vertex_handle(pointsnumber - 1);
				}
				else vh1 = newmesh.vertex_handle(startid - 1);
				vh2 = newmesh.vertex_handle(startid);
				for (int j = 1; j < offsetpntnum + 1; j++)
				{
					vh3 = newmesh.vertex_handle(startid + j);
					vh4 = newmesh.vertex_handle(pointsnumber + j - 1);
					facevhandle = { vh1, vh2, vh3, vh4 };
					newmesh.add_face(facevhandle);
					auto eh = newmesh.edge_handle(newmesh.find_halfedge(vh2, vh3));
					newmesh.data(eh).flag2 = true;
					vh1 = vh4;
					vh2 = vh3;
				}
				vh3 = newmesh.vertex_handle(endid);
				vh4 = newmesh.vertex_handle(endid+1);
				facevhandle = { vh1, vh2, vh3, vh4 };
				newmesh.add_face(facevhandle);
				auto eh = newmesh.edge_handle(newmesh.find_halfedge(vh2, vh3));
				newmesh.data(eh).flag2 = true;

				//add internal faces
				facevhandle.clear();
				for (auto f : aMesh.faces())
				{
					for (auto fv : aMesh.fv_range(f))
					{
						int id = fv.idx();
						auto p = aMesh.point(fv);
						if(id < beginid) facevhandle.push_back(newmesh.vertex_handle(id));
						else if (id >= pointsnumber - 2) facevhandle.push_back(newmesh.vertex_handle(id + offsetpntnum + 2));
						else if (endid < startid)
						{
							if (id >= startid + offsetpntnum - 1)
							{
								facevhandle.push_back(newmesh.vertex_handle(id + 2));
								auto p1 = newmesh.point(newmesh.vertex_handle(id + 2));
							}
							else if (beginid <= id && id <= startid - 2)
							{
								facevhandle.push_back(newmesh.vertex_handle(id + 1));
								auto p1 = newmesh.point(newmesh.vertex_handle(id + 1));
							}
							else
							{
								facevhandle.push_back(newmesh.vertex_handle(id - startid + pointsnumber + 1));
								auto p1 = newmesh.point(newmesh.vertex_handle(id - startid + pointsnumber + 1));
							}
						}
						else
						{
							if (id >= startid + offsetpntnum) facevhandle.push_back(newmesh.vertex_handle(id + 2));
							else if(id < startid) facevhandle.push_back(newmesh.vertex_handle(id));
							else facevhandle.push_back(newmesh.vertex_handle(id - startid + pointsnumber));
						}
					}
					newmesh.add_face(facevhandle);
					facevhandle.clear();
				}
			}
			
			//if (!OpenMesh::IO::write_mesh(newmesh, "two.obj"))
			//{
			//	std::cerr << "fail";
			//}

			startid = 0;
			for (int j = 0; j < wires.size(); j++)
			{
				beginid = startid;
				auto &edges = wires[j];
				for (int k = 0; k < edges.size(); k++)
				{
					auto &aedge = edgeshape[edges[k]];
					int cols = aedge.parameters.cols() - 1;
					if (aedge.if_C0)
					{
						if (k == edges.size() - 1) endid = beginid;
						else endid = startid + cols;
						for (int m = startid; m < startid + cols - 1; m++)
						{
							auto e = newmesh.edge_handle(newmesh.find_halfedge(newmesh.vertex_handle(m), newmesh.vertex_handle(m + 1)));
							newmesh.data(e).flag1 = true;
						}
						auto e = newmesh.edge_handle(newmesh.find_halfedge(newmesh.vertex_handle(startid + cols - 1), newmesh.vertex_handle(endid)));
						newmesh.data(e).flag1 = true;
					}
					startid += cols;
				}
			}
			
			double k1, k2;
			for (auto v : aMesh.vertices())
			{
				auto p = aMesh.point(v);
				Surface->PrincipalCurvature(p[0], p[1], k1, k2);
				if (isnan(k1*k2))
				{
					dprint("bug", v.idx());
				}
				aMesh.data(v).GaussCurvature = std::max(std::fabs(k1), std::fabs(k2));
			}

			for (auto tv : newmesh.vertices())
			{
				auto pos = newmesh.point(tv);
				auto v = asurface->Value(pos[0], pos[1]);
				newmesh.set_point(tv, Mesh::Point(v.X(), v.Y(), v.Z()));
			}
		}

		dprint("Piecewise PolyMesh Done!");
	}

	Matrix2Xd OccReader::Subdomain(Matrix2Xd &parameters, Matrix2Xd preedgepara, Matrix2Xd nextedgepara)
	{
		int pntsnum = parameters.cols();
		Matrix2Xd offsetline(2, pntsnum-2), dire(2, pntsnum), p1 = Matrix2Xd::Zero(2, 1), p2, p3, p4;
		p2 = p3 = p4 = p1;

		//calculate the offset step
		p1 = parameters.col(0);
		p2 = preedgepara - p1;
		p3 = parameters.col(1) - p1;
		double step = std::abs(p2(0)*p3(1) - p2(1)*p3(0))/p3.norm();
		p1 = parameters.col(pntsnum-1);
		p2 = nextedgepara - p1;
		p3 = parameters.col(pntsnum - 2) - p1;
		step = (step + std::abs(p2(0)*p3(1) - p2(1)*p3(0)) / p3.norm())*0.5;
		step = std::min(step, expected_edge_length * 0.4);

		//calculate the mid-angle line, i.e, the basic offset direction
		dire.col(0) = (preedgepara - parameters.col(0)).normalized();
		dire.col(pntsnum - 1) = (nextedgepara - parameters.col(pntsnum - 1)).normalized();
		for (int i = 1; i < pntsnum - 1; i++)
		{
			p1 = parameters.col(i);
			p2 = (parameters.col(i - 1) - p1).normalized();
			p3 = (parameters.col(i + 1) - p1).normalized();
			if ((p2 + p3).norm() < 0.001) p4 << -p3(1), p3(0);
			else if (p3(0)*p2(1) - p3(1)*p2(0) < 0) p4 = -((p2 + p3)*0.5).normalized();
			else p4 = ((p2 + p3)*0.5).normalized();
			dire.col(i) = p4;
		}
		//adjust the offset direction 
		Matrix2Xd dire_temp(2, pntsnum);
		dire_temp.col(0) = dire.col(0);
		dire_temp.col(pntsnum - 1) = dire.col(pntsnum - 1);
		for (int i = 1; i < pntsnum - 1; i++)
		{
			dire_temp.col(i) = ((dire.col(i - 1) + dire.col(i + 1))*0.5).normalized();
		}
		for (int i = 1; i < pntsnum - 1; i++)
		{
			dire.col(i) = ((dire_temp.col(i - 1) + dire_temp.col(i + 1))*0.5).normalized();
		}
		for (int i = 1; i < pntsnum - 1; i++)
		{
			p1 = parameters.col(i);
			offsetline.col(i-1) = p1 + step * dire.col(i);
		}
		return offsetline;
	}

	void OccReader::ComputeFaceAndEdge()
	{
		TopoDS_Shape &aShape = globalmodel.aShape;

		Vector3d ma(-DBL_MAX, -DBL_MAX, -DBL_MAX);
		Vector3d mi(DBL_MAX, DBL_MAX, DBL_MAX);
		for (TopExp_Explorer vertexExp(aShape, TopAbs_VERTEX); vertexExp.More(); vertexExp.Next())
		{
			gp_Pnt v = BRep_Tool::Pnt(TopoDS::Vertex(vertexExp.Current()));
			ma(0) = std::max(ma(0), v.X()); mi(0) = std::min(mi(0), v.X());
			ma(1) = std::max(ma(1), v.Y()); mi(1) = std::min(mi(1), v.Y());
			ma(2) = std::max(ma(2), v.Z()); mi(2) = std::min(mi(2), v.Z());
		}
		expected_edge_length = initialRate * (ma - mi).norm();
		epsratio *= (ma - mi).norm();//根据Boundingbox设定网格和曲面之间的误差

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
		double vertexThreshold = expected_edge_length * degeneratedRate;
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

		dprint("CAD Info:\nModel Size: [", mi(0), ma(0), "],[", mi(1), ma(1), "],[", mi(2), ma(2), "]"
			"\nSurface Number: ", faceshape.size(),
			"\nEdge Number: ", edgeshape.size(),
			"\n\nCompute Topology Done!");
	}

	void OccReader::narrow_surface()
	{
		//detect the narrow surface, consider only rectangle surface
		auto &faceshape = globalmodel.faceshape;
		auto &edgeshape = globalmodel.edgeshape;
		for (int i = 0; i < faceshape.size(); i++)
		{
			//if (i) continue;
			auto &wire = (faceshape[i]).wires;
			if (wire.empty() || wire.size() > 1) continue;
			auto &edges = wire.front();
			if (edges.size() != 4) continue;
			std::vector<gp_Pnt> corners;
			for (int j = 0; j < 4; j++)
			{
				auto &aedge = edgeshape[edges[j]].edge;
				if (aedge.Orientation() == TopAbs_FORWARD) corners.push_back(BRep_Tool::Pnt(TopExp::FirstVertex(aedge)));
				else corners.push_back(BRep_Tool::Pnt(TopExp::LastVertex(aedge)));
			}
			auto &p1 = corners[0], &p2 = corners[1], p3 = corners[2], p4 = corners[3];
			auto e1 = Point(p2.X() - p1.X(), p2.Y() - p1.Y(), p2.Z() - p1.Z());
			auto e2 = Point(p3.X() - p2.X(), p3.Y() - p2.Y(), p3.Z() - p2.Z());
			auto e3 = -Point(p4.X() - p3.X(), p4.Y() - p3.Y(), p4.Z() - p3.Z());
			auto e4 = -Point(p1.X() - p4.X(), p1.Y() - p4.Y(), p1.Z() - p4.Z());
			//dprint("jiuyg:", (e1.normalized().dot(e3.normalized())), (e2.normalized().dot(e4.normalized())));
			if ((e1.normalized().dot(e3.normalized())) < 0.8 || (e2.normalized().dot(e4.normalized())) < 0.8) continue;
			double length = edgeshape[edges[0]].length, width = edgeshape[edges[1]].length;
			//dprint("kiu:", length, width);
			int id0, id1;
			if (length > 1e4*width) { id0 = 0; id1 = 1; }
			else if (width > 1e4*length) { id0 = 1; id1 = 0; }
			else continue;

#if 1
			faceshape[i].if_exisited = false;
			ShapeEdge &m0 = edgeshape[edgeshape[edges[id0]].reversed_edge];
			ShapeEdge &m1 = edgeshape[edgeshape[edges[id0 + 2]].reversed_edge];
			m0.reversed_edge = m1.id;
			m1.reversed_edge = m0.id;
			m0.secondary_face = m1.main_face;
			m1.secondary_face = m0.main_face;

			edgeshape[edges[id1]].if_exisited = false;
			if (edgeshape[edges[id1]].reversed_edge != -1)
				edgeshape[edgeshape[edges[id1]].reversed_edge].if_exisited = false;
			edgeshape[edges[id1 + 2]].if_exisited = false;
			if (edgeshape[edges[id1 + 2]].reversed_edge != -1)
				edgeshape[edgeshape[edges[id1 + 2]].reversed_edge].if_exisited = false;
			edgeshape[edges[id0]].if_exisited = false;
			edgeshape[edges[id0 + 2]].if_exisited = false;

			int er0 = edgeshape[edges[id1]].reversed_edge;
			if (er0 != -1)
			{
				ShapeFace &f0 = faceshape[edgeshape[edges[id1]].secondary_face];
				for (auto &edges : f0.wires)
				{
					int foundid = -1;
					for (int j = 0; j < edges.size(); ++j)
					{
						if (edges[j] == er0)
							foundid = j;
					}
					if (foundid != -1)
					{
						edges.erase(edges.begin() + foundid);
						break;
					}
				}
			}
			int er1 = edgeshape[edges[id1 + 2]].reversed_edge;
			if (er1 != -1)
			{
				ShapeFace &f1 = faceshape[edgeshape[edges[id1 + 2]].secondary_face];
				for (auto &edges : f1.wires)
				{
					int foundid = -1;
					for (int j = 0; j < edges.size(); ++j)
					{
						if (edges[j] == er1)
							foundid = j;
					}
					if (foundid != -1)
					{
						edges.erase(edges.begin() + foundid);
						break;
					}
				}
			}
#endif
		}
	}

	void OccReader::Discrete_Edge()
	{
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;

		const double GL_c[8] = { 0.1012285363,0.2223810345,0.3137066459,0.3626837834,0.3626837834 ,0.3137066459,0.2223810345,0.1012285363 };
		const double GL_x[8] = { -0.9602898565,-0.7966664774,-0.5255354099,-0.1834346425,0.1834346425,0.5255354099,0.7966664774,0.9602898565 };

		for(auto edge = edgeshape.begin(); edge != edgeshape.end(); ++edge)
		{
			/*static int ti = 0;
			dprint(ti++);*/
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
			int segment_number = std::max(BRep_Tool::IsClosed(aedge) ? 4 : 4, int(curve_length / expected_edge_length));
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

	void OccReader::Face_type()
	{
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;

		for (int i = 0; i < faceshape.size(); i++)
		{
			if (!faceshape[i].if_exisited) continue;
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
					faceshape[i].Surface = new BezierSurface(geom_beziersurface->UDegree(), geom_beziersurface->VDegree(), U1, U2, V1, V2, cp );
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
					if(isvclosed) cp[r - 1].push_back(cp[r - 1][0]);
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
						if(isvclosed) w[r - 1].push_back(w[r - 1][0]);
					}
					if (isuclosed) w.push_back(w.front());
					faceshape[i].Surface = new BSplineSurface(geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), u, v, w, cp);
				}
				else
					faceshape[i].Surface = new BSplineSurface(geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), u, v, cp);
				//dprint(i, cp.size(), cp[0].size(), geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree());
				//gp_Vec u2, v2, uu, uv, vv;
				//gp_Pnt oringe;
				////double ut = u.back()-(u.back()-u.front())*0.5, vt = v.back()-(v.back()-v.front())*0.5;
				//double ut = 179.999998, vt = 0.999998;
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
			if (!edge.if_exisited) continue;
			auto &facer = faceshape[edge.main_face].Surface;
			auto &pnts = edge.parameters;
			Point4 UV = facer->Getbounds();
			double u1 = UV(0), u2 = UV(1), v1 = UV(2), v2 = UV(3);
			for (int j = 0; j < pnts.cols(); j++)
			{
				double &u = pnts(0, j);
				double &v = pnts(1, j);
				if (u <= u1) u = u1 + (u2 - u1)*10e-8;
				else if (u >= u2) u = u2 - (u2 - u1)*10e-8;
				if (v <= v1) v = v1 + (v2 - v1)*10e-8;
				else if (v >= v2) v = v2 - (v2 - v1)*10e-8;
			}
		}
		dprint("face type done!");
	}

	void OccReader::C0_Feature()
	{
		auto &edge = globalmodel.edgeshape;
		auto &face = globalmodel.faceshape;
		for (int i = 0; i < edge.size(); i++)
		{
			//if (i != 11) continue;
			auto &aedge = edge[i];
			if (!aedge.if_exisited) continue;
			//dprint(i, aedge.main_face, aedge.secondary_face);
			if (aedge.reversed_edge == -1)
			{
				aedge.if_C0 = true;
				continue;
			}	
			if (aedge.if_C0) continue;
			double C0 = 0.92;
			int single_flag = 0;
			auto &face1 = face[aedge.main_face].Surface;
			auto &face2 = face[aedge.secondary_face].Surface;
			auto &pnts1 = aedge.parameters;
			auto &pnts2 = edge[aedge.reversed_edge].parameters;
			for(int j = 0; j < pnts1.cols(); j++)
			{
				if (j - single_flag > pnts1.cols()*0.4) break;
				double u = pnts1(0, j), v = pnts1(1, j);
				Point ru = face1->PartialDerivativeU(u, v);
				Point rv = face1->PartialDerivativeV(u, v);
				Point n1 = (ru.cross(rv)).normalized();
				if ((face[aedge.main_face].face).Orientation() == TopAbs_REVERSED) n1 = -n1;

				u = pnts2(0, pnts1.cols() - 1 - j), v = pnts2(1, pnts1.cols() - 1 - j);
				ru = face2->PartialDerivativeU(u, v);
				rv = face2->PartialDerivativeV(u, v);
				Point n2 = (ru.cross(rv)).normalized();
				if ((face[aedge.secondary_face].face).Orientation() == TopAbs_REVERSED) n2 = -n2;
				if (n1.dot(n2) < C0) single_flag++;
			}
			if (single_flag > pnts1.cols()*0.6)
			{
				aedge.if_C0 = true;
				edge[aedge.reversed_edge].if_C0 = true;
			}
		}

		//edge[3].if_C0 = true;
		dprint("C0 feature done!");
	}

	void OccReader::curvature_feature()
	{
		auto &edge = globalmodel.edgeshape;
		auto &face = globalmodel.faceshape;
		double min_curv = 1/expected_edge_length ;
		for (int i = 0; i < edge.size(); i++)
		{
			//if (i != 3) continue;
			auto &aedge = edge[i];
			if (!aedge.if_exisited) continue;
			if (!aedge.if_curvature) continue;
			if ((aedge.parameters).cols() < 10)
			{
				aedge.if_curvature = false;
				continue;
			}
			if (aedge.if_C0 && aedge.reversed_edge != -1)
			{
				aedge.if_curvature = false;
				continue;
			}
			auto &facer = face[aedge.main_face].Surface;
			auto &pnts = aedge.parameters;
			double u1 = pnts(0, 0), v1 = pnts(1, 0), u2 = pnts(0, pnts.cols() - 1), v2 = pnts(1, pnts.cols() - 1);
			Eigen::Vector2d step(v1 - v2, u2 - u1);
			step /= (pnts.cols() - 1)*4;
			if ((face[aedge.main_face].face).Orientation() == TopAbs_REVERSED) step = -step;
			double base_stepu = step(0), base_stepv = step(1);
			Point4 UV = facer->Getbounds();
			double umin = UV(0), umax = UV(1), vmin = UV(2), vmax = UV(3);
			//dprint(min_curv, umin, umax, vmin, vmax);
			Eigen::Matrix2d Weingarten, value, vec;
			int single_flag = 0;
			for (int j = 1; j < pnts.cols()-1; j++)
			{			
				if (j - 1 - single_flag > (pnts.cols() - 2)*0.4) break;
				double stepu = base_stepu, stepv = base_stepv;
				int count = 0;
				while (true) 
				{
					u1 = pnts(0, j);
					v1 = pnts(1, j);
					for (int m = 0; m < 15; m++)
					{
						u1 += stepu;
						v1 += stepv;
						if (m==14 || u1 < umin || u1 > umax || v1 < vmin || v1 > vmax)
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
				stepu *= count*0.06;
				stepv *= count*0.06;
				u1 = pnts(0, j);
				v1 = pnts(1, j);
				//Point start_nomal = ((facer->PartialDerivativeU(u1, v1)).cross(facer->PartialDerivativeV(u1, v1))).normalized(), end_nomal;
				std::vector<double> curvature;
				double K, K1;
				for (int m = 0; m < 15; m++)
				{				
					facer->NormalCurvature(u1, v1, stepu, stepv, K);
					//dprint(j, u1, v1, K);
					curvature.push_back(K);
					//if (m == 14)
					//{
					//	end_nomal = ((facer->PartialDerivativeU(u1, v1)).cross(facer->PartialDerivativeV(u1, v1))).normalized();
					//}
					u1 += stepu;
					v1 += stepv;
				}
				K = curvature.front();
				std::sort(curvature.begin(), curvature.end());
				//dprint("dot", start_nomal.dot(end_nomal));
				if (curvature.back() < 1.5*min_curv || curvature[12] < 0.1*min_curv) continue;
				if (curvature.back() - curvature.front() > min_curv)
				{
					//dprint("jiuyy", j);
					single_flag++;
				}
			}
			//dprint(single_flag, pnts.cols() - 2);
			if (single_flag <= (pnts.cols()-2)*0.6)
			{
				aedge.if_curvature = false;
				if (aedge.reversed_edge != -1) edge[aedge.reversed_edge].if_curvature = false;			
			}
			else if(aedge.reversed_edge == -1) aedge.if_C0 = false;
		}
		dprint("curvature feature done!");
	}

	void OccReader::curvature_discrete()
	{
		auto &edge = globalmodel.edgeshape;
		auto &face = globalmodel.faceshape;
		for (int i = 0; i < edge.size(); i++)
		{
			auto &aedge = edge[i];
			if (!aedge.if_exisited) continue;
			if (!aedge.if_curvature) continue;
			auto f = face[aedge.main_face];
			dprint(i, aedge.main_face, f.wires.size(), f.wires[0].size());
		}
	}
}



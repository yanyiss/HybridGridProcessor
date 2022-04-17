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
			//dprint(i);
			//if (i != 0) continue;
			TriMesh &aMesh = Surface_TriMeshes[i];
			auto &wires = faceshape[i].wires;
			if (wires.empty())
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


			double ra = (y_step < epsilonerror ? 1 : x_step / y_step) * if_reverse;

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
			/*for (int pp = 0; pp < all_pnts.cols(); ++pp)
			{
				dprintwithprecision(15, pp, all_pnts(0, pp), all_pnts(1, pp));
			}*/
			//all_pnts(0, 871) -= 1e-10;
			all_pnts.block(1, 0, 1, all_pnts.cols()) *= ra;
			triangulate(all_pnts, bnd, sqrt(3) * x_step * x_step / 4 * mu, aMesh);
			double ra_inv = 1.0 / ra;
			for (auto tv : aMesh.vertices())
			{
				auto pos = aMesh.point(tv);
				aMesh.set_point(tv, TriMesh::Point(pos[0], pos[1] * ra_inv, 0));
			}

			//Remesh in domain
			//BSplineSurface* B = static_cast<BSplineSurface*>(faceshape[i].Surface);
			auto &Surface = faceshape[i].Surface;
			Riemannremesh Remesh(Surface, &aMesh);
			Remesh.remesh();
			double k1, k2;
			for (auto v : aMesh.vertices())
			{
				auto p = aMesh.point(v);
				Surface->PrincipalCurvature(p[0], p[1], k1, k2);
				aMesh.data(v).GaussCurvature = std::max(std::fabs(k1), std::fabs(k2));
			}

			for (auto tv : aMesh.vertices())
			{
				auto pos = aMesh.point(tv);
				auto v = asurface->Value(pos[0], pos[1]);
				aMesh.set_point(tv, TriMesh::Point(v.X(), v.Y(), v.Z()));
			}

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

#if 0
	void OccReader::Set_PolyMesh()
	{
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;

		Surface_PolyMeshes.resize(faceshape.size());
		for (int i = 0; i < faceshape.size(); i++)
		{
			dprint(i);
			Mesh &newmesh = Surface_PolyMeshes[i];
			TriMesh aMesh;
			auto &wires = faceshape[i].wires;
			if (wires.empty())
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

			Matrix2Xd all_pnts(2, pointsnumber), newall_pnts(2, pointsnumber);
			pointsnumber = 0;
			int s = 0;
			double if_reverse = aface.Orientation() ? -1.0 : 1.0;
			double ra = (y_step < epsilonerror ? 1 : x_step / y_step) * if_reverse;
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
			all_pnts.block(1, 0, 1, all_pnts.cols()) *= if_reverse;
			int startid = 0, endid = -1, beginid;
			for (int j = 0; j < wires.size(); j++)
			{
				beginid = startid;
				auto &edges = wires[j];
				for (int k = 0; k < edges.size(); k++)
				{
					auto &aedge = edgeshape[edges[k]];
					if (aedge.if_curvature)
					{
						if (k == edges.size() - 1) endid = beginid;
						else endid = startid + (aedge.parameters).cols()-1;

					}
				}
			}
			newall_pnts = Subdomain(all_pnts, bnd, pointsnumber);
			newall_pnts.block(1, 0, 1, newall_pnts.cols()) *= ra* if_reverse;
			triangulate(newall_pnts, bnd, sqrt(3) * x_step * x_step / 4 * mu, aMesh);
			double ra_inv = 1.0 / ra;
			for (auto v = aMesh.vertices_begin(); v != aMesh.vertices_end(); v++)
			{
				auto p = aMesh.point(*v);
				aMesh.set_point(*v, Mesh::Point{ p[0],p[1]* ra_inv,0 });
			}

			//Remesh in domain
			Riemannremesh Remesh(faceshape[i].Surface, &aMesh);
			Remesh.remesh();

			all_pnts.block(1, 0, 1, all_pnts.cols()) *= if_reverse;				
			
			Mesh::VertexHandle vh1, vh2, vh3, vh4, vend3, vend4;
			std::vector<Mesh::VertexHandle> facevhandle;
			Mesh::Point newp;

			//add boundary points
			for (int j = 0; j < pointsnumber; j++)
			{
				newp = Mesh::Point(all_pnts(0, j), all_pnts(1, j), 0);
				vh1 = newmesh.add_vertex(newp);
			}

			//add internal points
			for (auto v : aMesh.vertices())
			{
				newp = Mesh::Point(aMesh.point(v));
				vh1 = newmesh.add_vertex(newp);
			}

			//add boundary faces
			int count = 0;
			for (int j = 0; j < bnd.size(); j++)
			{
				vend3 = vh1 = newmesh.vertex_handle(count);
				vend4 = vh4 = newmesh.vertex_handle(count + pointsnumber);
				for (int k = 0; k < bnd[j].cols() - 1; k++)
				{
					vh2 = newmesh.vertex_handle(count + 1);
					vh3 = newmesh.vertex_handle(count + pointsnumber + 1);
					facevhandle = { vh1,vh2,vh3,vh4 };
					newmesh.add_face(facevhandle);
					facevhandle.clear();
					vh1 = vh2;
					vh4 = vh3;
					count++;
				}
				facevhandle = { vh1,vend3,vend4,vh4 };
				newmesh.add_face(facevhandle);
				facevhandle.clear();
				count++;
			}

			//add internal faces
			for (auto f : aMesh.faces())
			{
				for (auto fv = aMesh.fv_begin(f); fv.is_valid(); fv++)
				{
					vh1 = newmesh.vertex_handle((*fv).idx() + pointsnumber);
					facevhandle.push_back(vh1);
				}
				newmesh.add_face(facevhandle);
				facevhandle.clear();
			}
			//if (!OpenMesh::IO::write_mesh(newmesh, "one.obj"))
			//{
			//	std::cerr << "fail";
			//}

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
							auto e = newmesh.edge_handle(newmesh.find_halfedge(newmesh.vertex_handle(m), newmesh.vertex_handle(m + 1)));
							newmesh.data(e).flag2 = true;
						}
						auto e = newmesh.edge_handle(newmesh.find_halfedge(newmesh.vertex_handle(id + cols - 1), newmesh.vertex_handle(endid)));
						newmesh.data(e).flag2 = true;
					}
					id += cols;
				}
			}
			for (auto tv : newmesh.vertices())
			{
				auto pos = newmesh.point(tv);
				auto v = asurface->Value(pos[0], pos[1]);
				newmesh.set_point(tv, TriMesh::Point(v.X(), v.Y(), v.Z()));
			}
		}

		dprint("Piecewise PolyMesh Done!");
	}
#endif

	/*Matrix2Xd OccReader::Subdomain(Matrix2Xd &all_pnts, vector<Matrix2Xi> &bnd, int &pointsnumber)
	{
		Matrix2Xd newall_pnts(2, pointsnumber),edge;
		Matrix2Xi P;
		Matrix2Xd p1 = Matrix2Xd::Zero(2, 1), p2, p3, p4;
		p4 = p3 = p2 = p1;

		int cols;   //edges number
		double step, avelen=0;
		bool dire;
		std::vector<int>inflexion;

		for (int i = 0; i < bnd.size(); i++)
		{
			cols = bnd[i].cols();
			edge = Matrix2Xd::Zero(2, cols+1);
			P = Matrix2Xi::Zero(2, cols + 1);
			P << bnd[i].col(cols - 1), bnd[i];
			int startid = P(1, 0), id1, id2;

			for (int j = 1; j < cols+1; j++)
			{
				edge.col(j) = all_pnts.col(P(1, j)) - all_pnts.col(P(0, j));
				avelen += edge.col(j).norm();
			}
			edge.col(0) = edge.col(cols);
			avelen /= cols;
			step = avelen *1.0;

			//提取拐点并按角平分线设定步长
			Matrix2Xd midAngle(2, cols+2);
			std::vector<int> pntAngle(cols, 0);
			for (int j = 1; j < cols + 1; j++)
			{
				p1 = edge.col(j - 1);
				p2 = edge.col(j);
				double length1 = p1.norm(), length2 = p2.norm();
				p1.normalize();
				p2.normalize();
				if ((p2 - p1).norm() < 0.001)
				{
					p3 << -p1(1), p1(0);
					midAngle.col(j) = p3;
					pntAngle[j - 1] = PI;
					if (length1 < 0.75*avelen || length2 < 0.75*avelen) continue;
				}
				else
				{
					double angle = acos(-p1(0)*p2(0) - p1(1)*p2(1));	
					dire = p1(0)*p2(1) - p1(1)*p2(0) > 0;
					midAngle.col(j) = ((p2 - p1).normalized())* (dire ? 1 : -1);
					if (!dire) angle = 2 * PI - angle;
					pntAngle[j - 1] = angle;
				}
				inflexion.push_back(startid + j - 1);   //边界的拐点
			}
			midAngle.col(0) = midAngle.col(cols);
			midAngle.col(cols+1) = midAngle.col(1);
			for (int j = 0; j < inflexion.size(); j++)
			{
				id1 = inflexion[j]-startid;
				p1 = (edge.col(id1)).normalized();
				if (pntAngle[id1] < 1) p2 = midAngle.col(id1+1);//角比较小时，按自身角平分线
				else p2 = ((midAngle.col(id1) + midAngle.col(id1 + 2))*0.5).normalized(); 
				p2 *= step / std::max(std::sqrt(1 - std::pow(p1(0)*p2(0) + p1(1)*p2(1), 2.0)), 0.1);
				newall_pnts.col(id1+ startid) = all_pnts.col(id1+startid) + p2;
			}
			
			//判断自交
			std::vector<int>::iterator iter;
			Matrix2d A1, A2, A3, A4;
			bool global_intersect = true;
			while (global_intersect)
			{
				global_intersect = false;

				//角平分线自交
				for (int j = 0; j < inflexion.size() - 1; j++)
				{
					id1 = inflexion[j];
					p1 = all_pnts.col(id1);
					p2 = newall_pnts.col(id1);
					for (int k = j + 1; k < inflexion.size(); k++)
					{
						id2 = inflexion[k];
						p3 = all_pnts.col(id2);
						p4 = newall_pnts.col(id2);
						bool intersect = true, modified = false;
						while (intersect)
						{
							A2.col(0) = A1.col(0) = p4 - p3;
							A1.col(1) = p1 - p3;
							A2.col(1) = p2 - p3;
							A4.col(0) = A3.col(0) = p2 - p1;
							A3.col(1) = p3 - p1;
							A4.col(1) = p4 - p1;
							if (A1.determinant()*A2.determinant() < 0 && A3.determinant()*A4.determinant() < 0)
							{
								modified = true;
								p2 = 0.75*p2 + 0.25*p1;
								p4 = 0.75*p4 + 0.25*p3;
							}
							else intersect = false;
						}
						if (modified)
						{
							global_intersect = true;
							newall_pnts.col(id1) = p2;
							newall_pnts.col(id2) = p4;
						}
					}
				}
				Matrix2Xd p1_, p2_, p3_, p4_;

				//角平分线与平行线自交				
				for (int j = 0; j < inflexion.size(); j++)
				{
					id1 = inflexion[j];
					p1 = newall_pnts.col(id1);
					p1_ = all_pnts.col(id1);
					for (int k = 0; k < inflexion.size(); k++)
					{
						if (k == j || k == j - 1) continue;
						if (!j && k == inflexion.size() - 1) continue;
						id2 = inflexion[k];
						int id3 = (k < inflexion.size() - 1 ? inflexion[k + 1] : inflexion.front());
						p2 = newall_pnts.col(id2);
						p2_ = all_pnts.col(id2);
						p3 = newall_pnts.col(id3);
						p3_ = all_pnts.col(id3);
						bool intersect = true, modified = false;
						while (intersect)
						{
							A2.col(0) = A1.col(0) = p1 - p1_;
							A1.col(1) = p2 - p1_;
							A2.col(1) = p3 - p1_;
							A4.col(0) = A3.col(0) = p2 - p3;
							A3.col(1) = p1 - p3;
							A4.col(1) = p1_ - p3;
							if (A1.determinant()*A2.determinant() < 0 && A3.determinant()*A4.determinant() < 0)
							{								
								modified = true;
								p1 = 0.75*p1 + 0.25*p1_;
								p2 = 0.75*p2 + 0.25*p2_;
								p3 = 0.75*p3 + 0.25*p3_;
							}
							else intersect = false;
						}
						if (modified)
						{
							global_intersect = true;
							newall_pnts.col(id1) = p1;
							newall_pnts.col(id2) = p2;
							newall_pnts.col(id3) = p3;
						}
					}
				}

				//平行线与平行线自交
				for (int j = 0; j < inflexion.size()-1; j++)
				{
					//可以考虑用容器指针
					id1 = inflexion[j];
					p1 = newall_pnts.col(id1);
					p2 = newall_pnts.col(inflexion[j + 1]);
					p1_ = all_pnts.col(id1);
					p2_ = all_pnts.col(inflexion[j + 1]);
					for (int k = j + 2; k < inflexion.size(); k++)
					{
						if (!j && k == inflexion.size() - 1) continue;
						id2 = inflexion[k];
						p3 = newall_pnts.col(id2);
						p3_ = all_pnts.col(id2);
						int id3 = (k < inflexion.size() - 1 ? inflexion[k+1] : inflexion.front());
						p4 = newall_pnts.col(id3);
						p4_ = all_pnts.col(id3);
						bool intersect = true, modified = false;
						while (intersect)
						{
							A2.col(0) = A1.col(0) = p2 - p1;
							A1.col(1) = p3 - p1;
							A2.col(1) = p4 - p1;
							A4.col(0) = A3.col(0) = p4 - p3;
							A3.col(1) = p1 - p3;
							A4.col(1) = p2 - p3;
							if (A1.determinant()*A2.determinant() < 0 && A3.determinant()*A4.determinant() < 0)
							{
								modified = true;
								p1 = 0.75*p1 + 0.25*p1_;
								p2 = 0.75*p2 + 0.25*p2_;
								p3 = 0.75*p3 + 0.25*p3_;
								p4 = 0.75*p4 + 0.25*p4_;
							}
							else intersect = false;
						}
						//可以考虑变量引用，就无需再赋值
						if (modified)
						{
							global_intersect = true;
							newall_pnts.col(id1) = p1;
							newall_pnts.col(inflexion[j + 1]) = p2;
							newall_pnts.col(id2) = p3;
							newall_pnts.col(id3) = p4;
						}
					}
				}
			}
			
			//插值其余点
			int n;
			for (int j = 0; j < inflexion.size() - 1; j++)
			{
				id1 = inflexion[j];
				id2 = inflexion[j + 1];
				n = id2 - id1;  //n等分
				if (n == 1) continue;
				p1 = newall_pnts.col(id1);
				p2 = newall_pnts.col(id2);
				for (double k = 1; k < n; k++)
				{
					newall_pnts.col(id1 + k) = (1 - k / n)*p1 + k / n * p2;
				}
			}
			id1 = inflexion.back();
			id2 = inflexion.front();
			n = P(0, 0) - id1 + id2 - P(0, 1) + 1;  //n等分 
			if (n > 1)
			{
				p1 = newall_pnts.col(id1);
				p2 = newall_pnts.col(id2);
				for (double k = 1; k < n; k++)
				{
					if (id1 + k <= P(0, 0))
					{
						newall_pnts.col(id1 + k) = (1 - k / n)*p1 + k / n * p2;
					}
					else
					{
						newall_pnts.col(id1 + k - P(0, 0) + P(0, 1) - 1) = (1 - k / n)*p1 + k / n * p2;
					}
				}
			}
			inflexion.clear();
		}
		return newall_pnts;
	}*/

	/*Matrix2Xd OccReader::Subdomain(Matrix2Xd &all_pnts, vector<Matrix2Xi> &bnd, int &pointsnumber)
	{

	}*/

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
				dprint(i, "CylindricalSurface");
				opencascade::handle<Geom_CylindricalSurface> geom_cylindricalsurface = opencascade::handle<Geom_CylindricalSurface>::DownCast(geom_surface);;
				auto local_coordinate = &geom_cylindricalsurface->Position();
				gp_Pnt oringe = local_coordinate->Location();
				gp_Dir xdir = local_coordinate->XDirection();
				gp_Dir ydir = local_coordinate->YDirection();
				gp_Dir zdir = (geom_cylindricalsurface->Axis()).Direction();
				double r = geom_cylindricalsurface->Radius();
				//dprint("cur:", 1 / r);
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
				//gp_Vec u2, v2, uu, uv, vv;
				//gp_Pnt oringe;
				//geom_bsplinesurface->D2(u[1], v[1], oringe, u2, v2, uu, vv, uv);
				//auto u1 = faceshape[i].Surface->PartialDerivativeU(u[1], v[1]);
				//auto v1 = faceshape[i].Surface->PartialDerivativeV(u[1], v[1]);
				//auto uu1 = faceshape[i].Surface->PartialDerivativeUU(u[1], v[1]);
				//auto vv1 = faceshape[i].Surface->PartialDerivativeVV(u[1], v[1]);
				//auto uv1 = faceshape[i].Surface->PartialDerivativeUV(u[1], v[1]);
				//dprint("D1U:", u2.X(), u2.Y(), u2.Z(), u1(0), u1(1), u1(2));
				//dprint("D1V:", v2.X(), v2.Y(), v2.Z(), v1(0), v1(1), v1(2));
				//dprint("D2U:", uu.X(), uu.Y(), uu.Z(), uu1(0), uu1(1), uu1(2));
				//dprint("D2V:", vv.X(), vv.Y(), vv.Z(), vv1(0), vv1(1), vv1(2));
				//dprint("D2UV:", uv.X(), uv.Y(), uv.Z(), uv1(0), uv1(1), uv1(2));
				//delete faceshape[i].Surface;
				//faceshape[i].Surface = nullptr;
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
					dprint(weights->NbRows(), weights->NbColumns());
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

				TColgp_Array2OfPnt controlpoints = geom_bsplinesurface->Poles();
				vector<vector<Point>> cp(controlpoints.NbRows());
				for (int r = 1; r <= controlpoints.NbRows(); r++)
				{
					cp[r - 1].reserve(controlpoints.NbColumns());
					for (int c = 1; c <= controlpoints.NbColumns(); c++) {
						gp_Pnt pos = controlpoints.Value(r, c);
						cp[r - 1].emplace_back(pos.X(), pos.Y(), pos.Z());
					}
				}
				const TColStd_Array2OfReal* weights = geom_bsplinesurface->Weights();
				if (weights) {
					vector<vector<double>> w(weights->NbRows());
					for (int r = 1; r <= weights->NbRows(); r++)
					{
						w[r - 1].reserve(weights->NbColumns());
						for (int c = 1; c <= weights->NbColumns(); c++)
							w[r - 1].push_back(weights->Value(r, c));
					}
					faceshape[i].Surface = new BSplineSurface(geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), u, v, w, cp);
				}
				else
					faceshape[i].Surface = new BSplineSurface(geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), u, v, cp);
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
			auto &facer = faceshape[edge.main_face].Surface;
			auto &pnts = edge.parameters;
			Point4 UV = facer->Getbounds();
			double u1 = UV(0), u2 = UV(1), v1 = UV(2), v2 = UV(3);
			for (int j = 0; j < pnts.cols(); j++)
			{
				double &u = pnts(0, j);
				double &v = pnts(1, j);
				if (u < u1) u = u1;
				else if (u > u2) u = u2;
				if (v < v1) v = v1;
				else if (v > v2) v = v2;
			}
		}
	}

	void OccReader::Trim_Edge()  
	{
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		for (int i = 0; i < faceshape.size(); i++)
		{
			auto &face = faceshape[i];

			//获取该面片的参数域边界
			Point4 UV = face.Surface->Getbounds();
			double u1 = UV(0), u2 = UV(1), v1 = UV(2), v2 = UV(3);
			double epsilonu = (u2 - u1)*0.001, epsilonv = (v2 - v1)*0.001;
			//dprint("参数域端点：", u1, u2, v1, v2, epsilonu, epsilonv);

			//判断该边是否被裁剪
			auto &wires = face.wires;
			for (int j = 0; j < wires.size(); j++)
			{
				if (wires.empty()) continue;
				auto &edge = wires[j];
				for (int k = 0; k < edge.size(); k++)
				{
					auto &aedge = globalmodel.edgeshape[edge[k]];
					if (!aedge.if_trimmed) continue;
					auto &parameters = aedge.parameters;
					int line;
					double epsilon;
					double a1 = parameters(0, 0), a2 = parameters(1, 0), a3 = parameters(0, parameters.cols()-1), a4 = parameters(1, parameters.cols()-1);
					//dprint(a1, a2, a3, a4);
					if (std::abs(a1 - a3) > epsilonu)
					{
						if (std::abs(a2 - a4) > epsilonv)
						{
							//dprint("yy", k);
							continue;
						}
						line = 1;
						epsilon = epsilonv;
					}
					else
					{
						line = 0;
						epsilon = epsilonu;
					}
					double standard = parameters(line, 0);
					int nonparallel = 0;
					for (int m = 1; m < parameters.cols()-1; m++)
					{
						//dprint(m, parameters(line, m), standard, std::abs(parameters(line, m) - standard), epsilon);
						if (std::abs(parameters(line, m) - standard) > epsilon) nonparallel++;
						
					}
					if (nonparallel <= parameters.cols()*0.5)
					{
						aedge.if_trimmed = false;
						if(aedge.reversed_edge != -1)
							globalmodel.edgeshape[aedge.reversed_edge].if_trimmed = false;
						//dprint("kk", k);
					}
				}
			}
		}		
	}

	void OccReader::C0_Feature()
	{
		auto &edge = globalmodel.edgeshape;
		auto &face = globalmodel.faceshape;
		for (int i = 0; i < edge.size(); i++)
		{
			auto &aedge = edge[i];
			if (aedge.if_C0) continue;
			if (aedge.reversed_edge == -1)
			{
				aedge.if_C0 = true;
				continue;
			}	
			double C0;
			int single_flag = 0;
			if (aedge.if_trimmed) C0 = 0.95;
			else C0 = 0.9;
			//dprint(i, aedge.main_face, aedge.secondary_face);
			//dprint(C0);
			auto &face1 = face[aedge.main_face].Surface;
			auto &face2 = face[aedge.secondary_face].Surface;
			auto &pnts1 = aedge.parameters;
			auto &pnts2 = edge[aedge.reversed_edge].parameters;
			for(int j = 0; j < pnts1.cols(); j++)
			{
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
				//dprint("the dot:", n1[0], n1[1], n1[2], n2[0], n2[1], n2[2], n1.dot(n2));
				if (n1.dot(n2) < C0) single_flag++;
			}
			//dprint(pnts1.cols(), single_flag);
			if (single_flag > pnts1.cols()*0.6)
			{
				//dprint("non_trim:", i);
				aedge.if_C0 = true;
				edge[aedge.reversed_edge].if_C0 = true;
			}
		}
	}

	void OccReader::curvature_feature()
	{
		auto &edge = globalmodel.edgeshape;
		auto &face = globalmodel.faceshape;
		for (int i = 0; i < edge.size(); i++)
		{
			//dprint("edge", i);
			auto &aedge = edge[i];
			if (!aedge.if_curvature) continue;
			if (aedge.if_C0 )
			{
				aedge.if_curvature = false;
				continue;
			}
			auto &facer = face[aedge.main_face].Surface;
			//dprint("surface_id:", aedge.main_face, aedge.reversed_edge);
			auto &pnts = aedge.parameters;
			double u1 = pnts(0, 0), v1 = pnts(1, 0), u2 = pnts(0, pnts.cols() - 1), v2 = pnts(1, pnts.cols() - 1);
			Eigen::Vector2d step(v1 - v2, u2 - u1);
			step.normalize();
			double base_stepu = step(0), base_stepv = step(1);
			//dprint(u1, v1, u2, v2, base_stepu, base_stepv);
			Point4 UV = facer->Getbounds();
			double umin = UV(0), umax = UV(1), vmin = UV(2), vmax = UV(3);
			//dprint(umin, umax, vmin, vmax);
			Eigen::Matrix2d Weingarten, value, vec;
			int single_flag = 0;
			for (int j = 1; j < pnts.cols()-1; j++)
			{				
				//dprint("point", j);
				double stepu = base_stepu, stepv = base_stepv;
				int count = 0;
				while (true) 
				{
					u1 = pnts(0, j);
					v1 = pnts(1, j);
					for (int m = 0; m < 15; m++)
					{
						//dprint(u1, v1);
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
				std::vector<double> curvature;
				double K1, K2;
				for (int m = 0; m < 15; m++)
				{				
					facer->PrincipalCurvature(u1, v1, K1, K2);
					//dprint(j, u1, v1, K1);
					curvature.push_back(K1);
					u1 += stepu;
					v1 += stepv;
				}
				K1 = curvature.front();
				std::sort(curvature.begin(), curvature.end());
				//if (curvature.back() - curvature[curvature.size() - 2] > 4) curvature.pop_back();
				if (curvature.back() - curvature.front() > 4 && K1 > curvature[10])
				{
					single_flag++;
					//dprint("juidjfhrythg", j);
				}
			}
			if (single_flag <= (pnts.cols()-2)*0.6)
			{
				aedge.if_curvature = false;
				edge[aedge.reversed_edge].if_curvature = false;
			}
		}
	}
}



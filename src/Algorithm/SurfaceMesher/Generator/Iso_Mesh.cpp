#include "Iso_Mesh.h"
#include <fstream>
#include <iostream>
#include<algorithm>
#include "..\src\Algorithm\SurfaceMesher\Optimizer\TriangleMeshRemeshing.h"


namespace CADMesher
{
	Iso_Mesh::Iso_Mesh(QString & fileName)
	{
		occ_reader = new OccReader(fileName);
#if 1
		occ_reader->Set_TriMesh();
		//occ_reader->Surface_delete();
		MergeModel();
		ResetFeature();
		TriangleMeshRemeshing trm(&(globalmodel.initial_trimesh));
		trm.run();
		Write_Obj(globalmodel.initial_trimesh);
#else 
		occ_reader->Set_PolyMesh();
		MergeModel();
		ResetFeature1();
		TriangleMeshRemeshing trm(&(globalmodel.initial_polymesh));
		trm.run();
		//Write_Obj(globalmodel.initial_polymesh);
#endif
	}

	void Iso_Mesh::MergeModel()
	{
#if 1
		TriMesh &model_mesh = globalmodel.initial_trimesh;
		auto &surface_meshes = occ_reader->Surface_TriMeshes;
#else
		Mesh &model_mesh = globalmodel.initial_polymesh;
		auto &surface_meshes = occ_reader->Surface_PolyMeshes;
#endif
		vector<ShapeFace> faceshape = globalmodel.faceshape;
		vector<ShapeEdge> edgeshape = globalmodel.edgeshape;
		model_mesh.clear();
		for (auto &face : faceshape)
		{
			auto &wires = face.wires;
			for (auto &edges : wires)
			{
				int s = edges.size();
				for (int i = 0; i < s; ++i)
				{
					int ei = edges[i];
					if (edgeshape[ei].reversed_edge == -1)
					{
						continue;
					}
					auto &redge = edgeshape[edgeshape[ei].reversed_edge];
					if (!BRep_Tool::IsClosed(edgeshape[ei].edge) && !BRep_Tool::IsClosed(redge.edge))
					{
						continue;
					}

					int splitdump = 2;
					if (redge.if_splitted)
					{
						splitdump = edgeshape[ei].parameters.cols() - 3;
					}
					edges.insert(edges.begin() + i + 1, edgeshape.size());
					edgeshape[ei].prev_edge = edgeshape.size();
					edgeshape.emplace_back(edgeshape.size(), edgeshape[ei].edge);

					ShapeEdge &newedge = edgeshape.back();
					newedge.main_face = edgeshape[ei].main_face;
					newedge.secondary_face = edgeshape[ei].secondary_face;
					newedge.parameters = edgeshape[ei].parameters.block(0, splitdump, 2, edgeshape[ei].parameters.cols() - splitdump);
					newedge.reversed_edge = edgeshape[ei].reversed_edge;

					edgeshape[ei].parameters.conservativeResize(2, splitdump + 1);
					edgeshape[ei].if_splitted = true;
					newedge.if_splitted = true;
					if (redge.if_splitted)
					{
						edgeshape[ei].reversed_edge = redge.prev_edge;
						newedge.reversed_edge = redge.id;
						redge.reversed_edge = newedge.id;
						edgeshape[redge.prev_edge].reversed_edge = edgeshape[ei].id;
					}
					++i;
					++s;
				}
			}
		}
		dprint("Split Edges Done!");
		vector<unsigned> &triangle_surface_index = globalmodel.triangle_surface_index;
		int pointsnum = 0;
		for (auto &frac_mesh : surface_meshes)
		{
			pointsnum += frac_mesh.n_vertices();
		}
		triangle_surface_index.reserve(pointsnum);
		pointsnum = 0;
		int id = 0;
		for (int i = 0; i < surface_meshes.size(); i++)
		{
			auto &frac_mesh = surface_meshes[i];
			vector<Mesh::VertexHandle> vhandle;
			Mesh::VertexHandle vh;
			vhandle.reserve(frac_mesh.n_vertices());
			for (auto tv : frac_mesh.vertices())
			{
				auto v = frac_mesh.point(tv);
				vh = model_mesh.add_vertex(Mesh::Point(v[0], v[1], v[2]));
				model_mesh.data(vh).GaussCurvature = frac_mesh.data(tv).GaussCurvature;
				vhandle.push_back(vh);
			}
			for (auto tf : frac_mesh.faces())
			{
				vector<TriMesh::VertexHandle> pos;
				/*for (auto tfv = frac_mesh.cfv_begin(tf); tfv != frac_mesh.cfv_end(tf); tfv++)
				{
					pos.push_back(tfv->idx());
				}*/
				for (auto tfv : frac_mesh.fv_range(tf))
					pos.push_back(vhandle[tfv.idx()]);
				model_mesh.add_face(pos);
				triangle_surface_index.push_back(i);
			}
			for (auto e : frac_mesh.edges())
			{
				if (frac_mesh.data(e).flag1)
				{
					auto he = model_mesh.find_halfedge(model_mesh.vertex_handle((e.v0()).idx() + id), model_mesh.vertex_handle((e.v1()).idx() + id));
					model_mesh.data(model_mesh.edge_handle(he)).flag1 = true;
				}
				if (frac_mesh.data(e).flag2)
				{
					auto he = model_mesh.find_halfedge(model_mesh.vertex_handle((e.v0()).idx() + id), model_mesh.vertex_handle((e.v1()).idx() + id));
					model_mesh.data(model_mesh.edge_handle(he)).flag2 = true;
				}
			}
			id = model_mesh.n_vertices();

			auto &wires = faceshape[i].wires;
			for (auto &edges : wires)
			{
				int start_id = pointsnum;
				for (auto itr = edges.begin(); itr != edges.end(); itr++)
				{
					auto &aedge = edgeshape[*itr];
					aedge.begin_id = pointsnum;
					pointsnum += aedge.parameters.cols() - 1;
					aedge.end_id = itr != edges.end() - 1 ? pointsnum : start_id;
					aedge.prev_edge = itr != edges.begin() ? *(itr - 1) : edges.back();
				}
			}
			pointsnum = model_mesh.n_vertices();
		}

		for (auto te : model_mesh.edges())
		{
			if (model_mesh.is_boundary(te))
				model_mesh.data(te).set_edgeflag(true);
			else
				model_mesh.data(te).set_edgeflag(false);
		}
#if 1
		for (int i = 0; i < edgeshape.size(); i++)
		{
			auto &edge0 = edgeshape[i];
			if (edge0.if_merged || edge0.reversed_edge == -1) continue;
			auto &edge1 = edgeshape[edgeshape[i].reversed_edge];
			edge0.if_merged = true;
			edge1.if_merged = true;
			int m0 = edge0.begin_id;
			int m1 = edge1.begin_id;
			int n0 = EndId(edgeshape, edge0.id);
			int n1 = EndId(edgeshape, edge1.id);
			edgeshape[edge1.prev_edge].next_reversed_edge = edge0.id;
			edgeshape[edge0.prev_edge].next_reversed_edge = edge1.id;
			if (m0 != n1)
			{
				model_mesh.add_face(model_mesh.vertex_handle(m0), model_mesh.vertex_handle(n1), model_mesh.vertex_handle(m0 + 1));
			}
			if (m1 != n0)
			{
				model_mesh.add_face(model_mesh.vertex_handle(m1), model_mesh.vertex_handle(n0), model_mesh.vertex_handle(m1 + 1));
			}
			int length = edge0.parameters.cols();
			model_mesh.add_face(model_mesh.vertex_handle(m0 + 1), model_mesh.vertex_handle(n1), model_mesh.vertex_handle(m1 + length - 2));
			model_mesh.add_face(model_mesh.vertex_handle(m1 + 1), model_mesh.vertex_handle(n0), model_mesh.vertex_handle(m0 + length - 2));
			for (int j = 1; j < length - 2; j++)
			{
				model_mesh.add_face(model_mesh.vertex_handle(m0 + j),
					model_mesh.vertex_handle(m1 + length - j - 1), model_mesh.vertex_handle(m0 + j + 1));
				model_mesh.add_face(model_mesh.vertex_handle(m0 + j + 1)
					, model_mesh.vertex_handle(m1 + length - j - 1), model_mesh.vertex_handle(m1 + length - j - 2));
			}
			if (m0 != n1)
			{
				vector<OV> fv, fv1, fv2;
				for (auto fe : model_mesh.voh_range(model_mesh.vertex_handle(m0))) {
					if (model_mesh.data(model_mesh.edge_handle(fe)).get_edgeflag())
						fv.push_back(model_mesh.to_vertex_handle(fe));
					if (model_mesh.data(model_mesh.edge_handle(fe)).flag1)
						fv1.push_back(model_mesh.to_vertex_handle(fe));
					if (model_mesh.data(model_mesh.edge_handle(fe)).flag2)
						fv2.push_back(model_mesh.to_vertex_handle(fe));
				}
				model_mesh.data(model_mesh.vertex_handle(n1)).GaussCurvature += model_mesh.data(model_mesh.vertex_handle(m0)).GaussCurvature;
				model_mesh.data(model_mesh.vertex_handle(n1)).GaussCurvature *= 0.5;
				model_mesh.collapse(model_mesh.find_halfedge(model_mesh.vertex_handle(m0), model_mesh.vertex_handle(n1)));
				for (OV v : fv) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(n1));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).set_edgeflag(true);

				}
				for (OV v : fv1) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(n1));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).flag1 = true;

				}
				for (OV v : fv2) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(n1));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).flag2 = true;

				}
			}
			if (m1 != n0)
			{
				vector<OV> fv, fv1, fv2;
				for (auto fe : model_mesh.voh_range(model_mesh.vertex_handle(m1))) {
					if (model_mesh.data(model_mesh.edge_handle(fe)).get_edgeflag())
						fv.push_back(model_mesh.to_vertex_handle(fe));
					if (model_mesh.data(model_mesh.edge_handle(fe)).flag1)
						fv1.push_back(model_mesh.to_vertex_handle(fe));
					if (model_mesh.data(model_mesh.edge_handle(fe)).flag2)
						fv2.push_back(model_mesh.to_vertex_handle(fe));
				}
				model_mesh.data(model_mesh.vertex_handle(n0)).GaussCurvature += model_mesh.data(model_mesh.vertex_handle(m1)).GaussCurvature;
				model_mesh.data(model_mesh.vertex_handle(n0)).GaussCurvature *= 0.5;
				model_mesh.collapse(model_mesh.find_halfedge(model_mesh.vertex_handle(m1), model_mesh.vertex_handle(n0)));
				for (OV v : fv) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(n0));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).set_edgeflag(true);
				}
				for (OV v : fv1) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(n0));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).flag1 = true;
				}
				for (OV v : fv2) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(n0));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).flag2 = true;
				}
			}
			for (int j = 1; j < length - 1; j++)
			{
				vector<OV> fv, fv1, fv2;
				for (auto fe : model_mesh.voh_range(model_mesh.vertex_handle(m0 + j))) {
					if (model_mesh.data(model_mesh.edge_handle(fe)).get_edgeflag())
						fv.push_back(model_mesh.to_vertex_handle(fe));
					if (model_mesh.data(model_mesh.edge_handle(fe)).flag1)
						fv1.push_back(model_mesh.to_vertex_handle(fe));
					if (model_mesh.data(model_mesh.edge_handle(fe)).flag2)
						fv2.push_back(model_mesh.to_vertex_handle(fe));
				}
				model_mesh.data(model_mesh.vertex_handle(m1 + length - j - 1)).GaussCurvature += model_mesh.data(model_mesh.vertex_handle(m0 + j)).GaussCurvature;
				model_mesh.data(model_mesh.vertex_handle(m1 + length - j - 1)).GaussCurvature *= 0.5;
				model_mesh.collapse(model_mesh.find_halfedge(model_mesh.vertex_handle(m0 + j), model_mesh.vertex_handle(m1 + length - j - 1)));
				for (OV v : fv) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(m1 + length - j - 1));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).set_edgeflag(true);
				}
				for (OV v : fv1) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(m1 + length - j - 1));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).flag1 = true;
				}
				for (OV v : fv2) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(m1 + length - j - 1));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).flag2 = true;
				}
			}
		}
		model_mesh.garbage_collection();
#endif
		dprint("Merge Meshes Done!");
	}

	void Iso_Mesh::InitTree()
	{
		typedef ClosestPointSearch::AABBTree tree;
		vector<TriMesh> mt = occ_reader->Surface_TriMeshes;
		globalmodel.init_surfacemesh_tree = new tree*[mt.size()];
		for (int i = 0; i < mt.size(); i++)
		{
			globalmodel.init_surfacemesh_tree[i] = new tree(mt[i]);
		}
		globalmodel.init_trimesh_tree = new tree(globalmodel.initial_trimesh);
	}

	int Iso_Mesh::EndId(vector<ShapeEdge> &edgeshape, int edge_id)
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

	void Iso_Mesh::ResetFeature()
	{
		TriMesh &model_mesh = globalmodel.initial_trimesh;
		for (auto te : model_mesh.edges())
		{
			if (!model_mesh.data(te).flag1 && !model_mesh.data(te).flag2) continue;
			model_mesh.data(te.v0()).set_vertflag(true);
			model_mesh.data(te.v1()).set_vertflag(true);
		}
		dprint("Reset Feature Done!");
	}

	void Iso_Mesh::ResetFeature1()
	{
		Mesh &model_mesh = globalmodel.initial_polymesh;
		int offsetnum = 5;
		//double offsetlen = occ_reader->expected_edge_length*(2 - std::pow(0.5, offsetnum));
		double offsetlen = occ_reader->expected_edge_length*offsetnum*0.125;
		TriMesh standardmesh(model_mesh);
		auto aabbtree = new ClosestPointSearch::AABBTree(standardmesh);
		Mesh::VertexHandle v1, v2, v3, v4;
		Mesh::Point p1, p2, p3, p4;
		std::vector<Mesh::VertexHandle> facevhandle;

		//set vertex flag
		for (auto te : model_mesh.edges())
		{
			if (model_mesh.data(te).flag1)
			{
				model_mesh.data(te.v0()).set_vertflag(true);
				model_mesh.data(te.v1()).set_vertflag(true);
				model_mesh.data(te.v0()).flag1 = true;
				model_mesh.data(te.v1()).flag1 = true;
			}
			if (!model_mesh.data(te).flag2) continue;
			model_mesh.data(te.v0()).set_vertflag(true);
			model_mesh.data(te.v1()).set_vertflag(true);
			model_mesh.data(te.v0()).flag2 = true;
			model_mesh.data(te.v1()).flag2 = true;
			for (int i = 0; i < 2; i++)
			{
				if (i && model_mesh.is_boundary(te)) continue;
				auto he = model_mesh.halfedge_handle(te, i);
				v1 = model_mesh.to_vertex_handle(model_mesh.next_halfedge_handle(he));
				model_mesh.data(v1).set_vertflag(true);
				model_mesh.data(v1).flag2 = true;
				v1 = model_mesh.from_vertex_handle(model_mesh.prev_halfedge_handle(he));
				model_mesh.data(v1).set_vertflag(true);
				model_mesh.data(v1).flag2 = true;
			}
		}

		//collapse the input edge ended by vertex on quad and remove abnormal quad
		AdjustQuad(model_mesh, offsetlen);
		//auto v = model_mesh.vertex_handle(1937);
		//for (auto vf : model_mesh.vf_range(v))
		//{
		//	dprint(vf.idx());
		//}
		//auto f = model_mesh.face_handle(58139);
		//dprint("hh");
		//for (auto fv : model_mesh.fv_range(f))
		//{
		//	dprint(fv.idx());
		//}
		//f = model_mesh.face_handle(58140);
		//dprint("hh");
		//for (auto fv : model_mesh.fv_range(f))
		//{
		//	dprint(fv.idx());
		//}
		//f = model_mesh.face_handle(42740);
		//dprint("hh");
		//for (auto fv : model_mesh.fv_range(f))
		//{
		//	dprint(fv.idx());
		//}
		//f = model_mesh.face_handle(42832);
		//dprint("hh");
		//for (auto fv : model_mesh.fv_range(f))
		//{
		//	dprint(fv.idx());
		//}
		//f = model_mesh.face_handle(58467);
		//dprint("hh");
		//for (auto fv : model_mesh.fv_range(f))
		//{
		//	dprint(fv.idx());
		//}

		//set offset domain
		/*for (auto te : model_mesh.edges())
		{
			if (!model_mesh.data(te).flag2) continue;
			for (int i = 0; i < 2; i++)
			{
				if (i && model_mesh.is_boundary(te)) continue;
				auto he = model_mesh.halfedge_handle(te, i);
				auto he1 = model_mesh.next_halfedge_handle(model_mesh.next_halfedge_handle(he));
				v1 = model_mesh.from_vertex_handle(he1);
				if (!model_mesh.data(v1).is_adjusted)
				{
					p1 = model_mesh.point(model_mesh.to_vertex_handle(he));
					p2 = model_mesh.point(v1);
					p2 = aabbtree->closest_point(p1 + (p2 - p1).normalized()*offsetlen);
					p2 = aabbtree->closest_point(p1 + (p2 - p1).normalized()*offsetlen);
					model_mesh.set_point(v1, p2);
					model_mesh.data(v1).is_adjusted = true;
				}
				v1 = model_mesh.to_vertex_handle(he1);
				if (!model_mesh.data(v1).is_adjusted)
				{
					p1 = model_mesh.point(model_mesh.from_vertex_handle(he));
					p2 = model_mesh.point(v1);
					p2 = aabbtree->closest_point(p1 + (p2 - p1).normalized()*offsetlen);
					p2 = aabbtree->closest_point(p1 + (p2 - p1).normalized()*offsetlen);
					model_mesh.set_point(v1, p2);
					model_mesh.data(v1).is_adjusted = true;
				}
			}
		}*/
		delete aabbtree;
		aabbtree = nullptr;
		dprint("Reset Feature Done!");
	}

#pragma region

	void Iso_Mesh::AdjustQuad(Mesh &model_mesh, double &offsetlen)
	{
		Mesh::VertexHandle v1, v2, v3, v4;
		std::vector<Mesh::VertexHandle> facevhandle;
		for (int iter = 0;  iter < 2; iter++)
		{
			for (auto te : model_mesh.edges())
			{
				if (!model_mesh.data(te).flag2) continue;
				for (int i = 0; i < 2; i++)
				{
					if (i && model_mesh.is_boundary(te)) continue;
					auto he = model_mesh.halfedge_handle(te, i);
					auto he1 = model_mesh.next_halfedge_handle(model_mesh.next_halfedge_handle(he));
					v1 = model_mesh.from_vertex_handle(he1);
					auto offset_direction = (model_mesh.point(v1) - model_mesh.point(model_mesh.to_vertex_handle(he))).normalized();
					if (!model_mesh.data(v1).is_adjusted)
					{
						std::vector<OH> VIH;
						for (auto vih : model_mesh.vih_range(v1))
						{
							v2 = model_mesh.from_vertex_handle(vih);
							if (!model_mesh.is_collapse_ok(vih) || model_mesh.data(v2).flag2)
								continue;
							if (model_mesh.data(v2).flag1 && !model_mesh.data(v1).flag1)
								continue;
							if ((model_mesh.point(v2) - model_mesh.point(v1)).dot(offset_direction) < offsetlen)
								VIH.push_back(vih);
						}
						if (VIH.size())
							for (int i = 0; i < VIH.size(); i++) model_mesh.collapse(VIH[i]);
						model_mesh.data(v1).is_adjusted = true;
					}
					v1 = model_mesh.to_vertex_handle(he1);
					offset_direction = (model_mesh.point(v1) - model_mesh.point(model_mesh.to_vertex_handle(he))).normalized();
					if (!model_mesh.data(v1).is_adjusted)
					{
						std::vector<OH> VIH;
						for (auto vih : model_mesh.vih_range(v1))
						{
							v2 = model_mesh.from_vertex_handle(vih);
							if (!model_mesh.is_collapse_ok(vih) || model_mesh.data(v2).flag2)
								continue;
							if (model_mesh.data(v2).flag1 && !model_mesh.data(v1).flag1)
								continue;
							if ((model_mesh.point(v2) - model_mesh.point(v1)).dot(offset_direction) < offsetlen)
								VIH.push_back(vih);
						}
						for (int i = 0; i < VIH.size(); i++) model_mesh.collapse(VIH[i]);
						model_mesh.data(v1).is_adjusted = true;
					}
				}
			}
			model_mesh.garbage_collection();

			//remove abnomal quad
			bool is_Ok = false;
			while (is_Ok)
			{
				is_Ok = false;
				for (auto te : model_mesh.edges())
				{
					if (!model_mesh.data(te).flag2) continue;
					for (int i = 0; i < 2; i++)
					{
						if (i && model_mesh.is_boundary(te)) continue;
						auto he = model_mesh.halfedge_handle(te, i);
						he = model_mesh.opposite_halfedge_handle(model_mesh.next_halfedge_handle(model_mesh.next_halfedge_handle(he)));
						auto f1 = model_mesh.face_handle(he);
						v1 = model_mesh.to_vertex_handle(model_mesh.next_halfedge_handle(he));
						if (!model_mesh.data(v1).flag2) continue;
						auto he1 = model_mesh.next_halfedge_handle(he);
						auto f2 = model_mesh.opposite_face_handle(he1);
						bool is_remove = false;
						if (model_mesh.is_valid_handle(f2) && f2.valence() == 3)
						{
							for (auto fv : model_mesh.fv_range(f2))
							{
								if (!model_mesh.data(fv).flag2)
								{
									is_remove = true;
									break;
								}
							}
						}
						else
						{
							he1 = model_mesh.prev_halfedge_handle(he);
							f2 = model_mesh.opposite_face_handle(he1);
							if (model_mesh.is_valid_handle(f2) && f2.valence() == 3)
							{
								for (auto fv : model_mesh.fv_range(f2))
								{
									if (!model_mesh.data(fv).flag2)
									{
										is_remove = true;
										break;
									}
								}
							}
						}
						v1 = model_mesh.to_vertex_handle(model_mesh.next_halfedge_handle(he1));
						model_mesh.delete_face(f1);
						if (is_remove)
						{
							is_Ok = true;
							he1 = model_mesh.opposite_halfedge_handle(he1);
							v2 = model_mesh.to_vertex_handle(he1);
							v3 = model_mesh.from_vertex_handle(he1);
							v4 = model_mesh.to_vertex_handle(model_mesh.next_halfedge_handle(he1));
							model_mesh.delete_face(f2);
							facevhandle = { v2, v4, v1 };
							model_mesh.add_face(facevhandle);
							facevhandle = { v4, v3, v1 };
							model_mesh.add_face(facevhandle);
						}
					}
				}
				model_mesh.garbage_collection();
			}
			
			//inital is_ajusted
			for (auto v : model_mesh.vertices())
			{
				model_mesh.data(v).is_adjusted = false;
			}

		}
		
	}

	void Iso_Mesh::Open_File(std::ofstream &file_writer)
	{
		try
		{
			std::fstream fout("step_to_obj.obj", std::ios::out | std::ios::trunc);
		}
		catch (std::exception& e)
		{
			dprint("error happened:", e.what());
		}
		file_writer.open("step_to_obj.obj");
		if (file_writer.fail())
		{
			dprint("failed to open");
			exit(1);
		}
	}
}
#pragma endregion
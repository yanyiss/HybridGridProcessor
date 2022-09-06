#include "Iso_Mesh.h"
#include <fstream>
#include <iostream>
#include<algorithm>
#include "..\src\Algorithm\SurfaceMesher\Optimizer\TriangleMeshRemeshing.h"
#include "..\src\Algorithm\SurfaceMesher\Optimizer\AnisoMeshRemeshing.h"

#define USETRI

namespace CADMesher
{
	Iso_Mesh::Iso_Mesh(QString & fileName)
	{
		occ_reader = new OccReader(fileName);
#ifdef USETRI
		occ_reader->Set_TriMesh();
		MergeModel();
		ResetFeature();
#if 1
		TriangleMeshRemeshing trm(&(globalmodel.initial_trimesh));
		trm.run();
#else
		TriMesh temp(globalmodel.initial_trimesh);
		AnisoMeshRemeshing amr(&(globalmodel.initial_trimesh),&temp);
		amr.run(-1, 1.5);
#endif
		//Write_Obj(globalmodel.initial_trimesh);
#else 
		occ_reader->Set_PolyMesh();
		dprint(globalmodel.initial_polymesh.n_vertices());
		MergeModel();
		dprint(globalmodel.initial_polymesh.n_vertices());
		ResetFeature1();
		//TriangleMeshRemeshing trm(&(globalmodel.initial_polymesh));
		//trm.run();
		//Write_Obj(globalmodel.initial_polymesh);
#endif
	}

	void Iso_Mesh::MergeModel()
	{
#ifdef USETRI
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
			if (!face.if_exisited) continue;
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
			for (auto &tv : frac_mesh.vertices())
			{
				auto v = frac_mesh.point(tv);
				vh = model_mesh.add_vertex(Mesh::Point(v[0], v[1], v[2]));
				model_mesh.data(vh).GaussCurvature = frac_mesh.data(tv).GaussCurvature;
				vhandle.push_back(vh);
			}
			for (auto &tf : frac_mesh.faces())
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
			for (auto &e : frac_mesh.edges())
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

#ifdef USETRI
		globalmodel.init_trimesh_tree = new ClosestPointSearch::AABBTree(model_mesh);
#else
		TriMesh temp(model_mesh);
		globalmodel.init_trimesh_tree = new ClosestPointSearch::AABBTree(temp);
#endif // USETRI


		for (auto &te : model_mesh.edges())
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
			if (!edge0.if_exisited || edge0.if_merged || edge0.reversed_edge == -1) continue;
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

		for (auto te : model_mesh.edges())
		{
			if (model_mesh.is_boundary(te))
			{
				model_mesh.data(te).flag1 = true;
				model_mesh.data(te.v0()).set_vertflag(true);
				model_mesh.data(te.v1()).set_vertflag(true);
			}
		}

		dprint("Reset Feature Done!");
	}

	void Iso_Mesh::ResetFeature1()
	{
		Mesh &model_mesh = globalmodel.initial_polymesh;

		//set vertex flag
		for (auto te : model_mesh.edges())
		{
			if (model_mesh.data(te).flag1)
			{
				model_mesh.data(te.v0()).set_vertflag(true);
				model_mesh.data(te.v1()).set_vertflag(true);
			}
		}
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

#pragma region

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
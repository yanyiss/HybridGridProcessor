#include "Iso_Mesh.h"
#include <fstream>
#include <iostream>

namespace CADMesher
{
	Iso_Mesh::Iso_Mesh(QString & fileName)
	{
		occ_reader = new OccReader(fileName);

		occ_reader->Set_TriMesh();
#if 0
		globalmodel.initial_trimesh = occ_reader->Surface_TriMeshes[0];
#else
		MergeModel();
#endif
		Write_Obj(globalmodel.initial_trimesh);
		ResetFeature();

	}

	//void Iso_Mesh::MergeModel()
	//{
	//	using OpenMesh::Vec3d;
	//	TriMesh &model_mesh = globalmodel.initial_trimesh;
	//	model_mesh.clear();
	//	vector<ShapeFace> &faceshape = globalmodel.faceshape;
	//	vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;
	//	const auto &surface_meshes = occ_reader->Surface_TriMeshes;
	//	vector<unsigned> &triangle_surface_index = globalmodel.triangle_surface_index;
	//	int nv = 0;
	//	for (auto &frac_mesh : surface_meshes)
	//	{
	//		nv += frac_mesh.n_vertices();
	//	}
	//	triangle_surface_index.reserve(nv);
	//	int edgePointNum = 0;
	//	for (vector<ShapeEdge>::iterator edge = edgeshape.begin(); edge != edgeshape.end(); ++edge)
	//	{
	//		if (edge->reversed_edge = -1)
	//		{
	//			edgePointNum += 2 * (edge->parameters.rows() - 1);
	//		}
	//		else
	//		{
	//			edgePointNum += edge->parameters.rows() - 1;
	//		}
	//	}
	//	/*nv -= edgePointNum / 2;
	//	for (int i = 0; i < nv; ++i)
	//	{
	//		model_mesh.add_vertex(Vec3d());
	//	}*/
	//}

	void Iso_Mesh::MergeModel()
	{
		TriMesh &model_mesh = globalmodel.initial_trimesh;
		vector<ShapeFace> faceshape = globalmodel.faceshape;
		vector<ShapeEdge> edgeshape = globalmodel.edgeshape;
		model_mesh.clear();
		auto &surface_meshes = occ_reader->Surface_TriMeshes;
		for(auto &face : faceshape)
		{
			auto &wires = face.wires;
			for(auto &edges : wires)
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
		for (int i = 0; i < surface_meshes.size(); i++)
		{
			auto &frac_mesh = surface_meshes[i];
			vector<Mesh::VertexHandle> vhandle;
			vhandle.reserve(frac_mesh.n_vertices());
			for (auto tv : frac_mesh.vertices())
			{
				auto v = frac_mesh.point(tv);
				vhandle.push_back(model_mesh.add_vertex(Mesh::Point(v[0], v[1], v[2])));
			}
			for (auto tf : frac_mesh.faces())
			{
				vector<int> pos;
				for (auto tfv : frac_mesh.fv_range(tf))
					pos.push_back(tfv.idx());
				model_mesh.add_face(vhandle[pos[0]], vhandle[pos[1]], vhandle[pos[2]]);
				triangle_surface_index.push_back(i);
			}
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
		for (auto tv : model_mesh.vertices())
		{
			if (model_mesh.is_boundary(tv))
				model_mesh.data(tv).set_vertflag(true);
			else
				model_mesh.data(tv).set_vertflag(false);
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
				vector<OV> fv;
				for (auto fe : model_mesh.voh_range(model_mesh.vertex_handle(m0))) {
					if (model_mesh.data(model_mesh.edge_handle(fe)).get_edgeflag())
						fv.push_back(model_mesh.to_vertex_handle(fe));
				}
				model_mesh.collapse(model_mesh.find_halfedge(model_mesh.vertex_handle(m0), model_mesh.vertex_handle(n1)));
				for (OV v : fv) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(n1));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).set_edgeflag(true);
				}
			}
			if (m1 != n0)
			{
				vector<OV> fv;
				for (auto fe : model_mesh.voh_range(model_mesh.vertex_handle(m1))) {
					if (model_mesh.data(model_mesh.edge_handle(fe)).get_edgeflag())
						fv.push_back(model_mesh.to_vertex_handle(fe));
				}
				model_mesh.collapse(model_mesh.find_halfedge(model_mesh.vertex_handle(m1), model_mesh.vertex_handle(n0)));
				for (OV v : fv) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(n0));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).set_edgeflag(true);
				}
			}
			for (int j = 1; j < length - 1; j++)
			{
				vector<OV> fv;
				for (auto fe : model_mesh.voh_range(model_mesh.vertex_handle(m0 + j))) {
					if (model_mesh.data(model_mesh.edge_handle(fe)).get_edgeflag())
						fv.push_back(model_mesh.to_vertex_handle(fe));
				}
				model_mesh.collapse(model_mesh.find_halfedge(model_mesh.vertex_handle(m0 + j), model_mesh.vertex_handle(m1 + length - j - 1)));
				for (OV v : fv) {
					auto fe = model_mesh.find_halfedge(v, model_mesh.vertex_handle(m1 + length - j - 1));
					if (fe.is_valid())
						model_mesh.data(model_mesh.edge_handle(fe)).set_edgeflag(true);
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
		for (auto te : model_mesh.edges()) {
			if (!model_mesh.data(te).get_edgeflag()) continue;
			auto n0 = model_mesh.calc_face_normal(model_mesh.face_handle(te.h0()));
			auto n1 = model_mesh.calc_face_normal(model_mesh.face_handle(te.h1()));
			if (n0.dot(n1) > 0.8 && !model_mesh.is_boundary(te))
				model_mesh.data(te).set_edgeflag(false);
			else {
				model_mesh.data(te.v0()).set_vertflag(true);
				model_mesh.data(te.v1()).set_vertflag(true);
			}
		}
		for (auto tv : model_mesh.vertices()) {
			if (!model_mesh.data(tv).get_vertflag())
				continue;
			for (auto tve : model_mesh.ve_range(tv))
				if (model_mesh.data(tve).get_edgeflag())
					goto goto20210605;
			model_mesh.data(tv).set_vertflag(false);
		goto20210605:;
		}
		dprint("Reset Feature Done!");
	}

#pragma region
	std::string Iso_Mesh::Write_Obj(TriMesh &aMesh)
	{
		std::ofstream file_writer;
		Open_File(file_writer);

		for (auto tv = aMesh.vertices_begin(); tv != aMesh.vertices_end(); tv++)
		{
			auto pos = aMesh.point(*tv);
			file_writer << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
		}
		for (auto tf = aMesh.faces_begin(); tf != aMesh.faces_end(); tf++)
		{
			file_writer << "f";
			for (auto tfv = aMesh.fv_begin(*tf); tfv != aMesh.fv_end(*tf); tfv++)
			{
				file_writer << " " << tfv->idx() + 1;
			}
			file_writer << "\n";
		}
		file_writer.close();
		return "step_to_obj.obj";
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
#include "Iso_Mesh.h"
#include <fstream>
#include <iostream>
Iso_Mesh::Iso_Mesh(QString & fileName)
{
	occ_reader = new OccReader(fileName);
	//expected_area = pow(occ_reader->expected_edge_length, 2)*sqrt(3) / 4.0;

	occ_reader->Set_TriMesh();
	MergeModel(Tri);
	Write_Obj(globalmodel.initial_trimesh);
	ResetFeature();
	//InitTree();

	/*occ_reader->Set_PolyMesh();
	MergeModel(Quad);
	Write_Obj(globalmodel.initial_polymesh);*/

}

void Iso_Mesh::MergeModel(MergeType mt)
{
	switch (mt)
	{
	case Tri:
	{
		//{
	//	model_mesh.clear();
	//	auto &surface_meshes = occ_reader->Surface_Meshes;
	//	auto &faceshape = occ_reader->faceshape;
	//	auto &edgeshape = occ_reader->edgeshape;
	//	OpenMesh::VPropHandleT<std::pair<bool, OpenMesh::VertexHandle>> corresponding_mark;
	//	model_mesh.add_property(corresponding_mark);
	//
	//	auto faceshapeitr = faceshape.begin();
	//	for (auto &itr : surface_meshes) {
	//		int size = model_mesh.n_vertices();
	//		for(auto tv : itr.vertices()) {
	//			OpenMesh::VertexHandle v = model_mesh.add_vertex(itr.point(tv));
	//			model_mesh.property(corresponding_mark, v).first = false;
	//			//model_mesh.property(vertflag, v) = false;
	//		}
	//		for (auto tf : itr.faces()) {
	//			vector<int> id;
	//			for (auto tfv : itr.fv_range(tf)) id.push_back(tfv.idx() + size);
	//			OpenMesh::FaceHandle f = model_mesh.add_face(model_mesh.vertex_handle(id[0]),
	//				model_mesh.vertex_handle(id[1]), model_mesh.vertex_handle(id[2]));
	//			//for (auto tfe : model_mesh.fe_range(f)) model_mesh.property(edgeflag, tfe) = false;
	//		}
	//		for (auto &it : faceshapeitr->wires) {
	//			int start_id = size;
	//			for (auto &tr : it) {
	//				auto &aedge = edgeshape[tr];
	//				aedge.begin_id = size;
	//				size += aedge.parameters.rows() - 1;
	//				aedge.end_id = size;
	//			}
	//			edgeshape[*it.rbegin()].end_id = start_id;
	//			for (int id = start_id; id < size; id++) model_mesh.data(model_mesh.vertex_handle(id)).set_vertflag(true);
	//		}
	//		faceshapeitr++;
	//	}
	//	print("Vertex Number Before Merge:", model_mesh.n_vertices());
	//	typedef std::pair<bool, OpenMesh::VertexHandle> pbv;
	//	for (auto &itr : edgeshape) {
	//		if (itr.if_merged || itr.reversed_edge == -1) continue;
	//		
	//		auto &redge = edgeshape[itr.reversed_edge];
	//		itr.if_merged = true;
	//		redge.if_merged = true;
	//		int m0 = itr.begin_id;
	//		int n0 = itr.end_id;
	//		int m1 = redge.begin_id;
	//		int n1 = redge.end_id;
	//
	//		//保证每个拐点必被标记
	//		model_mesh.property(corresponding_mark, model_mesh.vertex_handle(m0++))
	//			= pbv(true, model_mesh.vertex_handle(n1));
	//		n1 = m1 + itr.parameters.rows() - 2;
	//		while (n1 > m1) {
	//			model_mesh.property(corresponding_mark, model_mesh.vertex_handle(m0++))
	//				= pbv(true, model_mesh.vertex_handle(n1--));
	//		}
	//		model_mesh.property(corresponding_mark, model_mesh.vertex_handle(m1))
	//			= pbv(true, model_mesh.vertex_handle(n0));
	//	}
	//
	//	for(auto tv : model_mesh.vertices()) {
	//		if (tv.idx() == 55680) {
	//			int p = 0;
	//		}
	//		//print("Merge Times:", tv.idx());
	//		if (!model_mesh.property(corresponding_mark, tv).first) continue;
	//		OpenMesh::VertexHandle cv = model_mesh.property(corresponding_mark, tv).second;
	//		if (model_mesh.property(corresponding_mark, cv).first) {
	//			//由tv引出的cv顺时针排列
	//			while (cv != tv) {
	//				vector<OpenMesh::VertexHandle> vertlist;
	//				for (auto tvv = model_mesh.cvv_begin(cv); tvv != model_mesh.cvv_end(cv); tvv++) vertlist.push_back(*tvv);
	//				model_mesh.property(corresponding_mark, cv).first = false;
	//				OpenMesh::VertexHandle tempvert = cv;
	//				cv = model_mesh.property(corresponding_mark, cv).second;
	//				for (int i = 1; i < vertlist.size(); i++) {
	//					model_mesh.delete_face(model_mesh.face_handle(model_mesh.find_halfedge(vertlist[i], vertlist[i - 1])));
	//					model_mesh.add_face(tv, vertlist[i], vertlist[i - 1]);
	//				}
	//				/*model_mesh.property(edgeflag, model_mesh.edge_handle(model_mesh.find_halfedge(tv, *vertlist.begin()))) = true;
	//				model_mesh.property(edgeflag, model_mesh.edge_handle(model_mesh.find_halfedge(tv, *vertlist.rbegin()))) = true;*/
	//				model_mesh.data(model_mesh.edge_handle(model_mesh.find_halfedge(tv, *vertlist.begin()))).set_edgeflag(true);
	//				model_mesh.data(model_mesh.edge_handle(model_mesh.find_halfedge(tv, *vertlist.rbegin()))).set_edgeflag(true);
	//				model_mesh.delete_vertex(tempvert);
	//				if (cv.idx() == -1) {
	//					model_mesh.property(corresponding_mark, tv).second = cv;
	//					break;
	//				}
	//			}
	//		}
	//		else {
	//			vector<OpenMesh::VertexHandle> vertlist;
	//			for (auto tvv = model_mesh.cvv_begin(cv); tvv != model_mesh.cvv_end(cv); tvv++) vertlist.push_back(*tvv);
	//			for (int i = 1; i < vertlist.size(); i++) {
	//				model_mesh.delete_face(model_mesh.face_handle(model_mesh.find_halfedge(vertlist[i], vertlist[i - 1])));
	//				model_mesh.add_face(tv, vertlist[i], vertlist[i - 1]);
	//			}
	//			/*model_mesh.property(edgeflag, model_mesh.edge_handle(model_mesh.find_halfedge(tv, *vertlist.begin()))) = true;
	//			model_mesh.property(edgeflag, model_mesh.edge_handle(model_mesh.find_halfedge(tv, *vertlist.rbegin()))) = true;*/
	//			model_mesh.data(model_mesh.edge_handle(model_mesh.find_halfedge(tv, *vertlist.begin()))).set_edgeflag(true);
	//			model_mesh.data(model_mesh.edge_handle(model_mesh.find_halfedge(tv, *vertlist.rbegin()))).set_edgeflag(true);
	//			model_mesh.delete_vertex(cv);
	//		}
	//	}
	//
	//
	//	model_mesh.remove_property(corresponding_mark);
	//	model_mesh.garbage_collection();
	//	print("Merge Done!\nVertex Number After Merge: ", model_mesh.n_vertices());
	//}
		TriMesh &model_mesh = globalmodel.initial_trimesh;
		vector<ShapeFace> faceshape = globalmodel.faceshape;
		vector<ShapeEdge> edgeshape = globalmodel.edgeshape;

		model_mesh.clear();
		auto &surface_meshes = occ_reader->Surface_TriMeshes;
		for (int i = 0; i < faceshape.size(); i++)
		{
			//dprint("f", i);
			auto &wires = faceshape[i].wires;
			for (auto it = wires.begin(); it != wires.end(); it++)
			{
				auto &edges = *it;
				for (auto itr_ = edges.begin(); itr_ != edges.end(); itr_++)
				{
					ShapeEdge &aedge = edgeshape[*itr_];
					if (aedge.reversed_edge == -1) continue;
					auto &redge = edgeshape[aedge.reversed_edge];
					if (!BRep_Tool::IsClosed(aedge.edge) && !BRep_Tool::IsClosed(redge.edge)) continue;

					int splitdump = 2;
					if (redge.if_splitted)
					{
						splitdump = aedge.parameters.rows() - 3;
					}
					itr_ = edges.insert(itr_ + 1, edgeshape.size());
					aedge.prev_edge = edgeshape.size();
					edgeshape.push_back(ShapeEdge(edgeshape.size(), aedge.edge));
					auto &newedge = edgeshape.back();
					newedge.main_face = aedge.main_face;
					newedge.secondary_face = aedge.secondary_face;
					newedge.parameters = aedge.parameters.block(splitdump, 0, aedge.parameters.rows() - splitdump, 2);
					newedge.reversed_edge = aedge.reversed_edge;
					aedge.parameters.conservativeResize(splitdump + 1, 2);
					//读M14_receiver.stp时此处会出问题 2021/09/10

					//MatrixXd t = aedge.parameters.block(0, 0, splitdump + 1, 2);
					////aedge.parameters = aedge.parameters.block(0, 0, splitdump + 1, 2);
					//if (aedge.parameters.rows() >= 47)
					//{
					//	print(aedge.parameters,"\n", t);
					//}
					//std::cout << splitdump << std::endl;
					//std::cout << aedge.parameters.rows() << std::endl;
					//aedge.parameters.resize(splitdump + 1, 2);
					//std::cout << "resize parameter" << std::endl;
					//aedge.parameters = t;
					//std::cout << "fdjksal";
					/*aedge.parameters.resize(48, 2);
					aedge.parameters.resize(47, 2);
					aedge.parameters.resize(3, 2);
					aedge.parameters.resize(40, 2);
					aedge.parameters.resize(splitdump + 1, 2);*/
					aedge.if_splitted = true;
					newedge.if_splitted = true;
					if (redge.if_splitted)
					{
						aedge.reversed_edge = redge.prev_edge;
						newedge.reversed_edge = redge.id;
						redge.reversed_edge = newedge.id;
						edgeshape[redge.prev_edge].reversed_edge = aedge.id;
					}
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
			//int i = 0;
		{
			//cout << "surfacemesh        " << i << endl;
			auto &frac_mesh = surface_meshes[i];
			vector<Mesh::VertexHandle> vhandle;
			vhandle.reserve(frac_mesh.n_vertices());
			//for (auto tv = frac_mesh.vertices_begin(); tv != frac_mesh.vertices_end(); tv++)
			for (auto tv : frac_mesh.vertices())
			{
				auto v = frac_mesh.point(tv);
				vhandle.push_back(model_mesh.add_vertex(Mesh::Point(v[0], v[1], v[2])));
				//model_mesh.data(vhandle.back()).set_vertflag(false);
			}
			//for (auto tf = frac_mesh.faces_begin(); tf != frac_mesh.faces_end(); tf++)
			for (auto tf : frac_mesh.faces())
			{
				vector<int> pos;
				/*for (auto tfv = frac_mesh.cfv_begin(tf); tfv != frac_mesh.cfv_end(tf); tfv++)
				{
					pos.push_back(tfv->idx());
				}*/
				for (auto tfv : frac_mesh.fv_range(tf))
					pos.push_back(tfv.idx());
				model_mesh.add_face(vhandle[pos[0]], vhandle[pos[1]], vhandle[pos[2]]);
				triangle_surface_index.push_back(i);
			}
			//cout << model_mesh.n_vertices() << endl;
			auto &wires = faceshape[i].wires;
			for (auto edges : wires)
			{
				int start_id = pointsnum;
				for (auto itr = edges.begin(); itr != edges.end(); itr++)
				{
					auto &aedge = edgeshape[*itr];
					aedge.begin_id = pointsnum;
					pointsnum += aedge.parameters.rows() - 1;
					aedge.end_id = itr != edges.end() - 1 ? pointsnum : start_id;
					aedge.prev_edge = itr != edges.begin() ? *(itr - 1) : edges.back();
				}
			}
			pointsnum = model_mesh.n_vertices();
		}

		for (OE te : model_mesh.edges())
		{
			if (model_mesh.is_boundary(te))
				model_mesh.data(te).set_edgeflag(true);
			else
				model_mesh.data(te).set_edgeflag(false);
		}
		for (OV tv : model_mesh.vertices())
		{
			if (model_mesh.is_boundary(tv))
				model_mesh.data(tv).set_vertflag(true);
			else
				model_mesh.data(tv).set_vertflag(false);
		}

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
			int length = edge0.parameters.rows();
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
				//auto tv = model_mesh.cvoh_begin(model_mesh.vertex_handle(m0));
				//for (; model_mesh.to_vertex_handle(*tv).idx() != n1; tv++);
				////model_mesh.property(vertflag, model_mesh.to_vertex_handle(tv)) = true;
				//model_mesh.collapse(model_mesh.halfedge_handle(tv->idx()));
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
				//auto tv = model_mesh.cvoh_begin(model_mesh.vertex_handle(m1));
				//for (; model_mesh.to_vertex_handle(*tv).idx() != n0; tv++);
				////model_mesh.property(vertflag, model_mesh.to_vertex_handle(tv)) = true;
				//model_mesh.collapse(model_mesh.halfedge_handle(tv->idx()));
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
				//auto tv = model_mesh.cvoh_begin(model_mesh.vertex_handle(m0 + j));
				//for (; model_mesh.to_vertex_handle(*tv).idx() != m1 + length - j - 1; tv++);
				////model_mesh.property(vertflag, model_mesh.to_vertex_handle(tv)) = true;
				//model_mesh.collapse(model_mesh.halfedge_handle(tv->idx()));
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
		dprint("Merge Meshes Done!");
	}
		break;
	case Quad:
	{

	}
		break;
	case Poly:
	{

	}
		break;
	default:
		break;
	}
	
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
		//if (!model_mesh.property(edgeflag, te)) continue;
		if (!model_mesh.data(te).get_edgeflag()) continue;
		auto n0 = model_mesh.calc_face_normal(model_mesh.face_handle(te.h0()));
		auto n1 = model_mesh.calc_face_normal(model_mesh.face_handle(te.h1()));
		if (n0.dot(n1) > 0.8 && !model_mesh.is_boundary(te))
			model_mesh.data(te).set_edgeflag(false);
		else {
			/*model_mesh.property(vertflag, te.v0()) = true;
			model_mesh.property(vertflag, te.v1()) = true;*/
			model_mesh.data(te.v0()).set_vertflag(true);
			model_mesh.data(te.v1()).set_vertflag(true);
		}
	}
	for (auto tv : model_mesh.vertices()) {
		/*if (!model_mesh.property(vertflag, tv)) continue;
		for (auto tve : model_mesh.ve_range(tv)) if (model_mesh.property(edgeflag, tve)) goto goto20210605;
		model_mesh.property(vertflag, tv) = false;*/
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
#pragma endregion
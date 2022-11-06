#ifndef OCCREADER_H
#define OCCREADER_H
#include "occheader.h"
#include "basic_def.h"
#include "DomainRemesh.h"
#include "..\src\Toolbox\Math\GeneralMathMethod.h"

namespace CADMesher
{
	class OccReader {
	public:
		explicit OccReader(QString &file);
		OccReader(const OccReader& or) = delete;
		~OccReader() {
			if (reader) 
			{ 
				delete reader; reader = nullptr; 
			}
		}

	protected:
		XSControl_Reader *reader;


	public:
		OpenMesh::Vec3d bbmin, bbmax;
		double expected_edge_length;
		double mu = 1.2;        //三角形面积允许扩张系数
		double epsratio = 0.005;//网格和曲面的误差系数
		int offset_quad_num = 9;
		double offset_initial_ratio = 0.02;
		double offset_increase_ratio = 1.5;
		vector<TriMesh> Surface_TriMeshes;
		vector<PolyMesh> Surface_PolyMeshes;

		void SetShape();
		void SetCADWireFrame();
		void ComputeFaceAndEdge();
		void Discrete_Edge();
		bool ProcessTangentialBoundary(int fid, int bid);
		void ClearBoundary(TriMesh &tm);
		void Face_type();
		void C0_Feature();
		void Curvature_Feature();
		void Set_TriMesh();
		void Offset_lines(Matrix2Xd &parameters, vector<Matrix2Xd> &offset_pnts, int begin, int pntnum, int quadnum);
		void Set_PolyMesh();
		void Set_Offset_Grid();
		void Re_discrete(ShapeEdge &edge, int id, int discrete_num, int quad_num, bool direction);
		bool If_decline(vector<double>& curvature);

		template<typename T0, typename T1>
		void MergeModel(T0 &model_mesh, T1 &surface_meshes)
		{
			vector<ShapeFace> faceshape = globalmodel.faceshape;
			vector<ShapeEdge> edgeshape = globalmodel.edgeshape;
			model_mesh.clear();
			for (auto& face : faceshape)
			{
				if (!face.if_exisited) continue;
				auto& wires = face.wires;
				for (auto& edges : wires)
				{
					int s = edges.size();
					for (int i = 0; i < s; ++i)
					{
						int ei = edges[i];
						if (edgeshape[ei].reversed_edge == -1)
						{
							continue;
						}
						auto& redge = edgeshape[edgeshape[ei].reversed_edge];
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

						ShapeEdge& newedge = edgeshape.back();
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
			vector<unsigned>& triangle_surface_index = globalmodel.triangle_surface_index;
			int pointsnum = 0;
			for (auto& frac_mesh : surface_meshes)
			{
				pointsnum += frac_mesh.n_vertices();
			}
			triangle_surface_index.reserve(pointsnum);
			pointsnum = 0;
			int id = 0;
			for (int i = 0; i < surface_meshes.size(); i++)
			{
				auto& frac_mesh = surface_meshes[i];
				vector<Mesh::VertexHandle> vhandle;
				Mesh::VertexHandle vh;
				vhandle.reserve(frac_mesh.n_vertices());
				for (auto& tv : frac_mesh.vertices())
				{
					auto v = frac_mesh.point(tv);
					vh = model_mesh.add_vertex(Mesh::Point(v[0], v[1], v[2]));
					model_mesh.data(vh).GaussCurvature = frac_mesh.data(tv).GaussCurvature;
					vhandle.push_back(vh);
				}
				for (auto& tf : frac_mesh.faces())
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
				for (auto& e : frac_mesh.edges())
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

				auto& wires = faceshape[i].wires;
				for (auto& edges : wires)
				{
					int start_id = pointsnum;
					for (auto itr = edges.begin(); itr != edges.end(); itr++)
					{
						auto& aedge = edgeshape[*itr];
						aedge.begin_id = pointsnum;
						pointsnum += aedge.parameters.cols() - 1;
						aedge.end_id = itr != edges.end() - 1 ? pointsnum : start_id;
						aedge.prev_edge = itr != edges.begin() ? *(itr - 1) : edges.back();
					}
				}
				pointsnum = model_mesh.n_vertices();
			}

			for (auto& te : model_mesh.edges())
			{
				if (model_mesh.is_boundary(te))
					model_mesh.data(te).set_edgeflag(true);
				else
					model_mesh.data(te).set_edgeflag(false);
			}
#if 1
			for (int i = 0; i < edgeshape.size(); i++)
			{
				auto& edge0 = edgeshape[i];
				if (!edge0.if_exisited || edge0.if_merged || edge0.reversed_edge == -1) continue;
				auto& edge1 = edgeshape[edgeshape[i].reversed_edge];
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
		int EndId(vector<ShapeEdge> &edgeshape, int edge_id);
		void SetTriFeature();
		void SetPolyFeature();

		template<typename T>
		void Set_Curvature(GeometryType * srf, T & aMesh)
		{
			double k1, k2;
			for (auto v : aMesh.vertices())
			{
				auto p = aMesh.point(v);
				if(srf->PrincipalCurvature(p[0], p[1], k1, k2))
					aMesh.data(v).GaussCurvature = std::max(std::fabs(k1), std::fabs(k2));		
				else aMesh.data(v).GaussCurvature = -1;
			}
			for (auto v : aMesh.vertices())
			{
				if (aMesh.data(v).GaussCurvature >= 0) continue;
				int count = 0;
				k1 = 0;
				for (auto vv : aMesh.vv_range(v))
				{
					if (aMesh.data(vv).GaussCurvature < 0) continue;
					k1 += aMesh.data(vv).GaussCurvature;
					count++;
				}
				if (!count)
				{
					dprint("Cannot compute the curvature of this surface");
					system("pause");
				}
				aMesh.data(v).GaussCurvature = k1 / count;
			}
		}


	//private:
		QString fileName;
		double initialRate = 0.004;
		//double degeneratedRate = 0.02;
	};

}
#endif // !OCCREADER_H


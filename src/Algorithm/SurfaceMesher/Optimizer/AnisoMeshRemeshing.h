#pragma once
#include "TriangleMeshRemeshing.h"
namespace CADMesher
{
	struct local_frame
	{
		local_frame()
			:e_x(OpenMesh::Vec3d(1, 0, 0)), e_y(OpenMesh::Vec3d(0, 1, 0)), n(OpenMesh::Vec3d(0, 0, 1))
		{}
		~local_frame() {}

		void find_e_x_y()
		{
			if (std::abs(n[2]) >= std::abs(n[1]) && std::abs(n[2]) >= std::abs(n[0]))
			{
				e_x[0] = 1.0; e_x[1] = 1.0; e_x[2] = (-n[0] - n[1]) / n[2];
			}
			else if (std::abs(n[1]) >= std::abs(n[2]) && std::abs(n[1]) >= std::abs(n[0]))
			{
				e_x[0] = 1.0; e_x[2] = 1.0; e_x[1] = (-n[0] - n[2]) / n[1];
			}
			else
			{
				e_x[1] = 1.0; e_x[2] = 1.0; e_x[0] = (-n[2] - n[1]) / n[0];
			}
			e_x.normalize();
			e_y = OpenMesh::cross(n, e_x);
		}

		OpenMesh::Vec3d e_x;
		OpenMesh::Vec3d e_y;
		OpenMesh::Vec3d n;
	};
	class AnisoMeshRemeshing
	{
	public:
		AnisoMeshRemeshing(TriMesh* mesh_, TriMesh* ref_mesh_) : mesh(mesh_), ref_mesh(ref_mesh_) { initRefHessian(); };
		~AnisoMeshRemeshing() {};
	public:
		void run(double len_, double rate_ = 1.5);
		double compute_src_mesh_ave_anisotropic_edge_length();

	private:
		void split(double cof);
		void collapse(double cof);
		void flip();
		void reposition(double step_length);

		void initRefHessian();
		void projectVertex(OV &tv, O3d &pos);
		void projectMesh();
		void projectMesh(std::vector<O3d> &pos);
		void localOptimize(int iter_num, double step_length);
		double calc_flip_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec3d& p3, const OpenMesh::Vec6d& M, bool use_area);

		ClosestPointSearch::AABBTree* aabbtree;
		TriMesh* mesh;
		TriMesh* ref_mesh;
		double len;
		double rate;
	};
}

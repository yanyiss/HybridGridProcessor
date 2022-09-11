#ifndef ANISOTROPICMESHREMESHING_H
#define ANISOTROPICMESHREMESHING_H

#include "TriangleMeshRemeshing.h"
#if 1
namespace CADMesher
{
	class AnisotropicMeshRemeshing
	{
	public:
		AnisotropicMeshRemeshing();
		~AnisotropicMeshRemeshing();

		void SetMesh(TriMesh* mesh)
		{
			reset_all_State();
			mesh_ = mesh;
		}
		void reset_all_State();

		void load_ref(TriMesh* ref_mesh);
		void sample_mesh_anisotropic_edge_length(double ref_edge_len = 1.0, double a = 1.5, bool add_flip = true);
		void do_remeshing(double ref_edge_len = 1.0, double a = 1.5);
		void build_AABB_tree_feature_edge_using_Ref();
		void calc_tri_quality();

		void LCOT_Optimize(int iter_num, double step_length);
		bool reposition_LCOT(double step_length);
		bool reposition_exp_LCOT(double step_length);
		bool flip_based_energy();
		void delete_boundary_small_tri();
		double get_ref_mesh_ave_anisotropic_edge_length() { return ref_mesh_ave_anisotropic_edge_length; }

	private:

		void project_on_reference_mesh_with_metric(TriMesh::VertexHandle vh, OpenMesh::Vec3d& p);
		void project_on_reference_edge_with_metric(TriMesh::VertexHandle vh, OpenMesh::Vec3d& p);
		void find_nearst_point_on_reference_mesh(OpenMesh::Vec3d& p, bool is_boundary);
		void project_on_reference();
		void project_on_reference_new_p(std::vector<OpenMesh::Vec3d>& np);

		bool split_one_edge(TriMesh::EdgeHandle& eh, OpenMesh::Vec3d& p);
		double calc_flip_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec3d& p3, const OpenMesh::Vec6d& M, bool use_area);
		void compute_src_mesh_ave_anisotropic_edge_length();

		TriMesh* mesh_; TriMesh* ref_mesh_;
		ClosestPointSearch::AABBTree* aabbtree;

		std::vector<int> below_30_tri;
		std::vector<double> below_30_tri_angle;
		double smallest_angle_th;

		double ref_mesh_ave_anisotropic_edge_length;

		double least_angle = PI / 90.0;
		double largest_angle = PI * 0.975;
	};
}
#endif

#if 0
namespace USTC_GCL_C
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

	class AnisotropicMeshRemeshing
	{
		//Q_OBJECT
	public:
		AnisotropicMeshRemeshing();
		~AnisotropicMeshRemeshing();

		void SetMesh(TriMesh* mesh)
		{
			reset_all_State();
			mesh_ = mesh;
		}
		void reset_all_State();

		void load_ref(GlobalGeometry* globalmodel);
		void sample_mesh_anisotropic_edge_length(double ref_edge_len = 1.0, double a = 1.5, bool add_flip = true);
		void do_remeshing(double ref_edge_len = 1.0, double a = 1.5);
		void build_AABB_tree_feature_edge_using_Ref();
		void calc_tri_quality();

		void LCOT_Optimize(int iter_num, double step_length);
		bool reposition_LCOT(double step_length);
		bool reposition_exp_LCOT(double step_length);
		bool flip_based_energy();

		//void uniform_optimize(int iter_num, double step_length);

		//void exp_mips_optimize(int iter_num, double area_angle_ratio, double energy_power);
		//void reposition_exp_mips(double area_angle_ratio, double energy_power);

		//void reposition_particle(double step_length);
		//void flip_based_particle_energy();

		void delete_boundary_small_tri();
		double get_ref_mesh_ave_anisotropic_edge_length() { return ref_mesh_ave_anisotropic_edge_length; }

	private:

		void project_on_reference_mesh_with_metric(TriMesh::VertexHandle vh, OpenMesh::Vec3d& p);
		void project_on_reference_edge_with_metric(TriMesh::VertexHandle vh, OpenMesh::Vec3d& p);
		void find_nearst_point_on_reference_mesh(OpenMesh::Vec3d& p, bool is_boundary);
		void project_on_reference();
		void project_on_reference_new_p(std::vector<OpenMesh::Vec3d>& np);
		void project_on_reference(OpenMesh::Vec3d& p, OpenMesh::Vec3d& sp, OpenMesh::Vec3d& dir, double& dis);

		bool split_one_edge(TriMesh::EdgeHandle& eh, OpenMesh::Vec3d& p);
		double calc_flip_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec3d& p3, const OpenMesh::Vec6d& M, bool use_area);
		//double calc_flip_particle_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec6d& M);

		TriMesh* mesh_;
		GlobalGeometry* referenceModel;
		/*std::vector<CGAL_3_Triangle> triangle_vectors;
		CGAL_AABB_Tree* AABB_tree;
		std::vector<CGAL_3_Segment> segment_vectors;
		std::vector<unsigned> segment_edge_id;
		CGAL_AABB_Segment_Tree* AABB_Segment_tree;*/
		//ClosestPointSearch::AABBTree* aabbtree;
		ParallelTools::AABBTree* aabbtree;

		std::vector<int> below_30_tri;
		std::vector<double> below_30_tri_angle;
		double smallest_angle_th;

		void compute_src_mesh_ave_anisotropic_edge_length();
		double ref_mesh_ave_anisotropic_edge_length;

		//double compute_exp_mips_area_energy(Mesh::VertexHandle vh, OpenMesh::Vec3d& np,
		//	const std::vector<OpenMesh::Vec6d>& vH, double area_angle_ratio, double energy_power);

		double least_angle = PI / 90.0;
		double largest_angle = PI * 0.975;
	};
}
#endif

#endif
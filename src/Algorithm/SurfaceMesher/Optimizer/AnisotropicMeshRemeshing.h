#ifndef ANISOTROPICMESHING_INTERFACE_H
#define ANISOTROPICMESHING_INTERFACE_H

//#include <QObject>
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

	class SegmentTree
	{
	public:
		SegmentTree(TriMesh* m)
		{
			for (auto te = m->edges_begin(); te != m->edges_end(); ++te)
			{
				if (m->data(te.handle()).get_edgeflag()) 
					hhs.push_back(m->halfedge_handle(te.handle(), 0));
			}
			mesh = m;
		}
		~SegmentTree() {};
	public:
		void closest_point_and_edgeid(O3d &in, O3d &out, int &eid)
		{
			double dis = DBL_MAX;
			eid = -1; out = in;
			O3d p01, pin;
			double proj, norm;
			for (auto &hh : hhs)
			{
				auto &p0 = mesh->point(mesh->from_vertex_handle(hh));
				auto &p1 = mesh->point(mesh->to_vertex_handle(hh));
				p01 = (p1 - p0).normalize();
				pin = in - p0;
				proj = pin.dot(p01);
				if (proj <= 0)
				{
					norm = pin.norm();
					if (norm < dis)
					{
						dis = norm;
						eid = hh.idx() / 2;
						out = p0;
					}
					continue;
				}
				pin = in - p1;
				proj = pin.dot(p01);
				if (proj >= 0)
				{
					norm = pin.norm();
					if (norm < dis)
					{
						dis = norm;
						eid = hh.idx() / 2;
						out = p1;
					}
					continue;
				}
				pin = p1 + proj * p01;
				norm = (pin - in).norm();
				if (norm < dis)
				{
					dis = norm;
					eid = hh.idx() / 2;
					out = pin;
				}
			}
		}
	private:
		std::vector<OH> hhs;
		TriMesh* mesh;
	};

#define PI 3.1415926535897932
	class AnisotropicMeshRemeshing// : public QObject
	{
		//Q_OBJECT
	public:
		AnisotropicMeshRemeshing();
		~AnisotropicMeshRemeshing();
		double MQE;
		double AQE;
		double MinAngle;
		double MaxAngle;
		double AveAngle;
		double MinArea;
		double MaxArea;
		double MinRER;
		double MaxRER;

		void SetMesh(TriMesh* mesh)
		{
			reset_all_State();
			mesh_ = mesh;
			initMeshStatusAndNormal(*mesh_);
		}
		void reset_all_State();

		void set_metric(OpenMesh::VertexHandle vh, OpenMesh::Vec6d& metric);
		void load_ref_mesh(TriMesh* aniso_ref_mesh, double tl, double model_size);
		void sample_mesh_anisotropic_edge_length(double ref_edge_len = 1.0, double a = 1.5, bool add_flip = true);
		void do_remeshing(double ref_edge_len = 1.0, double a = 1.5);
		void calc_tri_quality();

		void LCOT_Optimize(int iter_num, double step_length);
		bool reposition_LCOT(double step_length);
		bool reposition_exp_LCOT(double step_length);
		bool flip_based_energy();

		void uniform_optimize(int iter_num, double step_length);

		void exp_mips_optimize(int iter_num, double area_angle_ratio, double energy_power);
		void reposition_exp_mips(double area_angle_ratio, double energy_power);

		void reposition_particle(double step_length);
		void flip_based_particle_energy();

		void delete_boundary_small_tri();
		void set_draw_small_tri_ok(bool ok) { draw_small_tri_ok = ok; /*emit updateGL_Manual_signal();*/ };
		void set_smallest_angle_th(double th) { smallest_angle_th = th; /*emit updateGL_Manual_signal();*/ };
		double get_ref_mesh_ave_anisotropic_edge_length() { return ref_mesh_ave_anisotropic_edge_length; }


	private:
		bool draw_small_tri_ok;
		void project_on_reference_mesh_with_metric(Mesh::VertexHandle vh, OpenMesh::Vec3d& p);
		void project_on_reference_edge_with_metric(Mesh::VertexHandle vh, OpenMesh::Vec3d& p);
		void find_nearst_point_on_reference_mesh(OpenMesh::Vec3d& p, bool is_boundary);
		void project_on_reference();
		void project_on_reference_new_p(std::vector<OpenMesh::Vec3d>& np);
		void project_on_reference(OpenMesh::Vec3d& p, OpenMesh::Vec3d& sp, OpenMesh::Vec3d& dir, double& dis);

		bool split_one_edge(Mesh::EdgeHandle& eh, OpenMesh::Vec3d& p);
		void collapse_one_edge(OH hh);
		double calc_flip_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec3d& p3, const OpenMesh::Vec6d& M, bool use_area);
		double calc_flip_particle_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec6d& M);

		TriMesh* mesh_; TriMesh* ref_mesh_;
		ClosestPointSearch::AABBTree* aabbtree;
		SegmentTree* segtree;
		double min_cur;
		std::vector<int> below_30_tri;
		std::vector<double> below_30_tri_angle;
		double smallest_angle_th;

		void compute_src_mesh_ave_anisotropic_edge_length();
		double ref_mesh_ave_anisotropic_edge_length;

		double compute_exp_mips_area_energy(Mesh::VertexHandle vh, OpenMesh::Vec3d& np,
			const std::vector<OpenMesh::Vec6d>& vH, double area_angle_ratio, double energy_power);

		double least_angle = PI / 90.0;
		double largest_angle = PI * 0.975;
	};
}
#endif
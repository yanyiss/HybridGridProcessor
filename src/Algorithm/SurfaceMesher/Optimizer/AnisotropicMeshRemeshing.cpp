#include "AnisotropicMeshRemeshing.h"
#include <Eigen/Dense>
#if 1
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
namespace CADMesher
{
#define USE_FEATURE
#define USE_PROMOTION
#define USE_FEWERITERATION

	AnisotropicMeshRemeshing::AnisotropicMeshRemeshing()
	{
		//referenceModel = nullptr;
		aabbtree = nullptr;
		mesh_ = nullptr;
		smallest_angle_th = 1.0;
	}

	AnisotropicMeshRemeshing::~AnisotropicMeshRemeshing()
	{
	}

	void AnisotropicMeshRemeshing::reset_all_State()
	{
		below_30_tri.clear(); below_30_tri_angle.clear();
	}

	void AnisotropicMeshRemeshing::project_on_reference_mesh_with_metric(Mesh::VertexHandle vh, OpenMesh::Vec3d& p)
	{
		OpenMesh::Vec3d tp = mesh_->point(vh);
		auto voh_it = mesh_->voh_iter(vh);
		OpenMesh::Vec3d tq = mesh_->point(mesh_->to_vertex_handle(voh_it));
		OpenMesh::Vec3d p_ = tp;
		//auto iid = closestFaceAndPoint(p_, aabbtree);
		auto iid = aabbtree->closest_distance_and_face_handle(p_).second.idx();

#ifdef USE_PROMOTION
		double l_max = 0;
		for (auto tvv : mesh_->vv_range(vh))
			l_max = std::max(l_max, (mesh_->point(tvv) - tp).norm());
		if (l_max < (p_ - tp).norm() / 3.0)
			return;
#endif // USE_PROMOTION

		if ((p_ - tp).sqrnorm() > (tq - tp).sqrnorm())
		{
			printf("------------------------------------------------------------\n");
			printf("%d\n", vh.idx());
			printf("%f %f %f\n", p[0], p[1], p[2]);
			printf("%f %f %f\n", p_[0], p_[1], p_[2]);
			printf("%f %f %f\n", tp[0], tp[1], tp[2]);
			printf("------------------------------------------------------------\n");
			return;
		}
		mesh_->set_point(vh, p_);
		//auto ref_mesh_ = &(referenceModel->initial_trimesh);
		auto fv_it = ref_mesh_->cfv_iter(ref_mesh_->face_handle(iid));
		auto p0 = ref_mesh_->point(fv_it.handle());
		auto& h0 = ref_mesh_->data(fv_it).get_Hessian();
		auto p1 = ref_mesh_->point((++fv_it).handle());
		auto& h1 = ref_mesh_->data(fv_it).get_Hessian();
		auto p2 = ref_mesh_->point((++fv_it).handle());
		auto& h2 = ref_mesh_->data(fv_it).get_Hessian();

		OpenMesh::Vec3d bc;
		if (!baryCoord(p_, p0, p1, p2, bc))
			bc[0] = bc[1] = bc[2] = 1.0 / 3.0;

		OpenMesh::Vec6d h = bc[0] * h0 + bc[1] * h1 + bc[2] * h2;
		mesh_->data(vh).set_Hessian(h);
	}

	void AnisotropicMeshRemeshing::build_AABB_tree_feature_edge_using_Ref()
	{
		//if (ref_mesh_->n_edges() == 0) return;
//		std::vector<CGAL_double_3_Point> v_pos(ref_mesh_->n_vertices()); OpenMesh::Vec3d p;
//		for (Mesh::VertexIter v_it = ref_mesh_->vertices_begin(); v_it != ref_mesh_->vertices_end(); ++v_it)
//		{
//			int vertex_id = v_it.handle().idx();
//			p = ref_mesh_->point(v_it);
//			v_pos[vertex_id] = CGAL_double_3_Point(p[0], p[1], p[2]);
//		}
//
//		segment_vectors.clear(); segment_vectors.reserve(ref_mesh_->n_edges());
//		segment_edge_id.clear(); segment_edge_id.reserve(ref_mesh_->n_edges());
//		for (Mesh::EdgeIter e_it = ref_mesh_->edges_begin(); e_it != ref_mesh_->edges_end(); ++e_it)
//		{
//#ifdef USE_FEATURE
//			if (ref_mesh_->data(e_it).get_edgeflag())
//#else
//			if (ref_mesh_->is_boundary(e_it))
//#endif // USE_FEATURE
//			{
//				Mesh::HalfedgeHandle heh = ref_mesh_->halfedge_handle(e_it, 0);
//				int v0 = ref_mesh_->from_vertex_handle(heh).idx();
//				int v1 = ref_mesh_->to_vertex_handle(heh).idx();
//				segment_vectors.push_back(CGAL_3_Segment(v_pos[v0], v_pos[v1]));
//				segment_edge_id.push_back(e_it.handle().idx());
//			}
//		}
//
//		if (segment_vectors.size() > 0)
//		{
//			if (AABB_Segment_tree) delete AABB_Segment_tree;
//			AABB_Segment_tree = new CGAL_AABB_Segment_Tree(segment_vectors.begin(), segment_vectors.end());
//			AABB_Segment_tree->accelerate_distance_queries();
//
//			printf("------------------------------------------------------------\n");
//			printf("Build AABB Tree for feature edge.\n");
//		}

	}

	void AnisotropicMeshRemeshing::project_on_reference_edge_with_metric(Mesh::VertexHandle vh, OpenMesh::Vec3d& p)
	{
#if 0
		CGAL_AABB_Segment_Tree::Point_and_primitive_id point_primitive = AABB_Segment_tree->closest_point_and_primitive(CGAL_double_3_Point(p[0], p[1], p[2]));
		CGAL_double_3_Point pos = point_primitive.first;
		CGAL_Segment_Iterator it = point_primitive.second;
		unsigned edge_vector_id = std::distance(segment_vectors.begin(), it);
		unsigned edge_id = segment_edge_id[edge_vector_id];
		OpenMesh::Vec3d p_ = OpenMesh::Vec3d(pos.x(), pos.y(), pos.z());
		OpenMesh::Vec3d tp = mesh_->point(vh);
		Mesh::VertexOHalfedgeIter voh_it = mesh_->voh_iter(vh);
		OpenMesh::Vec3d tq = mesh_->point(mesh_->to_vertex_handle(voh_it));

#ifdef USE_PROMOTION
		double l_max = 0;
		for (auto tvv : mesh_->vv_range(vh))
			l_max = std::max(l_max, (mesh_->point(tvv) - tp).norm());
		if (l_max < (p_ - tp).norm() / 3.0)
			return;
#endif // USE_PROMOTION

		if ((p_ - tp).sqrnorm() > (tq - tp).sqrnorm())
		{
			printf("------------------------------------------------------------\n");
			printf("%d\n", vh.idx());
			printf("%f %f %f\n", p[0], p[1], p[2]);
			printf("%f %f %f\n", p_[0], p_[1], p_[2]);
			printf("%f %f %f\n", tp[0], tp[1], tp[2]);
			printf("------------------------------------------------------------\n");
		}
		mesh_->set_point(vh, p_);

		Mesh::EdgeHandle eh = ref_mesh_->edge_handle(edge_id);
		Mesh::HalfedgeHandle heh = ref_mesh_->halfedge_handle(eh, 0);
		Mesh::VertexHandle v0 = ref_mesh_->from_vertex_handle(heh);
		Mesh::VertexHandle v1 = ref_mesh_->to_vertex_handle(heh);
		Mesh::Point p0 = ref_mesh_->point(v0);
		OpenMesh::Vec6d& h0 = ref_mesh_->data(v0).get_Hessian();
		Mesh::Point p1 = ref_mesh_->point(v1);
		OpenMesh::Vec6d& h1 = ref_mesh_->data(v1).get_Hessian();

		double len = (p0 - p1).norm();
		double bc0 = (p_ - p1).norm() / len; double bc1 = 1.0 - bc0;
		mesh_->data(vh).set_Hessian(bc0*h0 + bc1 * h1);
#endif
	}

	void AnisotropicMeshRemeshing::find_nearst_point_on_reference_mesh(OpenMesh::Vec3d& p, bool is_boundary)
	{
		if (is_boundary)
		{
			/*CGAL_double_3_Point pos = AABB_Segment_tree->closest_point(CGAL_double_3_Point(p[0], p[1], p[2]));
			p = OpenMesh::Vec3d(pos.x(), pos.y(), pos.z());*/
		}
		else
		{
			//closestPoint(p, aabbtree);
			p = aabbtree->closest_point(p);
		}
	}

	void AnisotropicMeshRemeshing::project_on_reference()
	{
		unsigned nv = mesh_->n_vertices(); OpenMesh::Vec3d p;
		for (unsigned int i = 0; i < nv; ++i)
		{
			auto vh = mesh_->vertex_handle(i);
			p = mesh_->point(vh);

#ifdef USE_FEATURE
			if (mesh_->data(vh).get_vertflag())
#else
			if (mesh_->is_boundary(vh))
#endif // USE_FEATURE
			{
				project_on_reference_edge_with_metric(vh, p);
			}
			else
			{
				project_on_reference_mesh_with_metric(vh, p);
			}
		}
	}

	void AnisotropicMeshRemeshing::project_on_reference_new_p(std::vector<OpenMesh::Vec3d>& np)
	{
		unsigned nv = mesh_->n_vertices();
		for (unsigned int i = 0; i < nv; ++i)
		{
			auto vh = mesh_->vertex_handle(i);
#ifdef USE_FEATURE
			if (mesh_->data(vh).get_vertflag())
#else
			if (mesh_->is_boundary(vh))
#endif // USE_FEATURE
			{
				project_on_reference_edge_with_metric(vh, np[i]);
			}
			else
			{
				project_on_reference_mesh_with_metric(vh, np[i]);
			}
		}
	}

	void AnisotropicMeshRemeshing::load_ref(TriMesh* ref_mesh)
	{
		//aabbtree = globalmodel->initial_tree;
		aabbtree = globalmodel.init_trimesh_tree;
		//auto ref_mesh_ = &(globalmodel->initial_trimesh);
		ref_mesh_ = ref_mesh;
		ref_mesh_->request_vertex_status();
		ref_mesh_->request_edge_status();
		ref_mesh_->request_face_status();
		ref_mesh_->request_face_normals();
		ref_mesh_->request_vertex_normals();

		mesh_->request_vertex_status();
		mesh_->request_edge_status();
		mesh_->request_face_status();
		mesh_->request_face_normals();
		mesh_->request_vertex_normals();

		std::vector<double> K1, K2; std::vector<OpenMesh::Vec3d> D1, D2;
		compute_principal_curvature(ref_mesh_, K1, K2, D1, D2);

		int nv = ref_mesh_->n_vertices();
		Eigen::Matrix3d H; Eigen::Matrix3d D; D.setZero();
		std::vector<Eigen::Matrix3d> vH(nv); OpenMesh::Vec6d h;
		for (unsigned int i = 0; i < nv; ++i)
		{
			auto vh = ref_mesh_->vertex_handle(i);
			double k1 = K1[i]; k1 = std::abs(k1) < 1.0e-4 ? 1.0e-4 : k1;
			double k2 = K2[i]; k2 = std::abs(k2) < 1.0e-4 ? 1.0e-4 : k2;

			OpenMesh::Vec3d d1 = D1[i];
			OpenMesh::Vec3d d2 = D2[i];

			OpenMesh::Vec3d n = OpenMesh::cross(d1, d2).normalize();
			H(0, 0) = d1[0]; H(1, 0) = d1[1]; H(2, 0) = d1[2];
			H(0, 1) = d2[0]; H(1, 1) = d2[1]; H(2, 1) = d2[2];
			H(0, 2) = n[0]; H(1, 2) = n[1]; H(2, 2) = n[2];
			D(0, 0) = std::abs(k1); D(1, 1) = std::abs(k2);// D(2,2) = std::abs(k2) < std::abs(k1) ? std::abs(k2) : std::abs(k1);
			vH[i] = H * D * H.transpose();

			h[0] = vH[i](0, 0); h[1] = vH[i](0, 1); h[2] = vH[i](0, 2);
			h[3] = vH[i](1, 1); h[4] = vH[i](1, 2); h[5] = vH[i](2, 2);
			ref_mesh_->data(vh).set_Hessian(h);
		}
		build_AABB_tree_feature_edge_using_Ref();
		project_on_reference();
		calc_tri_quality();
		compute_src_mesh_ave_anisotropic_edge_length();
	}

	void AnisotropicMeshRemeshing::compute_src_mesh_ave_anisotropic_edge_length()
	{
		Eigen::Matrix3d H, H0, H1; Eigen::Vector3d P0, P1, P;
		Eigen::Matrix3d U; Eigen::Matrix3d V; Eigen::Vector3d sv; Eigen::Matrix3d diag_a; diag_a.setZero();
		Eigen::Matrix3d Q;
		unsigned ne = mesh_->n_edges(); double sum_edge_len = 0.0;
		for (unsigned i = 0; i < ne; ++i)
		{
			auto eh = mesh_->edge_handle(i);
			auto heh = mesh_->halfedge_handle(eh, 0);
			auto vh0_ = mesh_->from_vertex_handle(heh);
			auto vh1_ = mesh_->to_vertex_handle(heh);
			OpenMesh::Vec6d& h0 = mesh_->data(vh0_).get_Hessian();
			OpenMesh::Vec6d& h1 = mesh_->data(vh1_).get_Hessian();
			H0(0, 0) = h0[0]; H0(0, 1) = h0[1]; H0(0, 2) = h0[2];
			H0(1, 0) = h0[1]; H0(1, 1) = h0[3]; H0(1, 2) = h0[4];
			H0(2, 0) = h0[2]; H0(2, 1) = h0[4]; H0(2, 2) = h0[5];
			H1(0, 0) = h1[0]; H1(0, 1) = h1[1]; H1(0, 2) = h1[2];
			H1(1, 0) = h1[1]; H1(1, 1) = h1[3]; H1(1, 2) = h1[4];
			H1(2, 0) = h1[2]; H1(2, 1) = h1[4]; H1(2, 2) = h1[5];
			H = 0.5 * (H0 + H1);
			OpenMesh::Vec3d p0 = mesh_->point(vh0_);
			OpenMesh::Vec3d p1 = mesh_->point(vh1_);
			P0(0) = p0[0]; P0(1) = p0[1]; P0(2) = p0[2];
			P1(0) = p1[0]; P1(1) = p1[1]; P1(2) = p1[2];
			Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
			U = svd.matrixU(); V = svd.matrixV(); sv = svd.singularValues();
			diag_a(0, 0) = std::sqrt(sv(0)); diag_a(1, 1) = std::sqrt(sv(1)); diag_a(2, 2) = std::sqrt(sv(2));
			Q = U * diag_a*V.transpose();
			P0 = Q * P0; P1 = Q * P1;
			double edge_len = std::sqrt((P0(0) - P1(0))*(P0(0) - P1(0)) + (P0(1) - P1(1))*(P0(1) - P1(1)) + (P0(2) - P1(2))*(P0(2) - P1(2)));
			sum_edge_len += edge_len;
		}
		ref_mesh_ave_anisotropic_edge_length = sum_edge_len / ne;
		printf("------------------------------------------------------------\n");
		printf("Src mesh ave anisotropic edge length : %f\n", ref_mesh_ave_anisotropic_edge_length);
	}

	bool AnisotropicMeshRemeshing::split_one_edge(TriMesh::EdgeHandle& eh, OpenMesh::Vec3d& p)
	{
		auto heh0 = mesh_->halfedge_handle(eh, 0);
		auto heh1 = mesh_->halfedge_handle(eh, 1);
		auto vh0 = mesh_->to_vertex_handle(heh0); OpenMesh::Vec3d p0 = mesh_->point(vh0);
		auto vh1 = mesh_->to_vertex_handle(heh1); OpenMesh::Vec3d p1 = mesh_->point(vh1);

		std::vector<Mesh::VertexHandle> one_face(3);
		bool flag1 = mesh_->data(eh).flag1;
		bool flag2 = mesh_->data(eh).flag2;
		
		if (mesh_->is_boundary(eh))
		{
			if (TriMesh::InvalidFaceHandle != mesh_->face_handle(heh0))
			{
				auto vh2 = mesh_->to_vertex_handle(mesh_->next_halfedge_handle(heh0));
				OpenMesh::Vec3d p2 = mesh_->point(vh2);
				OpenMesh::Vec3d n = OpenMesh::cross(p1 - p2, p0 - p2).normalize();
				double a1 = OpenMesh::dot(n, OpenMesh::cross(p2 - p0, p - p0));
				double a2 = OpenMesh::dot(n, OpenMesh::cross(p1 - p2, p - p2));
				if (a1 < 1e-8 || a2 < 1e-8) return false;
				auto vh = mesh_->add_vertex(p);
				mesh_->split_edge(eh, vh);
				if (flag1 || flag2) {
					mesh_->data(vh).set_vertflag(true);
					auto hh = mesh_->edge_handle(mesh_->find_halfedge(vh, vh0));
					mesh_->data(hh).flag1 = flag1; mesh_->data(hh).flag2 = flag2;
					hh = mesh_->edge_handle(mesh_->find_halfedge(vh, vh1));
					mesh_->data(hh).flag1 = flag1; mesh_->data(hh).flag2 = flag2;
				}
			}
			else
			{
				auto vh3 = mesh_->to_vertex_handle(mesh_->next_halfedge_handle(heh1));
				OpenMesh::Vec3d p3 = mesh_->point(vh3);
				OpenMesh::Vec3d n = OpenMesh::cross(p0 - p3, p1 - p3).normalize();
				double a1 = OpenMesh::dot(n, OpenMesh::cross(p0 - p3, p - p3));
				double a2 = OpenMesh::dot(n, OpenMesh::cross(p3 - p1, p - p1));
				if (a1 < 1e-8 || a2 < 1e-8) return false;
				auto vh = mesh_->add_vertex(p);
				mesh_->split_edge(eh, vh);
				if (flag1 || flag2) {
					mesh_->data(vh).set_vertflag(true);
					auto hh = mesh_->edge_handle(mesh_->find_halfedge(vh, vh0));
					mesh_->data(hh).flag1 = flag1; mesh_->data(hh).flag2 = flag2;
					hh = mesh_->edge_handle(mesh_->find_halfedge(vh, vh1));
					mesh_->data(hh).flag1 = flag1; mesh_->data(hh).flag2 = flag2;
				}
			}
		}
		else
		{
			auto vh2 = mesh_->to_vertex_handle(mesh_->next_halfedge_handle(heh0)); OpenMesh::Vec3d p2 = mesh_->point(vh2);
			OpenMesh::Vec3d n1 = OpenMesh::cross(p1 - p2, p0 - p2).normalize();
			auto vh3 = mesh_->to_vertex_handle(mesh_->next_halfedge_handle(heh1)); OpenMesh::Vec3d p3 = mesh_->point(vh3);
			OpenMesh::Vec3d n2 = OpenMesh::cross(p0 - p3, p1 - p3).normalize();
			double a1 = OpenMesh::dot(n1, OpenMesh::cross(p2 - p0, p - p0));
			double a2 = OpenMesh::dot(n1, OpenMesh::cross(p1 - p2, p - p2));
			double a3 = OpenMesh::dot(n2, OpenMesh::cross(p0 - p3, p - p3));
			double a4 = OpenMesh::dot(n2, OpenMesh::cross(p3 - p1, p - p1));
			if (a1 < 1e-8 || a2 < 1e-8 || a3 < 1e-8 || a4 < 1e-8) return false;
			auto vh = mesh_->add_vertex(p);
			mesh_->split_edge(eh, vh);
			if (flag1 || flag2) {
				mesh_->data(vh).set_vertflag(true);
				auto hh = mesh_->edge_handle(mesh_->find_halfedge(vh, vh0));
				mesh_->data(hh).flag1 = flag1; mesh_->data(hh).flag2 = flag2;
				hh = mesh_->edge_handle(mesh_->find_halfedge(vh, vh1));
				mesh_->data(hh).flag1 = flag1; mesh_->data(hh).flag2 = flag2;
			}
		}
		return true;
	}

	void AnisotropicMeshRemeshing::sample_mesh_anisotropic_edge_length(double ref_edge_len, double a, bool add_flip)
	{
		//auto ref_mesh_ = &(referenceModel->initial_trimesh);
		if (ref_mesh_ == NULL) return;

		unsigned nv = mesh_->n_vertices();
		project_on_reference();

		unsigned ne = mesh_->n_edges();
		Eigen::Matrix3d H, H0, H1; Eigen::Vector3d P0, P1, P;
		Eigen::Matrix3d U; Eigen::Matrix3d V; Eigen::Vector3d sv; Eigen::Matrix3d diag_a; diag_a.setZero();
		Eigen::Matrix3d Q; std::vector<double> h(9);
		int iter_num = 0; int split_count = 0; int collapse_count = 0;
		int no_split_collapse = 0;

		std::cout << "Sampling time" << std::endl;

		for (unsigned ii = 0; ii < 5; ++ii)
		{
			no_split_collapse = 0;

			if (add_flip)
			{
				flip_based_energy();
			}

			split_count = 0; ne = mesh_->n_edges(); nv = mesh_->n_vertices();
			for (unsigned i = 0; i < ne; )
			{
				auto eh = mesh_->edge_handle(i);
				auto heh = mesh_->halfedge_handle(eh, 0);
				auto vh0_ = mesh_->from_vertex_handle(heh);
				auto vh1_ = mesh_->to_vertex_handle(heh);
				OpenMesh::Vec6d& h0 = mesh_->data(vh0_).get_Hessian();
				OpenMesh::Vec6d& h1 = mesh_->data(vh1_).get_Hessian();
				H0(0, 0) = h0[0]; H0(0, 1) = h0[1]; H0(0, 2) = h0[2];
				H0(1, 0) = h0[1]; H0(1, 1) = h0[3]; H0(1, 2) = h0[4];
				H0(2, 0) = h0[2]; H0(2, 1) = h0[4]; H0(2, 2) = h0[5];
				H1(0, 0) = h1[0]; H1(0, 1) = h1[1]; H1(0, 2) = h1[2];
				H1(1, 0) = h1[1]; H1(1, 1) = h1[3]; H1(1, 2) = h1[4];
				H1(2, 0) = h1[2]; H1(2, 1) = h1[4]; H1(2, 2) = h1[5];
				H = 0.5 * (H0 + H1);
				OpenMesh::Vec3d p0 = mesh_->point(vh0_);
				OpenMesh::Vec3d p1 = mesh_->point(vh1_);
				P0(0) = p0[0]; P0(1) = p0[1]; P0(2) = p0[2];
				P1(0) = p1[0]; P1(1) = p1[1]; P1(2) = p1[2];

				double edge_len = std::sqrt((P0 - P1).transpose()*H*(P0 - P1));
				if (edge_len > ref_edge_len * a)
				{
#ifdef USE_FEATURE
					bool is_boundary_edge = mesh_->data(eh).get_edgeflag();
#else
					bool is_boundary_edge = mesh_->is_boundary(eh);
#endif // USE_FEATURE
					OpenMesh::Vec3d p = 0.5*(p0 + p1);
					find_nearst_point_on_reference_mesh(p, is_boundary_edge);
					if (!split_one_edge(eh, p)) { ++i; continue; };
					if (is_boundary_edge)//feature and not boundary
					{
						auto vh_ = mesh_->vertex_handle(nv);
						project_on_reference_edge_with_metric(vh_, mesh_->point(vh_));
					}
					else
					{
						auto vh_ = mesh_->vertex_handle(nv);
						project_on_reference_mesh_with_metric(vh_, mesh_->point(vh_));
					}

					nv = mesh_->n_vertices(); ne = mesh_->n_edges();
					++split_count;
				}
				else
				{
					++i;
				}
			}

			if (split_count == 0) { no_split_collapse += 1; }
			else { printf("Split 0 : %d\n", split_count); }
			ne = mesh_->n_edges(); printf("Split Edges : %d\n", ne);

			//if(ii == 0) break;
			if (add_flip)
			{
				flip_based_energy();
			}
			iter_num = 0;
			while (iter_num < 10)
			{
				collapse_count = 0; ne = mesh_->n_edges();
				for (unsigned i = 0; i < ne; )
				{
					auto eh = mesh_->edge_handle(i);
					auto heh = mesh_->halfedge_handle(eh, 0);
					auto vh0_ = mesh_->from_vertex_handle(heh);
					auto vh1_ = mesh_->to_vertex_handle(heh);
					OpenMesh::Vec6d& h0 = mesh_->data(vh0_).get_Hessian();
					OpenMesh::Vec6d& h1 = mesh_->data(vh1_).get_Hessian();
					H0(0, 0) = h0[0]; H0(0, 1) = h0[1]; H0(0, 2) = h0[2];
					H0(1, 0) = h0[1]; H0(1, 1) = h0[3]; H0(1, 2) = h0[4];
					H0(2, 0) = h0[2]; H0(2, 1) = h0[4]; H0(2, 2) = h0[5];
					H1(0, 0) = h1[0]; H1(0, 1) = h1[1]; H1(0, 2) = h1[2];
					H1(1, 0) = h1[1]; H1(1, 1) = h1[3]; H1(1, 2) = h1[4];
					H1(2, 0) = h1[2]; H1(2, 1) = h1[4]; H1(2, 2) = h1[5];
					H = 0.5 * (H0 + H1);
					OpenMesh::Vec3d p0 = mesh_->point(vh0_);
					OpenMesh::Vec3d p1 = mesh_->point(vh1_);
					P0(0) = p0[0]; P0(1) = p0[1]; P0(2) = p0[2];
					P1(0) = p1[0]; P1(1) = p1[1]; P1(2) = p1[2];
					double edge_len = std::sqrt((P0 - P1).transpose()*H*(P0 - P1));
					/*Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV );
					U = svd.matrixU(); V = svd.matrixV(); sv = svd.singularValues();
					diag_a(0,0) = std::sqrt(sv(0)); diag_a(1,1) = std::sqrt(sv(1)); diag_a(2,2) = std::sqrt(sv(2));
					Q = U*diag_a*V.transpose();
					P0 = Q*P0; P1 = Q*P1;
					double edge_len_ = std::sqrt( (P0(0) - P1(0))*(P0(0) - P1(0)) + (P0(1) - P1(1))*(P0(1) - P1(1)) + (P0(2) - P1(2))*(P0(2) - P1(2)) );
					printf("%f %f\n", edge_len, edge_len_);*/
					if (edge_len < ref_edge_len / a)
					{
#ifdef USE_FEATURE
						if (!mesh_->data(eh).get_edgeflag())
#else
						if (!mesh_->is_boundary(eh))
#endif // USE_FEATURE
						{
#ifdef USE_FEATURE
							if (mesh_->is_collapse_ok(heh) && !mesh_->data(vh0_).get_vertflag())
#else
							if (mesh_->is_collapse_ok(heh) && !mesh_->is_boundary(vh0_)) //from vh0 to vh1
#endif // USE_FEATURE
							{
#ifdef USE_FEATURE
								if (!mesh_->data(vh1_).get_vertflag())
#else
								if (!mesh_->is_boundary(vh1_))
#endif // USE_FEATURE
								{
									mesh_->set_point(vh1_, 0.5*(p0 + p1));
									project_on_reference_mesh_with_metric(vh1_, mesh_->point(vh1_));
								}
								mesh_->collapse(heh); mesh_->garbage_collection();
								ne = mesh_->n_edges(); ++collapse_count;
							}
#ifdef USE_FEATURE
							else if (mesh_->is_collapse_ok(mesh_->opposite_halfedge_handle(heh)) && !mesh_->data(vh1_).get_vertflag())
#else
							else if (mesh_->is_collapse_ok(mesh_->opposite_halfedge_handle(heh)) && !mesh_->is_boundary(vh1_))
#endif // USE_FEATURE
							{
								heh = mesh_->opposite_halfedge_handle(heh);
#ifdef USE_FEATURE
								if (!mesh_->data(vh0_).get_vertflag())
#else
								if (!mesh_->is_boundary(vh0_))
#endif // USE_FEATURE
								{
									mesh_->set_point(vh0_, 0.5*(p0 + p1));
									project_on_reference_mesh_with_metric(vh0_, mesh_->point(vh0_));
								}
								mesh_->collapse(heh); mesh_->garbage_collection();

								ne = mesh_->n_edges(); ++collapse_count;
							}
							else
							{
								++i;
							}
						}
						else
						{
							if (mesh_->is_collapse_ok(heh))
							{
								mesh_->set_point(vh1_, 0.5*(p0 + p1));
								project_on_reference_edge_with_metric(vh1_, mesh_->point(vh1_));

								mesh_->collapse(heh); mesh_->garbage_collection();
								ne = mesh_->n_edges(); ++collapse_count;
							}
							else if (mesh_->is_collapse_ok(mesh_->opposite_halfedge_handle(heh)))
							{
								heh = mesh_->opposite_halfedge_handle(heh);
								mesh_->set_point(vh0_, 0.5*(p0 + p1));
								project_on_reference_mesh_with_metric(vh0_, mesh_->point(vh0_));
								mesh_->collapse(heh); mesh_->garbage_collection();
								ne = mesh_->n_edges(); ++collapse_count;
							}
							else
							{
								++i;
							}
						}
					}
					else
					{
						++i;
					}
				}
				if (collapse_count == 0) break;
				else { printf("Collapse : %d : %d\n", iter_num, collapse_count); }
				++iter_num;
			}
			ne = mesh_->n_edges(); printf("Collapse Edges : %d\n", ne);
			if (iter_num == 0) no_split_collapse += 1;

			if (no_split_collapse == 2) break;
		}
		calc_tri_quality();
	}

	void AnisotropicMeshRemeshing::do_remeshing(double ref_edge_len /* = 1.0 */, double a /* = 1.5 */)
	{
		//auto ref_mesh_ = &(referenceModel->initial_trimesh);
		if (!ref_mesh_)
		{
			dprint("empty reference mesh!");
			return;
		}
		unsigned nv = mesh_->n_vertices();
		//project_on_reference();

		Eigen::Matrix3d H, H0, H1; Eigen::Vector3d P0, P1, P;
		Eigen::Matrix3d U; Eigen::Matrix3d V; Eigen::Vector3d sv; Eigen::Matrix3d diag_a; diag_a.setZero();
		Eigen::Matrix3d Q; std::vector<double> h(9);
		int iter_num = 0; int split_count = 0; int collapse_count = 0;

		std::cout << "Samplint tp " << std::endl;

		for (unsigned ii = 0; ii < 5; ++ii)
		{
			clock_t iter_start = clock();
			int no_split_collapse = 0;
			flip_based_energy();

			double step_length = 0.5 - 0.04*ii;
			reposition_LCOT(step_length);

			split_count = 0; nv = mesh_->n_vertices();
			unsigned ne = mesh_->n_edges();
			double cof = std::min(1.0, 0.1*ii + 0.7);
			for (unsigned i = 0; i < mesh_->n_edges(); ++i)
			{
				if (i % 10000 == 0) std::cout << i << "," << mesh_->n_edges() << std::endl;
				if (mesh_->n_edges() > ne*1.2) break;
				auto eh = mesh_->edge_handle(i);
				auto heh = mesh_->halfedge_handle(eh, 0);
				auto vh0_ = mesh_->from_vertex_handle(heh);
				auto vh1_ = mesh_->to_vertex_handle(heh);
				OpenMesh::Vec6d& h0 = mesh_->data(vh0_).get_Hessian();
				OpenMesh::Vec6d& h1 = mesh_->data(vh1_).get_Hessian();
				H0(0, 0) = h0[0]; H0(0, 1) = h0[1]; H0(0, 2) = h0[2];
				H0(1, 0) = h0[1]; H0(1, 1) = h0[3]; H0(1, 2) = h0[4];
				H0(2, 0) = h0[2]; H0(2, 1) = h0[4]; H0(2, 2) = h0[5];
				H1(0, 0) = h1[0]; H1(0, 1) = h1[1]; H1(0, 2) = h1[2];
				H1(1, 0) = h1[1]; H1(1, 1) = h1[3]; H1(1, 2) = h1[4];
				H1(2, 0) = h1[2]; H1(2, 1) = h1[4]; H1(2, 2) = h1[5];
				H = 0.5 * (H0 + H1);
				OpenMesh::Vec3d p0 = mesh_->point(vh0_);
				OpenMesh::Vec3d p1 = mesh_->point(vh1_);
				P0(0) = p0[0]; P0(1) = p0[1]; P0(2) = p0[2];
				P1(0) = p1[0]; P1(1) = p1[1]; P1(2) = p1[2];
				Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
				U = svd.matrixU(); V = svd.matrixV(); sv = svd.singularValues();
				diag_a(0, 0) = std::sqrt(sv(0)); diag_a(1, 1) = std::sqrt(sv(1)); diag_a(2, 2) = std::sqrt(sv(2));
				Q = U * diag_a*V.transpose();
				P0 = Q * P0; P1 = Q * P1;
				double edge_len = std::sqrt((P0(0) - P1(0))*(P0(0) - P1(0)) + (P0(1) - P1(1))*(P0(1) - P1(1)) + (P0(2) - P1(2))*(P0(2) - P1(2)));
				double max_len = 0;
				using std::max;
				if (mesh_->face_handle(heh).is_valid()) {
					max_len = max(mesh_->calc_edge_length(heh),
						max(mesh_->calc_edge_length(mesh_->next_halfedge_handle(heh)), mesh_->calc_edge_length(mesh_->prev_halfedge_handle(heh))));
				}
				OpenMesh::HalfedgeHandle oppoheh = mesh_->opposite_halfedge_handle(heh);
				if (mesh_->face_handle(oppoheh).is_valid()) {
					max_len = max(max_len, max(mesh_->calc_edge_length(oppoheh),
						max(mesh_->calc_edge_length(mesh_->next_halfedge_handle(oppoheh)), mesh_->calc_edge_length(mesh_->prev_halfedge_handle(oppoheh)))));
				}
				if (edge_len > ref_edge_len * a / cof && mesh_->calc_edge_length(eh) > max_len / (a))
					//if (edge_len > ref_edge_len * a / cof)
				{
#ifdef USE_FEATURE
					bool is_boundary_edge = mesh_->data(eh).get_edgeflag();
#else
					bool is_boundary_edge = mesh_->is_boundary(eh);
#endif // USE_FEATURE

					OpenMesh::Vec3d p = 0.5*(p0 + p1);
					find_nearst_point_on_reference_mesh(p, is_boundary_edge);
					if (!split_one_edge(eh, p)) { continue; }

					if (is_boundary_edge)//feature and not boundary
					{
						auto vh_ = mesh_->vertex_handle(nv);
						project_on_reference_edge_with_metric(vh_, mesh_->point(vh_));
					}
					else
					{
						auto vh_ = mesh_->vertex_handle(nv);
						project_on_reference_mesh_with_metric(vh_, mesh_->point(vh_));
					}

					nv = mesh_->n_vertices();
					++split_count;
				}
			}
			mesh_->garbage_collection();
			std::cout << split_count << ", split count\n";
			if (split_count == 0)
				no_split_collapse += 1;

			flip_based_energy();

			reposition_LCOT(0.5);

			collapse_count = 0;
			//continue;
			//if (ii == 0) break;
			for (auto te = mesh_->edges_sbegin(); te != mesh_->edges_end(); te++)
			{
				//if (i % 10000 == 0) std::cout << i << std::endl;
				if (te->idx() % 10000 == 0) std::cout << te->idx() << std::endl;
				auto eh = *te;
				auto heh = mesh_->halfedge_handle(eh, 0);
				auto vh0_ = mesh_->from_vertex_handle(heh);
				auto vh1_ = mesh_->to_vertex_handle(heh);
				OpenMesh::Vec6d& h0 = mesh_->data(vh0_).get_Hessian();
				OpenMesh::Vec6d& h1 = mesh_->data(vh1_).get_Hessian();
				H0(0, 0) = h0[0]; H0(0, 1) = h0[1]; H0(0, 2) = h0[2];
				H0(1, 0) = h0[1]; H0(1, 1) = h0[3]; H0(1, 2) = h0[4];
				H0(2, 0) = h0[2]; H0(2, 1) = h0[4]; H0(2, 2) = h0[5];
				H1(0, 0) = h1[0]; H1(0, 1) = h1[1]; H1(0, 2) = h1[2];
				H1(1, 0) = h1[1]; H1(1, 1) = h1[3]; H1(1, 2) = h1[4];
				H1(2, 0) = h1[2]; H1(2, 1) = h1[4]; H1(2, 2) = h1[5];
				H = 0.5 * (H0 + H1);
				OpenMesh::Vec3d p0 = mesh_->point(vh0_);
				OpenMesh::Vec3d p1 = mesh_->point(vh1_);
				P0(0) = p0[0]; P0(1) = p0[1]; P0(2) = p0[2];
				P1(0) = p1[0]; P1(1) = p1[1]; P1(2) = p1[2];
				Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
				U = svd.matrixU(); V = svd.matrixV(); sv = svd.singularValues();
				diag_a(0, 0) = std::sqrt(sv(0)); diag_a(1, 1) = std::sqrt(sv(1)); diag_a(2, 2) = std::sqrt(sv(2));
				Q = U * diag_a*V.transpose();
				P0 = Q * P0; P1 = Q * P1;
				double edge_len = std::sqrt((P0(0) - P1(0))*(P0(0) - P1(0)) + (P0(1) - P1(1))*(P0(1) - P1(1)) + (P0(2) - P1(2))*(P0(2) - P1(2)));
				//double temp_e = (p0 - p1).norm();
				double min_len = mesh_->calc_edge_length(heh);
				using std::min;
				if (mesh_->face_handle(heh).is_valid()) {
					min_len = min(min_len, min(mesh_->calc_edge_length(mesh_->next_halfedge_handle(heh)), mesh_->calc_edge_length(mesh_->prev_halfedge_handle(heh))));
				}
				OpenMesh::HalfedgeHandle oppoheh = mesh_->opposite_halfedge_handle(heh);
				if (mesh_->face_handle(oppoheh).is_valid()) {
					min_len = min(min_len, min(mesh_->calc_edge_length(oppoheh),
						min(mesh_->calc_edge_length(mesh_->next_halfedge_handle(oppoheh)), mesh_->calc_edge_length(mesh_->prev_halfedge_handle(oppoheh)))));
				}

				if (edge_len < ref_edge_len * cof / a && mesh_->calc_edge_length(eh) < min_len*a)
					//if (edge_len < ref_edge_len * cof / a)
				{
#ifdef USE_FEATURE
					if (!mesh_->data(eh).get_edgeflag())
#else
					if (!mesh_->is_boundary(eh))
#endif // USE_FEATURE
					{
#ifdef USE_FEATURE
						if (mesh_->is_collapse_ok(heh) && !mesh_->data(vh0_).get_vertflag())
#else
						if (mesh_->is_collapse_ok(heh) && !mesh_->is_boundary(vh0_)) //from vh0 to vh1
#endif // USE_FEATURE
						{
#ifdef USE_FEATURE
							if (!mesh_->data(vh1_).get_vertflag())
#else
							if (!mesh_->is_boundary(vh1_))
#endif // USE_FEATURE
							{
								OpenMesh::Vec3d pri = mesh_->point(vh1_);
								mesh_->set_point(vh1_, 0.5*(p0 + p1));
								project_on_reference_mesh_with_metric(vh1_, mesh_->point(vh1_));
								OpenMesh::HalfedgeHandle dh = mesh_->opposite_halfedge_handle(mesh_->next_halfedge_handle(mesh_->opposite_halfedge_handle(heh)));
								if (!mesh_->is_boundary(dh) &&
									(mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(dh)) + mesh_->calc_sector_angle(mesh_->opposite_halfedge_handle(dh))) > PI)
								{
									mesh_->set_point(vh1_, pri);
									continue;
								}
								OpenMesh::HalfedgeHandle th = mesh_->opposite_halfedge_handle(mesh_->prev_halfedge_handle(heh));
								if (!mesh_->is_boundary(th) &&
									(mesh_->calc_sector_angle(th) + mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(mesh_->opposite_halfedge_handle(th)))) > PI) {
									mesh_->set_point(vh1_, pri);
									continue;
								}
							}
							mesh_->collapse(heh);
							++collapse_count;
						}
#ifdef USE_FEATURE
						else if (mesh_->is_collapse_ok(mesh_->opposite_halfedge_handle(heh)) && !mesh_->data(vh1_).get_vertflag())
#else
						else if (mesh_->is_collapse_ok(mesh_->opposite_halfedge_handle(heh)) && !mesh_->is_boundary(vh1_))
#endif // USE_FEATURE
						{
							heh = mesh_->opposite_halfedge_handle(heh);
#ifdef USE_FEATURE
							if (!mesh_->data(vh0_).get_vertflag())
#else
							if (!mesh_->is_boundary(vh0_))
#endif // USE_FEATURE
							{
								OpenMesh::Vec3d pri = mesh_->point(vh0_);
								mesh_->set_point(vh0_, 0.5*(p0 + p1));
								project_on_reference_mesh_with_metric(vh0_, mesh_->point(vh0_));
								OpenMesh::HalfedgeHandle dh = mesh_->opposite_halfedge_handle(mesh_->next_halfedge_handle(mesh_->opposite_halfedge_handle(heh)));
								if (!mesh_->is_boundary(dh) &&
									(mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(dh)) + mesh_->calc_sector_angle(mesh_->opposite_halfedge_handle(dh))) > PI)
								{
									mesh_->set_point(vh0_, pri);
									continue;
								}
								OpenMesh::HalfedgeHandle th = mesh_->opposite_halfedge_handle(mesh_->prev_halfedge_handle(heh));
								if (!mesh_->is_boundary(th) &&
									(mesh_->calc_sector_angle(th) + mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(mesh_->opposite_halfedge_handle(th)))) > PI) 
								{
									mesh_->set_point(vh0_, pri);
									continue;
								}
							}
							mesh_->collapse(heh);
							++collapse_count;
						}
					}
					else
					{
						if (mesh_->is_collapse_ok(heh))
						{
							OpenMesh::Vec3d pri = mesh_->point(vh1_);
							mesh_->set_point(vh1_, 0.5*(p0 + p1));
							project_on_reference_edge_with_metric(vh1_, mesh_->point(vh1_));
							OpenMesh::HalfedgeHandle dh = mesh_->opposite_halfedge_handle(mesh_->next_halfedge_handle(mesh_->opposite_halfedge_handle(heh)));
							if (!mesh_->is_boundary(dh) &&
								(mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(dh)) + mesh_->calc_sector_angle(mesh_->opposite_halfedge_handle(dh))) > PI)
							{
								mesh_->set_point(vh1_, pri);
								continue;
							}
							OpenMesh::HalfedgeHandle th = mesh_->opposite_halfedge_handle(mesh_->prev_halfedge_handle(heh));
							if (!mesh_->is_boundary(th) &&
								(mesh_->calc_sector_angle(th) + mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(mesh_->opposite_halfedge_handle(th)))) > PI) {
								mesh_->set_point(vh1_, pri);
								continue;
							}

							mesh_->collapse(heh);
							++collapse_count;
						}
						else if (mesh_->is_collapse_ok(mesh_->opposite_halfedge_handle(heh)))
						{
							heh = mesh_->opposite_halfedge_handle(heh);
							OpenMesh::Vec3d pri = mesh_->point(vh0_);
							mesh_->set_point(vh0_, 0.5*(p0 + p1));
							project_on_reference_mesh_with_metric(vh0_, mesh_->point(vh0_));
							OpenMesh::HalfedgeHandle dh = mesh_->opposite_halfedge_handle(mesh_->next_halfedge_handle(mesh_->opposite_halfedge_handle(heh)));
							if (!mesh_->is_boundary(dh) &&
								(mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(dh)) + mesh_->calc_sector_angle(mesh_->opposite_halfedge_handle(dh))) > PI)
							{
								mesh_->set_point(vh0_, pri);
								continue;
							}
							OpenMesh::HalfedgeHandle th = mesh_->opposite_halfedge_handle(mesh_->prev_halfedge_handle(heh));
							if (!mesh_->is_boundary(th) &&
								(mesh_->calc_sector_angle(th) + mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(mesh_->opposite_halfedge_handle(th)))) > PI) {
								mesh_->set_point(vh0_, pri);
								continue;
							}
							mesh_->collapse(heh);
							++collapse_count;
						}
					}
				}
			}

			mesh_->garbage_collection();
			std::cout << collapse_count << ", collapse cout\n";
			if (collapse_count == 0)
				no_split_collapse += 1;
			//if (ii == 0) break;
			if (no_split_collapse == 2)
				break;

			std::cout << ii << " iter is OK " << std::endl;
		}

#ifdef USE_FEWERITERATION
		LCOT_Optimize(10, 1.0);
		LCOT_Optimize(5, 0.25);
#else
		LCOT_Optimize(50, 1.0);
		LCOT_Optimize(25, 0.25);
#endif // USE_FEWERITERATION

		calc_tri_quality();
	}

	void AnisotropicMeshRemeshing::calc_tri_quality()
	{
		if (mesh_->n_vertices() == 0) return;

		std::vector<Eigen::Matrix3d> v_h(mesh_->n_vertices()); Eigen::Matrix3d H;
		for (auto v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			int vertex_id = v_it.handle().idx();
			OpenMesh::Vec6d& h = mesh_->data(v_it).get_Hessian();
			H(0, 0) = h[0]; H(1, 0) = h[1]; H(2, 0) = h[2];
			H(0, 1) = h[1]; H(1, 1) = h[3]; H(2, 1) = h[4];
			H(0, 2) = h[2]; H(1, 2) = h[4]; H(2, 2) = h[5];
			v_h[vertex_id] = H;
		}

		unsigned nf = mesh_->n_faces();
		std::vector<double > quality_error(nf, 0.0); std::vector<double > theta;

		Eigen::Matrix3d U; Eigen::Matrix3d V; Eigen::Vector3d a; Eigen::Matrix3d diag_a; diag_a.setZero();
		Eigen::Matrix3d Q; Eigen::Vector3d d;
		std::vector<OpenMesh::Vec3d> p(3); Eigen::Vector3d p_;
		std::vector<double> angle(3); double area; double l_e, s_e; double h_p;
		double min_quality_error = 1e30; double ave_quality_error = 0.0; double max_quality_error = 0.0;
		int min_quality_error_face, max_quality_error_face;
		double max_angle = 0.0; double min_angle = 360; int below_30_count = 0; double ave_angle = 0.0; int above_90_count = 0;
		double min_angle_vertex = -1; double min_non_boundary_angle = 360; below_30_tri.clear(); below_30_tri_angle.clear();
		double ave_transfrom_area = 0.0; double max_transfrom_area_ratio = 0.0; double min_transfrom_area_ratio = 1.0e30; int max_transfrom_area_id, min_transfrom_area_id;
		std::vector<double> transform_area(nf); std::vector<double> face_smalleset_angle(nf);
		std::vector<double> radius_edge_ratio(nf, 0.0); OpenMesh::Vec3d c;
		double smallest_radius_edge_ratio = 1e30; double largest_radius_edge_ratio = 0.0;
		for (auto f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
		{
			int face_id = f_it->idx();
			H.setZero();  int v_count = 0;
			for (auto fv_it = mesh_->fv_iter(f_it); fv_it; ++fv_it)
			{
				int fv_id = fv_it.handle().idx();
				H += v_h[fv_id] / 3.0;
				p[v_count] = mesh_->point(fv_it);
				++v_count;
			}
			Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
			U = svd.matrixU(); V = svd.matrixV(); a = svd.singularValues();
			diag_a(0, 0) = std::sqrt(a(0)); diag_a(1, 1) = std::sqrt(a(1)); diag_a(2, 2) = std::sqrt(a(2));
			Q = U * diag_a*V.transpose();

			for (unsigned i = 0; i < 3; ++i)
			{
				d(0) = p[i][0]; d(1) = p[i][1]; d(2) = p[i][2];
				p_ = Q * d;
				p[i][0] = p_(0); p[i][1] = p_(1); p[i][2] = p_(2);
			}

			area = OpenMesh::cross(p[1] - p[0], p[2] - p[0]).norm() * 0.5;

			ave_transfrom_area += area;
			if (max_transfrom_area_ratio < area) { max_transfrom_area_ratio = area; max_transfrom_area_id = face_id; }
			if (min_transfrom_area_ratio > area) { min_transfrom_area_ratio = area; min_transfrom_area_id = face_id; }
			transform_area[face_id] = area;

			l_e = 0.0; s_e = 1e30; h_p = 0.0;
			double sin_val = (OpenMesh::cross((p[2] - p[1]), (p[0] - p[1]))).norm();
			for (unsigned ii = 0; ii < 3; ++ii)
			{
				unsigned jj = (ii == 2 ? 0 : ii + 1);
				unsigned kk = (ii == 0 ? 2 : ii - 1);
				double temp_l = (p[jj] - p[ii]).norm();
				h_p += temp_l;
				if (temp_l > l_e) l_e = temp_l;
				if (temp_l < s_e) s_e = temp_l;
				double cos_val = OpenMesh::dot(p[jj] - p[ii], p[kk] - p[ii]);
				angle[ii] = std::atan2(sin_val, cos_val) * 180 / M_PI;
			}

			h_p *= 0.5;
			quality_error[face_id] = 2.0*sqrt(3.0) * area / (h_p*l_e);
			if (quality_error[face_id] < min_quality_error) { min_quality_error = quality_error[face_id]; min_quality_error_face = face_id; }
			if (quality_error[face_id] > max_quality_error) { max_quality_error = quality_error[face_id]; max_quality_error_face = face_id; }
			ave_quality_error += quality_error[face_id];

			double local_min_angle = 360.0;
			for (unsigned i = 0; i < angle.size(); ++i)
			{
				theta.push_back(angle[i]);

				if (max_angle < angle[i]) max_angle = angle[i];
				if (local_min_angle > angle[i]) { local_min_angle = angle[i]; }
				if (angle[i] > 90)
				{
					++above_90_count;
				}
			}

			if (min_angle > local_min_angle)
			{
				min_angle = local_min_angle;
			}
			if (local_min_angle < 30)
			{
				++below_30_count;
				below_30_tri.push_back(face_id);
				below_30_tri_angle.push_back(local_min_angle);
			}
			face_smalleset_angle[face_id] = local_min_angle;
			ave_angle += local_min_angle;

			//radius edge ratio
			OpenMesh::Vec3d n = OpenMesh::cross(p[2] - p[1], p[0] - p[1]);
			Eigen::Matrix3d A; Eigen::Vector3d b;
			A(0, 0) = n[0]; A(0, 1) = n[1]; A(0, 2) = n[2];
			b(0) = OpenMesh::dot(n, p[0]);
			b(1) = p[1].sqrnorm() - p[0].sqrnorm();
			b(2) = p[2].sqrnorm() - p[0].sqrnorm();

			A(1, 0) = 2.0 * (p[1][0] - p[0][0]); A(1, 1) = 2.0 * (p[1][1] - p[0][1]); A(1, 2) = 2.0 * (p[1][2] - p[0][2]);
			A(2, 0) = 2.0 * (p[2][0] - p[0][0]); A(2, 1) = 2.0 * (p[2][1] - p[0][1]); A(2, 2) = 2.0 * (p[2][2] - p[0][2]);

			Eigen::Vector3d x = A.fullPivHouseholderQr().solve(b);
			c[0] = x(0); c[1] = x(1); c[2] = x(2);

			double radius = (c - p[0]).norm();
			radius_edge_ratio[face_id] = radius / s_e;
			if (radius_edge_ratio[face_id] < smallest_radius_edge_ratio)
			{
				smallest_radius_edge_ratio = radius_edge_ratio[face_id];
			}
			if (radius_edge_ratio[face_id] > largest_radius_edge_ratio)
			{
				largest_radius_edge_ratio = radius_edge_ratio[face_id];
			}
		}

		//edge length
		std::vector<double> edge_length(mesh_->n_edges());
		double smallest_edge_length = 1e30; double largest_edge_length = 0.0; double avg_edge_length = 0.0;
		int smallest_edge_id, largest_edge_id = -1;
		for (auto e_it = mesh_->edges_begin(); e_it != mesh_->edges_end(); ++e_it)
		{
			H.setZero(); int v_count = 0;
			auto heh = mesh_->halfedge_handle(e_it.handle(), 0);
			auto vh = mesh_->from_vertex_handle(heh);
			OpenMesh::Vec3d p_v = mesh_->point(vh);
			H += v_h[vh.idx()] * 0.5;
			vh = mesh_->to_vertex_handle(heh);
			p_v -= mesh_->point(vh);
			H += v_h[vh.idx()] * 0.5;

			Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
			U = svd.matrixU(); V = svd.matrixV(); a = svd.singularValues();
			diag_a(0, 0) = std::sqrt(a(0)); diag_a(1, 1) = std::sqrt(a(1)); diag_a(2, 2) = std::sqrt(a(2));
			Q = U * diag_a*V.transpose();
			d(0) = p_v[0]; d(1) = p_v[1]; d(2) = p_v[2];
			p_ = Q * d; p_v[0] = p_(0); p_v[1] = p_(1); p_v[2] = p_(2);
			edge_length[e_it.handle().idx()] = p_v.norm();
			avg_edge_length += p_v.norm();
			if (edge_length[e_it.handle().idx()] < smallest_edge_length)
			{
				smallest_edge_length = edge_length[e_it.handle().idx()];
				smallest_edge_id = e_it.handle().idx();
			}
			if (edge_length[e_it.handle().idx()] > largest_edge_length)
			{
				largest_edge_length = edge_length[e_it.handle().idx()];
				largest_edge_id = e_it.handle().idx();
			}
		}

		ave_quality_error /= mesh_->n_faces();
		printf("\n[Quality Error] Min : %f, Ave : %f\n", min_quality_error, ave_quality_error);

		double min_theta = min_angle; double ave_theta = ave_angle / (double)(mesh_->n_faces());
		double theta_below_30 = (double)(below_30_count) / (double)(mesh_->n_faces());
		double theta_above_90 = (double)(above_90_count) / (double)(mesh_->n_faces());
		printf("[Angle] Min : %f, Max : %f, Ave : %f\n<30 : %f%%, %d; >90 : %f%%, %d\n", min_theta, max_angle, ave_theta, theta_below_30*100.0, below_30_count, theta_above_90, above_90_count);

		ave_transfrom_area /= mesh_->n_faces(); double below_05_transform_area = 0; double above_2_transform_area = 0;
		for (unsigned i = 0; i < transform_area.size(); ++i)
		{
			transform_area[i] /= ave_transfrom_area;
			if (transform_area[i] < 0.5) below_05_transform_area += 1.0;
			else if (transform_area[i] > 2.0) above_2_transform_area += 1.0;
		}
		min_transfrom_area_ratio /= ave_transfrom_area; max_transfrom_area_ratio /= ave_transfrom_area;
		printf("[Area Quality Error] Min : %f(%d), Max : %f(%d)\n", min_transfrom_area_ratio, min_transfrom_area_id, max_transfrom_area_ratio, max_transfrom_area_id);
		printf("Below 0.5 : %f%%, %d, Above 2.0 : %f%%, %d\n", below_05_transform_area*100.0 / mesh_->n_faces(), (int)below_05_transform_area,
			above_2_transform_area*100.0 / mesh_->n_faces(), (int)above_2_transform_area);

		avg_edge_length /= mesh_->n_edges();
		printf("Edge Length : Min : %f(%d), Max : %f(%d); Avg : %f\n", smallest_edge_length, smallest_edge_id, largest_edge_length, largest_edge_id, avg_edge_length);
		printf("Radius Edge Ratio : Min : %f, Max : %f\n", smallest_radius_edge_ratio, largest_radius_edge_ratio);
		printf("Mesh vertices : %d\n\n", mesh_->n_vertices());
	}

	double AnisotropicMeshRemeshing::calc_flip_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec3d& p3, const OpenMesh::Vec6d& M, bool use_area)
	{
		double m0 = M[0]; double m1 = M[1]; double m2 = M[2]; double m3 = M[3]; double m4 = M[4]; double m5 = M[5];
		double x1 = p1[0]; double y1 = p1[1]; double z1 = p1[2];
		double x2 = p2[0]; double y2 = p2[1]; double z2 = p2[2];
		double x3 = p3[0]; double y3 = p3[1]; double z3 = p3[2];

		double e = (m0*x1*x1) / 24 - (m0*x1*x2) / 24 - (m0*x1*x3) / 24 + (m1*x1*y1) / 12 - (m1*x1*y2) / 24 - (m1*x1*y3) / 24 + (m2*x1*z1) / 12
			- (m2*x1*z2) / 24 - (m2*x1*z3) / 24 + (m0*x2*x2) / 24 - (m0*x2*x3) / 24 - (m1*x2*y1) / 24 + (m1*x2*y2) / 12 - (m1*x2*y3) / 24
			- (m2*x2*z1) / 24 + (m2*x2*z2) / 12 - (m2*x2*z3) / 24 + (m0*x3*x3) / 24 - (m1*x3*y1) / 24 - (m1*x3*y2) / 24 + (m1*x3*y3) / 12
			- (m2*x3*z1) / 24 - (m2*x3*z2) / 24 + (m2*x3*z3) / 12 + (m3*y1*y1) / 24 - (m3*y1*y2) / 24 - (m3*y1*y3) / 24 + (m4*y1*z1) / 12
			- (m4*y1*z2) / 24 - (m4*y1*z3) / 24 + (m3*y2*y2) / 24 - (m3*y2*y3) / 24 - (m4*y2*z1) / 24 + (m4*y2*z2) / 12 - (m4*y2*z3) / 24
			+ (m3*y3*y3) / 24 - (m4*y3*z1) / 24 - (m4*y3*z2) / 24 + (m4*y3*z3) / 12 + (m5*z1*z1) / 24 - (m5*z1*z2) / 24 - (m5*z1*z3) / 24
			+ (m5*z2*z2) / 24 - (m5*z2*z3) / 24 + (m5*z3*z3) / 24;
		if (use_area)
		{
			double a1 = (y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2);
			double a2 = (x1 - x3)*(z1 - z2) - (x1 - x2)*(z1 - z3);
			double a3 = (x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2);
			return e * std::sqrt(a1*a1 + a2 * a2 + a3 * a3);//Area
			//return std::exp(e * std::sqrt(a1*a1 + a2*a2 + a3*a3) );//Area
		}
		else
		{
			return e;
		}
	}

	bool AnisotropicMeshRemeshing::flip_based_energy()
	{
		mesh_->update_face_normals();
		int nv = mesh_->n_vertices();
		Eigen::Matrix3d H; Eigen::Matrix3d D; D.setZero();
		Eigen::Matrix3d vH; TriMesh::VertexHandle vh;
		std::vector<OpenMesh::Vec6d> vM(nv);

		for (unsigned i = 0; i < nv; ++i)
		{
			auto vh = mesh_->vertex_handle(i);
			vM[i] = mesh_->data(vh).get_Hessian();
			//vM[i] = OpenMesh::Vec6d(1, 0, 0, 1, 0, 1);
		}

		OpenMesh::Vec3d p0; OpenMesh::Vec3d p1;
		OpenMesh::Vec3d p2; OpenMesh::Vec3d p3;
		double f0; double f1; double f2; double f3;
		double energy0; double energy1; double area;
		int iter_count = 0; OpenMesh::Vec3d n;
		std::vector<OpenMesh::Vec3d> tri_p(3); std::vector<double> tri_f_val(3);
		OpenMesh::Vec3d v1; OpenMesh::Vec3d v2;
		double angle_th = 0.01745240643728351281941897851632;
		while (iter_count < 10)
		{
			int flip_count = 0;
			for (auto e_it = mesh_->edges_begin(); e_it != mesh_->edges_end(); ++e_it)
			{
#ifdef USE_FEATURE
				//if (mesh_->data(e_it).get_edgeflag() || !is_flip_ok_openmesh(e_it.handle(), *mesh_))
				if (mesh_->data(e_it).get_edgeflag() || !mesh_->is_flip_ok(e_it.handle()))
#else
				if (mesh_->is_boundary(e_it) || !is_flip_ok_openmesh(e_it.handle(), *mesh_))
#endif // USE_FEATURE
					continue;

				OpenMesh::HalfedgeHandle a0 = mesh_->halfedge_handle(e_it, 0);
				OpenMesh::FaceHandle fh0 = mesh_->face_handle(a0); n = mesh_->normal(fh0);
				OpenMesh::HalfedgeHandle b0 = mesh_->halfedge_handle(e_it, 1);
				OpenMesh::FaceHandle fh1 = mesh_->face_handle(b0); n += mesh_->normal(fh1);
				n.normalize();

#ifdef USE_PROMOTION
				//一定程度上防止翻折
				if (mesh_->calc_sector_angle(a0) + mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(b0)) > PI
					|| mesh_->calc_sector_angle(b0) + mesh_->calc_sector_angle(mesh_->prev_halfedge_handle(a0)) > PI)
					continue;
#endif // USE_PROMOTION


				vh = mesh_->to_vertex_handle(a0); int v0_id = vh.idx();
				p0 = mesh_->point(vh);
				vh = mesh_->to_vertex_handle(b0); int v1_id = vh.idx();
				p1 = mesh_->point(vh);
				vh = mesh_->to_vertex_handle(mesh_->next_halfedge_handle(a0)); int v2_id = vh.idx();
				p2 = mesh_->point(vh);
				vh = mesh_->to_vertex_handle(mesh_->next_halfedge_handle(b0)); int v3_id = vh.idx();
				p3 = mesh_->point(vh);

				//#ifdef USE_PROMOTION
				//			//防止极端角出现
				//			auto ang0 = [&](OpenMesh::Vec3d p0, OpenMesh::Vec3d p1, OpenMesh::Vec3d p2) {return acos((p0 - p1).dot(p2 - p1)); };
				//			auto ang1 = [&](OpenMesh::Vec3d p0, OpenMesh::Vec3d p1, OpenMesh::Vec3d p2) {return ang0(p0, p1, p2) > least_angle && ang0(p0, p1, p2) < largest_angle
				//				&& ang0(p1,p2,p0) > least_angle && ang0(p1,p2,p0) < largest_angle && ang0(p2,p0,p1) > least_angle && ang0(p2,p0,p1) < largest_angle; };
				//			if (!ang1(p0, p2, p3) || !ang1(p1, p2, p3)) continue;
				//#endif // USE_PROMOTION


				energy0 = calc_flip_energy(p0, p2, p1, (vM[v0_id] + vM[v2_id] + vM[v1_id]) / 3.0, true);
				energy0 += calc_flip_energy(p0, p1, p3, (vM[v0_id] + vM[v1_id] + vM[v3_id]) / 3.0, true);

				energy1 = calc_flip_energy(p0, p2, p3, (vM[v0_id] + vM[v2_id] + vM[v3_id]) / 3.0, true);
				energy1 += calc_flip_energy(p2, p1, p3, (vM[v2_id] + vM[v1_id] + vM[v3_id]) / 3.0, true);

				double z_flag10 = OpenMesh::cross(p2 - p0, p3 - p0).norm(); //area
				double z_flag12 = z_flag10 / ((p0 - p2).norm() * (p3 - p2).norm());
				double z_flag13 = z_flag10 / ((p0 - p3).norm() * (p2 - p3).norm());
				z_flag10 /= ((p2 - p0).norm() * (p3 - p0).norm());
				double z_flag22 = OpenMesh::cross(p1 - p2, p3 - p2).norm();
				double z_flag21 = z_flag22 / ((p2 - p1).norm() * (p3 - p1).norm());
				double z_flag23 = z_flag22 / ((p2 - p3).norm() * (p1 - p3).norm());
				z_flag22 /= ((p1 - p2).norm() * (p3 - p2).norm());

				if (energy0 > energy1
					&& z_flag10 > angle_th && z_flag12 > angle_th && z_flag13 > angle_th
					&& z_flag22 > angle_th && z_flag21 > angle_th && z_flag23 > angle_th)
				{
					//flip_openmesh(e_it.handle(), *mesh_); ++flip_count;
					mesh_->flip(e_it.handle()); ++flip_count;
				}
			}
			//printf("%d Flip Count : %d\n",iter_count, flip_count);
			if (flip_count == 0) { break; };
			++iter_count;
		}
		return true;
	}

	bool AnisotropicMeshRemeshing::reposition_LCOT(double step_length)
	{
		unsigned nv = mesh_->n_vertices(); if (nv == 0) return false;
		std::vector<OpenMesh::Vec6d> vH(nv);
		for (unsigned i = 0; i < nv; ++i)
		{
			auto vh = mesh_->vertex_handle(i);
			vH[i] = mesh_->data(vh).get_Hessian();
		}
		unsigned nf = mesh_->n_faces();
		std::vector<OpenMesh::Vec6d> fH(nf); std::vector<double> face_area(nf); std::vector<OpenMesh::Vec3d> face_normal(nf);
		for (unsigned i = 0; i < nf; ++i)
		{
			auto fv_it = mesh_->fv_iter(mesh_->face_handle(i));
			OpenMesh::Vec3d v0 = mesh_->point(fv_it);
			fH[i] = vH[fv_it.handle().idx()]; ++fv_it;
			OpenMesh::Vec3d v1 = mesh_->point(fv_it);
			fH[i] += vH[fv_it.handle().idx()]; ++fv_it;
			OpenMesh::Vec3d v2 = mesh_->point(fv_it);
			fH[i] += vH[fv_it.handle().idx()]; fH[i] /= 3.0;

			face_normal[i] = OpenMesh::cross(v1 - v0, v2 - v0); // n * area * 2
			face_area[i] = face_normal[i].norm() * 0.5;
		}

		std::vector<local_frame> vertex_f(nv);
		for (auto v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			int vertex_id = v_it.handle().idx(); OpenMesh::Vec3d n(0, 0, 0);
			for (auto vf_it = mesh_->vf_iter(v_it); vf_it; ++vf_it)
			{
				int face_id = vf_it.handle().idx();
				n += face_normal[face_id];
			}
			vertex_f[vertex_id].n = n.normalize();
			vertex_f[vertex_id].find_e_x_y();
		}

		//adjust the position
		std::vector<OpenMesh::Vec3d> new_pos(nv); OpenMesh::Vec3d p;  OpenMesh::Vec3d tp;
		Eigen::Matrix3d Trans, H_3D; Eigen::Matrix3d vertex_h;
		double v_h1; double v_h2; double v_h3; double v_h4; double v_h5; double v_h6;
		double lamda = 0.05;
		for (auto v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			int vertex_id = v_it->idx(); p = mesh_->point(v_it);
#ifdef USE_FEATURE
			if (mesh_->data(v_it).get_vertflag())
#else
			if (mesh_->is_boundary(v_it))
#endif // USE_FEATURE
			{
				new_pos[vertex_id] = p;
				continue;
			}

			double x1 = p[0]; double y1 = p[1]; double z1 = p[2];

			double min_radius = 1.0e30;
			double gx = 0.0; double gy = 0.0; double gz = 0.0; OpenMesh::Vec3d g;
			v_h1 = v_h2 = v_h3 = v_h4 = v_h5 = v_h6 = 0.0;
			for (auto voh_it = mesh_->voh_iter(v_it); voh_it; ++voh_it)
			{
				auto vh = mesh_->to_vertex_handle(voh_it); tp = mesh_->point(vh);
				double x2 = tp[0]; double y2 = tp[1]; double z2 = tp[2];
				double mid_edge_len = (tp - p).norm() * 0.5;
				if (mid_edge_len < min_radius) min_radius = mid_edge_len;

				auto fh = mesh_->face_handle(voh_it);
				if (fh != TriMesh::InvalidFaceHandle)
				{
					int face_id = fh.idx();
					auto heh = mesh_->next_halfedge_handle(voh_it);
					vh = mesh_->to_vertex_handle(heh); tp = mesh_->point(vh);
					double x3 = tp[0]; double y3 = tp[1]; double z3 = tp[2];
					double h1 = fH[face_id][0]; double h2 = fH[face_id][1]; double h3 = fH[face_id][2];
					double h4 = fH[face_id][3]; double h5 = fH[face_id][4]; double h6 = fH[face_id][5];
					double f = (h1*x1*x1) / 24 - (h1*x1*x2) / 24 - (h1*x1*x3) / 24 + (h2*x1*y1) / 12 - (h2*x1*y2) / 24 - (h2*x1*y3) / 24 + (h3*x1*z1) / 12
						- (h3*x1*z2) / 24 - (h3*x1*z3) / 24 + (h1*x2*x2) / 24 - (h1*x2*x3) / 24 - (h2*x2*y1) / 24 + (h2*x2*y2) / 12 - (h2*x2*y3) / 24
						- (h3*x2*z1) / 24 + (h3*x2*z2) / 12 - (h3*x2*z3) / 24 + (h1*x3*x3) / 24 - (h2*x3*y1) / 24 - (h2*x3*y2) / 24 + (h2*x3*y3) / 12
						- (h3*x3*z1) / 24 - (h3*x3*z2) / 24 + (h3*x3*z3) / 12 + (h4*y1*y1) / 24 - (h4*y1*y2) / 24 - (h4*y1*y3) / 24 + (h5*y1*z1) / 12
						- (h5*y1*z2) / 24 - (h5*y1*z3) / 24 + (h4*y2*y2) / 24 - (h4*y2*y3) / 24 - (h5*y2*z1) / 24 + (h5*y2*z2) / 12 - (h5*y2*z3) / 24
						+ (h4*y3*y3) / 24 - (h5*y3*z1) / 24 - (h5*y3*z2) / 24 + (h5*y3*z3) / 12 + (h6*z1*z1) / 24 - (h6*z1*z2) / 24 - (h6*z1*z3) / 24
						+ (h6*z2*z2) / 24 - (h6*z2*z3) / 24 + (h6*z3*z3) / 24;

					double a = face_area[face_id];

					double dfx = (h1*x1) / 12 - (h1*x2) / 24 - (h1*x3) / 24 + (h2*y1) / 12 - (h2*y2) / 24 - (h2*y3) / 24 + (h3*z1) / 12 - (h3*z2) / 24 - (h3*z3) / 24;
					double dfy = (h2*x1) / 12 - (h2*x2) / 24 - (h2*x3) / 24 + (h4*y1) / 12 - (h4*y2) / 24 - (h4*y3) / 24 + (h5*z1) / 12 - (h5*z2) / 24 - (h5*z3) / 24;
					double dfz = (h3*x1) / 12 - (h3*x2) / 24 - (h3*x3) / 24 + (h5*y1) / 12 - (h5*y2) / 24 - (h5*y3) / 24 + (h6*z1) / 12 - (h6*z2) / 24 - (h6*z3) / 24;

					double dax = (2 * (y2 - y3)*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) + 2 * (z2 - z3)*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))) / (4 * std::sqrt(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) + ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2)) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))));
					double day = -(2 * (x2 - x3)*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) - 2 * (z2 - z3)*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))) / (4 * std::sqrt(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) + ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2)) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))));
					double daz = -(2 * (x2 - x3)*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2)) + 2 * (y2 - y3)*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))) / (4 * std::sqrt(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) + ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2)) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))));

					if (a < 1e-8) { dax = 0.0; day = 0.0; daz = 0.0; }

					g[0] = (f * dax + dfx * a); g[1] = (f * day + dfy * a); g[2] = (f * daz + dfz * a);
					gx -= OpenMesh::dot(g, vertex_f[vertex_id].e_x); gy -= OpenMesh::dot(g, vertex_f[vertex_id].e_y);

					double h = a * 2.0 / std::sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2) + (z3 - z2)*(z3 - z2));
					if (h < min_radius) min_radius = h;

					v_h1 += a * h1 / 12.0; v_h2 += a * h2 / 12.0; v_h3 += a * h3 / 12.0; v_h4 += a * h4 / 12.0; v_h5 += a * h5 / 12.0; v_h6 += a * h6 / 12.0;
				}
			}
			vertex_h(0, 0) = v_h1; vertex_h(0, 1) = v_h2; vertex_h(0, 2) = v_h3;
			vertex_h(1, 0) = v_h2; vertex_h(1, 1) = v_h4; vertex_h(1, 2) = v_h5;
			vertex_h(2, 0) = v_h3; vertex_h(2, 1) = v_h5; vertex_h(2, 2) = v_h6;

			Trans(0, 0) = vertex_f[vertex_id].e_x[0]; Trans(0, 1) = vertex_f[vertex_id].e_x[1]; Trans(0, 2) = vertex_f[vertex_id].e_x[2];
			Trans(1, 0) = vertex_f[vertex_id].e_y[0]; Trans(1, 1) = vertex_f[vertex_id].e_y[1]; Trans(1, 2) = vertex_f[vertex_id].e_y[2];
			Trans(2, 0) = vertex_f[vertex_id].n[0];   Trans(2, 1) = vertex_f[vertex_id].n[1];   Trans(2, 2) = vertex_f[vertex_id].n[2];
			H_3D = Trans * vertex_h * Trans.transpose(); //H_3D.setIdentity();
			double inv_det = 1.0 / (H_3D(0, 0) * H_3D(1, 1) - H_3D(0, 1) * H_3D(1, 0));
			double dx0 = inv_det * (gx * H_3D(1, 1) - H_3D(0, 1)* gy);
			double dx1 = inv_det * (gy * H_3D(0, 0) - H_3D(1, 0)* gx);
			double move_dis = std::sqrt(dx0*dx0 + dx1 * dx1);
			if (move_dis > min_radius)
			{
				dx0 = dx0 * min_radius / move_dis;
				dx1 = dx1 * min_radius / move_dis;
			}

			new_pos[vertex_id] = p + (dx0 * vertex_f[vertex_id].e_x + dx1 * vertex_f[vertex_id].e_y) * step_length;
		}

		project_on_reference_new_p(new_pos);
		return false;
	}

	bool AnisotropicMeshRemeshing::reposition_exp_LCOT(double step_length)
	{
		unsigned nv = mesh_->n_vertices(); if (nv == 0) return false;
		std::vector<OpenMesh::Vec6d> vH(nv);
		for (unsigned i = 0; i < nv; ++i)
		{
			auto vh = mesh_->vertex_handle(i);
			vH[i] = mesh_->data(vh).get_Hessian();
		}
		unsigned nf = mesh_->n_faces();
		std::vector<OpenMesh::Vec6d> fH(nf); std::vector<double> face_area(nf); std::vector<OpenMesh::Vec3d> face_normal(nf);
		for (unsigned i = 0; i < nf; ++i)
		{
			auto fv_it = mesh_->fv_iter(mesh_->face_handle(i));
			OpenMesh::Vec3d v0 = mesh_->point(fv_it);
			fH[i] = vH[fv_it.handle().idx()]; ++fv_it;
			OpenMesh::Vec3d v1 = mesh_->point(fv_it);
			fH[i] += vH[fv_it.handle().idx()]; ++fv_it;
			OpenMesh::Vec3d v2 = mesh_->point(fv_it);
			fH[i] += vH[fv_it.handle().idx()]; fH[i] /= 3.0;

			face_normal[i] = OpenMesh::cross(v1 - v0, v2 - v0); // n * area * 2
			face_area[i] = face_normal[i].norm() * 0.5;
		}

		std::vector<local_frame> vertex_f(nv);
		for (auto v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			int vertex_id = v_it.handle().idx(); OpenMesh::Vec3d n(0, 0, 0);
			for (auto vf_it = mesh_->vf_iter(v_it); vf_it; ++vf_it)
			{
				int face_id = vf_it.handle().idx();
				n += face_normal[face_id];
			}
			vertex_f[vertex_id].n = n.normalize();
			vertex_f[vertex_id].find_e_x_y();
		}

		//adjust the position
		std::vector<OpenMesh::Vec3d> new_pos(nv); OpenMesh::Vec3d p;  OpenMesh::Vec3d tp;
		Eigen::Matrix3d Trans, H_3D; Eigen::Matrix3d vertex_h;
		double v_h1; double v_h2; double v_h3; double v_h4; double v_h5; double v_h6;
		double lamda = 0.05;
		for (auto v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			int vertex_id = v_it->idx(); p = mesh_->point(v_it);
			double x1 = p[0]; double y1 = p[1]; double z1 = p[2];

			double min_radius = 1.0e30;
			double gx = 0.0; double gy = 0.0; double gz = 0.0; OpenMesh::Vec3d g;
			v_h1 = v_h2 = v_h3 = v_h4 = v_h5 = v_h6 = 0.0;
			for (auto voh_it = mesh_->voh_iter(v_it); voh_it; ++voh_it)
			{
				auto vh = mesh_->to_vertex_handle(voh_it); tp = mesh_->point(vh);
				double x2 = tp[0]; double y2 = tp[1]; double z2 = tp[2];
				double mid_edge_len = (tp - p).norm() * 0.5;
				if (mid_edge_len < min_radius) min_radius = mid_edge_len;

				auto fh = mesh_->face_handle(voh_it);
				if (fh != Mesh::InvalidFaceHandle)
				{
					int face_id = fh.idx();
					auto heh = mesh_->next_halfedge_handle(voh_it);
					vh = mesh_->to_vertex_handle(heh); tp = mesh_->point(vh);
					double x3 = tp[0]; double y3 = tp[1]; double z3 = tp[2];
					double h1 = fH[face_id][0]; double h2 = fH[face_id][1]; double h3 = fH[face_id][2];
					double h4 = fH[face_id][3]; double h5 = fH[face_id][4]; double h6 = fH[face_id][5];

					double f = (h1*x1*x1) / 24 - (h1*x1*x2) / 24 - (h1*x1*x3) / 24 + (h2*x1*y1) / 12 - (h2*x1*y2) / 24 - (h2*x1*y3) / 24 + (h3*x1*z1) / 12
						- (h3*x1*z2) / 24 - (h3*x1*z3) / 24 + (h1*x2*x2) / 24 - (h1*x2*x3) / 24 - (h2*x2*y1) / 24 + (h2*x2*y2) / 12 - (h2*x2*y3) / 24
						- (h3*x2*z1) / 24 + (h3*x2*z2) / 12 - (h3*x2*z3) / 24 + (h1*x3*x3) / 24 - (h2*x3*y1) / 24 - (h2*x3*y2) / 24 + (h2*x3*y3) / 12
						- (h3*x3*z1) / 24 - (h3*x3*z2) / 24 + (h3*x3*z3) / 12 + (h4*y1*y1) / 24 - (h4*y1*y2) / 24 - (h4*y1*y3) / 24 + (h5*y1*z1) / 12
						- (h5*y1*z2) / 24 - (h5*y1*z3) / 24 + (h4*y2*y2) / 24 - (h4*y2*y3) / 24 - (h5*y2*z1) / 24 + (h5*y2*z2) / 12 - (h5*y2*z3) / 24
						+ (h4*y3*y3) / 24 - (h5*y3*z1) / 24 - (h5*y3*z2) / 24 + (h5*y3*z3) / 12 + (h6*z1*z1) / 24 - (h6*z1*z2) / 24 - (h6*z1*z3) / 24
						+ (h6*z2*z2) / 24 - (h6*z2*z3) / 24 + (h6*z3*z3) / 24;

					double a = face_area[face_id];

					double dfx = (h1*x1) / 12 - (h1*x2) / 24 - (h1*x3) / 24 + (h2*y1) / 12 - (h2*y2) / 24 - (h2*y3) / 24 + (h3*z1) / 12 - (h3*z2) / 24 - (h3*z3) / 24;
					double dfy = (h2*x1) / 12 - (h2*x2) / 24 - (h2*x3) / 24 + (h4*y1) / 12 - (h4*y2) / 24 - (h4*y3) / 24 + (h5*z1) / 12 - (h5*z2) / 24 - (h5*z3) / 24;
					double dfz = (h3*x1) / 12 - (h3*x2) / 24 - (h3*x3) / 24 + (h5*y1) / 12 - (h5*y2) / 24 - (h5*y3) / 24 + (h6*z1) / 12 - (h6*z2) / 24 - (h6*z3) / 24;

					double dax = (2 * (y2 - y3)*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) + 2 * (z2 - z3)*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))) / (4 * std::sqrt(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) + ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2)) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))));
					double day = -(2 * (x2 - x3)*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) - 2 * (z2 - z3)*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))) / (4 * std::sqrt(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) + ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2)) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))));
					double daz = -(2 * (x2 - x3)*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2)) + 2 * (y2 - y3)*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))) / (4 * std::sqrt(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2)) + ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2)) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))));

					if (a < 1e-8) { dax = 0.0; day = 0.0; daz = 0.0; }
					double s = 10;
					double k = s * f*a;
					if (k > 120) k = 120;
					double exp_e = std::exp(k);
					g[0] = (f * dax + dfx * a); g[1] = (f * day + dfy * a); g[2] = (f * daz + dfz * a);
					gx -= OpenMesh::dot(s*exp_e*g, vertex_f[vertex_id].e_x); gy -= OpenMesh::dot(exp_e*g, vertex_f[vertex_id].e_y);

					double h = a * 2.0 / std::sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2) + (z3 - z2)*(z3 - z2));
					if (h < min_radius) min_radius = h;

					v_h1 += (s*a * h1 / 12.0 + s * s*g[0] * g[0])*exp_e; //xx
					v_h2 += (s*a * h2 / 12.0 + s * s*g[0] * g[1])*exp_e;
					v_h3 += (s*a * h3 / 12.0 + s * s*g[0] * g[2])*exp_e;
					v_h4 += (s*a * h4 / 12.0 + s * s*g[1] * g[1])*exp_e;
					v_h5 += (s*a * h5 / 12.0 + s * s*g[1] * g[2])*exp_e;
					v_h6 += (s*a * h6 / 12.0 + s * s*g[2] * g[2])*exp_e;
				}
			}
			vertex_h(0, 0) = v_h1; vertex_h(0, 1) = v_h2; vertex_h(0, 2) = v_h3;
			vertex_h(1, 0) = v_h2; vertex_h(1, 1) = v_h4; vertex_h(1, 2) = v_h5;
			vertex_h(2, 0) = v_h3; vertex_h(2, 1) = v_h5; vertex_h(2, 2) = v_h6;

			Trans(0, 0) = vertex_f[vertex_id].e_x[0]; Trans(0, 1) = vertex_f[vertex_id].e_x[1]; Trans(0, 2) = vertex_f[vertex_id].e_x[2];
			Trans(1, 0) = vertex_f[vertex_id].e_y[0]; Trans(1, 1) = vertex_f[vertex_id].e_y[1]; Trans(1, 2) = vertex_f[vertex_id].e_y[2];
			Trans(2, 0) = vertex_f[vertex_id].n[0];   Trans(2, 1) = vertex_f[vertex_id].n[1];   Trans(2, 2) = vertex_f[vertex_id].n[2];
			H_3D = Trans * vertex_h * Trans.transpose(); //H_3D.setIdentity();
			double inv_det = 1.0 / (H_3D(0, 0) * H_3D(1, 1) - H_3D(0, 1) * H_3D(1, 0));
			double dx0 = inv_det * (gx * H_3D(1, 1) - H_3D(0, 1)* gy);
			double dx1 = inv_det * (gy * H_3D(0, 0) - H_3D(1, 0)* gx);
			double move_dis = std::sqrt(dx0*dx0 + dx1 * dx1);
			if (move_dis > min_radius)
			{
				dx0 = dx0 * min_radius / move_dis;
				dx1 = dx1 * min_radius / move_dis;
			}

			new_pos[vertex_id] = p + (dx0 * vertex_f[vertex_id].e_x + dx1 * vertex_f[vertex_id].e_y) * step_length;
		}

		project_on_reference_new_p(new_pos);
		return false;
	}

	void AnisotropicMeshRemeshing::LCOT_Optimize(int iter_num, double step_length)
	{
		int iter_count = 0;
		//project_on_reference();

		std::cout << "LCOT Optimize " << std::endl;

		while (iter_count < iter_num)
		{
			flip_based_energy();
			reposition_LCOT(step_length);
			++iter_count;
			printf("%d ", iter_count);
		}
		calc_tri_quality();
	}

	void AnisotropicMeshRemeshing::delete_boundary_small_tri()
	{
		unsigned nf = mesh_->n_faces();
		for (unsigned i = 0; i < below_30_tri.size(); ++i)
		{
			if (below_30_tri_angle[i] < smallest_angle_th)
			{
				auto fh = mesh_->face_handle(below_30_tri[i]);
				for (auto fhe_it = mesh_->fh_iter(fh); fhe_it; ++fhe_it)
				{
					if (mesh_->is_boundary(mesh_->edge_handle(fhe_it)))
					{
						mesh_->delete_face(fh); break;
					}
				}
			}
		}
		mesh_->garbage_collection();
		calc_tri_quality();
	}
}

#endif
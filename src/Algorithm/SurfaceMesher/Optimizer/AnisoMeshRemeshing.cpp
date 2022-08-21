#include "AnisoMeshRemeshing.h"
namespace CADMesher
{
	void AnisoMeshRemeshing::run(double len_, double rate_)
	{
		len = len_;
		rate = rate_;
		aabbtree = globalmodel.init_trimesh_tree;
		projectMesh();
		for (int i = 0; i < 5; ++i)
		{
			double cof = std::min(1.0, 0.1*i + 0.7);
			flip();
			double step_length = 0.5 - 0.04*i;
			reposition(step_length);
			split(cof);
			flip();
			reposition(0.5);
			collapse(cof);
		}
		localOptimize(10, 1.0);
		localOptimize(5, 0.25);
	}

	void AnisoMeshRemeshing::split(double cof)
	{
		int ne = mesh->n_edges();
		Matrix3d H; OpenMesh::Vec6d h;
		Eigen::Vector3d p0, p1, sv;
		Matrix3d diag; diag.setZero();
		Matrix3d Q;
		auto openmesh2eigen = [&](O3d &op, Eigen::Vector3d &ep)
		{
			ep(0) = op[0]; ep(1) = op[1]; ep(2) = op[2];
		};
		auto eEnd = mesh->edges_end();
		for (auto te = mesh->edges_begin(); te != eEnd; ++te)
		{
			if (mesh->n_edges() > 1.2*ne)
				break;
			auto eh = te.handle();
			auto the = mesh->halfedge_handle(eh, 0);
			auto v0 = mesh->from_vertex_handle(the);
			auto v1 = mesh->to_vertex_handle(the);
			h = 0.5*(mesh->data(v0).get_Hessian() + mesh->data(v1).get_Hessian());
			H << h[0], h[1], h[2], h[1], h[3], h[4], h[2], h[4], h[5];
			Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV); 
			sv = svd.singularValues();
			diag(0, 0) = std::sqrt(sv(0)); diag(1, 1) = std::sqrt(sv(1)); diag(2, 2) = std::sqrt(sv(2));
			openmesh2eigen(mesh->point(v0), p0);
			openmesh2eigen(mesh->point(v1), p1);
			Q = svd.matrixU() * diag * svd.matrixV().transpose();
			p0 = Q * p0; p1 = Q * p1;
			double edge_len = std::sqrt((p0(0) - p1(0))*(p0(0) - p1(0)) + (p0(1) - p1(1))*(p0(1) - p1(1)) + (p0(2) - p1(2))*(p0(2) - p1(2)));

			if (edge_len > len*rate)
			{
				OpenMesh::Vec3d p = 0.5*(mesh->point(v0) + mesh->point(v1));
				OV newvert;
				if (mesh->data(eh).get_edgeflag())
				{
					bool f1 = mesh->data(eh).flag1;
					bool f2 = mesh->data(eh).flag2;
					mesh->split_edge(eh, newvert);
					if (f1)
					{
						mesh->data(mesh->edge_handle(mesh->find_halfedge(v0, newvert))).flag1 = true;
						mesh->data(mesh->edge_handle(mesh->find_halfedge(v1, newvert))).flag1 = true;
					}
					if (f2)
					{
						mesh->data(mesh->edge_handle(mesh->find_halfedge(v0, newvert))).flag2 = true;
						mesh->data(mesh->edge_handle(mesh->find_halfedge(v1, newvert))).flag2 = true;
					}
					mesh->data(newvert).set_vertflag(true);
				}
				else
				{
					newvert = mesh->add_vertex(p);
					mesh->split_edge(eh, newvert);
					projectVertex(newvert, p);
				}
			}

		}
		mesh->garbage_collection();
		dprint("split done");
	}

	void AnisoMeshRemeshing::collapse(double cof)
	{
		int ne = mesh->n_edges();
		Matrix3d H; OpenMesh::Vec6d h;
		Eigen::Vector3d p0, p1, sv;
		Matrix3d diag; diag.setZero();
		Matrix3d Q;
		auto openmesh2eigen = [&](O3d &op, Eigen::Vector3d &ep)
		{
			ep(0) = op[0]; ep(1) = op[1]; ep(2) = op[2];
		};
		for (auto te = mesh->edges_begin(); te != mesh->edges_end(); ++te)
		{
			if (!te->is_valid())
				continue;
			auto eh = te.handle();
			auto the = mesh->halfedge_handle(eh, 0);
			if (!mesh->is_collapse_ok(the))
				continue;
			auto v0 = mesh->from_vertex_handle(the);
			auto v1 = mesh->to_vertex_handle(the);
			if (mesh->data(v0).get_vertflag())
			{
				if (mesh->data(v1).get_vertflag())
					continue;
				the = mesh->opposite_halfedge_handle(the);
				std::swap(v0, v1);
			}
			h = 0.5*(mesh->data(v0).get_Hessian() + mesh->data(v1).get_Hessian());
			H << h[0], h[1], h[2], h[1], h[3], h[4], h[2], h[4], h[5];
			Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
			sv = svd.singularValues();
			diag(0, 0) = std::sqrt(sv(0)); diag(1, 1) = std::sqrt(sv(1)); diag(2, 2) = std::sqrt(sv(2));
			openmesh2eigen(mesh->point(v0), p0);
			openmesh2eigen(mesh->point(v1), p1);
			Q = svd.matrixU() * diag * svd.matrixV().transpose();
			p0 = Q * p0; p1 = Q * p1;
			double edge_len = std::sqrt((p0(0) - p1(0))*(p0(0) - p1(0)) + (p0(1) - p1(1))*(p0(1) - p1(1)) + (p0(2) - p1(2))*(p0(2) - p1(2)));


			if (edge_len < len * cof / rate)// && mesh_->calc_edge_length(eh) < min_len*a)
			{
				if (mesh->data(v1).get_vertflag())
				{
					mesh->collapse(the);
				}
				else
				{
					projectVertex(v1, 0.5*(mesh->point(v0) + mesh->point(v1)));
					mesh->collapse(the);
				}
			}
		}
		mesh->garbage_collection();
		dprint("collapse done");
	}

	void AnisoMeshRemeshing::flip()
	{
		mesh->update_face_normals();
		int nv = mesh->n_vertices();
		Eigen::Matrix3d H; Eigen::Matrix3d D; D.setZero();
		Eigen::Matrix3d vH; OV vh;
		std::vector<OpenMesh::Vec6d> vM(nv);

		for (unsigned i = 0; i < nv; ++i)
		{
			vh = mesh->vertex_handle(i);
			vM[i] = mesh->data(vh).get_Hessian();
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
		double oneThird = 1.0 / 3.0;
		while (iter_count < 10)
		{
			dprint(iter_count);
			int flip_count = 0;
			for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
			{
				if (!e_it->is_valid() || mesh->data(e_it).get_edgeflag() || !mesh->is_flip_ok(*e_it))
					continue;
				
				int edge_id = e_it->idx();
				OpenMesh::HalfedgeHandle a0 = mesh->halfedge_handle(e_it, 0);
				OpenMesh::FaceHandle fh0 = mesh->face_handle(a0); n = mesh->normal(fh0);
				OpenMesh::HalfedgeHandle b0 = mesh->halfedge_handle(e_it, 1);
				OpenMesh::FaceHandle fh1 = mesh->face_handle(b0); n += mesh->normal(fh1);
				n.normalize();

				//一定程度上防止翻折
				if (mesh->calc_sector_angle(a0) + mesh->calc_sector_angle(mesh->prev_halfedge_handle(b0)) > PI
					|| mesh->calc_sector_angle(b0) + mesh->calc_sector_angle(mesh->prev_halfedge_handle(a0)) > PI)
					continue;


				vh = mesh->to_vertex_handle(a0); int v0_id = vh.idx();
				p0 = mesh->point(vh);
				vh = mesh->to_vertex_handle(b0); int v1_id = vh.idx();
				p1 = mesh->point(vh);
				vh = mesh->to_vertex_handle(mesh->next_halfedge_handle(a0)); int v2_id = vh.idx();
				p2 = mesh->point(vh);
				vh = mesh->to_vertex_handle(mesh->next_halfedge_handle(b0)); int v3_id = vh.idx();
				p3 = mesh->point(vh);

				//#ifdef USE_PROMOTION
				//			//防止极端角出现
				//			auto ang0 = [&](OpenMesh::Vec3d p0, OpenMesh::Vec3d p1, OpenMesh::Vec3d p2) {return acos((p0 - p1).dot(p2 - p1)); };
				//			auto ang1 = [&](OpenMesh::Vec3d p0, OpenMesh::Vec3d p1, OpenMesh::Vec3d p2) {return ang0(p0, p1, p2) > least_angle && ang0(p0, p1, p2) < largest_angle
				//				&& ang0(p1,p2,p0) > least_angle && ang0(p1,p2,p0) < largest_angle && ang0(p2,p0,p1) > least_angle && ang0(p2,p0,p1) < largest_angle; };
				//			if (!ang1(p0, p2, p3) || !ang1(p1, p2, p3)) continue;
				//#endif // USE_PROMOTION


				energy0 = calc_flip_energy(p0, p2, p1, (vM[v0_id] + vM[v2_id] + vM[v1_id]) * oneThird, true);
				energy0 += calc_flip_energy(p0, p1, p3, (vM[v0_id] + vM[v1_id] + vM[v3_id]) * oneThird, true);

				energy1 = calc_flip_energy(p0, p2, p3, (vM[v0_id] + vM[v2_id] + vM[v3_id]) * oneThird, true);
				energy1 += calc_flip_energy(p2, p1, p3, (vM[v2_id] + vM[v1_id] + vM[v3_id]) * oneThird, true);

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
					mesh->flip(*e_it);
					++flip_count;
				}
			}
			if (flip_count == 0) { break; };
			++iter_count;
		}
		dprint("flip done");
		//return true;
	}

	void AnisoMeshRemeshing::reposition(double step_length)
	{
		unsigned nv = mesh->n_vertices(); //if (nv == 0) return false;
		std::vector<OpenMesh::Vec6d> vH(nv);
		for (unsigned i = 0; i < nv; ++i)
		{
			auto vh = mesh->vertex_handle(i);
			vH[i] = mesh->data(vh).get_Hessian();
		}
		unsigned nf = mesh->n_faces();
		std::vector<OpenMesh::Vec6d> fH(nf); std::vector<double> face_area(nf); std::vector<OpenMesh::Vec3d> face_normal(nf);
		double oneThird = 1.0 / 3.0;
		for (unsigned i = 0; i < nf; ++i)
		{
			auto fv_it = mesh->fv_iter(mesh->face_handle(i));
			OpenMesh::Vec3d v0 = mesh->point(fv_it);
			fH[i] = vH[fv_it.handle().idx()]; ++fv_it;
			OpenMesh::Vec3d v1 = mesh->point(fv_it);
			fH[i] += vH[fv_it.handle().idx()]; ++fv_it;
			OpenMesh::Vec3d v2 = mesh->point(fv_it);
			fH[i] += vH[fv_it.handle().idx()]; fH[i] *= oneThird;

			face_normal[i] = OpenMesh::cross(v1 - v0, v2 - v0); // n * area * 2
			face_area[i] = face_normal[i].norm() * 0.5;
		}

		std::vector<local_frame> vertex_f(nv);
		for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
		{
			int vertex_id = v_it.handle().idx(); OpenMesh::Vec3d n(0, 0, 0);
			for (auto vf_it = mesh->vf_iter(v_it); vf_it; ++vf_it)
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
		for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
		{
			int vertex_id = v_it->idx(); p = mesh->point(v_it);
			if (mesh->data(v_it).get_vertflag())
			{
				new_pos[vertex_id] = p;
				continue;
			}

			double x1 = p[0]; double y1 = p[1]; double z1 = p[2];

			double min_radius = 1.0e30;
			double gx = 0.0; double gy = 0.0; double gz = 0.0; OpenMesh::Vec3d g;
			v_h1 = v_h2 = v_h3 = v_h4 = v_h5 = v_h6 = 0.0;
			for (auto voh_it = mesh->voh_iter(v_it); voh_it; ++voh_it)
			{
				auto vh = mesh->to_vertex_handle(voh_it); tp = mesh->point(vh);
				double x2 = tp[0]; double y2 = tp[1]; double z2 = tp[2];
				double mid_edge_len = (tp - p).norm() * 0.5;
				if (mid_edge_len < min_radius) min_radius = mid_edge_len;

				auto fh = mesh->face_handle(voh_it);
				if (fh.is_valid())
				{
					int face_id = fh.idx();
					auto heh = mesh->next_halfedge_handle(voh_it);
					vh = mesh->to_vertex_handle(heh); tp = mesh->point(vh);
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
		projectMesh(new_pos);
		dprint("reposition done");
		//project_on_reference_new_p(new_pos);
		//return false;
	}

	void AnisoMeshRemeshing::initRefHessian()
	{
		initMeshStatusAndNormal(*mesh);
		initMeshStatusAndNormal(*ref_mesh);
		std::vector<double> K1, K2; std::vector<OpenMesh::Vec3d> D1, D2;
		compute_principal_curvature(ref_mesh, K1, K2, D1, D2);

		int nv = ref_mesh->n_vertices();
		Eigen::Matrix3d H; Eigen::Matrix3d D; D.setZero();
		std::vector<Eigen::Matrix3d> vH(nv); OpenMesh::Vec6d h;
		for (unsigned int i = 0; i < nv; ++i)
		{
			Mesh::VertexHandle vh = ref_mesh->vertex_handle(i);
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
			ref_mesh->data(vh).set_Hessian(h);
		}
	}

	void AnisoMeshRemeshing::projectVertex(OV &tv, O3d &pos)
	{
		auto posAndId = aabbtree->closest_point_and_face_handle(pos);
		mesh->set_point(tv, posAndId.first);

		auto fv_it = ref_mesh->cfv_iter(posAndId.second);
		auto p0 = ref_mesh->point(fv_it.handle());
		OpenMesh::Vec6d& h0 = ref_mesh->data(fv_it).get_Hessian();
		auto p1 = ref_mesh->point((++fv_it).handle());
		OpenMesh::Vec6d& h1 = ref_mesh->data(fv_it).get_Hessian();
		auto p2 = ref_mesh->point((++fv_it).handle());
		OpenMesh::Vec6d& h2 = ref_mesh->data(fv_it).get_Hessian();

		OpenMesh::Vec3d bc;
		if (!baryCoord(posAndId.first, p0, p1, p2, bc))
			bc[0] = bc[1] = bc[2] = 0.3333333333333;

		OpenMesh::Vec6d h = bc[0] * h0 + bc[1] * h1 + bc[2] * h2;
		mesh->data(tv).set_Hessian(h);
	}

	void AnisoMeshRemeshing::projectMesh()
	{
		for (auto tv = mesh->vertices_begin(); tv != mesh->vertices_end(); ++tv)
		{
			if (mesh->data(tv.handle()).get_vertflag())
				continue;
			projectVertex(tv.handle(), mesh->point(*tv));
		}
	}

	void AnisoMeshRemeshing::projectMesh(std::vector<O3d> &pos)
	{
		for (auto tv = mesh->vertices_begin(); tv != mesh->vertices_end(); ++tv)
		{
			if (mesh->data(tv.handle()).get_vertflag())
				continue;
			projectVertex(tv.handle(),pos[tv->idx()]);
		}
	}

	void AnisoMeshRemeshing::localOptimize(int iter_num, double step_length)
	{
		int iter_count = 0;
		projectMesh();

		std::cout << "LCOT Optimize " << std::endl;

		while (iter_count < iter_num)
		{
			flip();
			reposition(step_length);
			++iter_count;
			printf("%d ", iter_count);
		}
	}

	double AnisoMeshRemeshing::calc_flip_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec3d& p3, const OpenMesh::Vec6d& M, bool use_area)
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

	double AnisoMeshRemeshing::compute_src_mesh_ave_anisotropic_edge_length()
	{
		Eigen::Matrix3d H, H0, H1; Eigen::Vector3d P0, P1, P;
		Eigen::Matrix3d U; Eigen::Matrix3d V; Eigen::Vector3d sv; Eigen::Matrix3d diag_a; diag_a.setZero();
		Eigen::Matrix3d Q;
		unsigned ne = mesh->n_edges(); double sum_edge_len = 0.0;
		for (unsigned i = 0; i < ne; ++i)
		{
			Mesh::EdgeHandle eh = mesh->edge_handle(i);
			Mesh::HalfedgeHandle heh = mesh->halfedge_handle(eh, 0);
			Mesh::VertexHandle vh0_ = mesh->from_vertex_handle(heh);
			Mesh::VertexHandle vh1_ = mesh->to_vertex_handle(heh);
			OpenMesh::Vec6d& h0 = mesh->data(vh0_).get_Hessian();
			OpenMesh::Vec6d& h1 = mesh->data(vh1_).get_Hessian();
			H0(0, 0) = h0[0]; H0(0, 1) = h0[1]; H0(0, 2) = h0[2];
			H0(1, 0) = h0[1]; H0(1, 1) = h0[3]; H0(1, 2) = h0[4];
			H0(2, 0) = h0[2]; H0(2, 1) = h0[4]; H0(2, 2) = h0[5];
			H1(0, 0) = h1[0]; H1(0, 1) = h1[1]; H1(0, 2) = h1[2];
			H1(1, 0) = h1[1]; H1(1, 1) = h1[3]; H1(1, 2) = h1[4];
			H1(2, 0) = h1[2]; H1(2, 1) = h1[4]; H1(2, 2) = h1[5];
			H = 0.5 * (H0 + H1);
			OpenMesh::Vec3d p0 = mesh->point(vh0_);
			OpenMesh::Vec3d p1 = mesh->point(vh1_);
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
		double ref_mesh_ave_anisotropic_edge_length = sum_edge_len / ne;
		dprint("Src mesh ave anisotropic edge length :", ref_mesh_ave_anisotropic_edge_length);
		return ref_mesh_ave_anisotropic_edge_length;
	}
}
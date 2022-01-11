#include "TriangleMeshRemeshing.h"

namespace CADMesher
{
	double TriangleMeshRemeshing::expected_length = 0;

	TriangleMeshRemeshing::~TriangleMeshRemeshing()
	{
		if (aabbtree) { delete aabbtree; aabbtree = nullptr; }
		//if (AABB_Segment_tree) { delete AABB_Segment_tree; AABB_Segment_tree = nullptr; }
	}

	void TriangleMeshRemeshing::RemeshingMethod()
	{
		for (auto tv : mesh->vertices())
			mesh->data(tv).set_targetlength(expected_length);

		for (int i = 0; i < 20; i++)
		{
			dprint("\niteration times:", i + 1);
			split();
			collapse();
			equalize_valence();
			if (i > 10 && minAngle() < lowerAngleBound)
			{
				adjustTargetLength();
				processAngle();
			}
			tangential_relaxation();
			dprint("mesh vertices number:", mesh->n_vertices());
		}
		dprint("Remeshing Done!");

		//compute angle quality
		double min_angle = 2 * PI;
		double max_angle = 0;
		for (auto tv : mesh->vertices())
			mesh->data(tv).set_vertflag(false);

		for (auto te : mesh->edges())
			mesh->data(te).set_edgeflag(false);

		for (auto the : mesh->halfedges()) {
			if (mesh->is_boundary(the))
				continue;
			double ang = mesh->calc_sector_angle(the);
			if (ang < lowerAngleBound) {
				mesh->data(mesh->edge_handle(the)).set_edgeflag(true);
				mesh->data(mesh->from_vertex_handle(the)).set_vertflag(true);
				mesh->data(mesh->to_vertex_handle(the)).set_vertflag(true);
			}
			min_angle = std::min(ang, min_angle);
			max_angle = std::max(ang, max_angle);
		}
		std::ofstream file_writer;
		file_writer.open("tes.txt");
		double qr = 0;
		double avg = 0;
		for (auto tf : mesh->faces())
		{
			double sum = 0;
			double l = 0;
			for (auto tfe : mesh->fe_range(tf))
			{
				sum += mesh->calc_edge_length(tfe);
				l = std::max(l, mesh->calc_edge_length(tfe));
			}
			file_writer << 12 * mesh->calc_face_area(tf) / (sqrt(3)*sum*l) << std::endl;
			qr = std::min(qr, 12 * mesh->calc_face_area(tf) / (sqrt(3)*sum*l));
			avg += 12 * mesh->calc_face_area(tf) / (sqrt(3)*sum*l);
		}
		file_writer.close();
		dprintwithprecision(15, "the min angle:", min_angle, "\nthe max angle:", max_angle, qr, avg / mesh->n_faces());
	}

	void TriangleMeshRemeshing::split()
	{
		//for (auto te = mesh->edges_begin(), tee = mesh->edges_end(); te != tee; te++) {
		//	if (mesh->calc_edge_length(*te) < high || mesh->is_boundary(*te)) continue;
		//	auto newvert = mesh->add_vertex(mesh->calc_centroid(*te));
		//	bool flag = mesh->property(edgeflag, *te);
		//	pair<OpenMesh::VertexHandle, OpenMesh::VertexHandle> vert(te->v0(), te->v1());
		//	mesh->split_edge(*te, newvert);
		//	//if (mesh->property(vertflag, vert.first) && mesh->property(vertflag, vert.second)) mesh->property(vertflag, newvert) = true;
		//	if (flag) {
		//		mesh->property(vertflag, newvert) = true;
		//		mesh->property(edgeflag, mesh->edge_handle(mesh->find_halfedge(vert.first, newvert))) = true;
		//		mesh->property(edgeflag, mesh->edge_handle(mesh->find_halfedge(vert.second, newvert))) = true;
		//		continue;
		//	}
		//	OpenMesh::Vec3d pos(0.0, 0.0, 0.0);
		//	/*for (auto tvv : mesh->vv_range(newvert)) pos += mesh->point(tvv);
		//	mesh->set_point(newvert, pos / mesh->valence(newvert));*/
		//	mesh->set_point(newvert, aabbtree->closest_point(GravityPos(newvert)));
		//} 

		for (auto te : mesh->edges()) {
			//dprint(te.idx());
			double t0 = mesh->data(te.v0()).get_targetlength();
			double t1 = mesh->data(te.v1()).get_targetlength();
			if (mesh->calc_edge_length(te) < 4.0 / 3.0 * std::min(t0, t1) || mesh->is_boundary(te))
				continue;
			OV newvert = mesh->add_vertex(mesh->calc_edge_midpoint(te));
			mesh->data(newvert).set_targetlength(0.5*(t0 + t1));
			bool flag = mesh->data(te).get_edgeflag();
			std::pair<OV, OV> vert(te.v0(), te.v1());
#ifdef  OPENMESH_TRIMESH_ARRAY_KERNEL_HH
			mesh->split_edge(te, newvert);
#else
			split_one_edge(te, mesh->point(newvert));
#endif

			if (flag) {
				mesh->data(newvert).set_vertflag(true);
				mesh->data(mesh->edge_handle(mesh->find_halfedge(vert.first, newvert))).set_edgeflag(true);
				mesh->data(mesh->edge_handle(mesh->find_halfedge(vert.second, newvert))).set_edgeflag(true);
				continue;
			}
			mesh->set_point(newvert, aabbtree->closest_point(GravityPos(newvert)));
			//mesh->set_point(newvert, aabbtree->closest_point(mesh->point(newvert)));
		}
		mesh->garbage_collection();
		mesh->update_normals();
		dprint("split done");
	}

	void TriangleMeshRemeshing::collapse()
	{
		//for (auto the = mesh->halfedges_begin(); the != mesh->halfedges_end(); the++) {
		//	if (mesh->calc_edge_length(*the) >= low || !mesh->is_collapse_ok(*the)) continue;
		//	auto fromvert = mesh->from_vertex_handle(*the);
		//	auto tovert = mesh->to_vertex_handle(*the);
		//	double x = mesh->calc_edge_length(*the);
		//	if (mesh->property(edgeflag, mesh->edge_handle(*the))) {
		//		if (x >= high * 0.5) continue;
		//		int count = 0;
		//		OpenMesh::VertexHandle v = fromvert;
		//		for (auto tvoh : mesh->voh_range(fromvert)) {
		//			if (!mesh->property(edgeflag, mesh->edge_handle(tvoh))) continue;
		//			if (mesh->to_vertex_handle(tvoh)==tovert) continue;
		//			count++;
		//			v = mesh->to_vertex_handle(tvoh);
		//		}
		//		if (count == 1) {
		//			double y = (mesh->point(v) - mesh->point(tovert)).norm();
		//			double z = (mesh->point(v) - mesh->point(fromvert)).norm();
		//			double a = sqrt((x + y + z)*(x + y - z)*(x + z - y)*(y + z - x)) / 4.0;
		//			//if (400 * a > expected_length*y) continue;//要保证由x,y,z构成的三角形在y上的高，小于err=expected_length*0.05
		//			//这里为了进一步固定特征，要求更严格
		//			if (z > x) mesh->set_point(fromvert, aabbtree->closest_point((mesh->point(v) + mesh->point(tovert))*0.5));
		//			//if (z > x) mesh->set_point(fromvert, (mesh->point(v) + mesh->point(tovert))*0.5);
		//			else {
		//				mesh->collapse(*the);
		//				mesh->property(edgeflag, mesh->edge_handle(mesh->find_halfedge(v, tovert))) = true;
		//			}
		//			continue;
		//		}
		//	}
		//	//if (mesh->property(vertflag, fromvert) && mesh->property(vertflag, tovert)) {
		//	//	vector<OpenMesh::VertexHandle> movepos;
		//	//	OpenMesh::VertexHandle vflag = fromvert;
		//	//	auto tvv = mesh->vv_begin(fromvert);
		//	//	for (; tvv != mesh->vv_end(fromvert); tvv++) {
		//	//		if (mesh->property(vertflag, *tvv)) {
		//	//			auto ve0 = (mesh->point(fromvert) - mesh->point(*tvv)).normalize();
		//	//			auto ve1 = mesh->calc_edge_vector(*the).normalize();
		//	//			movepos.push_back(*tvv);
		//	//			//cout << acos(max(-1.0, min(1.0, ve0.dot(ve1)))) << endl;
		//	//			if (acos(max(-1.0, min(1.0, ve0.dot(ve1)))) < 0.1&&*tvv!=tovert) {
		//	//				vflag = *tvv;
		//	//			}
		//	//		}
		//	//	}
		//	//	if (movepos.size() == 2 && vflag != fromvert) {
		//	//		mesh->set_point(fromvert, 0.5*(mesh->point(vflag) + mesh->point(tovert)));
		//	//	}
		//	//	continue;
		//	//}
		//	if (mesh->property(vertflag, fromvert)) {
		//		/*if (mesh->property(vertflag, tovert)) {
		//			if (x <= high * 0.1) {
		//				mesh->collapse(*the);
		//			}
		//		}*/
		//		continue;
		//	}
		//	auto pos = mesh->property(vertflag, tovert) ? mesh->point(tovert) : mesh->calc_centroid(*the);
		//	//stop if collapsing results in long edges
		//	for (auto thev = mesh->vv_begin(fromvert); thev != mesh->vv_end(fromvert); thev++)
		//		if ((pos - mesh->point(*thev)).norm() > high) goto goto20210523;
		//	for (auto thev = mesh->vv_begin(tovert); thev != mesh->vv_end(tovert); thev++)
		//		if ((pos - mesh->point(*thev)).norm() > high) goto goto20210523;
		//	//stop if collapsing results in small angles
		//	auto t_f = mesh->opposite_halfedge_handle(*the);
		//	auto he = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(*the)));
		//	while (he != t_f) {
		//		auto p1 = mesh->point(mesh->from_vertex_handle(he));
		//		auto p2 = mesh->point(mesh->opposite_vh(he));
		//		auto p01 = (p1 - pos).normalize();
		//		auto p12 = (p2 - p1).normalize();
		//		auto p20 = (pos - p2).normalize();
		//		if (acos(-p01.dot(p12)) < 0.1) goto goto20210523;
		//		if (acos(-p12.dot(p20)) < 0.1) goto goto20210523;
		//		if (acos(-p20.dot(p01)) < 0.1) goto goto20210523;
		//		he = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(he));
		//	}
		//	//stop if collapsing results in intersection
		//	mesh->set_point(tovert, pos);
		//	mesh->collapse(*the);
		//	//mesh->set_point(tovert, aabbtree->closest_point(GravityPos(tovert)));
		//goto20210523:;
		//}

		for (auto the : mesh->halfedges()) {
			if (!mesh->is_collapse_ok(the)) continue;

			OV fromvert = mesh->from_vertex_handle(the);
			OV tovert = mesh->to_vertex_handle(the);
			double t0 = mesh->data(fromvert).get_targetlength();
			double t1 = mesh->data(tovert).get_targetlength();
			double x = mesh->calc_edge_length(the);
			double min_of_t0_t1 = std::min(t0, t1);

			if (x >= 4.0 / 5.0 * min_of_t0_t1) continue;

			if (mesh->data(mesh->edge_handle(the)).get_edgeflag()) {
				if (x >= 2.0 / 3.0 * min_of_t0_t1) continue;
				int count = 0;
				OV v = fromvert;
				for (auto tvoh : mesh->voh_range(fromvert)) {
					if (!mesh->data(mesh->edge_handle(tvoh)).get_edgeflag()) continue;
					if (mesh->to_vertex_handle(tvoh) == tovert) continue;
					count++;
					v = mesh->to_vertex_handle(tvoh);
				}
				if (count == 1) {
					double y = (mesh->point(v) - mesh->point(tovert)).norm();
					double z = (mesh->point(v) - mesh->point(fromvert)).norm();
					double a = sqrt((x + y + z)*(x + y - z)*(x + z - y)*(y + z - x)) / 4.0;
					//if (400 * a > expected_length*y) continue;//要保证由x,y,z构成的三角形在y上的高，小于err=expected_length*0.05
					//这里为了进一步固定特征，要求更严格
					if (z > x) mesh->set_point(fromvert, aabbtree->closest_point((mesh->point(v) + mesh->point(tovert))*0.5));
					else {
						mesh->collapse(the);
						//mesh->property(edgeflag, mesh->edge_handle(mesh->find_halfedge(v, tovert))) = true;
						mesh->data(mesh->edge_handle(mesh->find_halfedge(v, tovert))).set_edgeflag(true);
					}
				}
				continue;
			}
			if (mesh->data(mesh->edge_handle(the)).get_edgeflag()) continue;
			if (mesh->data(fromvert).get_vertflag()) {
				if (mesh->data(tovert).get_vertflag() && x < min_of_t0_t1*0.2) mesh->collapse(the);
				continue;
			}
			O3d pos = mesh->data(tovert).get_vertflag() ? mesh->point(tovert) : mesh->calc_centroid(the);

			//stop if collapsing results in long edges
			for (OV thev : mesh->vv_range(fromvert))
				if ((pos - mesh->point(thev)).norm() > 4.0 / 3.0 * min_of_t0_t1) goto goto20210523;
			for (OV thev : mesh->vv_range(tovert))
				if ((pos - mesh->point(thev)).norm() > 4.0 / 3.0 * min_of_t0_t1) goto goto20210523;

			//stop if collapsing results in small angles
			OH t_f = mesh->opposite_halfedge_handle(the);
			OH he = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(the)));
			while (he != t_f) {
				O3d p1 = mesh->point(mesh->from_vertex_handle(he));
#ifdef  OPENMESH_TRIMESH_ARRAY_KERNEL_HH
				O3d p2 = mesh->point(mesh->opposite_vh(he));
#else 
				O3d p2 = mesh->point(mesh->to_vertex_handle(mesh->next_halfedge_handle(he)));
#endif
				O3d p01 = (p1 - pos).normalize();
				O3d p12 = (p2 - p1).normalize();
				O3d p20 = (pos - p2).normalize();
				if (acos(-p01.dot(p12)) < 0.1) goto goto20210523;
				if (acos(-p12.dot(p20)) < 0.1) goto goto20210523;
				if (acos(-p20.dot(p01)) < 0.1) goto goto20210523;
				he = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(he));
			}

			//stop if collapsing results in selfintersection
			//to be completed

			mesh->set_point(tovert, pos);
			mesh->collapse(the);
		goto20210523:;
		}
		mesh->garbage_collection();
		dprint("collapse done");
	}

	void TriangleMeshRemeshing::equalize_valence()
	{

		for (auto te : mesh->edges()) {
			if (
#ifdef  OPENMESH_TRIMESH_ARRAY_KERNEL_HH
				!mesh->is_flip_ok(te) ||
#endif
				mesh->data(te).get_edgeflag()
				)
				continue;
			auto va = [&](OV v) {return mesh->valence(v) + (mesh->is_boundary(v) ? 2 : 0); };
			OH h0 = te.h0();
			OH h1 = te.h1();
			int v0 = va(te.v0());
			int v1 = va(te.v1());
#ifdef  OPENMESH_TRIMESH_ARRAY_KERNEL_HH
			int u0 = va(mesh->opposite_vh(h0));
			int u1 = va(mesh->opposite_vh(h1));
#else
			int u0 = va(mesh->to_vertex_handle(mesh->next_halfedge_handle(h0)));
			int u1 = va(mesh->to_vertex_handle(mesh->next_halfedge_handle(h1)));
#endif
			if (fabs(v0 - 6) + fabs(v1 - 6) + fabs(u0 - 6) + fabs(u1 - 6) <= fabs(v0 - 7) + fabs(v1 - 7) + fabs(u0 - 5) + fabs(u1 - 5)) continue;

			//在一定程度上可防止翻折
			double alpha0, alpha1;
			alpha0 = acos(-mesh->calc_edge_vector(h0).dot(mesh->calc_edge_vector(mesh->next_halfedge_handle(h0))));
			alpha1 = acos(-mesh->calc_edge_vector(h1).dot(mesh->calc_edge_vector(mesh->prev_halfedge_handle(h1))));
			if (alpha0 + alpha1 > PI) continue;
			alpha0 = acos(-mesh->calc_edge_vector(h0).dot(mesh->calc_edge_vector(mesh->prev_halfedge_handle(h0))));
			alpha1 = acos(-mesh->calc_edge_vector(h1).dot(mesh->calc_edge_vector(mesh->next_halfedge_handle(h1))));
			if (alpha0 + alpha1 > PI) continue;

			//检查二面角
			auto n0 = mesh->calc_face_normal(mesh->face_handle(h0));
			auto n1 = mesh->calc_face_normal(mesh->face_handle(h1));
			if (n0.dot(n1) < 0.8)
			{
				if (mesh->data(te.v0()).get_vertflag() && mesh->data(te.v1()).get_vertflag())
				{
					mesh->data(te).set_edgeflag(true);
					continue;
				}
			}


			////防止出现狭长三角形
			//auto V0 = mesh->point(te->v0());
			//auto V1 = mesh->point(te->v1());
			//auto U0 = mesh->point(mesh->opposite_vh(h0));
			//auto U1 = mesh->point(mesh->opposite_vh(h1));
			//if (((U0 - V1).norm() + (U1 - V1).norm()) / (U0 - U1).norm() < 1.1) continue;
			//if (((U0 - V0).norm() + (U1 - V0).norm()) / (U0 - U1).norm() < 1.1) continue;

			//假设flip，检查局部网格角度是否被优化
			auto opt = [&](OpenMesh::HalfedgeHandle he) {
				return fabs(mesh->calc_sector_angle(he) - PI / 3.0) + fabs(mesh->calc_sector_angle(mesh->prev_halfedge_handle(he)) - PI / 3.0) +
					fabs(mesh->calc_sector_angle(mesh->next_halfedge_handle(he)) - PI / 3.0);
			};
			double opt_before = opt(h0) + opt(h1);
#ifdef  OPENMESH_TRIMESH_ARRAY_KERNEL_HH
			mesh->flip(te);
#else
			flip_openmesh(te, *mesh);
#endif
			//若局部网格角度未被优化，则再次flip回到初始状态
#ifdef  OPENMESH_TRIMESH_ARRAY_KERNEL_HH
			if (opt_before < opt(te.h0()) + opt(te.h1()))
				mesh->flip(te);
#else
			if (opt_before < opt(te.h0()) + opt(te.h1()))
				flip_openmesh(te, *mesh);
#endif
		}
		dprint("equalize done");
	}

	void TriangleMeshRemeshing::adjustTargetLength()
	{
		double maxL = 0;
		double minL = 0;
		double sum = 0;
		double threshold = 0.1*expected_length;
		for (auto tv : mesh->vertices())
		{
			maxL = 0;
			minL = DBL_MAX;
			sum = 0;
			for (auto tve : mesh->ve_range(tv))
			{
				double l = mesh->calc_edge_length(tve);
				maxL = std::max(maxL, l);
				minL = std::min(minL, l);
				sum += l;
			}
			//mesh->data(tv).set_targetlength(mesh->data(tv).get_targetlength()*0.5 + 0.2*maxL + 0.3*minL);
			//mesh->data(tv).set_targetlength(sum / mesh->valence(tv));
			mesh->data(tv).set_targetlength(std::min(expected_length, std::max(threshold, 1.2*sum / mesh->valence(tv))));
		}
	}

	void TriangleMeshRemeshing::processAngle()
	{
#if 1
		int id = 0;
		double threshold = 0.1*expected_length;
		double minAngle = 4.0;
		for (auto tf : mesh->faces())
		{
			if (!tf.is_valid()) continue;
			id = 0;
			minAngle = 4.0;
			int i = 0;
			for (auto tfh : mesh->fh_range(tf))
			{
				double angle = mesh->calc_sector_angle(tfh);
				if (angle < minAngle)
				{
					minAngle = angle;
					id = i;
				}
				++i;
			}
			if (minAngle < lowerAngleBound)
			{
				auto h_iter = mesh->fh_begin(tf);
				for (i = 0; i < id; ++i) ++h_iter;
				auto th = mesh->prev_halfedge_handle(*h_iter);
				auto te = mesh->edge_handle(th);
				if (mesh->is_collapse_ok(th) && mesh->calc_edge_length(te) < threshold)
				{
					if (mesh->data(mesh->from_vertex_handle(th)).get_vertflag())
					{
						mesh->data(mesh->to_vertex_handle(th)).set_vertflag(true);
					}
					if (mesh->data(mesh->edge_handle(mesh->prev_halfedge_handle(th))).get_edgeflag())
					{
						mesh->data(mesh->edge_handle(*h_iter)).set_edgeflag(true);
					}
					mesh->collapse(th);
				}
				else
				{
					if (mesh->calc_sector_angle(th) < mesh->calc_sector_angle(mesh->prev_halfedge_handle(th)))
					{
						th = mesh->prev_halfedge_handle(th);
					}
					auto flagvert = mesh->to_vertex_handle(th);
					auto ph = mesh->prev_halfedge_handle(th);
					auto ne = mesh->edge_handle(mesh->next_halfedge_handle(th));
					te = mesh->edge_handle(th);
					auto pe = mesh->edge_handle(ph);

					if (!mesh->is_flip_ok(pe))
						continue;
					if (mesh->data(pe).get_edgeflag())
					{
						mesh->data(flagvert).set_vertflag(true);
						mesh->data(te).set_edgeflag(true);
						mesh->data(ne).set_edgeflag(true);
						mesh->data(pe).set_edgeflag(false);
						if (!mesh->data(flagvert).get_vertflag())
						{
							mesh->set_point(flagvert, mesh->calc_edge_midpoint(pe));
						}
					}
					mesh->flip(pe);
				}
			}
		}
		mesh->garbage_collection();
#endif
	}

	void TriangleMeshRemeshing::tangential_relaxation()
	{
		for (auto tv : mesh->vertices()) {
			if (mesh->data(tv).get_vertflag() || mesh->is_boundary(tv))
				continue;
			//mesh->set_point(tv, GravityPos(tv));
			mesh->set_point(tv, aabbtree->closest_point(GravityPos(tv)));
		}
		dprint("project done");
	}

	O3d TriangleMeshRemeshing::GravityPos(const OV &v)
	{
		OpenMesh::Vec3d sum(0, 0, 0);
		for (auto tv : mesh->vv_range(v))
		{
			sum += mesh->point(tv);
		}
		return sum / mesh->valence(v);
		//	//该函数貌似有问题 2021/11/24 by yanyisheshou
		//	vector<double> area;
		//	area.reserve(mesh->n_faces());
		//#ifdef  OPENMESH_TRIMESH_ARRAY_KERNEL_HH
		//	for (auto tf : mesh->vf_range(v))
		//		area.push_back(mesh->calc_face_area(tf));
		//#else
		//	for (auto tf : mesh->vf_range(v))
		//	{
		//		std::vector<double> el;
		//		el.reserve(3);
		//		for (auto tfe : mesh->fe_range(tf))
		//			el.push_back(mesh->calc_edge_length(tfe));
		//		double c = 0.5*(el[0] + el[1] + el[2]);
		//		dprint(c*(c - el[0])*(c - el[1])*(c - el[2]));
		//		area.push_back(std::sqrt(c*(c - el[0])*(c - el[1])*(c - el[2])));
		//	}
		//#endif
		//	double sum = accumulate(area.begin(), area.end(), 0.);
		//	area.insert(area.begin(), *area.rbegin());
		//	O3d g(0.0, 0.0, 0.0);
		//	auto itr = area.begin();
		//	for (auto tvoh : mesh->voh_range(v)) {
		//		//dprint(mesh->point(mesh->to_vertex_handle(tvoh)));
		//		g += (*itr + *(itr + 1))*mesh->point(mesh->to_vertex_handle(tvoh));
		//		itr++;
		//	}
		//	return g / (sum * 2.0);
	}

	double TriangleMeshRemeshing::minAngle()
	{
		double angle = 4;
		for (auto th : mesh->halfedges())
		{
			angle = std::min(angle, mesh->calc_sector_angle(th));
		}
		return angle;
	}

	bool TriangleMeshRemeshing::split_one_edge(Mesh::EdgeHandle& eh, OpenMesh::Vec3d& p)
	{
		Mesh::HalfedgeHandle heh0 = mesh->halfedge_handle(eh, 0);
		Mesh::HalfedgeHandle heh1 = mesh->halfedge_handle(eh, 1);
		Mesh::VertexHandle vh0 = mesh->to_vertex_handle(heh0); OpenMesh::Vec3d p0 = mesh->point(vh0);
		Mesh::VertexHandle vh1 = mesh->to_vertex_handle(heh1); OpenMesh::Vec3d p1 = mesh->point(vh1);

		std::vector<Mesh::VertexHandle> one_face(3);
		bool flag = mesh->data(eh).get_edgeflag();
		if (mesh->is_boundary(eh))
		{
			if (Mesh::InvalidFaceHandle != mesh->face_handle(heh0))
			{
				Mesh::VertexHandle vh2 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh0));
				OpenMesh::Vec3d p2 = mesh->point(vh2);
				OpenMesh::Vec3d n = OpenMesh::cross(p1 - p2, p0 - p2).normalize();
				double a1 = OpenMesh::dot(n, OpenMesh::cross(p2 - p0, p - p0));
				double a2 = OpenMesh::dot(n, OpenMesh::cross(p1 - p2, p - p2));
				if (a1 < 1e-8 || a2 < 1e-8) return false;
				Mesh::VertexHandle vh = mesh->add_vertex(p);
				mesh->delete_edge(eh, false); mesh->garbage_collection();
				one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh->add_face(one_face);
				one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh->add_face(one_face);
			}
			else
			{
				Mesh::VertexHandle vh3 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh1));
				OpenMesh::Vec3d p3 = mesh->point(vh3);
				OpenMesh::Vec3d n = OpenMesh::cross(p0 - p3, p1 - p3).normalize();
				double a1 = OpenMesh::dot(n, OpenMesh::cross(p0 - p3, p - p3));
				double a2 = OpenMesh::dot(n, OpenMesh::cross(p3 - p1, p - p1));
				if (a1 < 1e-8 || a2 < 1e-8) return false;
				Mesh::VertexHandle vh = mesh->add_vertex(p);
				mesh->delete_edge(eh, false); mesh->garbage_collection();
				one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh->add_face(one_face);
				one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh->add_face(one_face);
			}
		}
		else
		{
			Mesh::VertexHandle vh2 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh0)); OpenMesh::Vec3d p2 = mesh->point(vh2);
			OpenMesh::Vec3d n1 = OpenMesh::cross(p1 - p2, p0 - p2).normalize();
			double a1 = OpenMesh::dot(n1, OpenMesh::cross(p2 - p0, p - p0));
			double a2 = OpenMesh::dot(n1, OpenMesh::cross(p1 - p2, p - p2));
			Mesh::VertexHandle vh3 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh1)); OpenMesh::Vec3d p3 = mesh->point(vh3);
			OpenMesh::Vec3d n2 = OpenMesh::cross(p0 - p3, p1 - p3).normalize();
			double a3 = OpenMesh::dot(n2, OpenMesh::cross(p0 - p3, p - p3));
			double a4 = OpenMesh::dot(n2, OpenMesh::cross(p3 - p1, p - p1));
			if (a1 < 1e-8 || a2 < 1e-8 || a3 < 1e-8 || a4 < 1e-8) return false;
			Mesh::VertexHandle vh = mesh->add_vertex(p);
			mesh->delete_edge(eh, false); mesh->garbage_collection();
			one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh->add_face(one_face);
			one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh->add_face(one_face);
			one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh->add_face(one_face);
			one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh->add_face(one_face);
		}
		return true;
	}
}
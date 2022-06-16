#include "TriangleMeshRemeshing.h"

namespace CADMesher
{
	void TriangleMeshRemeshing::run()
	{
#if 0
		for (auto &tv : mesh->vertices())
		{
			if (tv.valence() <= 2 && !tv.is_boundary() || isnan(mesh->data(tv).GaussCurvature))
			{
				int p = 0;
			}
			int p = 0;
			mesh->data(tv).set_targetlength(mesh->data(tv).GaussCurvature > 1.0e-4 ?
				expected_length * 2 / (6 + Log10(mesh->data(tv).GaussCurvature)) : expected_length);
		}
#else
		initTargetLength();
#endif
		
		TriMeshQualityHelper tmqh(mesh);
		tmqh.print();
		dprint("\nMesh vertices number in initialization:", mesh->n_vertices());
		//mesh->delete_isolated_vertices();
		//mesh->garbage_collection();

		tr.refresh();

		int itertimes = 0;
		/*for (; itertimes < iterPeriod[0]; ++itertimes)
		{
			tr.mark();
			split();
			collapse();

			dprint("mesh vertices number:", mesh->n_vertices());
#ifdef printRemeshingInfo
			tmqh.update();
			tmqh.print();
			tr.pastMark("the " + std::to_string(itertimes + 1) + "-th iteration time:");
			dprint();
#endif
		}*/
		//return;

		dprint("global remeshing");
		itertimes = 0;
		for (; itertimes < 5; ++itertimes)
		{
			tr.mark();
			split();
			collapse();
			equalize_valence();
			tangential_relaxation();

			dprint("mesh vertices number:", mesh->n_vertices());
#ifdef printRemeshingInfo
			tmqh.update();
			tmqh.print();
			tr.pastMark("the " + std::to_string(itertimes + 1) + "-th iteration time:");
			dprint();
#endif
		}
		globalProject();

		dprint("remeshing while eliminating small angle");
		itertimes = 0;
		for (; itertimes < 5; ++itertimes)
		{
			tr.mark();
			if (processFeatureConstraintAngle() < 20)
			{
				tmqh.update();
				if(tmqh.getAvgQuality() > 0.88)
					break;
			}
			adjustTargetLength();
			split();
			collapse(true);
			equalize_valence(true);
			tangential_relaxation();

			dprint("mesh vertices number:", mesh->n_vertices());
#ifdef printRemeshingInfo
			tmqh.update();
			tmqh.print();
			tr.pastMark("the " + std::to_string(itertimes + 1) + "-th iteration time:");
			dprint();
#endif
		}

		dprint("aggressively eliminating small angle");
		itertimes = 0;
		while (processFeatureConstraintAngle(true) && ++itertimes < 10);

		for (auto &tv : mesh->vertices())
		{
			if (tv.valence() <= 2 && !tv.is_boundary())
			{
				int p = 0;
				dprint(tv.idx(), tv.valence());
			}
		}

		tmqh.update();
		tmqh.print();
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
		if (polymeshInput)
		{
			assembleMesh();
		}
#endif
		dprint("Remeshing Done!\n");
	}

	void TriangleMeshRemeshing::split()
	{
		int splitnumber = 0;
		for (auto &te : mesh->edges()) {
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
			if (polymeshInput && te.v0().idx() < boundaryNum && te.v1().idx() < boundaryNum)
			{
				continue;
			}
#endif
			if (split_one_edge(te))
				++splitnumber;
		}
		mesh->garbage_collection();
		//mesh->update_normals();
		//dprint("split done");
#ifdef printRemeshingInfo
		tr.out("split:");
		dprint("split number:", splitnumber);
#endif
	}

	bool TriangleMeshRemeshing::split_one_edge(const OpenMesh::SmartEdgeHandle &te, bool ifRelaxCondition)
	{
		//dprint(te.idx());
		SOV vert[2] = { te.v0(), te.v1() };
		double t0 = mesh->data(vert[0]).get_targetlength();
		double t1 = mesh->data(vert[1]).get_targetlength();
		if (!ifRelaxCondition && mesh->calc_edge_length(te) < 1.333 * std::min(t0, t1))
			return false;
		//dprint(te.idx(), te.v0().idx(), te.v1().idx());
		OV newvert = mesh->add_vertex(mesh->calc_edge_midpoint(te));
		mesh->data(newvert).set_targetlength(std::min(0.5*(t0 + t1), 1.5 * std::min(t0, t1)));
		//bool flag = mesh->data(te).get_edgeflag();
		bool flag1 = mesh->data(te).flag1;
		bool flag2 = mesh->data(te).flag2;
		mesh->split_edge(te, newvert);

		if (flag1)
		{
			mesh->data(newvert).set_vertflag(true);
			mesh->data(mesh->edge_handle(mesh->find_halfedge(vert[0], newvert))).flag1 = true;
			mesh->data(mesh->edge_handle(mesh->find_halfedge(vert[1], newvert))).flag1 = true;
		}
		else if (flag2)
		{
			mesh->data(newvert).set_vertflag(true);
			mesh->data(mesh->edge_handle(mesh->find_halfedge(vert[0], newvert))).flag2 = true;
			mesh->data(mesh->edge_handle(mesh->find_halfedge(vert[1], newvert))).flag2 = true;
		}
		return true;
		/*else
		{
			mesh->set_point(newvert, aabbtree->closest_point(GravityPos(newvert)));
		}*/
	}

	void TriangleMeshRemeshing::collapse(bool ifEnhanced)
	{
		int collapsenumber = 0;
		auto edgeFlag = [&](OpenMesh::SmartEdgeHandle &e)
		{
			return mesh->data(e).flag1 || mesh->data(e).flag2;
		};
		for (auto the : mesh->halfedges()) {
			if (!mesh->is_collapse_ok(the) || mesh->is_boundary(the.from()) && !mesh->is_boundary(the.to()))
			{
				continue;
			}
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
			if (polymeshInput && the.from().idx() < boundaryNum)
			{
				continue;
			}
#endif
			/*if (mesh->data(the.edge()).flag1 || mesh->data(the.edge()).flag2)
				continue;
			if (mesh->is_boundary(the))
				continue;*/

			OV fromvert = mesh->from_vertex_handle(the);
			if (mesh->data(fromvert).get_vertflag())
				continue;

			OV tovert = mesh->to_vertex_handle(the);
			double t0 = mesh->data(fromvert).get_targetlength();
			double t1 = mesh->data(tovert).get_targetlength();
			double x = mesh->calc_edge_length(the);
			double min_of_t0_t1 = std::min(t0, t1);

			if (x >= 0.8 * min_of_t0_t1)
				continue;

#if 1
			//特征线上移动点
			//if (mesh->data(the.edge()).get_edgeflag())
			if (edgeFlag(the.edge()))
			{
				if (x >= 0.666 * min_of_t0_t1) continue;
				int count = 0;
				OV v = fromvert;
				for (auto &tvoh : mesh->voh_range(fromvert)) 
				{
					//if (!mesh->data(tvoh.edge()).get_edgeflag()) continue;
					if (!edgeFlag(tvoh.edge())) continue;
					if (tvoh.to() == tovert) continue;
					++count;
					v = tvoh.to();
				}
				if (count == 1) 
				{
					double y = (mesh->point(v) - mesh->point(tovert)).norm();
					double z = (mesh->point(v) - mesh->point(fromvert)).norm();
					//double a = sqrt((x + y + z)*(x + y - z)*(x + z - y)*(y + z - x)) * 0.25;
					//if (400 * a > expected_length*y) continue;//要保证由x,y,z构成的三角形在y上的高，小于err=expected_length*0.05
					//这里为了进一步固定特征，要求更严格
					/*if (z > x)
					{
						mesh->set_point(fromvert, aabbtree->closest_point((mesh->point(v) + mesh->point(tovert))*0.5));
					}
					else {
						mesh->collapse(the);
						mesh->data(mesh->edge_handle(mesh->find_halfedge(v, tovert))).set_edgeflag(true);
					}*/
					if (mesh->calc_sector_angle(the) < lowerAngleBound)
					{
						if (mesh->is_collapse_ok(the.prev())) 
						{
							mesh->set_point(fromvert, mesh->point(v));
							mesh->collapse(the.prev());
						}
					}
					else if (mesh->calc_sector_angle(the.next()) < lowerAngleBound)
					{
						mesh->set_point(fromvert, mesh->point(tovert));
						mesh->collapse(the);
					}
				}
				continue;
			}
			//if (mesh->data(mesh->edge_handle(the)).get_edgeflag()) continue;
			//移除短边
			/*if (mesh->data(fromvert).get_vertflag()) 
			{
				if (mesh->data(tovert).get_vertflag() && x < min_of_t0_t1*0.1)
				{
					mesh->collapse(the);
				}
				continue;
			}*/
#else
			if (mesh->data(the.edge()).flag1 || mesh->data(the.edge()).flag2)
				continue;
#endif
			if (mesh->is_boundary(the))
				continue;
			O3d pos = mesh->data(tovert).get_vertflag() ? mesh->point(tovert) : mesh->calc_centroid(the);

			if (ifEnhanced)
			{
				//stop if collapsing results in long edges
				for (OV thev : mesh->vv_range(fromvert))
					if ((pos - mesh->point(thev)).norm() > 1.33 * min_of_t0_t1)
						goto goto20210523;
				for (OV thev : mesh->vv_range(tovert))
					if ((pos - mesh->point(thev)).norm() > 1.33 * min_of_t0_t1)
						goto goto20210523;

				//stop if collapsing results in small angles
				auto t_f = the.opp();
				auto he = the.prev().opp().prev();
				while (he != t_f) {
					O3d p1 = mesh->point(he.from());
					O3d p2 = mesh->point(mesh->opposite_vh(he));
					O3d p01 = (p1 - pos).normalize();
					O3d p12 = (p2 - p1).normalize();
					O3d p20 = (pos - p2).normalize();
					if (acos(-p01.dot(p12)) < lowerAngleBound) goto goto20210523;
					if (acos(-p12.dot(p20)) < lowerAngleBound) goto goto20210523;
					if (acos(-p20.dot(p01)) < lowerAngleBound) goto goto20210523;
					he = he.opp().prev();
				}
			}
			//stop if collapsing results in selfintersection
			//to be completed

			mesh->set_point(tovert, pos);
			mesh->collapse(the);
			++collapsenumber;
		goto20210523:;
		}
		mesh->garbage_collection();
		//dprint("collapse done");
#ifdef printRemeshingInfo
		tr.out("collapse:");
		dprint("collapse number:", collapsenumber);
#endif
	}

	void TriangleMeshRemeshing::equalize_valence(bool ifEnhanced)
	{
		auto va = [&](OV v) {return mesh->valence(v) + (mesh->is_boundary(v) ? 2 : 0); };
		double cosLowerAngleBound = cos(lowerAngleBound);

		/*double PiDividedByThree = PI * 0.3333;
		double temp = PiDividedByThree - lowerAngleBound; temp *= temp;
		auto penaltyFunction = [&](double a)
		{
			double t = a - PiDividedByThree; t *= t;
			return t > temp ? t * t : t;
		};
		auto opt = [&](OpenMesh::SmartHalfedgeHandle &he) {
			return penaltyFunction(mesh->calc_sector_angle(he)) + penaltyFunction(mesh->calc_sector_angle(he.prev()))
				+ penaltyFunction(mesh->calc_sector_angle(he.next()));
		};*/

//		//除去度为3的顶点
//		for (auto &tv : mesh->vertices())
//		{
//			if (mesh->valence(tv) == 3)
//			{
//				for (auto &tve : mesh->ve_range(tv))
//				{
//					if (tve.is_boundary())
//						continue;
//#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
//					if (polymeshInput && tve.v0().idx() < boundaryNum && tve.v1().idx() < boundaryNum)
//					{
//						continue;
//					}
//#endif
//					split_one_edge(tve, true);
//					break;
//				}
//			}
//		}
//		mesh->garbage_collection();

		int equalizenumber = 0;
		for (auto te : mesh->edges()) {
			if (!mesh->is_flip_ok(te) || mesh->data(te).flag1 || mesh->data(te).flag2)
				continue;
			auto h0 = te.h0();
			auto h1 = te.h1();
			int v0 = va(h0.from());
			int v1 = va(h0.to());
			int u0 = va(mesh->opposite_vh(h0));
			int u1 = va(mesh->opposite_vh(h1));
			if (fabs(v0 - 6) + fabs(v1 - 6) + fabs(u0 - 6) + fabs(u1 - 6) <= fabs(v0 - 7) + fabs(v1 - 7) + fabs(u0 - 5) + fabs(u1 - 5))
				continue;

			//在一定程度上可防止翻折
			double alpha0, alpha1;
			alpha0 = acos(-mesh->calc_edge_vector(h0).dot(mesh->calc_edge_vector(/*mesh->next_halfedge_handle(h0)*/h0.next())));
			alpha1 = acos(-mesh->calc_edge_vector(h1).dot(mesh->calc_edge_vector(/*mesh->prev_halfedge_handle(h1)*/h1.prev())));
			if (alpha0 + alpha1 > PI)
				continue;
			alpha0 = acos(-mesh->calc_edge_vector(h0).dot(mesh->calc_edge_vector(/*mesh->prev_halfedge_handle(h0)*/h0.prev())));
			alpha1 = acos(-mesh->calc_edge_vector(h1).dot(mesh->calc_edge_vector(/*mesh->next_halfedge_handle(h1)*/h1.next())));
			if (alpha0 + alpha1 > PI)
				continue;

			////检查二面角
			//auto n0 = mesh->calc_face_normal(/*mesh->face_handle(h0)*/h0.face());
			//auto n1 = mesh->calc_face_normal(/*mesh->face_handle(h1)*/h1.face());
			//if (n0.dot(n1) < 0.8)
			//{
			//	if (mesh->data(te.v0()).get_vertflag() && mesh->data(te.v1()).get_vertflag())
			//	{
			//		mesh->data(te).set_edgeflag(true);
			//		continue;
			//	}
			//}
			////尽量防止出现小角
			//if (mesh->calc_sector_angle(h0.next()) < lowerAngleBound)
			//	continue;
			//if (mesh->calc_sector_angle(h1.next()) < lowerAngleBound)
			//	continue;
			//假设flip，检查局部网格角度是否被优化
			//double opt_before = opt(h0) + opt(h1);


			////防止出现狭长三角形
			//auto V0 = mesh->point(te->v0());
			//auto V1 = mesh->point(te->v1());
			//auto U0 = mesh->point(mesh->opposite_vh(h0));
			//auto U1 = mesh->point(mesh->opposite_vh(h1));
			//if (((U0 - V1).norm() + (U1 - V1).norm()) / (U0 - U1).norm() < 1.1) continue;
			//if (((U0 - V0).norm() + (U1 - V0).norm()) / (U0 - U1).norm() < 1.1) continue;

			//检测flip是否会导致小角产生
			if (ifEnhanced)
			{
				auto U01 = (mesh->point(h1.next().to()) - mesh->point(h0.next().to())).normalized();
				auto temp = (mesh->point(h0.from()) - mesh->point(h0.next().to())).normalized();
				if (U01.dot(temp) > cosLowerAngleBound)
					continue;
				temp = (mesh->point(h0.to())- mesh->point(h0.next().to())).normalized();
				if (U01.dot(temp) > cosLowerAngleBound)
					continue;
				temp = (mesh->point(h1.next().to()) - mesh->point(h0.from())).normalized();
				if (U01.dot(temp) > cosLowerAngleBound)
					continue;
				temp = (mesh->point(h1.next().to()) - mesh->point(h0.to())).normalized();
				if (U01.dot(temp) > cosLowerAngleBound)
					continue;
			}

			mesh->flip(te);
			++equalizenumber;
			//若局部网格角度未被优化，则再次flip回到初始状态
			//if (opt_before < opt(te.h0()) + opt(te.h1()))
				//mesh->flip(te);
		}
		mesh->garbage_collection();

		//dprint("equalize done");
#ifdef printRemeshingInfo
		tr.out("equalize:");
		dprint("equalize number:", equalizenumber);
#endif
	}

	void TriangleMeshRemeshing::adjustTargetLength()
	{
		//double maxL = 0;
		//double minL = expected_length;
		//double sum = 0;
		//double threshold = 0.15 * expected_length;
		//for (auto tv : mesh->vertices())
		//{
		//	if (!mesh->data(tv).get_vertflag())
		//		continue;
		//	mesh->data(tv).set_targetlength(mesh->data(tv).get_targetlength()+)
		//	/*sum = 0;
		//	for (auto tve : mesh->ve_range(tv))
		//	{
		//		sum += mesh->calc_edge_length(tve);
		//	}
		//	mesh->data(tv).set_targetlength(std::min(expected_length, std::max(threshold, 1.2*sum / mesh->valence(tv))));*/
		//}
		std::vector<double> tl; tl.reserve(mesh->n_vertices());
		for (auto &tv : mesh->vertices())
		{
			double l = 0;
			for (auto &tvv : mesh->vv_range(tv))
			{
				l += mesh->data(tvv).get_targetlength();
			}
			tl.push_back(0.5 * (l / mesh->valence(tv) + mesh->data(tv).get_targetlength()));
		}
		for (auto &tv : mesh->vertices())
			mesh->data(tv).set_targetlength(tl[tv.idx()]);

		//for (auto &te : mesh->edges())
		//{
		//	/*if (!mesh->data(te).flag1 && !mesh->data(te).flag2)
		//		continue;
		//	double l = mesh->calc_edge_length(te);
		//	mesh->data(te.v0()).set_targetlength(l);
		//	mesh->data(te.v1()).set_targetlength(l);*/
		//	
		//}
#ifdef printRemeshingInfo
		tr.out("adjust target length:");
#endif
	}

	int TriangleMeshRemeshing::processFeatureConstraintAngle(bool ifEnhanced)
	{
		auto edgeFlag = [&](OpenMesh::SmartEdgeHandle &e)
		{
			return mesh->data(e).flag1 || mesh->data(e).flag2;
		};
		auto setAvgLength = [&](OV &v)
		{
			double l = 0;
			for (auto &tvv : mesh->vv_range(v))
			{
				l += mesh->data(tvv).get_targetlength();
			}
			mesh->data(v).set_targetlength(0.5 * (l / mesh->valence(v) + mesh->data(v).get_targetlength()));
		};

		int processNumber = 0;
		int findNumber = 0;
		//处理有两个角小于阈值的情况
		for (auto &th : mesh->halfedges())
		{
			if (th.is_boundary() || !th.is_valid() || mesh->calc_sector_angle(th) > lowerAngleBound)
				continue;
			++findNumber;
			double angle0 = mesh->calc_sector_angle(th.next());
			double angle1 = mesh->calc_sector_angle(th.prev());
			SOH CB = th;
			bool isTwoSmallAngles = false;
			if (angle0 < lowerAngleBound)
				isTwoSmallAngles = true;
			else if (angle1 < lowerAngleBound)
			{
				CB = th.prev();
				isTwoSmallAngles = true;
			}
			if (isTwoSmallAngles)
			{
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
				if (polymeshInput && CB.to().idx() < boundaryNum && CB.next().to().idx() < boundaryNum)
				{
					continue;
				}
#endif
				OE BA(CB.next().edge());
				OE AC(CB.prev().edge());
				if (mesh->data(BA).flag1)
				{
					mesh->data(BA).flag1 = false;
					mesh->data(AC).flag1 = true;
					mesh->data(CB.edge()).flag1 = true;
				}
				if (mesh->data(BA).flag2)
				{
					mesh->data(BA).flag2 = false;
					mesh->data(AC).flag2 = true;
					mesh->data(CB.edge()).flag2 = true;
				}
				if (edgeFlag(CB.edge()) || edgeFlag(CB.next().edge()) || edgeFlag(CB.prev().edge()))
				{
					mesh->data(CB.from()).set_vertflag(true);
				}
				if (mesh->is_flip_ok(BA))
					mesh->flip(BA);
				else
				{
					split_one_edge(CB.next().edge(), true);
					OV newvert = mesh->vertex_handle(mesh->n_vertices() - 1);
					mesh->set_point(newvert, aabbtree->closest_point(GravityPos(newvert)));
					//setAvgLength(newvert);
					//setAvgLength(CB.from());
				}
				++processNumber;
			}
		}
		mesh->garbage_collection();
		if (!findNumber)
		{
#ifdef printRemeshingInfo
			tr.out("processFeatureAngle:");
			dprint("processNumber:", processNumber, "\tfindNunber:", findNumber);
#endif
			return findNumber;
		}
		//return findNumber;
		//处理有一个角小于阈值的情况
		for (auto &th : mesh->halfedges())
		{
			if (th.is_boundary() || !th.is_valid() || mesh->calc_sector_angle(th) > lowerAngleBound)
				continue;
			SOH CB = th;
			if (!mesh->is_collapse_ok(th.prev()))
			{
				//dprint(th.to().idx(), th.prev().opp().prev().from().idx(), th.to().valence(), th.prev().opp().prev().from().valence(), th.from());
				//dprint(mesh->calc_sector_angle(th), mesh->calc_sector_angle(th.next()), mesh->calc_sector_angle(th.prev()));
				if (ifEnhanced)
				{
					if (th.to().valence() == 3)
						CB = th.next().opp().prev();
					else if (th.prev().opp().prev().from().valence() == 3)
						CB = th.prev().opp().prev().opp().prev();
					if (mesh->is_collapse_ok(CB.prev()) && !mesh->is_boundary(CB.prev().edge()))
						goto goto20220511;
					else
						continue;
				}
				else
					continue;
			}
			++findNumber;
			double angle0 = mesh->calc_sector_angle(th.next());
			double angle1 = mesh->calc_sector_angle(th.prev());
			if (angle0 > lowerAngleBound && angle1 > lowerAngleBound)
			{
			goto20220511:;
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
				if (polymeshInput && CB.from().idx() < boundaryNum && CB.prev().from().idx() < boundaryNum)
				{
					continue;
				}
#endif
				if (mesh->data(CB.next().edge()).flag1)
					mesh->data(CB.edge()).flag1 = true;
				if (mesh->data(CB.next().edge()).flag2)
					mesh->data(CB.edge()).flag2 = true;
				auto temp = CB.prev().opp().prev();
				if (mesh->data(temp.prev().edge()).flag1)
					mesh->data(temp.edge()).flag1 = true;
				if (mesh->data(temp.prev().edge()).flag2)
					mesh->data(temp.edge()).flag2 = true;
				if (edgeFlag(CB.edge()) || edgeFlag(CB.next().edge()) || edgeFlag(CB.prev().edge()) || edgeFlag(temp.edge()) || edgeFlag(temp.prev().edge()))
					mesh->data(CB.from()).set_vertflag(true);
				mesh->collapse(CB.prev());
				++processNumber;
			}
		}
		mesh->garbage_collection();
#ifdef printRemeshingInfo
		tr.out("processFeatureAngle:");
		dprint("processNumber:", processNumber, "\tfindNunber:", findNumber);
#endif
		return findNumber;
	}

	void TriangleMeshRemeshing::processAngle()
	{
		int id = 0;
		/*double upperAngleBound = PI - lowerAngleBound * 2;
		for (auto &th : mesh->halfedges())
		{
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
			if (polymeshInput && th.from().idx() < boundaryNum && th.to().idx() < boundaryNum)
			{
				continue;
			}
#endif
			if (th.is_boundary() || !th.is_valid())
				continue;
			if (mesh->calc_sector_angle(th) < upperAngleBound)
				continue;
			split_one_edge(th.prev().edge(), id);
		}
		mesh->garbage_collection();*/

		for (auto &th : mesh->halfedges())
		{
			if (th.is_boundary() || !th.is_valid())
				continue;
			if (mesh->calc_sector_angle(th) > lowerAngleBound)
				continue;
			double angle0 = mesh->calc_sector_angle(th.next());
			double angle1 = mesh->calc_sector_angle(th.prev());

			if (angle0 < lowerAngleBound)
			{
				if (mesh->data(th.next().edge()).flag1 || mesh->data(th.next().edge()).flag2 || mesh->is_boundary(th.from()))
					continue;
				mesh->data(th.next().edge()).flag1 = mesh->data(th.edge()).flag1 || mesh->data(th.prev().edge()).flag1;
				mesh->data(th.next().edge()).flag2 = mesh->data(th.edge()).flag2 || mesh->data(th.prev().edge()).flag2;
				mesh->data(th.edge()).flag1 = false;
				mesh->data(th.edge()).flag2 = false;
				mesh->data(th.prev().edge()).flag1 = false;
				mesh->data(th.prev().edge()).flag2 = false;
				mesh->data(th.from()).set_vertflag(false);
				mesh->set_point(th.from(), aabbtree->closest_point(GravityPos(th.from())));
			}
			else if (angle1 < lowerAngleBound)
			{
				if (mesh->data(th.edge()).flag1 || mesh->data(th.edge()).flag2 || mesh->is_boundary(th.prev().from()))
					continue;
				mesh->data(th.edge()).flag1 = mesh->data(th.next().edge()).flag1 || mesh->data(th.prev().edge()).flag1;
				mesh->data(th.edge()).flag2 = mesh->data(th.next().edge()).flag2 || mesh->data(th.prev().edge()).flag2;
				mesh->data(th.next().edge()).flag1 = false;
				mesh->data(th.next().edge()).flag2 = false;
				mesh->data(th.prev().edge()).flag1 = false;
				mesh->data(th.prev().edge()).flag2 = false;
				mesh->data(th.prev().from()).set_vertflag(false);
				mesh->set_point(th.prev().from(), aabbtree->closest_point(GravityPos(th.prev().from())));
			}
#if 1
			else if (mesh->is_collapse_ok(th.prev()) && !mesh->is_boundary(th.from()) && !mesh->is_boundary(th.prev().from()))
			{
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
				if (polymeshInput && th.from().idx() < boundaryNum)
				{
					continue;
				}

#endif
				mesh->data(th.edge()).flag1 = mesh->data(th.edge()).flag1 || mesh->data(th.next().edge()).flag1;
				mesh->data(th.edge()).flag2 = mesh->data(th.edge()).flag2 || mesh->data(th.next().edge()).flag2;
				mesh->data(th.next().edge()).flag1 = false;
				mesh->data(th.next().edge()).flag2 = false;
				mesh->data(th.prev().edge()).flag1 = false;
				mesh->data(th.prev().edge()).flag2 = false;
				mesh->data(th.from()).set_vertflag(mesh->data(th.prev().from()).get_vertflag() || mesh->data(th.from()).get_vertflag());
				mesh->data(th.prev().from()).set_vertflag(false);

				auto tv = th.from();
				{
					mesh->collapse(th.prev());
					mesh->set_point(tv, aabbtree->closest_point(GravityPos(tv)));

				}
			}
#else
			else 
			{
				mesh->data(th.from()).set_targetlength(mesh->data(th.from()).get_targetlength()*0.5);
				mesh->data(th.prev().from()).set_targetlength(mesh->data(th.prev().from()).get_targetlength()*0.5);
			}
#endif
		}
		mesh->garbage_collection();

		//#if 1
		//		for (auto &tf : mesh->faces())
		//		{
		//			if (!tf.is_valid()) 
		//				continue;
		//			minAngle = 4.0;
		//			for (auto &tfh : mesh->fh_range(tf))
		//			{
		//				double angle = mesh->calc_sector_angle(tfh);
		//				if (angle < minAngle)
		//				{
		//					minAngle = angle;
		//					th = tfh;
		//				}
		//			}
		//			if (minAngle < lowerAngleBound)
		//			{
		//				double angle0 = mesh->calc_sector_angle(th.next());
		//				double angle1 = mesh->calc_sector_angle(th.prev());
		//#if 1
		//				if (angle0 < lowerAngleBound)
		//				{
		//					split_one_edge(th.next().edge(), id);
		//				}
		//				else if (angle1 < lowerAngleBound)
		//				{
		//					split_one_edge(th.edge(), id);
		//				}
		//				else
		//				{
		//					mesh->data(th.edge()).flag1 = mesh->data(th.edge()).flag1 || mesh->data(th.next().edge()).flag1;
		//					mesh->data(th.edge()).flag2 = mesh->data(th.edge()).flag2 || mesh->data(th.next().edge()).flag2;
		//					mesh->data(th.next().edge()).flag1 = false;
		//					mesh->data(th.next().edge()).flag2 = false;
		//					mesh->data(th.prev().edge()).flag1 = false;
		//					mesh->data(th.prev().edge()).flag2 = false;
		//					mesh->data(th.prev().from()).set_vertflag(false);
		//				}
		//
		//#else
		//				//三角形三点共线
		//				if (angle0 < lowerAngleBound)
		//				{
		//					if (mesh->is_flip_ok(th.next().edge()))
		//					{
		//						mesh->data(th.next().edge()).flag1 = false;
		//						mesh->data(th.next().edge()).flag2 = false;
		//						mesh->flip(th.next().edge());
		//					}
		//					//split_one_edge(th.next().edge(), id);
		//				}
		//				else if (angle1 < lowerAngleBound)
		//				{
		//					if (mesh->is_flip_ok(th.edge()))
		//					{
		//						//mesh->data(th.edge()).flag1 = false;
		//						mesh->data(th.edge()).flag2 = false;
		//						mesh->flip(th.edge());
		//					}
		//					//th = th.prev();
		//					//split_one_edge(th.next().edge(), id);
		//				}
		//				//th = mesh->find_halfedge(mesh->vertex_handle(mesh->n_vertices() - 1), mesh->vertex_handle(th.from().idx())).next();
		//				//只有一个顶点的角很小
		//				else
		//				{
		//					if (!mesh->is_collapse_ok(th.prev()))
		//					{
		//						int p = 0;
		//						//auto tv = th.from().valence() == 3 ? th.from() : (th.to().valence() == 3 ? th.to() : th.next().to());
		//
		//						split_one_edge(th.prev().edge(), p);
		//						split_one_edge(th.next().opp().prev().edge(), p);
		//					}
		//					else
		//					{
		//
		//#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
		//						if (polymeshInput)
		//						{
		//							if (th.prev().from().idx() >= boundaryNum)
		//							{
		//								mesh->collapse(th.prev());
		//							}
		//							else if (th.prev().to().idx() >= boundaryNum)
		//							{
		//								mesh->collapse(th.prev().opp());
		//							}
		//							continue;
		//						}
		//#endif
		//						mesh->collapse(th.prev());
		//					}
		//				}
		//#endif
		//			}
		//		}
		//#else
		//		for (auto tf : mesh->faces())
		//		{
		//			if (!tf.is_valid()) continue;
		//			id = 0;
		//			minAngle = 4.0;
		//			int i = 0;
		//			for (auto tfh : mesh->fh_range(tf))
		//			{
		//				double angle = mesh->calc_sector_angle(tfh);
		//				if (angle < minAngle)
		//				{
		//					minAngle = angle;
		//					id = i;
		//				}
		//				++i;
		//			}
		//			if (minAngle < lowerAngleBound)
		//			{
		//				auto h_iter = mesh->fh_begin(tf);
		//				for (i = 0; i < id; ++i)
		//				{
		//					++h_iter;
		//				}
		//				auto th = mesh->prev_halfedge_handle(*h_iter);
		//				auto te = mesh->edge_handle(th);
		//				if (mesh->is_boundary(th.from()))
		//				{
		//					continue;
		//				}
		//				if (mesh->is_collapse_ok(th) && mesh->calc_edge_length(te) < threshold)
		//				{
		//					if (mesh->data(th.from()).get_vertflag())
		//					{
		//						mesh->data(th.to()).set_vertflag(true);
		//					}
		//					/*if (mesh->data(th.prev().edge()).get_edgeflag())
		//					{
		//						mesh->data(h_iter->edge()).set_edgeflag(true);
		//					}*/
		//					mesh->data(h_iter->edge()).flag1 = mesh->data(th.prev().edge()).flag1;
		//					mesh->data(h_iter->edge()).flag2 = mesh->data(th.prev().edge()).flag2;
		//
		//					mesh->collapse(th);
		//				}
		//				else
		//				{
		//					if (mesh->calc_sector_angle(th) < mesh->calc_sector_angle(th.prev()))
		//					{
		//						th = th.prev();
		//					}
		//					auto flagvert = th.to();
		//					auto ph = th.prev();
		//					auto ne = th.next().edge();
		//					te = th.edge();
		//					auto pe = ph.edge();
		//
		//					if (!mesh->is_flip_ok(pe))
		//						continue;
		//					//if (mesh->data(pe).get_edgeflag())
		//					//{
		//					//	mesh->data(flagvert).set_vertflag(true);
		//					//	mesh->data(te).set_edgeflag(true);
		//					//	mesh->data(ne).set_edgeflag(true);
		//					//	mesh->data(pe).set_edgeflag(false);
		//					//	/*if (!mesh->data(flagvert).get_vertflag())
		//					//	{
		//					//		mesh->set_point(flagvert, mesh->calc_edge_midpoint(pe));
		//					//	}*/
		//					//}
		//					if (mesh->data(pe).flag1)
		//					{
		//						mesh->data(flagvert).set_vertflag(true);
		//						mesh->data(te).flag1 = true;
		//						mesh->data(ne).flag1 = true;
		//						mesh->data(pe).flag1 = false;
		//					}
		//					if (mesh->data(pe).flag2)
		//					{
		//						mesh->data(flagvert).set_vertflag(true);
		//						mesh->data(te).flag2 = true;
		//						mesh->data(ne).flag2 = true;
		//						mesh->data(pe).flag2 = false;
		//					}
		//					mesh->flip(pe);
		//				}
		//			}
		//		}
		//#endif


#ifdef printRemeshingInfo
		tr.out("process small angle:");
#endif
	}

	void TriangleMeshRemeshing::tangential_relaxation()
	{
		std::vector<OpenMesh::Vec3d> normal;
		normal.reserve(mesh->n_faces());
		std::vector<double> area;
		area.reserve(mesh->n_faces());
		for (auto &tf : mesh->faces())
		{
			normal.push_back(mesh->calc_face_normal(tf));
			area.push_back(mesh->calc_face_area(tf));
		}
		auto updateVFRingNormals = [&](SOV &tv)
		{
			for (auto &tvf : mesh->vf_range(tv))
			{
				normal[tvf.idx()] = mesh->calc_face_normal(tvf);
			}
		};
		auto updateVFRingAreas = [&](SOV &tv)
		{
			for (auto &tvf : mesh->vf_range(tv))
			{
				area[tvf.idx()] = mesh->calc_face_area(tvf);
			}
		};

		for (auto tv : mesh->vertices()) {
			if (mesh->data(tv).get_vertflag() || mesh->is_boundary(tv))
				continue;
			mesh->set_point(tv, aabbtree->closest_point(GravityPos(tv, normal, area)));
			updateVFRingAreas(tv);
			updateVFRingNormals(tv);
		}
#ifdef printRemeshingInfo
		tr.out("project:");
#endif
	}

	void TriangleMeshRemeshing::globalProject()
	{
		//vector<unsigned> &triangle_surface_index = globalmodel.triangle_surface_index;
		//vector<vector<unsigned>> vertex_surface_index(globalmodel.faceshape.size());
		//for (auto tv : mesh_->vertices()) {
		//	OpenMesh::Vec3d p = mesh_->point(tv);
		//	CGAL_AABB_Tree::Point_and_primitive_id point_primitive = AABB_tree->closest_point_and_primitive(CGAL_double_3_Point(p[0], p[1], p[2]));
		//	CGAL_Triangle_Iterator it = point_primitive.second;
		//	unsigned face_id = std::distance(triangle_vectors.begin(), it);
		//	vertex_surface_index[triangle_surface_index[face_id]].push_back(tv.idx());
		//}
		///*print(vertex_surface_index[67].size());
		//for (int i = 0; i < vertex_surface_index[67].size(); i++)
		//	std::cout << vertex_surface_index[67][i] << ",";*/
		//MeshProjectToSurface(mesh_, vertex_surface_index, &globalmodel);
		vector<unsigned> &triangle_surface_index = globalmodel.triangle_surface_index;
		vector<vector<unsigned>> vertex_surface_index(globalmodel.faceshape.size());
		for (auto &tv : mesh->vertices())
		{
			auto cdf = aabbtree->closest_distance_and_face_handle(mesh->point(tv));
			vertex_surface_index[triangle_surface_index[cdf.second.idx()]].push_back(tv.idx());

		}
		//MeshProjectToSurface(mesh, vertex_surface_index, &globalmodel);
	}

	O3d TriangleMeshRemeshing::GravityPos(const OV &v)
	{
		OpenMesh::Vec3d posSum(0, 0, 0);
		for (auto &tvoh : mesh->voh_range(v))
		{
			posSum += mesh->point(tvoh.to());
		}
		posSum /= mesh->valence(v);

#if 0
		return posSum;
#else
		OpenMesh::Vec3d n(0, 0, 0);
		for (auto &tvf : mesh->vf_range(v))
		{
			n += mesh->calc_face_area(tvf)*mesh->calc_face_normal(tvf);
		}
		if (n.sqrnorm() > epsilonerror)
			n.normalize();
		return posSum - (posSum - mesh->point(v)).dot(n)*n;
#endif
	}

	O3d TriangleMeshRemeshing::GravityPos(const OV &v, const std::vector<OpenMesh::Vec3d> &normal, const std::vector<double> &area)
	{
		OpenMesh::Vec3d posSum(0, 0, 0);
		for (auto &tvoh : mesh->voh_range(v))
		{
			posSum += mesh->point(tvoh.to());
		}
		posSum /= mesh->valence(v);

#if 0
		return posSum;
#else
		OpenMesh::Vec3d n(0, 0, 0);
		for (auto &tvf : mesh->vf_range(v))
		{
			n += area[tvf.idx()] * normal[tvf.idx()];
		}
		if (n.sqrnorm() > epsilonerror)
			n.normalize();
		return posSum - (posSum - mesh->point(v)).dot(n)*n;
#endif
	}

	void TriangleMeshRemeshing::initTargetLength()
	{
		//double kMin = 1.0e-2;
		//double kMax = 1.0e1;
		////expected_length *= 2;
		//double error = (6 / kMin + std::sqrt(36 / (kMin*kMin) - 12 * expected_length*expected_length)) / 6.0;
		//double lowerLengthBound = std::sqrt(3 * error*(2 / kMax - error));
		////double kMin = sqrt(3) / expected_length;

		double error = 0.01*expected_length;
		double minL = 0.1*expected_length;
		double maxL = expected_length;

		double kMin = 6 * error / (maxL*maxL + 3 * error*error);
		double kMax = 6 * error / (minL*minL + 3 * error*error);


		for (auto &tv : mesh->vertices())
		{
			if (tv.valence() <= 2 && !tv.is_boundary() || isnan(mesh->data(tv).GaussCurvature))
			{
				int p = 0;
			}
			int p = 0;
			/*mesh->data(tv).set_targetlength(mesh->data(tv).GaussCurvature > 1.0e-4 ?
				expected_length * 2 / (6 + Log10(mesh->data(tv).GaussCurvature)) : expected_length);*/
			//mesh->data(tv).set_targetlength()
			double gc = mesh->data(tv).GaussCurvature;
			if (gc <= kMin * 1.001)
				mesh->data(tv).set_targetlength(maxL);
			else if (gc > kMax)
				mesh->data(tv).set_targetlength(minL);
			else
				mesh->data(tv).set_targetlength(std::sqrt(3 * error*(2 / gc - error)));
			//dprint(gc, mesh->data(tv).get_targetlength());
		}
		int p = 0;
	}

#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
	TriangleMeshRemeshing::TriangleMeshRemeshing(PolyMesh *mesh_, double target_length)
		:expected_length(target_length), polymeshInput(true)
	{
		initMeshStatusAndNormal(*mesh_);
		polymesh = new Mesh();
		polymesh->reserve(mesh_->n_vertices(), mesh_->n_edges(), mesh_->n_faces());

		enum vertexType { tri, poly, mixed };
		OpenMesh::VPropHandleT<int> idMap;
		OpenMesh::VPropHandleT<vertexType> typeMap;
		mesh_->add_property(idMap);
		mesh_->add_property(typeMap);
		int i = 0, j = 0, k = 0;
		{
			bool allTri = true, allPoly = true;
			int triangleNum = 0;
			for (auto &tv : mesh_->vertices())
			{
				allTri = allPoly = true;
				for (auto &tvf : mesh_->vf_range(tv))
				{
					if (tvf.valence() == 3)
					{
						allPoly = false;
						++triangleNum;
					}
					else
					{
						allTri = false;
					}
				}
				if (allTri)
				{
					mesh_->property(idMap, tv) = i++;
					mesh_->property(typeMap, tv) = tri;
				}
				else if (allPoly)
				{
					mesh_->property(idMap, tv) = j++;
					mesh_->property(typeMap, tv) = poly;
				}
				else
				{
					mesh_->property(idMap, tv) = k++;
					mesh_->property(typeMap, tv) = mixed;
				}
			}
			triangleNum /= 3;
			for (auto &tv : mesh_->vertices())
			{
				if (mesh_->property(typeMap, tv) != mixed)
				{
					mesh_->property(idMap, tv) += k;
				}
			}
			mesh = new TriMesh();
			mesh->reserve(i + k, i + k + triangleNum, triangleNum);//V+F-E=2-2g-b  =>  E~=V+F
		}

		for (auto &tv : mesh_->vertices())
		{
			if (mesh_->property(typeMap, tv) == mixed)
			{
				auto &av = mesh->add_vertex(mesh_->point(tv));
				auto &pav = polymesh->add_vertex(mesh_->point(tv));
				if (mesh_->data(tv).get_vertflag())
				{
					mesh->data(av).set_vertflag(true);
					polymesh->data(pav).set_vertflag(true);
				}
			}
		}
		for (auto &tv : mesh_->vertices())
		{
			switch (mesh_->property(typeMap, tv))
			{
			case tri:
			{
				auto &av = mesh->add_vertex(mesh_->point(tv));
				if (mesh_->data(tv).get_vertflag())
				{
					mesh->data(av).set_vertflag(true);
				}
				break;
			}
			case poly:
			{
				auto &av = polymesh->add_vertex(mesh_->point(tv));
				if (mesh_->data(tv).get_vertflag())
				{
					polymesh->data(av).set_vertflag(true);
				}
				break;
			}
			default:
				break;
			}
		}

		int id[4];
		int ii = 0;
		for (auto &tf : mesh_->faces())
		{
			ii = 0;
			if (tf.valence() == 3)
			{
				for (auto &tfv : mesh_->fv_range(tf))
				{
					id[ii++] = mesh_->property(idMap, tfv);
				}
				mesh->add_face(mesh->vertex_handle(id[0]), mesh->vertex_handle(id[1]), mesh->vertex_handle(id[2]));
			}
			else
			{
				for (auto &tfv : mesh_->fv_range(tf))
				{
					id[ii++] = mesh_->property(idMap, tfv);
				}
				polymesh->add_face(polymesh->vertex_handle(id[0]), polymesh->vertex_handle(id[1]),
					polymesh->vertex_handle(id[2]), polymesh->vertex_handle(id[3]));
			}
		}
		for (auto &te : mesh_->edges())
		{
			vertexType vt0 = mesh_->property(typeMap, te.v0());
			vertexType vt1 = mesh_->property(typeMap, te.v1());

			if (mesh_->data(te).flag1)
			{
				if (vt0 == tri || vt1 == tri)
				{
					mesh->data(mesh->edge_handle(mesh->find_halfedge(mesh->vertex_handle(mesh_->property(idMap, te.v0())),
						mesh->vertex_handle(mesh_->property(idMap, te.v1()))))).flag1 = true;
				}
				else if (vt0 == poly || vt1 == poly)
				{
					polymesh->data(polymesh->edge_handle(polymesh->find_halfedge(polymesh->vertex_handle(mesh_->property(idMap, te.v0())),
						polymesh->vertex_handle(mesh_->property(idMap, te.v1()))))).flag1 = true;
				}
				else
				{
					mesh->data(mesh->edge_handle(mesh->find_halfedge(mesh->vertex_handle(mesh_->property(idMap, te.v0())),
						mesh->vertex_handle(mesh_->property(idMap, te.v1()))))).flag1 = true;
					polymesh->data(polymesh->edge_handle(polymesh->find_halfedge(polymesh->vertex_handle(mesh_->property(idMap, te.v0())),
						polymesh->vertex_handle(mesh_->property(idMap, te.v1()))))).flag1 = true;
				}
			}
			if (mesh_->data(te).flag2)
			{
				if (vt0 == tri || vt1 == tri)
				{
					mesh->data(mesh->edge_handle(mesh->find_halfedge(mesh->vertex_handle(mesh_->property(idMap, te.v0())),
						mesh->vertex_handle(mesh_->property(idMap, te.v1()))))).flag2 = true;
				}
				else if (vt0 == poly || vt1 == poly)
				{
					polymesh->data(polymesh->edge_handle(polymesh->find_halfedge(polymesh->vertex_handle(mesh_->property(idMap, te.v0())),
						polymesh->vertex_handle(mesh_->property(idMap, te.v1()))))).flag2 = true;
				}
				else
				{
					mesh->data(mesh->edge_handle(mesh->find_halfedge(mesh->vertex_handle(mesh_->property(idMap, te.v0())),
						mesh->vertex_handle(mesh_->property(idMap, te.v1()))))).flag2 = true;
					polymesh->data(polymesh->edge_handle(polymesh->find_halfedge(polymesh->vertex_handle(mesh_->property(idMap, te.v0())),
						polymesh->vertex_handle(mesh_->property(idMap, te.v1()))))).flag2 = true;
				}
		}
		}
		boundaryNum = k;
		mesh_->remove_property(idMap);
		mesh_->remove_property(typeMap);

		*mesh_ = *polymesh;
		delete polymesh;
		polymesh = mesh_;
		initMeshStatusAndNormal(*mesh);

		//*polymesh = Mesh(*mesh);

		if (expected_length <= 0)
		{
			expected_length = meshAverageLength(*mesh);
		}
		aabbtree = new ClosestPointSearch::AABBTree(*mesh);
	}

	void TriangleMeshRemeshing::assembleMesh()
	{
#if 1
		int nv = polymesh->n_vertices();
		auto vItr = mesh->vertices_begin();
		for (int i = 0; i < boundaryNum; ++i, ++vItr);
		for (; vItr != mesh->vertices_end(); ++vItr)
		{
			auto &av = polymesh->add_vertex(mesh->point(*vItr));
			if (mesh->data(*vItr).get_vertflag())
			{
				polymesh->data(av).set_vertflag(true);
			}
		}

		int id[3];
		int ii = 0;
		for (auto &tf : mesh->faces())
		{
			ii = 0;
			for (auto &tfv : mesh->fv_range(tf))
			{
				id[ii++] = tfv.idx() < boundaryNum ? tfv.idx() : (tfv.idx() + nv - boundaryNum);
			}
			polymesh->add_face(polymesh->vertex_handle(id[0]), polymesh->vertex_handle(id[1]), polymesh->vertex_handle(id[2]));
		}
		for (auto &te : mesh->edges())
		{
			auto ph = polymesh->edge_handle(polymesh->find_halfedge(
				polymesh->vertex_handle(te.v0().idx() < boundaryNum ? te.v0().idx() : te.v0().idx() + nv - boundaryNum),
				polymesh->vertex_handle(te.v1().idx() < boundaryNum ? te.v1().idx() : te.v1().idx() + nv - boundaryNum)));
			polymesh->data(ph).flag1 = mesh->data(te).flag1;
			polymesh->data(ph).flag2 = mesh->data(te).flag2;
		}
#else
		*polymesh = Mesh(*mesh);
#endif
		initMeshStatusAndNormal(*polymesh);

		delete mesh;
		mesh = nullptr;
	}
#endif

	/*bool TriangleMeshRemeshing::split_one_edge(Mesh::EdgeHandle& eh, OpenMesh::Vec3d& p)
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
	}*/
}
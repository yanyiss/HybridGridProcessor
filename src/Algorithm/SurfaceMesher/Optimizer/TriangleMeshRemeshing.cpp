#include "TriangleMeshRemeshing.h"

namespace CADMesher
{
	TriangleMeshRemeshing::TriangleMeshRemeshing(TriMesh* mesh_, double target_length)
		:mesh(mesh_), expected_length(target_length)
	{
		if (expected_length <= 0)
		{
			expected_length = meshAverageLength(*mesh);
		}
		aabbtree = globalmodel.init_trimesh_tree;
		boundaryNum = 0;
	}

	void TriangleMeshRemeshing::run()
	{
		//tmqh用来监控网格的质量
		TriMeshQualityHelper tmqh(mesh);
		tmqh.print();
		dprint("\nMesh vertices number in initialization:", mesh->n_vertices());
		//tr用来监控算法时间
		tr.refresh();

		//remesh算法分为三个阶段：
		//1.提高网格平均质量
		//2.消除大部分小角，并且继续优化网格平均质量
		//3.消除所有小角

		initTargetLength();//根据曲率设置目标边长
		//提高网格平均质量
		dprint("global remeshing");
		int itertimes = 0;
		for (; itertimes < 5; ++itertimes)
		{
			tr.mark();
			adjustTargetLength();
			split();//将长边分割为短边
			collapse();//去除、删除短边
			equalize_valence();//平衡度数（一个顶点的边数）
			tangential_relaxation();//把顶点移动到周围点的中心
			dprint("mesh vertices number:", mesh->n_vertices());
#ifdef printRemeshingInfo
			tmqh.update();
			tmqh.print();
			tr.pastMark("the " + std::to_string(itertimes + 1) + "-th iteration time:");
			dprint();
#endif
		}
		tr.mark();
		globalProject();//点到曲面的投影
		tr.pastMark("project to the origin surface time:");

		//消除大部分小角，并且继续优化网格质量
		dprint("remeshing while eliminating small angle");
		itertimes = 0;
		for (; itertimes < 5; ++itertimes)
		{
			tr.mark();
			if (processFeatureConstraintAngle() < 20)//处理最小角
			{
				tmqh.update();
				if(tmqh.getAvgQuality() > 0.88)
					break;
			}

			adjustTargetLength();
			split();
			collapse(true);
			equalize_valence(true);
			tangential_relaxation();//把顶点移动到周围点的中心

			dprint("mesh vertices number:", mesh->n_vertices());
#ifdef printRemeshingInfo
			tmqh.update();
			tmqh.print();
			tr.pastMark("the " + std::to_string(itertimes + 1) + "-th iteration time:");
			dprint();
#endif
		}

		tr.mark();
		globalProject();
		tr.pastMark("project to the origin surface time:");

		//return;
		//消除所有小角
		dprint("aggressively eliminating small angle");
		itertimes = 0;
		while (processFeatureConstraintAngle(true) && ++itertimes < 5);

		tmqh.update();
		tmqh.print();
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
		if (polymeshInput)
		{
			assembleMesh();//将非三角形装配入网格
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
	}

	void TriangleMeshRemeshing::collapse(bool ifEnhanced)
	{
		int collapsenumber = 0;
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
			if (mesh->data(the.edge()).get_edgeflag())
			{
				if (x >= 0.666 * min_of_t0_t1) continue;
				int count = 0;
				OV v = fromvert;
				for (auto &tvoh : mesh->voh_range(fromvert)) 
				{
					if (!mesh->data(tvoh.edge()).get_edgeflag()) continue;
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
#else
			if (mesh->data(the.edge()).flag1 || mesh->data(the.edge()).flag2)
				continue;
#endif
			if (mesh->is_boundary(the))
				continue;
			O3d pos = mesh->data(tovert).get_vertflag() ? mesh->point(tovert) : mesh->calc_centroid(the);

			if (ifEnhanced)
			{
				//若collapse后出现长边，则退出
				for (OV thev : mesh->vv_range(fromvert))
					if ((pos - mesh->point(thev)).norm() > 1.33 * min_of_t0_t1)
						goto goto20210523;
				for (OV thev : mesh->vv_range(tovert))
					if ((pos - mesh->point(thev)).norm() > 1.33 * min_of_t0_t1)
						goto goto20210523;

				//若collapse后出现小角，则退出
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

			//若collapse后出现翻折三角形，则退出
			/*bool ifDetect = false;
			for (auto tvv : mesh->vv_range(the.from()))
				if (mesh->data(tvv).get_vertflag())
					ifDetect = true;
			if (ifDetect)*/
			{
				auto& fromPos = mesh->point(the.from());
				//auto& pos = mesh->point(the.to());
				OH flipH = the.prev().opp();
				int endId = the.opp().next().idx();
				while (flipH.idx() != endId)
				{
					auto& p0 = mesh->point(mesh->to_vertex_handle(flipH));
					auto& p1 = mesh->point(mesh->to_vertex_handle(mesh->next_halfedge_handle(flipH)));
					if ((p0 - fromPos).cross(p1 - fromPos).dot((p0 - pos).cross(p1 - pos)) < 0)
						goto goto20210523;
					flipH = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(flipH));
				}
				//若tovert无标记，则collapse后的点坐标pos是边的中点，此时tovert周围可能出现翻折三角形
				if (!mesh->data(tovert).get_vertflag())
				{
					auto& toPos = mesh->point(the.to()); 
					flipH = the.opp().prev();
					endId = the.next().opp().idx();
					while (flipH.idx() != endId)
					{
						auto& p0 = mesh->point(mesh->from_vertex_handle(flipH));
						auto& p1 = mesh->point(mesh->from_vertex_handle(mesh->prev_halfedge_handle(flipH)));
						if ((p0 - toPos).cross(p1 - toPos).dot((p0 - pos).cross(p1 - pos)) < 0)
						{
							//pos = toPos;
							//break;
							goto goto20210523;
						}
						flipH = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(flipH));
					}
				}
			}

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
			bool if_equal = false;
			//判断flip是否会优化度数分布
			/*int before = (v0 - 6)*(v0 - 6) + (v1 - 6)*(v1 - 6) + (u0 - 6)*(u0 - 6) + (u1 - 6)*(u1 - 6);
			int later = (v0 - 7)*(v0 - 7) + (v1 - 7)*(v1 - 7) + (u0 - 5)*(u0 - 5) + (u1 - 5)*(u1 - 5);
			if (before < later)
			{
				continue;
			}*/
			if (fabs(v0 - 6) + fabs(v1 - 6) + fabs(u0 - 6) + fabs(u1 - 6) < fabs(v0 - 7) + fabs(v1 - 7) + fabs(u0 - 5) + fabs(u1 - 5))
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
			auto opt = [&](OpenMesh::HalfedgeHandle he) {
				return fabs(mesh->calc_sector_angle(he) - PI / 3.0) + fabs(mesh->calc_sector_angle(mesh->prev_halfedge_handle(he)) - PI / 3.0) +
					fabs(mesh->calc_sector_angle(mesh->next_halfedge_handle(he)) - PI / 3.0);
			};
			double opt_before = opt(h0) + opt(h1);
			mesh->flip(te);
			if (opt_before < opt(te.h0()) + opt(te.h1())) 
				mesh->flip(te);
			++equalizenumber;
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
		//对目标边长做光滑处理
		std::vector<double> tl; tl.reserve(mesh->n_vertices());
		for (auto& tv : mesh->vertices())
		{
			double l = 0;
			if (mesh->is_boundary(tv))
			//if (mesh->data(tv).get_vertflag())
			{
				tl.push_back(mesh->data(tv).get_targetlength());
			}
			else
			{
				for (auto& tvv : mesh->vv_range(tv))
				{
					l += mesh->data(tvv).get_targetlength();
				}
				tl.push_back(0.5 * (l / mesh->valence(tv) + mesh->data(tv).get_targetlength()));
			}
		}
		for (auto &tv : mesh->vertices())
			mesh->data(tv).set_targetlength(tl[tv.idx()]);
#ifdef printRemeshingInfo
		tr.out("adjust target length:");
#endif
	}

	int TriangleMeshRemeshing::processFeatureConstraintAngle(bool ifEnhanced)
	{
		//由于在之前remesh的过程中，我们始终保持特征线不动，因此一些小角难以消除。该函数在最大程度上保持特征线的同时，可以消除小角。
		//但由于网格边界始终被锁定，以及网格尺寸在局部不合适，小角仍可能无法完全消除
		auto setAvgLength = [&](OV &v)
		{
			double l = 0;
			for (auto &tvv : mesh->vv_range(v))
			{
				l += mesh->data(tvv).get_targetlength();
			}
			mesh->data(v).set_targetlength(0.5 * (l / mesh->valence(v) + mesh->data(v).get_targetlength()));
		};
		//三角形中有小角的情况分为两种，即有两个小角和只有一个小角，分别进行处理
		int processNumber = 0;
		int findNumber = 0;
		//处理有两个角小于阈值的情况
		/*
		如下图的三角形ABC，这里A,B皆为要消除的小角
		消除小角的方法是翻转(flip)AB或者分割(split)AB
		               C
		               *
		            *       *
		         *               *
		      *                      *
	      B *  *  *  *  *  *  *  *  *  * A
		*/
		for (auto &th : mesh->halfedges())
		{
			if (th.edge().is_boundary() || !th.is_valid() || mesh->calc_sector_angle(th) > lowerAngleBound)
				continue;
			++findNumber;
			//calc_sector_angle(h)是计算三角形在顶点h.to()上的角的弧度值
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
			if (isTwoSmallAngles && !CB.next().edge().is_boundary())
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
				//if (edgeFlag(CB.edge()) || edgeFlag(CB.next().edge()) || edgeFlag(CB.prev().edge()))
				if (mesh->data(CB.edge()).get_edgeflag() || mesh->data(CB.next().edge()).get_edgeflag() 
					|| mesh->data(CB.prev().edge()).get_edgeflag())
				{
					mesh->data(CB.from()).set_vertflag(true);
				}
				if (mesh->is_flip_ok(BA))
					mesh->flip(BA);
				else
				{
					split_one_edge(CB.next().edge(), true);
					OV newvert = mesh->vertex_handle(mesh->n_vertices() - 1);
					/*auto query = GravityPos(newvert);
					closestPoint(query, aabbtree);
					mesh->set_point(newvert, query);*/
					//mesh->set_point(newvert, aabbtree->closest_point(GravityPos(newvert)));
					//auto cp = aabbtree->closest_point(ParallelTools::Vec3d(GravityPos(newvert).data()));
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

		//处理有一个角小于阈值的情况
		/*
		如下图的三角形ABC，这里B是要消除的小角
		消除小角的基本方法是坍缩(collapse)AC
		若AC不可坍缩，这一般意味着B的度数为3，可以做特殊处理
		 A
		  *
		  *       *
		  *              *
		  *                     *
		  *                            *
		  *                                   * B
		  *                      *
		  *           *
		  *
		 C
		*/
		for (auto &th : mesh->halfedges())
		{
			if (th.edge().is_boundary() || !th.is_valid() || mesh->calc_sector_angle(th) > lowerAngleBound)
				continue;
			SOH CB = th;
			if (!mesh->is_collapse_ok(th.prev()))
			{
				if (ifEnhanced)
				{
					if (th.to().valence() == 3)
						CB = th.next().opp().prev();
					else if (th.prev().opp().prev().from().valence() == 3)
						CB = th.prev().opp().prev().opp().prev();
					if (mesh->is_collapse_ok(CB.prev()) && !mesh->is_boundary(CB.prev().edge()) && !CB.prev().from().is_boundary())
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
			if (angle0 > lowerAngleBound && angle1 > lowerAngleBound && !CB.prev().edge().is_boundary() && !CB.prev().from().is_boundary())
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
				//if (edgeFlag(CB.edge()) || edgeFlag(CB.next().edge()) || edgeFlag(CB.prev().edge()) || edgeFlag(temp.edge()) || edgeFlag(temp.prev().edge()))
				if (mesh->data(CB.edge()).get_edgeflag() || mesh->data(CB.next().edge()).get_edgeflag() || 
					mesh->data(CB.prev().edge()).get_edgeflag() || mesh->data(temp.edge()).get_edgeflag() ||
					mesh->data(temp.prev().edge()).get_edgeflag())
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
		auto ifFlip = [&](SOV & tv)
		{
			for (auto& tvoh : mesh->voh_range(tv))
			{
				//if (tvoh.is_boundary() || tvoh.opp().is_boundary())
				//	continue;
				if (normal[tvoh.face().idx()].dot(normal[tvoh.opp().face().idx()]) < 0)
					return true;
			}
			return false;
		};

		for (auto tv : mesh->vertices()) {
			if (mesh->data(tv).get_vertflag() || mesh->is_boundary(tv))
				continue;
			O3d oldPos = mesh->point(tv);
			mesh->set_point(tv, aabbtree->closest_point(GravityPos(tv, normal, area)));
			updateVFRingAreas(tv);
			updateVFRingNormals(tv);
			//若移动点坐标导致三角形翻折，则退回操作
			if (ifFlip(tv))
			{
				mesh->set_point(tv, oldPos);
				updateVFRingAreas(tv);
				updateVFRingNormals(tv);
			}

#if 0
			O3d oldPos = mesh->point(tv);
			auto query = GravityPos(tv, normal, area);
			int fid = closestFaceAndPoint(query, aabbtree);
			//将点向初始网格投影得到面索引fid后，即可查找面fid所在的曲面索引，从而查找背景网格在该点的步长
			mesh->data(tv).set_targetlength(referenceModel->bgshape[referenceModel->triangle_surface_index[fid]].spacing(query.data()));
			//dprint(mesh->data(tv).get_targetlength());
			if (mesh->data(tv).get_vertflag())
				continue;
			mesh->set_point(tv, query);
			updateVFRingAreas(tv);
			updateVFRingNormals(tv);
			//若移动点坐标导致三角形翻折，则退回操作
			if (ifFlip(tv))
			{
				mesh->set_point(tv, oldPos);
				updateVFRingAreas(tv);
				updateVFRingNormals(tv);
			}
#endif
		}
#ifdef printRemeshingInfo
		tr.out("project:");
#endif
	}

	void TriangleMeshRemeshing::globalProject()
	{
#if 0
		vector<unsigned>& triangle_surface_index = referenceModel->triangle_surface_index;
		for (auto& tv : mesh->vertices())
		{
			if (mesh->data(tv).get_vertflag())
				continue;
			auto fid = closestFace(mesh->point(tv), aabbtree);
			mesh->set_point(tv, project_pnt_to_surface(triangle_surface_index[fid], mesh->point(tv)));
		}
#endif
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
		//已知每个顶点附近的最大法曲率，并且给定误差，根据参考材料中的公式可换算得到该顶点处的目标边长
		double error = 0.01*expected_length;
		double minL = 0.1*expected_length;
		double maxL = expected_length;

		double kMin = 6 * error / (maxL*maxL + 3 * error*error);
		double kMax = 6 * error / (minL*minL + 3 * error*error);

		for (auto &tv : mesh->vertices())
		{
			double gc = mesh->data(tv).GaussCurvature;//最大主曲率的绝对值
			if (gc <= kMin * 1.001)
				mesh->data(tv).set_targetlength(maxL);
			else if (gc > kMax)
				mesh->data(tv).set_targetlength(minL);
			else
				mesh->data(tv).set_targetlength(std::sqrt(3 * error*(2 / gc - error)));

			//mesh->data(tv).set_targetlength(0.001);
			//dprint(gc, mesh->data(tv).get_targetlength());
		}
	}

#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
	TriangleMeshRemeshing::TriangleMeshRemeshing(PolyMesh *mesh_, double target_length)
		:expected_length(target_length), polymeshInput(true)
	{
		//若输入的是三角形与四边形的混合网格，则现将四边形全部删除，优化完成后再重新加入四边形
		initMeshStatusAndNormal(*mesh_);
		polymesh = new Mesh();
		polymesh->reserve(mesh_->n_vertices(), mesh_->n_edges(), mesh_->n_faces());

		//标记顶点的类型，即周围全是三角形（tri），全是四边形（poly）以及两者都存在（mixed）
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
		//将三角形和四边形分别存入mesh和polymesh
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
		//更新特征标记
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
	}

	void TriangleMeshRemeshing::assembleMesh()
	{
		//将优化前删除的四边形重新装入原网格中
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
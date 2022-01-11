#pragma once

#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>
#include <numeric>
#include "..\src\Dependency\CPS\CPS_AABBTree.h"
#include "..\src\Algorithm\SurfaceMesher\Generator\basic_def.h"

namespace CADMesher
{
	class TriangleMeshRemeshing
	{
	public:
		TriangleMeshRemeshing(TriMesh *mesh_)
			:mesh(mesh_)
		{
			high = 4.0 / 3.0*expected_length;
			low = 4.0 / 5.0*expected_length;
			aabbtree = new ClosestPointSearch::AABBTree(*mesh);
			/*std::vector<CGAL_3_Segment> segment_vectors;
			for (auto te : mesh_->edges()) {
				if (!mesh_->data(te).get_edgeflag()) continue;
				O3d p0 = mesh_->point(te.v0());
				O3d p1 = mesh_->point(te.v1());
				segment_vectors.emplace_back(CGAL_double_3_Point(p0[0], p0[1], p0[2]), CGAL_double_3_Point(p1[0], p1[1], p1[2]));
			}
			AABB_Segment_tree = new CGAL_AABB_Segment_Tree(segment_vectors.begin(), segment_vectors.end());
			AABB_Segment_tree->accelerate_distance_queries();*/
		};
		~TriangleMeshRemeshing();

	public:
		void RemeshingMethod();

		//private:
		void split();
		void collapse();
		void equalize_valence();
		void adjustTargetLength();
		void processAngle();
		void tangential_relaxation();

		double minAngle();

	private:
		O3d GravityPos(const OV &v);

		bool split_one_edge(Mesh::EdgeHandle& eh, OpenMesh::Vec3d& p);

	private:
		double high;
		double low;

		double lowerAngleBound = 0.05;
		/*double beta_min = 7.0 * PI / 36.0;
		double beta_max = 17.0 * PI / 36.0;*/
		TriMesh *mesh = nullptr;
		ClosestPointSearch::AABBTree *aabbtree = nullptr;
		//CGAL_AABB_Segment_Tree* AABB_Segment_tree = nullptr;

	public:
		static double expected_length;

	};
}


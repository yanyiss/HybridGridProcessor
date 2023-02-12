#pragma once

#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>
#include <numeric>
#include "..\src\Algorithm\SurfaceMesher\Generator\basic_def.h"
//#define printRemeshingInfo

namespace CADMesher
{
	class TriangleMeshRemeshing
	{
	public:
		explicit TriangleMeshRemeshing(TriMesh *mesh_, double target_length = -1, bool proj = false);
		TriangleMeshRemeshing(const TriangleMeshRemeshing &tmr) = delete;
		~TriangleMeshRemeshing() { 
			//if (aabbtree) { delete aabbtree; aabbtree = nullptr; } 
		}

	public:
		void run();

	private:
		//main step
		void split();
		bool split_one_edge(const OpenMesh::SmartEdgeHandle &te, bool ifRelaxCondition = false);
		void collapse(bool ifEnhanced = false);
		void equalize_valence(bool ifEnhanced = false);
		void tangential_relaxation();
		//auxiliary step
		void adjustTargetLength();
		int processFeatureConstraintAngle(bool ifEnhanced = false);
		void globalProject();
		//geometry support
		void initTargetLength();
		O3d GravityPos(const OV &v);
		O3d GravityPos(const OV &v, const std::vector<OpenMesh::Vec3d> &normal, const std::vector<double> &area);


	private:
		double expected_length;
		bool if_proj = false;

		timeRecorder tr;

		double coerciveAngleBound = 0.06;
		double lowerAngleBound = 0.1;
		std::vector<double> initial_FaceTargetLength;

		TriMesh *mesh = nullptr;
		ClosestPointSearch::AABBTree *aabbtree = nullptr;

#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
	public:
		explicit TriangleMeshRemeshing(PolyMesh *mesh_, double target_length = -1, bool proj = false);//, polymeshInput(true)
	private:
		bool polymeshInput = false;
		int boundaryNum;
		PolyMesh *polymesh;
	private:
		void assembleMesh();
#endif
	};
}


#pragma once
/*
-- An algorithm designed for optimizing general quad mesh into those with orthogonal diagonals
-- Inspired from the paper: Checkerboard Patterns with Black Rectangles
-- By yanyisheshou at GCL in USTC
-- Email: duweiyou@mail.ustc.edu.cn
-- 2021/12/23
*/
#include <Eigen\Core>
#include "..\MeshViewer\MeshDefinition.h"
#include "..\Toolbox\dprint.h"
#include "..\CPS\CPS_AABBTree.h"

class CheckBoardGenerator
{
public:
	CheckBoardGenerator(Mesh &trimesh_) :trimesh(&trimesh_) {
		init();
	}
	~CheckBoardGenerator() {
		if (aabbtree) { delete aabbtree; aabbtree = nullptr; }
	}

public:
	typedef int ui;

	typedef Eigen::Matrix<double, 3, -1> m3xd;
	typedef Eigen::Matrix<ui, -1, 1> vxu;
	typedef Eigen::VectorXd vxd;
	typedef Eigen::Vector3d v3d;
	typedef Eigen::Matrix3d m3d;
	

	void init();
	void run();

	void getCheckBoard(m3xd &V, std::vector<vxu> &F);
	void getPrimal(m3xd &V, std::vector<vxu> &F);
	void getDual(m3xd &V, std::vector<vxu> &F);
	void getMesh(Mesh &m, m3xd &V, std::vector<vxu> &F);

	//private:
	struct triple {
		int x, y, z;
		triple(int x_,int y_,int z_):x(x_),y(y_),z(z_){}
	};

	ui num[3];                                         //number of vertices, face and edges of polymesh
	Mesh* trimesh = nullptr;
	ClosestPointSearch::AABBTree* aabbtree = nullptr;
	m3xd cof;                                           //coefficent of plane function
	vxd ct;                                             //constant term of place function

	m3xd v[2];                                          //vertex position
	std::vector<std::vector<triple>> adj[2];            //connectivity
	//OpenMesh::VPropHandleT<bool> diagonalMeshIndex;

	void updateMesh(Mesh& m);
	void printCurrentInfo();

	//void setDiagonalMeshIndex();//divide original mesh vertices into two diagonal mesh classes and mark with true & false

//protected:

};


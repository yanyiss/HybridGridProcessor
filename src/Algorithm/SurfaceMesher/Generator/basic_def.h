#pragma once
#define PI 3.1415926535897932
#define epsilonerror 1.1e-15

#include <TopoDS_Face.hxx>
#include <vector>
using std::vector;
struct ShapeFace
{
	int id;
	TopoDS_Face face;
	vector<vector<int>> wires;
	ShapeFace(int id_, TopoDS_Face face_)
	{
		id = id_;
		face = face_;
	}
};

#include <TopoDS_Edge.hxx>
#include <Eigen/Dense>
using namespace Eigen;
struct ShapeEdge
{
	int id;
	bool if_merged;
	bool if_splitted;
	TopoDS_Edge edge;
	int main_face;
	int secondary_face;
	int reversed_edge;
	int prev_edge;
	int next_reversed_edge;
	MatrixX2d parameters;
	int begin_id;
	int end_id;

	int next_edge;
	bool if_visited;
	ShapeEdge(int id_, TopoDS_Edge edge_)
	{
		id = id_;
		edge = edge_;
		if_merged = false;
		if_splitted = false;
		next_reversed_edge = -1;
		main_face = -1;
		secondary_face = -1;
		reversed_edge = -1;
		begin_id = -1;

		next_edge = -1;
		if_visited = false;
	}
};

#include "..\src\MeshViewer\MeshDefinition.h"
#include "..\src\Dependency\CPS\CPS_AABBTree.h"
struct GlobalGeometry {
	TopoDS_Shape aShape;
	vector<ShapeFace> faceshape;
	vector<ShapeEdge> edgeshape;
	//vector < particle_state_t, Eigen::aligned_allocator<particle_state_t>
	//vector<ShapeEdge, Eigen::aligned_allocator<ShapeEdge>> faceshape;
	//vector<ShapeFace, Eigen::aligned_allocator<ShapeFace>> edgeshape;
	vector<unsigned> triangle_surface_index;
	TriMesh initial_trimesh;
	TriMesh isotropic_trimesh;
	TriMesh anisotropic_trimesh;
	ClosestPointSearch::AABBTree* init_trimesh_tree = nullptr;
	ClosestPointSearch::AABBTree** init_surfacemesh_tree = nullptr;

	PolyMesh initial_polymesh;
	GlobalGeometry() {}
	void clear() {
		faceshape.clear();
		edgeshape.clear();
		triangle_surface_index.clear();
		initial_trimesh.clear();
		isotropic_trimesh.clear();
		anisotropic_trimesh.clear();

		if(init_trimesh_tree) 
			init_trimesh_tree->clear();
		if (init_surfacemesh_tree)
		{
			delete[] init_surfacemesh_tree;
			init_surfacemesh_tree = nullptr;
		}
		initial_polymesh.clear();
	}
};

extern GlobalGeometry globalmodel;

void MeshProjectToSurface(Mesh* mesh, vector<vector<unsigned>> &vertex_surface_index, GlobalGeometry* model);

#include "..\src\Toolbox\dprinter\dprint.h"

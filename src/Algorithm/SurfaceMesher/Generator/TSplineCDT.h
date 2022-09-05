#pragma once
#include <iostream>
#include <fstream>
#include "DelaunayTriangulation.h"

namespace PDF3D
{
namespace CDT
{

class TSplineCDT
{
public:
	TSplineCDT() {}
	~TSplineCDT() {}

public:
	
	void delaunayTriangle(const std::vector<DT_Point>& Fixed_Vertices);

	
	void constraintedDelaunayTriangle(const std::vector<DT_Point>& Fixed_Vertices, const std::vector<DT_Edge>& Fixed_Edge);


	
	void constraintedDelaunayTriangle_with_inner_boundary(const std::vector<DT_Point>& Fixed_Vertices, const std::vector<DT_Edge>& Fixed_Edge);


	//output the vertices_ and triangles_ to V and F which represent a mesh.
	void output2io(std::vector<double>& V, std::vector<int>& F);

public:

	void load_info(const std::string& _filename_1, const std::string& _filename_2, std::vector<DT_Point>& Fixed_Vertices, std::vector<DT_Edge>& Fixed_Edge);

	void print_info(const std::string& _filename_1, const std::string& _filename_2);

	void print_edge(const std::vector<DT_Edge>& Fixed_Edge, const std::string& _filename_);

	void print_vertex(std::vector<DT_Point>& Fixed_Vertices, const std::string& _filename_);

	void addPoint(const std::vector<DT_Edge>& Fixed_Edge, std::vector<bool>& vertices_flag, int v_index);

	void addEdge(const std::vector<DT_Edge>& Fixed_Edge, std::vector<bool>& vertices_flag, DT_Edge target_edge);


private:
	double sign(const DT_Point& p, const DT_Point& p0, const DT_Point& p1);

	bool isPointInTriangle(const DT_Point& pt, const DT_Point& v1, const DT_Point& v2, const DT_Point& v3);

	bool isPointOnFixEdge(const DT_Point& P1, const DT_Point& P2, const DT_Point& tar_p);

	bool isTwoVertexInFixEdge(int v_index1, int v_index2);

	bool isTriangleContainFixEdge(const DT_Edge& fixedge);

	int find_opposite_target(int target, int neighbor_1, int neighbor_2, int& face_index);

	bool judge_Circumcircle_1(const DT_Point& vA, const DT_Point& vB, const DT_Point& vC, const DT_Point& vD);

	bool judge_Circumcircle_2(const DT_Point& vA, const DT_Point& vB, const DT_Point& vC, const DT_Point& vD);

	bool isTowSegmentsInsersect(const DT_Point& e00, const DT_Point& e01, const DT_Point& e10, const DT_Point& e11);

	DT_Point calcIntersection(const DT_Point& e00, const DT_Point& e01, const DT_Point& e10, const DT_Point& e11);

	bool judgeCrossingInDTMesh(const DT_Edge& Edge_1, const DT_Edge& Edge_2);

private:
	void splitTriangle(int v_index, int face_index);

	void swapTriangle(int target, int neighbor_1, int neighbor_2, int opposite, int f_index1, int f_index2);

	void addTriangleInOrder(int vertex_A, int vertex_B, int vertex_C);

	void delete_face(int face_index);

private:
	int find_TruncateFace_StartingA(const DT_Edge& target_edge);

	int opposedTriangle(int face_index, int v_index, int& v_seg, int& v_above, int& v_below, double k_AB, int Target_A);

	void triangulatePseudopolygonDelaunay(std::vector<int> PU, const DT_Edge& edge);

private:
	void deleteAuxiliaryTriangle();

	void delete_Outer_Face(const std::vector<DT_Edge>& Fixed_Edge);

	void delete_Inner_Face(const std::vector<DT_Edge>& Fixed_Edge);

	bool isPointInRegion(const std::vector<DT_Edge>& boundary_edge, const DT_Point& pt);

	void findBoundaryFixedEdge(std::vector<DT_Edge> Fixed_Edge, std::vector<DT_Edge>& Fixed_Boundary, std::vector<DT_Edge>& inner_Fixed_Edge);

	void find_outer_and_inner_Boundary(std::vector<DT_Edge> Fixed_Edge, std::vector<DT_Edge>& outer_Fixed_Boundary, std::vector<std::vector<DT_Edge>>& inner_Fixed_Boundary_list);

	void splitCrossingFixedEdge(std::vector<DT_Edge>& fixed_edge);

	void delete_noninner_Fixed_Edge(std::vector<DT_Edge>& Fixed_Boundary, std::vector<DT_Edge>& inner_Fixed_Edge);

private:
	std::vector<DT_Point> vertices_;
	std::vector<DT_Triangle> triangle_;

	DT_Point bbmin_;
	DT_Point bbmax_;

	std::vector< std::vector<int> > vert_fix_edge_;

public:
	const std::vector<DT_Point>& getVertices() { return vertices_; }
	const std::vector<DT_Triangle>& getTriangle() { return triangle_; }

private:
	void computerBBox();

	void computerVertEdge(const std::vector<DT_Edge>& fixed_edge);

private:
	void orientation_reserve();
};

}// namespace CDT 
}// namespace PDF3D
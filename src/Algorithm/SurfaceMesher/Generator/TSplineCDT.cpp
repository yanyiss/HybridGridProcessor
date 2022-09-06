#include "TSplineCDT.h"

namespace PDF3D
{
namespace CDT
{
void TSplineCDT::delaunayTriangle(const std::vector<DT_Point>& Fixed_Vertices)
{
	vertices_ = Fixed_Vertices;
	std::vector<DT_Edge> fixed_edge;

	int num_vertices = vertices_.size();
	int num_fixed_edges = fixed_edge.size();

	bool setBBox = false;
	if (!setBBox)
	{
		computerBBox();
	}

	vertices_.push_back({ 2 * bbmin_[0] - bbmax_[0] - 0.01, 2 * bbmin_[1] - bbmax_[1] });
	vertices_.push_back({ 2 * bbmax_[0] - bbmin_[0] + 0.01, 2 * bbmin_[1] - bbmax_[1] });
	vertices_.push_back({ (bbmax_[0] + bbmin_[0]) / 2, 2 * bbmax_[1] - bbmin_[1] + 0.01 });

	fixed_edge.push_back({ num_vertices,num_vertices + 1 });
	fixed_edge.push_back({ num_vertices + 1,num_vertices + 2 });
	fixed_edge.push_back({ num_vertices,num_vertices + 2 });

	computerVertEdge(fixed_edge);

	triangle_.clear();
	triangle_.push_back(DT_Triangle(num_vertices, num_vertices + 1, num_vertices + 2));

	std::vector<bool> vertices_flag(num_vertices, false);

	for (int i = 0; i < num_vertices; i++)
	{
		addPoint(fixed_edge, vertices_flag, i);
	}

	deleteAuxiliaryTriangle();

	fixed_edge.clear();
	vertices_.pop_back(); vertices_.pop_back(); vertices_.pop_back();
}

void TSplineCDT::constraintedDelaunayTriangle(const std::vector<DT_Point>& Fixed_Vertices, const std::vector<DT_Edge>& Fixed_Edge)
{
	vertices_ = Fixed_Vertices;
	std::vector<DT_Edge> fixed_edge = Fixed_Edge;


	std::vector<DT_Edge> Fixed_Boundary;
	std::vector<DT_Edge> Inner_Fixed_Edge;
	findBoundaryFixedEdge(fixed_edge, Fixed_Boundary, Inner_Fixed_Edge);
	splitCrossingFixedEdge(Inner_Fixed_Edge);
	delete_noninner_Fixed_Edge(Fixed_Boundary, Inner_Fixed_Edge);

	int num_vertices = vertices_.size();
	int num_fixed_edges = fixed_edge.size();

	bool setBBox = false;
	if (!setBBox)
	{
		computerBBox();
	}

	vertices_.push_back({ 2 * bbmin_[0] - bbmax_[0], 2 * bbmin_[1] - bbmax_[1] });
	vertices_.push_back({ 2 * bbmax_[0] - bbmin_[0], 2 * bbmin_[1] - bbmax_[1] });
	vertices_.push_back({ (bbmax_[0] + bbmin_[0]) / 2, 2.5 * bbmax_[1] - bbmin_[1] +1 });

	fixed_edge.push_back({ num_vertices,num_vertices + 1 });
	fixed_edge.push_back({ num_vertices + 1,num_vertices + 2 });
	fixed_edge.push_back({ num_vertices,num_vertices + 2 });

	computerVertEdge(fixed_edge);

	triangle_.clear();
	triangle_.push_back(DT_Triangle(num_vertices, num_vertices + 1, num_vertices + 2));


	std::vector<bool> vertices_flag(num_vertices, false);


	for (int i = 0; i < num_fixed_edges; i++)
	{
		const DT_Edge& target_edge = fixed_edge[i];
		addEdge(fixed_edge, vertices_flag, target_edge);
	}


	for (int i = 0; i < num_vertices; i++)
	{
		addPoint(fixed_edge, vertices_flag, i);
	}

	deleteAuxiliaryTriangle();

	fixed_edge.pop_back(); fixed_edge.pop_back(); fixed_edge.pop_back();
	vertices_.pop_back(); vertices_.pop_back(); vertices_.pop_back();


	delete_Outer_Face(Fixed_Boundary);




}


void TSplineCDT::constraintedDelaunayTriangle_with_inner_boundary(const std::vector<DT_Point>& Fixed_Vertices, const std::vector<DT_Edge>& Fixed_Edge)
{
	vertices_ = Fixed_Vertices;
	std::vector<DT_Edge> fixed_edge = Fixed_Edge;


	std::vector<DT_Edge> outer_Fixed_Boundary; outer_Fixed_Boundary.clear();
	std::vector<std::vector<DT_Edge>> inner_Fixed_Boundary_list; inner_Fixed_Boundary_list.clear();

	find_outer_and_inner_Boundary(Fixed_Edge, outer_Fixed_Boundary, inner_Fixed_Boundary_list);


	for (unsigned i = 0; i < inner_Fixed_Boundary_list.size(); i++)
	{
		delete_noninner_Fixed_Edge(outer_Fixed_Boundary, inner_Fixed_Boundary_list[i]);
	}

	int num_vertices = vertices_.size();
	int num_fixed_edges = fixed_edge.size();

	bool setBBox = false;
	if (!setBBox)
	{
		computerBBox();
	}



	vertices_.push_back({ 2 * bbmin_[0] - bbmax_[0], 2 * bbmin_[1] - bbmax_[1] });
	vertices_.push_back({ 2 * bbmax_[0] - bbmin_[0], 2 * bbmin_[1] - bbmax_[1] });
	vertices_.push_back({ (bbmax_[0] + bbmin_[0]) / 2, 2.5 * bbmax_[1] - bbmin_[1]+0.1 });

	/*vertices_.push_back({ 2 * bbmin_[0] - bbmax_[0], 2 * bbmin_[1] - bbmax_[1] - 0.1 });
	vertices_.push_back({ bbmax_[0] + 100.0, 2 * bbmin_[1] - bbmax_[1] });
	vertices_.push_back({ 2 * bbmax_[0] - bbmin_[0], bbmax_[1] + 100.0 });*/

	fixed_edge.push_back({ num_vertices,num_vertices + 1 });
	fixed_edge.push_back({ num_vertices + 1,num_vertices + 2 });
	fixed_edge.push_back({ num_vertices,num_vertices + 2 });

	computerVertEdge(fixed_edge);

	triangle_.clear();
	triangle_.push_back(DT_Triangle(num_vertices, num_vertices + 1, num_vertices + 2));


	std::vector<bool> vertices_flag(num_vertices, false);


	for (int i = 0; i < num_fixed_edges; i++)
	{
		const DT_Edge& target_edge = fixed_edge[i];
		addEdge(fixed_edge, vertices_flag, target_edge);
	}


	for (int i = 0; i < num_vertices; i++)
	{
		addPoint(fixed_edge, vertices_flag, i);
	}

	deleteAuxiliaryTriangle();

	fixed_edge.pop_back(); fixed_edge.pop_back(); fixed_edge.pop_back();
	vertices_.pop_back(); vertices_.pop_back(); vertices_.pop_back();


	delete_Outer_Face(outer_Fixed_Boundary);


	for (unsigned i = 0; i < inner_Fixed_Boundary_list.size(); i++)
	{
		delete_Inner_Face(inner_Fixed_Boundary_list[i]);
	}


	orientation_reserve();
}

void TSplineCDT::addPoint(const std::vector<DT_Edge>& Fixed_Edge, std::vector<bool>& vertices_flag, int v_index)
{
	if (vertices_flag[v_index] == true)
		return;

	const DT_Point& tar_p = vertices_[v_index];

	int num_face = triangle_.size();
	int face_index = -1;
	for (int i = 0; i < num_face; i++)
	{
		const DT_Triangle& T = triangle_[i];
		const DT_Point& A = vertices_[T[0]];
		const DT_Point& B = vertices_[T[1]];
		const DT_Point& C = vertices_[T[2]];
		if (isPointInTriangle(tar_p, A, B, C))
		{
			face_index = i;
			break;
		}
	}

	if (face_index == -1)
	{
		for (int i = 0; i < num_face; i++)
		{
			const DT_Triangle& T = triangle_[i];
			const DT_Point& A = vertices_[T[0]];
			const DT_Point& B = vertices_[T[1]];
			const DT_Point& C = vertices_[T[2]];
			if (isPointInTriangle(tar_p, A, B, C))
			{
				face_index = i;
				break;
			}
		}
	}


	splitTriangle(v_index, face_index);

	vertices_flag[v_index] = true;

	std::vector<size_t> new_face;
	new_face.push_back(face_index);
	new_face.push_back(triangle_.size() - 2);
	new_face.push_back(triangle_.size() - 1);


	while (new_face.size() > 0)
	{
		int face_curr_index = new_face[new_face.size() - 1];
		const DT_Triangle& face_curr = triangle_[face_curr_index];


		int neighbor_1_index = face_curr[1];
		int neighbor_2_index = face_curr[2];
		if (face_curr[1] == v_index)
		{
			neighbor_1_index = face_curr[2];
			neighbor_2_index = face_curr[0];
		}
		else if (face_curr[2] == v_index)
		{
			neighbor_1_index = face_curr[0];
			neighbor_2_index = face_curr[1];
		}
		if (isTwoVertexInFixEdge(neighbor_1_index, neighbor_2_index))
		{
			new_face.pop_back();
			continue;
		}


		int face_target_index;
		int opposite_index = find_opposite_target(v_index, neighbor_1_index, neighbor_2_index, face_target_index);
		if (judge_Circumcircle_1(vertices_[v_index], vertices_[neighbor_1_index], vertices_[neighbor_2_index], vertices_[opposite_index]))
		{
			new_face.pop_back();
			continue;
		}
		else
		{
			swapTriangle(v_index, neighbor_1_index, neighbor_2_index, opposite_index, face_curr_index, face_target_index);
			new_face.pop_back();
			new_face.push_back(face_curr_index);
			new_face.push_back(face_target_index);
		}
	}
	return;
}

void TSplineCDT::addEdge(const std::vector<DT_Edge>& Fixed_Edge, std::vector<bool>& vertices_flag, DT_Edge target_edge)
{
	int target_A = target_edge[0];
	int target_B = target_edge[1];
	addPoint(Fixed_Edge, vertices_flag, target_A);
	addPoint(Fixed_Edge, vertices_flag, target_B);
	if (isTriangleContainFixEdge(target_edge))
	{
		return;
	}


	int Tri = find_TruncateFace_StartingA(target_edge);

	if (Tri < 0)
	{
		//std::cout << "error " << std::endl;
		return;
	}


	int v_curr = target_A;
	std::vector<int> PU, PL;

	int t_seg, v_seg, v_above, v_below;
	double k_AB = (vertices_[target_B][1] - vertices_[target_A][1]) / (vertices_[target_B][0] - vertices_[target_A][0]);
	while (1)
	{
		t_seg = opposedTriangle(Tri, v_curr, v_seg, v_above, v_below, k_AB, target_A);
		if (v_seg == target_B)
		{
			PU.push_back(v_above);
			PL.push_back(v_below);

			delete_face(Tri);

			if (t_seg > Tri)
			{
				t_seg = t_seg - 1;
			}
			delete_face(t_seg);
			break;
		}
		if ((vertices_[v_seg][1] - vertices_[target_A][1]) - k_AB * (vertices_[v_seg][0] - vertices_[target_A][0]) > 0)
		{
			PU.push_back(v_above);
			v_curr = v_above;
		}
		else
		{
			PL.push_back(v_below);
			v_curr = v_below;
		}

		delete_face(Tri);
		if (t_seg > Tri)
		{
			t_seg = t_seg - 1;
		}
		Tri = t_seg;
	}


	triangulatePseudopolygonDelaunay(PU, target_edge);
	triangulatePseudopolygonDelaunay(PL, target_edge);
}




double TSplineCDT::sign(const DT_Point& p, const DT_Point& p0, const DT_Point& p1)
{
	return ((p[0] - p0[0]) * (p1[1] - p0[1]))
		- ((p1[0] - p0[0]) * (p[1] - p0[1]));
}

bool TSplineCDT::isPointInTriangle(const DT_Point& pt, const DT_Point& v1, const DT_Point& v2, const DT_Point& v3)
{
	double d1 = sign(pt, v1, v2);
	double d2 = sign(pt, v2, v3);
	double d3 = sign(pt, v3, v1);

	bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
	bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

	return !(has_neg && has_pos);
}

bool TSplineCDT::isTwoVertexInFixEdge(int v_index1, int v_index2)
{
	std::vector<int>& v1_edge = vert_fix_edge_[v_index1];
	for (unsigned i = 0; i < v1_edge.size(); i++)
	{
		if (v1_edge[i] == v_index2)
			return true;
	}
	return false;
}

bool TSplineCDT::isTriangleContainFixEdge(const DT_Edge& fixedge)
{
	int num_face = triangle_.size();
	int Edge_A = fixedge[0];
	int Edge_B = fixedge[1];
	for (int i = 0; i < num_face; i++)
	{
		const DT_Triangle& face = triangle_[i];

		if (face[0] == Edge_A || face[1] == Edge_A || face[2] == Edge_A)
		{
			if (face[0] == Edge_B || face[1] == Edge_B || face[2] == Edge_B)
				return true;
		}
	}
	return false;
}


int TSplineCDT::find_opposite_target(int target, int neighbor_1, int neighbor_2, int& face_index)
{
	int num_face = triangle_.size();

	int A[2] = {};
	int face_index_set[2] = {};
	int index = 0;
	for (int i = 0; i < num_face; i++)
	{
		const DT_Triangle& face = triangle_[i];
		if (!face.isContainVertex(neighbor_1) || !face.isContainVertex(neighbor_2))
			continue;

		if (face.isContainVertex(target))
			continue;

		if ((neighbor_1 == face[0] && neighbor_2 == face[1]) || (neighbor_1 == face[1] && neighbor_2 == face[0]))
		{
			face_index = i;
			return face[2];
		}
		if ((neighbor_1 == face[1] && neighbor_2 == face[2]) || (neighbor_1 == face[2] && neighbor_2 == face[1]))
		{
			face_index = i;
			return face[0];
		}
		if ((neighbor_1 == face[0] && neighbor_2 == face[2]) || (neighbor_1 == face[2] && neighbor_2 == face[0]))
		{
			face_index = i;
			return face[1];
		}
	}

	return -1;
}

bool TSplineCDT::judge_Circumcircle_1(const DT_Point& vA, const DT_Point& vB, const DT_Point& vC, const DT_Point& vD)
{
	double BA_u = vA[0] - vB[0];	double BA_v = vA[1] - vB[1];
	double CA_u = vA[0] - vC[0];	double CA_v = vA[1] - vC[1];
	double BD_u = vD[0] - vB[0];	double BD_v = vD[1] - vB[1];
	double CD_u = vD[0] - vC[0];	double CD_v = vD[1] - vC[1];

	double cos_theta_1 = (BA_u * CA_u + BA_v * CA_v) / (sqrt(BA_u * BA_u + BA_v * BA_v) * sqrt(CA_u * CA_u + CA_v * CA_v));

	double cos_theta_2 = (BD_u * CD_u + BD_v * CD_v) / (sqrt(BD_u * BD_u + BD_v * BD_v) * sqrt(CD_u * CD_u + CD_v * CD_v));


	if (cos_theta_1 + cos_theta_2 >= 0)
		return true;
	else
		return false;
}

bool TSplineCDT::judge_Circumcircle_2(const DT_Point& vA, const DT_Point& vB, const DT_Point& vC, const DT_Point& vD)
{
	double BA_u = vA[0] - vB[0];	double BA_v = vA[1] - vB[1];
	double CA_u = vA[0] - vC[0];	double CA_v = vA[1] - vC[1];
	double BD_u = vD[0] - vB[0];	double BD_v = vD[1] - vB[1];
	double CD_u = vD[0] - vC[0];	double CD_v = vD[1] - vC[1];

	double cos_theta_1 = (BA_u * CA_u + BA_v * CA_v) / (sqrt(BA_u * BA_u + BA_v * BA_v) * sqrt(CA_u * CA_u + CA_v * CA_v));
	double cos_theta_2 = (BD_u * CD_u + BD_v * CD_v) / (sqrt(BD_u * BD_u + BD_v * BD_v) * sqrt(CD_u * CD_u + CD_v * CD_v));

	if (cos_theta_1 >= cos_theta_2)
		return true;
	else
		return false;
}

bool TSplineCDT::isTowSegmentsInsersect(const DT_Point& e00, const DT_Point& e01, const DT_Point& e10, const DT_Point& e11)
{
	double sign1 = sign(e00, e10, e11);
	double sign2 = sign(e01, e10, e11);
	double sign3 = sign(e10, e00, e01);
	double sign4 = sign(e11, e00, e01);

	if (sign1 * sign2 < 0 && sign3 * sign4 < 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

DT_Point TSplineCDT::calcIntersection(const DT_Point& e00, const DT_Point& e01, const DT_Point& e10, const DT_Point& e11)
{
	double ax = e00[0], ay = e00[1];
	double bx = e01[0], by = e01[1];
	double cx = e10[0], cy = e10[1];
	double dx = e11[0], dy = e11[1];

	double num = (ay - cy) * (dx - cx) - (ax - cx) * (dy - cy);
	double den = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx);

	double r = num / den;
	return DT_Point(ax + r * (bx - ax), ay + r * (by - ay));
}

bool TSplineCDT::judgeCrossingInDTMesh(const DT_Edge& Edge_1, const DT_Edge& Edge_2)
{
	int target_A = Edge_1[0];	int target_B = Edge_1[1];
	int target_C = Edge_2[0];	int target_D = Edge_2[1];


	if (target_A == target_C || target_A == target_D || target_B == target_C || target_B == target_D)
	{
		return false;
	}

	return isTowSegmentsInsersect(vertices_[target_A], vertices_[target_B], vertices_[target_C], vertices_[target_D]);
}






void TSplineCDT::splitTriangle(int v_index, int face_index)
{
	DT_Triangle& face = triangle_[face_index];

	int index_a = face[0];
	int index_b = face[1];
	int index_c = face[2];

	face[0] = v_index; face[1] = index_a; face[2] = index_b;

	triangle_.push_back(DT_Triangle(index_b, index_c, v_index));
	triangle_.push_back(DT_Triangle(index_c, index_a, v_index));
	return;
}

void TSplineCDT::swapTriangle(int vt, int e1, int e2, int vo, int f_index1, int f_index2)
{
	DT_Triangle& tri1 = triangle_[f_index1];
	for (int i = 0; i < 3; i++)
	{
		if (tri1[i] == e1)
		{
			tri1[i] = vo;
		}
	}
	DT_Triangle& tri2 = triangle_[f_index2];
	for (int i = 0; i < 3; i++)
	{
		if (tri2[i] == e2)
		{
			tri2[i] = vt;
		}
	}
}

void TSplineCDT::addTriangleInOrder(int vertex_A, int vertex_B, int vertex_C)
{
	const DT_Point& pA = vertices_[vertex_A];
	const DT_Point& pB = vertices_[vertex_B];
	const DT_Point& pC = vertices_[vertex_C];

	if (sign(pA, pB, pC) <= 0)
	{
		triangle_.push_back(DT_Triangle(vertex_A, vertex_B, vertex_C));
	}
	else
	{
		triangle_.push_back(DT_Triangle(vertex_A, vertex_C, vertex_B));
	}
}

void TSplineCDT::delete_face(int face_index)
{
	std::vector<DT_Triangle>::iterator it = triangle_.begin();
	it = it + face_index;
	it = triangle_.erase(it);
}


int TSplineCDT::find_TruncateFace_StartingA(const DT_Edge& target_edge)
{
	int num_face = triangle_.size();
	int target_A = target_edge[0];
	int target_B = target_edge[1];

	for (int i = 0; i < num_face; i++)
	{
		const DT_Triangle& face = triangle_[i];
		DT_Edge Edge_2(-1, -1);
		if (face[0] == target_A)
		{
			Edge_2 = DT_Edge(face[1], face[2]);
		}
		else if (face[1] == target_A)
		{
			Edge_2 = DT_Edge(face[2], face[0]);
		}
		else if (face[2] == target_A)
		{
			Edge_2 = DT_Edge(face[0], face[1]);
		}

		if (Edge_2[0] > -1)
		{
			if (judgeCrossingInDTMesh(target_edge, Edge_2))
			{
				return i;
			}
		}
	}
	return -1;
}

int TSplineCDT::opposedTriangle(int face_index, int v_index, int& v_seg, int& v_above, int& v_below, double k_AB, int target_A)
{
	DT_Triangle face = triangle_[face_index];

	int neighbor_1 = -1, neighbor_2 = -1;
	if (face[0] == v_index)
	{
		neighbor_1 = face[1];
		neighbor_2 = face[2];
	}
	else if (face[1] == v_index)
	{
		neighbor_1 = face[2];
		neighbor_2 = face[0];
	}
	else if (face[2] == v_index)
	{
		neighbor_1 = face[0];
		neighbor_2 = face[1];
	}
	else
	{
		//std::cout << "TsplineCDT: oppesed triangle: the face not contain vert" << std::endl;
	}

	const DT_Point& vA = vertices_[target_A];
	const DT_Point& vn1 = vertices_[neighbor_1];
	if ((vn1[1] - vA[1]) - k_AB * (vn1[0] - vA[0]) > 0)
	{
		v_above = neighbor_1;
		v_below = neighbor_2;
	}
	else
	{
		v_above = neighbor_2;
		v_below = neighbor_1;
	}

	int target_face;
	v_seg = find_opposite_target(v_index, neighbor_1, neighbor_2, target_face);
	return target_face;
}

void TSplineCDT::triangulatePseudopolygonDelaunay(std::vector<int> P, const DT_Edge& edge)
{
	int P_size = P.size();

	if (P_size > 1)
	{
		int c = 0;
		for (int i = 1; i < P_size; i++)
		{
			if (judge_Circumcircle_2(vertices_[P[i]], vertices_[edge[0]], vertices_[edge[1]], vertices_[P[c]]))
				continue;
			else
				c = i;
		}

		std::vector<int> P1, P2;
		for (int i = 0; i < P_size; i++)
		{
			if (i < c)
				P1.push_back(P[i]);
			else if (i > c)
				P2.push_back(P[i]);
		}
		DT_Edge edge1(edge[0], P[c]);
		DT_Edge edge2(P[c], edge[1]);
		triangulatePseudopolygonDelaunay(P1, edge1);
		triangulatePseudopolygonDelaunay(P2, edge2);
		addTriangleInOrder(edge[0], edge[1], P[c]);
	}
	else if (P_size == 1)
	{
		addTriangleInOrder(edge[0], edge[1], P[0]);
	}
	else
	{
		return;
	}
}





void TSplineCDT::delete_Outer_Face(const std::vector<DT_Edge>& Fixed_Edge)
{
	int num_face = triangle_.size();
	for (int i = num_face - 1; i >= 0; i--)
	{
		const DT_Triangle& face = triangle_[i];

		const DT_Point& va = vertices_[face[0]];
		const DT_Point& vb = vertices_[face[1]];
		const DT_Point& vc = vertices_[face[2]];

		DT_Point pt((va[0] + vb[0] + vc[0]) / 3, (va[1] + vb[1] + vc[1]) / 3);

		if (!isPointInRegion(Fixed_Edge, pt))
		{
			delete_face(i);
		}
	}
	return;

}

void TSplineCDT::delete_Inner_Face(const std::vector<DT_Edge>& Fixed_Edge)
{
	int num_face = triangle_.size();
	for (int i = num_face - 1; i >= 0; i--)
	{
		const DT_Triangle& face = triangle_[i];

		const DT_Point& va = vertices_[face[0]];
		const DT_Point& vb = vertices_[face[1]];
		const DT_Point& vc = vertices_[face[2]];

		DT_Point pt((va[0] + vb[0] + vc[0]) / 3, (va[1] + vb[1] + vc[1]) / 3);

		if (isPointInRegion(Fixed_Edge, pt))
		{
			delete_face(i);
		}
	}
	return;
}

bool TSplineCDT::isPointInRegion(const std::vector<DT_Edge>& boundary_edge, const DT_Point& pt)
{
	int nCross = 0;
	for (unsigned i = 0; i < boundary_edge.size(); i++)
	{
		const DT_Edge& edge = boundary_edge[i];
		DT_Point& p1 = vertices_[edge[0]];
		DT_Point& p2 = vertices_[edge[1]];


		if (abs(p1[1] - p2[1]) < M_EPSILON)
		{
			continue;
		}

		if ((p1[1] <= pt[1] && p2[1] > pt[1])
			|| (p1[1] >= pt[1] && p2[1] < pt[1]))
		{
			double x = (pt[1] - p1[1]) * (p2[0] - p1[0]) / (p2[1] - p1[1]) + p1[0];
			if (x > pt[0])
			{
				nCross++;
			}
		}














	}

	if (nCross % 2 == 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}





void TSplineCDT::computerBBox()
{
	bbmin_ = DT_Point(DBL_MAX, DBL_MAX);
	bbmax_ = DT_Point(-DBL_MAX, -DBL_MAX);

	for (unsigned i = 0; i < vertices_.size(); i++)
	{
		DT_Point& p = vertices_[i];
		if (bbmin_[0] > p[0])
			bbmin_[0] = p[0];
		if (bbmin_[1] > p[1])
			bbmin_[1] = p[1];

		if (bbmax_[0] < p[0])
			bbmax_[0] = p[0];
		if (bbmax_[1] < p[1])
			bbmax_[1] = p[1];
	}
}

void TSplineCDT::computerVertEdge(const std::vector<DT_Edge>& fixed_edge)
{
	vert_fix_edge_.clear();
	vert_fix_edge_.resize(vertices_.size());
	for (unsigned i = 0; i < fixed_edge.size(); i++)
	{
		const DT_Edge& e = fixed_edge[i];
		vert_fix_edge_[e[0]].push_back(e[1]);
		vert_fix_edge_[e[1]].push_back(e[0]);
	}

}

void TSplineCDT::deleteAuxiliaryTriangle()
{
	int num_vertices = vertices_.size();
	int num_face = triangle_.size();
	int a = num_vertices - 1;
	int b = num_vertices - 2;
	int c = num_vertices - 3;
	for (int i = num_face - 1; i >= 0; i--)
	{
		const DT_Triangle& face = triangle_[i];
		if (face[0] == a || face[1] == a || face[2] == a)
		{
			delete_face(i);
			continue;
		}
		if (face[0] == b || face[1] == b || face[2] == b)
		{
			delete_face(i);
			continue;
		}
		if (face[0] == c || face[1] == c || face[2] == c)
		{
			delete_face(i);
			continue;
		}
	}
}


void TSplineCDT::findBoundaryFixedEdge(std::vector<DT_Edge> Fixed_Edge, std::vector<DT_Edge>& Fixed_Boundary, std::vector<DT_Edge>& inner_Fixed_Edge)
{
	int begin = Fixed_Edge[0][0];
	int index = 0;
	for (unsigned i = 0; i < Fixed_Edge.size(); i++)
	{
		Fixed_Boundary.push_back(Fixed_Edge[i]);
		if (Fixed_Edge[i][1] == begin)
		{
			index = i + 1;
			break;
		}
	}
	for (unsigned i = index; i < Fixed_Edge.size(); i++)
	{
		inner_Fixed_Edge.push_back(Fixed_Edge[i]);
	}
	return;
}

void TSplineCDT::find_outer_and_inner_Boundary(std::vector<DT_Edge> Fixed_Edge, std::vector<DT_Edge>& outer_Fixed_Boundary, std::vector<std::vector<DT_Edge>>& inner_Fixed_Boundary_list)
{
	int begin = Fixed_Edge[0][0];
	unsigned index = 0;


	for (unsigned i = 0; i < Fixed_Edge.size(); i++)
	{
		outer_Fixed_Boundary.push_back(Fixed_Edge[i]);
		if (Fixed_Edge[i][1] == begin)
		{
			index = i + 1;
			break;
		}
	}

	while (index < Fixed_Edge.size())
	{
		begin = Fixed_Edge[index][0];


		std::vector<DT_Edge> inner_Fixed_Boundary; inner_Fixed_Boundary.clear();
		for (unsigned i = index; i < Fixed_Edge.size(); i++)
		{
			inner_Fixed_Boundary.push_back(Fixed_Edge[i]);
			if (Fixed_Edge[i][1] == begin)
			{
				index = i + 1;
				break;
			}
		}
		inner_Fixed_Boundary_list.push_back(inner_Fixed_Boundary);
	}

	return;
}

void TSplineCDT::splitCrossingFixedEdge(std::vector<DT_Edge>& fixed_edge)
{
	std::vector<DT_Edge> new_fixed_Edge; new_fixed_Edge.clear();

	bool is_insert;
	while (fixed_edge.size() > 0)
	{
		is_insert = false;
		DT_Edge Edge_1 = fixed_edge[0];
		for (unsigned i = 1; i < fixed_edge.size(); i++)
		{
			DT_Edge Edge_2 = fixed_edge[i];
			if (judgeCrossingInDTMesh(Edge_1, Edge_2) == true)
			{
				is_insert = true;
				DT_Point inter_Pt = calcIntersection(
					vertices_[Edge_1[0]], vertices_[Edge_1[1]], vertices_[Edge_2[0]], vertices_[Edge_2[1]]);

				vertices_.push_back(inter_Pt);


				std::vector<DT_Edge>::iterator it = fixed_edge.begin();
				fixed_edge.erase(it + i);
				fixed_edge.erase(it);
				int fix_v_num = vertices_.size() - 1;
				fixed_edge.push_back(DT_Edge(fix_v_num, Edge_1[0]));
				fixed_edge.push_back(DT_Edge(fix_v_num, Edge_1[1]));
				fixed_edge.push_back(DT_Edge(fix_v_num, Edge_2[0]));
				fixed_edge.push_back(DT_Edge(fix_v_num, Edge_2[1]));
				break;
			}
		}

		if (is_insert == 0)
		{
			std::vector<DT_Edge>::iterator it = fixed_edge.begin();
			fixed_edge.erase(it);
			new_fixed_Edge.push_back(Edge_1);
		}
	}
	fixed_edge = new_fixed_Edge;
}

void TSplineCDT::delete_noninner_Fixed_Edge(std::vector<DT_Edge>& Fixed_Boundary, std::vector<DT_Edge>& inner_Fixed_Edge)
{
	std::vector<DT_Edge> old_fix_edge = inner_Fixed_Edge;
	inner_Fixed_Edge.clear();

	bool isCross;
	for (unsigned i = 0; i < old_fix_edge.size(); i++)
	{
		isCross = false;
		const DT_Edge& Edge_1 = old_fix_edge[i];
		for (unsigned j = 0; j < Fixed_Boundary.size(); j++)
		{
			const DT_Edge& Edge_2 = Fixed_Boundary[j];
			if (judgeCrossingInDTMesh(Edge_1, Edge_2))
			{
				isCross = true;
				break;
			}
		}

		if (!isCross)
		{
			inner_Fixed_Edge.push_back(Edge_1);
		}
	}
}

void TSplineCDT::orientation_reserve()
{
	for (unsigned i = 0; i < triangle_.size(); i++)
	{
		const DT_Triangle& face = triangle_[i];

		DT_Point p1 = vertices_[face[0]];
		DT_Point p2 = vertices_[face[1]];
		DT_Point p3 = vertices_[face[2]];
		if (sign(p1, p2, p3) < 0)
		{
			triangle_[i] = DT_Triangle(face[1], face[0], face[2]);
		}
		else
		{
			continue;
		}
	}
	return;
}


void TSplineCDT::load_info(const std::string& _filename_1, const std::string& _filename_2, std::vector<DT_Point>& Fixed_Vertices, std::vector<DT_Edge>& Fixed_Edge)
{


	int sum_v = Fixed_Vertices.size();

	std::ifstream fin1(_filename_1, std::ios::in);
	if (!fin1.is_open())
	{
		//std::cout << "δ�ɹ����ļ�" << std::endl;
	}
	double data1, data2, data3;
	while (!fin1.eof())
	{
		fin1 >> data1 >> data2 >> data3;
		Fixed_Vertices.push_back({ data1,data2 });
	}
	fin1.close();
	Fixed_Vertices.pop_back();


	std::ifstream fin2(_filename_2);
	if (!fin2.is_open())
	{
		//std::cout << "δ�ɹ����ļ�" << std::endl;
	}
	int idata1, idata2;
	while (fin2)
	{
		fin2 >> idata1 >> idata2;
		Fixed_Edge.push_back({ idata1 + sum_v,idata2 + sum_v });
	}
	fin2.close();
	Fixed_Edge.pop_back();

}

void TSplineCDT::print_info(const std::string& _filename_1, const std::string& _filename_2)
{


	std::ofstream fout_1(_filename_1, std::ios::trunc);
	for (unsigned i = 0; i < vertices_.size(); i++)
	{
		DT_Point p = vertices_[i];
		fout_1 << p[0] << "\t" << p[1] << "\t" << 0.0 << std::endl;
	}
	fout_1.close();


	std::ofstream fout_2(_filename_2, std::ios::trunc);
	for (unsigned i = 0; i < triangle_.size(); i++)
	{
		DT_Triangle t = triangle_[i];
		fout_2 << t[0] << "\t" << t[1] << "\t" << t[2] << std::endl;
	}
	fout_2.close();
}

void TSplineCDT::output2io(std::vector<double>& V, std::vector<int>& F)
{
	V.resize(vertices_.size() * 2);
	F.resize(triangle_.size() * 3);

	for (unsigned i = 0; i < vertices_.size(); i++)
	{
		V[2 * i] = vertices_[i][0];
		V[2 * i + 1] = vertices_[i][1];
	}

	for (unsigned i = 0; i < triangle_.size(); i++)
	{
		F[3 * i] = triangle_[i][0];
		F[3 * i + 1] = triangle_[i][1];
		F[3 * i + 2] = triangle_[i][2];
	}

}


void TSplineCDT::print_edge(const std::vector<DT_Edge>& Fixed_Edge, const std::string& _filename_)
{
	std::ofstream fout(_filename_, std::ios::trunc);
	for (unsigned i = 0; i < Fixed_Edge.size(); i++)
	{
		fout << Fixed_Edge[i][0] << "\t" << Fixed_Edge[i][1] << std::endl;
	}
	fout.close();
}
}
}
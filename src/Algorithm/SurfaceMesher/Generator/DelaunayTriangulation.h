#pragma once
#include <math.h>
#include <vector>
#include <cfloat>
#include "MathConstant.h"

namespace PDF3D
{
namespace CDT
{
struct DT_Edge
{
private:
	int p_[2];
public:
	int& operator[](int i) { return p_[i]; }
	int operator[](int i) const { return p_[i]; }

	DT_Edge(int a, int b) { p_[0] = a; p_[1] = b; }
	DT_Edge() { p_[0] = p_[1] = 0; }
};

struct DT_Triangle
{
private:
	int p_[3];
public:
	int& operator[](int i) { return p_[i]; }
	int operator[](int i) const { return p_[i]; }

	DT_Triangle(int a, int b, int c) { p_[0] = a; p_[1] = b; p_[2] = c; }
	DT_Triangle() { p_[0] = p_[1] = p_[2] = 0; }

	bool isContainVertex(int vert) const
	{
		return (vert == p_[0] || vert == p_[1] || vert == p_[2]);
	}
};

struct DT_Point
{
private:
	double p_[2];
public:
	double& operator[](int i) { return p_[i]; }
	double operator[](int i) const { return p_[i]; }
	bool operator==(DT_Point& point)const
	{
		return((p_[0] == point[0]) && (p_[1] == point[1]));
	}

	DT_Point(double a, double b) { p_[0] = a; p_[1] = b; }
	DT_Point() { p_[0] = p_[1] = 0; }
};

class DelaunayTriangulation
{};
}// namespace CDT
}// namespace PDF3D
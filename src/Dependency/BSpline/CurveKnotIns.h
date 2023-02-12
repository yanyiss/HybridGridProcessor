#pragma once
#include <vector>
#include <Eigen\Dense>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include "BezierCurve.h"
#include "BSplineCurve.h"

typedef Eigen::Vector4d Vec4d;
typedef Eigen::Vector4d Point4;

class CurveKnotIns
{
public:
	CurveKnotIns();
	~CurveKnotIns();

	// Insert a knot.
	CurveKnotIns(int &np, int &p, std::vector<double> &UP, std::vector<Point4> &Pw, double &u, int &k, int &s, int &r, int &nq, std::vector<double> &UQ, std::vector<Point4> &Qw);

	// Insert many knots.
	CurveKnotIns(int &n, int &p, std::vector<double> &U, std::vector<Point4> &Pw, std::vector<double> &X,  int &r,std::vector<double> &Ubar, std::vector<Point4> &Qw);

	void FindSpanMult(double u, int &knot_span, int &multi);

	// Split Bspline to Bezier curves.
	void BsplineToBeziers_Curve(int &n, int &p,std::vector<double> U, std::vector<Point4> Pw, int &nb, std::vector<std::vector<Point4>> &Qw);

	// Split Bezier to Bezier curves.
	void BezierDecomposition_Curve(int p, std::vector<Point4> &Pw, double u_min, double u_max, double u, std::vector<Point4> &left_Qw, std::vector<Point4> &right_Qw);

private:
	BSplineCurve	curve_;

};

/* 残留问题 */
/* 细分之后每一段的参数范围 */



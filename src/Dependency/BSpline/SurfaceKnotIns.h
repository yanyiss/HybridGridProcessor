#pragma once
#include <vector>
#include <Eigen\Eigen>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include "BezierSurface.h"
#include "BSplineSurface.h"

typedef Eigen::Vector4d Vec4d;
typedef Eigen::Vector4d Point4;

class SurfaceKnotIns
{
public:
	SurfaceKnotIns();
	~SurfaceKnotIns();

	// Insert a knot.
	SurfaceKnotIns(int &np, int &p, std::vector<double> &UP, int &mp, int &q, std::vector<double> &VP, std::vector<std::vector<Point4>> &Pw, DIRECTION dir, int &uv, int &k, int &s, int &r, int &nq, std::vector<double> &UQ, int &mq, std::vector<double> &VQ, std::vector<std::vector<Point4>> &Qw);

	// Insert many knots.
	SurfaceKnotIns(int &n, int &p, std::vector<double> &U, int &m, int &q, std::vector<double> &V, std::vector<std::vector<Point4>> &Pw, std::vector<double> &X, int &r, DIRECTION dir,  std::vector<double> &Ubar, std::vector<double> &Vbar, std::vector<std::vector<Point4>> &Qw);

	// Split Bspline to Bezier surfaces.
	void  BsplineToBeziers_Surface(int &n, int &p, std::vector<double> &U, int &m, int &q, std::vector<double> &V, std::vector<std::vector<Point4>> &Pw, DIRECTION &dir, int &nb, std::vector<std::vector<std::vector<Point4>>>Qw);

	// Split Bezier to Bezier suefaces.
	void BezierDecomposition_Surface(int &p, int &q, std::vector<std::vector<Point4>> Pw, double u_min, double u_max, double v_max, double v_min, double uv, DIRECTION &dir, std::vector<std::vector<Point4>> &left_Qw, std::vector<std::vector<Point4>> &right_Qw);
	/* 细分之后每一段的参数范围 */



	void FindSpanMult(double uv, DIRECTION dir, int &knot_span, int &multi);

private:
	BSplineSurface surface_;

};



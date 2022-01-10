#pragma once
#ifndef NEWTON_RAPHSON
#define NEWTON_RAPHSON

#include "BSplineCurve.h"
#include "BSplineSurface.h"
class NewtonRaphson
{
public:
	NewtonRaphson();
	~NewtonRaphson();
	NewtonRaphson(Point P, BezierCurve &curve, double &t);
	NewtonRaphson(Point P, BezierSurface &surface, double &u, double &v);
	NewtonRaphson(Point P, BSplineCurve &curve, double &t);
	NewtonRaphson(Point P, BSplineSurface &surface, double &u, double &v);

	void UpdateParameter(double &t, double dt, double t_min, double t_max);
	void UpdateParameter(double &u, double &v, double du, double dv, double u_min, double u_max, double v_min, double v_max);

private:	
	double	eps = 1e-14;
	double	ERROR_BOUND = 1.0E-8;
	double	COSINE_BOUND = cos(M_PI * 89 / 180);	// cos 89бу= 0.01745
};

#endif
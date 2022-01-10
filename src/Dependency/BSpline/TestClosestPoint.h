#pragma once
#ifndef TEST_CLOSEST_POINT
#define TEST_CLOSEST_POINT

#include "SquaredDistanceFunc_Curve.h"
#include "SquaredDistanceFunc_Surface.h"
#include <deque>
#include "NewtonRaphson.h"

//double Binomial(int n, int k);
//std::vector<double> Binomial(int n);

struct ProjectionPointToCurve
{
	double t;
	double dist;
};

struct ProjectionPointToSurface
{
	double u;
	double v;
	double dist;
};

class TestClosestPoint
{
public:
	TestClosestPoint(BezierCurve &curve, Point &testpoint, double &dist, std::vector<ProjectionPointToCurve> &closest_pts);
	TestClosestPoint(BezierSurface &surface, Point &testpoint, std::vector<ProjectionPointToSurface> &closest_pts);
	TestClosestPoint(BSplineCurve &curve, Point &testpoint, double &dist, std::vector<ProjectionPointToCurve> &closest_pts);
	TestClosestPoint(BSplineSurface &surface, Point &testpoint, std::vector<ProjectionPointToSurface> &closest_pts);

	std::vector<ProjectionPointToCurve> GetProjrctionToCurve() { return closest_pts_c; }
	std::vector<ProjectionPointToSurface> GetProjrctionToSurface() { return closest_pts_s; }
	double GetDistance() { return curr_min_dist; }

	~TestClosestPoint();

	bool TestEndPointOfCurve(BezierCurve &curve, Point P);
	bool TestIntUniquenessOfCurve(BezierCurve &curve, Point P);
	
	bool TestCornerPointOfSurface(BezierSurface &surface, Point P);
	bool TestBoudaryCurveOfSurface(BezierSurface &surface, Point P);
	bool TestIntUniquenessOfSurface(BezierSurface &surface, Point P);
	bool TestIntUniquenessOfSurface(BezierSurface &surface, Point P, bool);
	
	//void UpdateInitParamOfCurve();
	//void UpdateLowerBoundOfDist();
	void UpdateCurveMessage(double t, double dist);
	void UpdateSurfaceMessage(double u, double v, double dist);

private:
	Point P;
	double curr_min_dist;
	std::deque<BezierCurve> RemainingCurveSegments;
	std::deque<BezierSurface> RemainingSurfacePatches;
	std::vector<ProjectionPointToCurve> closest_pts_c;
	std::vector<ProjectionPointToSurface> closest_pts_s;
};

#endif
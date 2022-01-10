#pragma once
#ifndef SQUARED_DISTANCE_FUN_C
#define SQUARED_DISTANCE_FUN_C
#include "BezierCurve.h"

class SquaredDistanceFunc_Curve
{
public:
	SquaredDistanceFunc_Curve();
	~SquaredDistanceFunc_Curve();
	SquaredDistanceFunc_Curve(BezierCurve &bezier_, Point P);
	SquaredDistanceFunc_Curve(int deg, double tmin, double tmax, std::vector<double> coefficients);
	SquaredDistanceFunc_Curve(int deg, double tmin, double tmax, std::vector<double> w, std::vector<double> coefficients);
	
	void SetDegree(const int deg) { degree = deg; }
	void SetParameterSpan(const double min, const double max) { t_min = min; t_max = max; }
	void SetRational(const bool is_rational) { isRational = is_rational; }
	void SetWeights(const std::vector<double> w) { weights.clear(); weights = w; }
	void SetCoeffs(const std::vector<double> & coeffs_) { coeffs.clear(); coeffs = coeffs_; }

	const int GetDegree(void) const { return degree; }
	const void GetParameterSpan(double &min, double &max) const { min = t_min; max = t_max; }
	const bool GetRational(void) const { return isRational; }
	const std::vector<double> GetWeights(void) const { return weights; }
	const std::vector<double> GetCoeffs(void) const { return coeffs; }

	double DeCasteljau(const std::vector<double> &coefficients, int deg, const double t)const;
	double Derivative(const std::vector<double> &coefficients, const double t, const int k) const;

	double operator()(const double t) const { return DeCasteljau(t); }
	double DeCasteljau(const double t)const;
	double Derivative(const double t)const;
	double Derivative2(const double t)const;

private:
	int degree;
	double t_min;
	double t_max;
	bool isRational;
	std::vector<double> coeffs;
	std::vector<double> weights;
};

double Binomial(int n, int k);

std::vector<double> Binomial(int n);

#endif
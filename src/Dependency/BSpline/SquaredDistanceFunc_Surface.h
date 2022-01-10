#pragma once
#ifndef SQUARED_DISTANCE_FUN_S
#define SQUARED_DISTANCE_FUN_S

#include "BezierSurface.h"
double Binomial(int n, int k);
std::vector<double> Binomial(int n);
class SquaredDistanceFunc_Surface
{
public:
	SquaredDistanceFunc_Surface();
	SquaredDistanceFunc_Surface(int u_deg, int v_deg, double umin, double umax, double vmin, double vmax, std::vector<std::vector<double>> w, std::vector<std::vector<double>> coefficients);
	SquaredDistanceFunc_Surface(int u_deg, int v_deg, double umin, double umax, double vmin, double vmax, std::vector<std::vector<double>> coefficients);
	SquaredDistanceFunc_Surface(BezierSurface &bezier_, Point P);
	~SquaredDistanceFunc_Surface();

	void SetDegree(const int degu, const int degv) { u_degree = degu; v_degree = degv; }
	void SetParameterSpan(const double umin, const double umax, const double vmin, const double vmax)
	{
		u_min = umin; u_max = umax; v_min = vmin, v_max = vmax;
	}
	void SetRational(const bool is_rational) { isRational = is_rational; }
	void SetWeights(const std::vector<std::vector<double>> w) { weights.clear(); weights = w; }
	void SetCoeffs(const std::vector<std::vector<double>> & coeffs_) { coeffs.clear(); coeffs = coeffs_; }

	const int GetDegreeU(void) const { return u_degree; }
	const int GetDegreeV(void) const { return v_degree; }
	const void GetParameterSpan(double &umin, double &umax, double &vmin, double &vmax) const
	{
		umin = u_min; umax = u_max; vmin = v_min; vmax = v_max;
	}
	const bool GetRational(void) const { return isRational; }
	const std::vector<std::vector<double>> GetWeights(void) const { return weights; }
	const std::vector<std::vector<double>> GetCoeffs(void) const { return coeffs; }

	double DeCasteljau(const std::vector<std::vector<double>> &coefficients, int degu, int degv, const double u, const double v) const;
	double Derivative(const std::vector<std::vector<double>> &coefficients, const double u, const double v, const int k, const int l) const;

	double operator()(const double u, const double v) const { return DeCasteljau(u, v); }
	double DeCasteljau(const double u, const double v) const;
	double PartialDerivativeU(const double u, const double v) const;
	double PartialDerivativeV(const double u, const double v)const;
	double PartialDerivativeUU(const double u, const double v) const;
	double PartialDerivativeUV(const double u, const double v)const;
	double PartialDerivativeVV(const double u, const double v) const;

private:
	size_t u_degree;
	size_t v_degree;
	double u_max;
	double u_min;
	double v_max;
	double v_min;
	bool isRational;
	std::vector<std::vector<double>> weights;
	std::vector<std::vector<double>> coeffs;
};

#endif
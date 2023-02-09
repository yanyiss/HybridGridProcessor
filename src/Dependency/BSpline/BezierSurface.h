#pragma once
#include <vector>
#include <iostream>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include"GeometryType.h"

enum DIRECTION { U_DIRECTION, V_DIRECTION };

class BezierSurface : public GeometryType
{
public:
	BezierSurface();
	BezierSurface(int u_deg, int v_deg, double umin, double umax, double vmin, double vmax, std::vector<std::vector<double>> w, std::vector<std::vector<Point>> controlpoints);
	BezierSurface(int u_deg, int v_deg, double umin, double umax, double vmin, double vmax, std::vector<std::vector<Point>> controlpoints);
	~BezierSurface();

	void SetDegree(const int degu, const int degv) { u_degree = degu; v_degree = degv; }
	void SetParameterSpan(const double umin, const double umax, const double vmin, const double vmax)
		{ u_min = umin; u_max = umax; v_min = vmin, v_max = vmax; }
	void SetRational(const bool is_rational) { isRational = is_rational; }
	void SetWeights(const std::vector<std::vector<double>> w) 
	{
		weights.clear(); 
		weights = w; 
		isRational = true;
	}
	void SetControlPoints(const std::vector<std::vector<Point>> & points) 
		{ ctrlpoints.clear(); ctrlpoints = points; }

	const int GetDegreeU(void) const { return u_degree; }
	const int GetDegreeV(void) const { return v_degree; }
	const void GetParameterSpan(double &umin, double &umax, double &vmin, double &vmax) const
		{ umin = u_min; umax = u_max; vmin = v_min; vmax = v_max; }
	const bool GetRational(void) const { return isRational; }
	const std::vector<std::vector<double>> GetWeights(void) const { return weights; }
	const std::vector<std::vector<Point>> GetControlPoints(void) const { return ctrlpoints; }

	bool NormalValue(const double u, const double v, Point &n) const
	{
		n = PartialDerivativeU(u, v).cross(PartialDerivativeV(u, v));
		if (n.norm() < 1e-10) return false;
		n.normalize();
		return true;
	}

	bool PrincipalCurvature(const double u, const double v, double &k1, double &k2) const  //计算主曲率
	{
		Eigen::Matrix2d Weingarten, value;
		Point ru = PartialDerivativeU(u, v);
		Point rv = PartialDerivativeV(u, v);
		double E = ru.dot(ru);
		double F = ru.dot(rv);
		double G = rv.dot(rv);
		Point n = ru.cross(rv);
		if (n.norm() < 1e-10) return false;
		n.normalize();
		double L = PartialDerivativeUU(u, v).dot(n);
		double M = PartialDerivativeUV(u, v).dot(n);
		double N = PartialDerivativeVV(u, v).dot(n);
		Weingarten << L * G - M * F, M * E - L * F,
			M * G - N * F, N * E - M * F;
		Eigen::EigenSolver<Eigen::Matrix2d> es(Weingarten);
		value = (es.pseudoEigenvalueMatrix()) / (E*G - pow(F, 2));
		k1 = std::max(std::abs(value(0, 0)), std::abs(value(1, 1)));
		k2 = std::min(std::abs(value(0, 0)), std::abs(value(1, 1)));
		return true;
	}

	bool NormalCurvature(const double u, const double v, const double x, const double y, double &k) const
	{
		Eigen::Matrix2d Weingarten, value, eigenvector;
		Point ru = PartialDerivativeU(u, v);
		Point rv = PartialDerivativeV(u, v);
		double E = ru.dot(ru);
		double F = ru.dot(rv);
		double G = rv.dot(rv);
		Point n = ru.cross(rv);
		if (n.norm() < 1e-10) return false;
		n.normalize();
		double L = PartialDerivativeUU(u, v).dot(n);
		double M = PartialDerivativeUV(u, v).dot(n);
		double N = PartialDerivativeVV(u, v).dot(n);
		k = std::abs((L*x*x + 2 * M*x*y + N * y*y) / (E*x*x + 2 * F*x*y + G * y*y));
		return true;
	}

	template<typename T>
	T DeCasteljau(const std::vector<std::vector<T>> controlpoints, const double u, const double v) const
	{
		if (u < u_min || u > u_max) return T();
		if (v < v_min || v > v_max) return T();

		double u_ = (u - u_min) / (u_max - u_min);
		double v_ = (v - v_min) / (v_max - v_min);

		std::vector<std::vector<T>> d(controlpoints);
		for (int j = 0; j <= v_degree; j++) // P_(0,vj)^(u_degree)
		{
			for (int k = 1; k <= u_degree; k++)
			{
				for (int i = 0; i <= u_degree - k; i++)
				{
					d[i][j] = (1 - u_) * d[i][j] + u_ * d[i + 1][j];
				}
			}
		}

		for (int l = 1; l <= v_degree; l++)
		{
			for (int j = 0; j <= v_degree - l; j++)
			{
				d[0][j] = (1 - v_) * d[0][j] + v_ * d[0][j + 1];
			}
		}

		return d[0][0];
	}

	template<typename T>
	T Derivative(std::vector<std::vector<T>> controlpoints, const double u, const double v, const int k, const int l) const
	{
		if (k > u_degree || l > v_degree)
			return T();

		if (u < u_min || u > u_max) return T();
		if (v < v_min || v > v_max) return T();

		std::vector<std::vector<T>> d(controlpoints);
		for (int a = 1; a <= k; a++)
		{
			for (int i = 0; i <= u_degree - a; i++)
			{
				for (int j = 0; j <= v_degree; j++)
				{
					d[i][j] = (u_degree - a + 1) / (u_max - u_min) * (d[i + 1][j] - d[i][j]);
				}
			}
		}

		for (int b = 1; b <= l; b++)
		{
			for (int i = 0; i <= u_degree - k; i++)
			{
				for (int j = 0; j <= v_degree - b; j++)
				{
					d[i][j] = (v_degree - b + 1) / (v_max - v_min) * (d[i][j + 1] - d[i][j]);
				}
			}
		}

		std::vector<std::vector<T>> dev_ctrlpts(u_degree + 1 - k, std::vector<T>(v_degree + 1 - l));
		for (int i = 0; i <= u_degree - k; i++)
		{
			for (int j = 0; j <= v_degree - l; j++)
			{
				dev_ctrlpts[i][j] = d[i][j];
			}
		}

		BezierSurface dev_surface;
		dev_surface.SetDegree(u_degree - k, v_degree - l);
		dev_surface.SetParameterSpan(u_min, u_max, v_min, v_max);

		return dev_surface.DeCasteljau(dev_ctrlpts, u, v);
	}

	Point operator()(const double u, const double v) const { return DeCasteljau(u, v); }
	Point DeCasteljau(const double u, const double v) const;
	Point PartialDerivativeU(const double u, const double v) const;
	Point PartialDerivativeV(const double u, const double v)const;
	Point PartialDerivativeUU(const double u, const double v) const;
	Point PartialDerivativeUV(const double u, const double v)const;
	Point PartialDerivativeVV(const double u, const double v) const;

	void Split(double uv, DIRECTION direction, BezierSurface &left_s, BezierSurface &right_s);

	bool SaveControlPoints(const std::string & filename) const;
	bool LoadControlPoints(const std::string & filename);

private:
	size_t u_degree;
	size_t v_degree;
	double u_max;
	double u_min;
	double v_max;
	double v_min;
	bool isRational;
	std::vector<std::vector<double>> weights;
	std::vector<std::vector<Point>> ctrlpoints;
};

/* 参数区间不是[0,1]，求值求导怎么办*/

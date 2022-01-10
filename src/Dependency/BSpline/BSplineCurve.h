#pragma once
#ifndef BSPLINE_C
#define BSPLINE_C
#include "BezierCurve.h"
//#include <Eigen/Eigen>

typedef Eigen::Vector3d Point;
typedef Eigen::Vector4d Point4;

class BSplineCurve
{
public:
	BSplineCurve();
	BSplineCurve(int deg, std::vector<double> knots, std::vector<double> weights, std::vector<Point> ctrlpoints);
	BSplineCurve(int deg, std::vector<double> knots, std::vector<Point> ctrlpoints);
	~BSplineCurve();
	
	void SetDegree(const int deg) { degree = deg; }
	void SetRational(const bool rational) { isRational = rational; }
	void SetKnots(const std::vector<double> knots_) { knots.clear(); knots = knots_; }
	void SetWeights(const std::vector<double> weights_) 
	{
		weights.clear();
		weights = weights_; 
		isRational = true;
	}
	void SetControlPoints(const std::vector<Point> &points) { ctrlpoints.clear(); ctrlpoints = points; }
	void SetDataPoints(const std::vector<Point> &points) { datapoints.clear(); datapoints = points; }

	const int GetDegree(void) { return degree; }
	const bool GetRational(void) { return isRational; }
	const int GetNumOfKnots(void) const { return knots.size(); }
	const int GetNumOfCtrlpts(void) const { return ctrlpoints.size(); }
	const std::vector<double> GetKnots(void) const { return knots; }
	const std::vector<double> GetWeights(void) const { return weights; }
	const std::vector<Point> GetControlPoints(void) const { return ctrlpoints; }
	const int FindSpan(int n, int p, double u, std::vector<double> K) const ;
	
	template<typename T>
	T DeBoor(std::vector<T> controlpoints, const double t) const
	{
		if (controlpoints.empty()) return T();
		if (t < knots.front() || t > knots.back()) return T();

		int n = controlpoints.size() - 1;
		int r = FindSpan(n, degree, t, knots);

		std::vector<T> d(controlpoints.begin() + r - degree, controlpoints.begin() + r + 1);
		for (int j = 1; j <= degree; j++)		// compute d^j_i backwards
		{
			for (int i = degree; i >= j; i--)
			{
				double alpha = (t - knots[r - degree + i]) / (knots[i + 1 + r - j] - knots[i + r - degree]);
				d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
			}
		}

		return d[degree];
	}

	template<typename T>
	T Derivative(std::vector<T> controlpoints, const double t, const int k) const
	{
		if (controlpoints.empty()) return T();
		if (k > degree) return T();
		if (t < knots.front() || t > knots.back()) return T();

		int n = controlpoints.size() - 1;
		std::vector<T> d(controlpoints);
		for (int a = 1; a <= k; a++)
		{
			for (int i = 0; i <= n - a; i++)
			{
				d[i] = (d[i + 1] - d[i]) * (degree - a + 1) / (knots[i + degree + 1] - knots[i + a]);
			}
		}

		std::vector<T> dev_ctrlpts(d.begin(), d.end() - k);
		std::vector<double> dt_knots(knots.begin() + k, knots.end() - k);

		BSplineCurve dev_curve;
		dev_curve.SetDegree(degree - k);
		dev_curve.SetKnots(dt_knots);

		return dev_curve.DeBoor(dev_ctrlpts, t);
	}

	Point operator()(const double t) const { return DeBoor(t); }
	Point DeBoor(const double t) const;
	Point Derivative(const double t) const;
	Point Derivative2(const double t) const;

	template<typename T>
	void KnotInsertion(double t, int k, const std::vector<double> &U, const std::vector<T> &Pw, std::vector<double> &Unew, std::vector<T> &Qw)
	{
		int m = U.size() - 1;
		int n = Pw.size() - 1;

		int mult = 0;
		int r = FindSpan(n, degree, t, knots);
		for (int i = r; i >= degree; i--)
		{
			if (abs(t - U[i]) < DBL_EPSILON)
				mult++;
		}
		// if (k + r > degree)	k = degree - mult;

		Unew.resize(m + 1 + k);
		Qw.resize(n + 1 + k);

		/* Load new knot vector */
		for (int i = 0; i <= r; i++) Unew[i] = U[i];
		for (int i = 1; i <= k; i++) Unew[r + i] = t;
		for (int i = r + 1; i <= m; i++) Unew[i + k] = U[i];

		/* Save unaltered control points */
		for (int i = 0; i <= r - degree; i++)	Qw[i] = Pw[i];
		for (int i = r - mult; i <= n; i++)		Qw[i + k] = Pw[i];
		std::vector<T> Rw(Pw.begin() + r - degree, Pw.begin() + r - mult + 1);

		/* Insert the knot k times */
		int L;
		for (int j = 1; j <= k; j++)
		{
			L = r - degree + j;
			for (int i = 0; i <= degree - mult - j; i++)
			{
				double alpha = (t - U[i + L]) / (U[i + r + 1] - U[i + L]);
				Rw[i] = alpha * Rw[i + 1] + (1.0 - alpha) * Rw[i];
			}
			Qw[L] = Rw[0];
			Qw[r - mult + k - j] = Rw[degree - mult - j];
		}

		/* Load remaining control points */
		for (int i = 1; i < degree - mult - k; i++)
		{
			Qw[L + i] = Rw[i];
		}

	}

	template<typename T>
	void BSplineToBezier(int p, std::vector<double> &U, std::vector<T> &Pw, int &nb, std::vector<double> &Unew, std::vector<std::vector<T>> &Qw)
	{
		// int n = Pw.size() - 1;
		int m = U.size() - 1;
		int a = p, b = p + 1;
		nb = 0;	// number of bezier segments, count from 0;
		Qw.clear();
		Unew.clear();
		Unew.push_back(U[a]);

		std::vector<T> Qw_0(p + 1);
		for (int i = 0; i <= p; i++) Qw_0[i] = Pw[i];
		Qw.push_back(Qw_0);

		while (b < m)
		{
			int multi;
			int i = b;
			while (b < m && ((U[b + 1] - U[b]) < DBL_EPSILON))
				b++;
			multi = b - i + 1;

			std::vector<T> Qw_(p + 1);
			if (multi < p)
			{
				/* Compute alphas */
				std::vector<double> alphas(p - multi);
				double numerator = U[b] - U[a];
				for (int j = 0; j < p - multi; j++)
					alphas[j] = numerator / (U[a + multi + 1 + j] - U[a]);

				/* Insert U[b] r times */
				int r = p - multi;
				for (int j = 1; j <= r; j++)
				{
					int save = r - j;
					int s = multi + j;
					for (int k = p; k >= s; k--)
					{
						double alpha = alphas[k - s];
						Qw[nb][k] = alpha * Qw[nb][k] + (1.0 - alpha) * Qw[nb][k - 1];
					}
					if (b < m)	// control points of next segment
						Qw_[save] = Qw[nb][p];		//Qw[nb + 1][save] = Qw[nb][p]; 				
				}
			}
			Unew.push_back(U[b]);

			/* Initialize for next segment */
			if (b < m)
			{
				for (int i = p - multi; i <= p; i++)
					Qw_[i] = Pw[b - p + i];		//Qw[nb][i] = Pw[b - p + i];
				a = b;
				b = b + 1;
				Qw.push_back(Qw_);
				nb++;
			}
		}
	}


	void KnotInsertion(double t, int k, BSplineCurve &new_curve);
	void Split(double t, BSplineCurve &left_c, BSplineCurve &right_c);
	void BSplineToBezier(std::vector<BezierCurve> &bezier_curves);

		
	void Parameterization(std::vector<double>& t);
	const std::vector<Point> & DataPoints(void) const { return datapoints; }
	void UpdateKnots(void);
	void Regression(void);

private:
	int degree;
	bool isRational;
	std::vector<double> knots;
	std::vector<double> weights;
	std::vector<Point> ctrlpoints;
	
	std::vector<Point> datapoints;
};
#endif

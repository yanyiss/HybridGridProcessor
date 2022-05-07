#include "BSplineSurface.h"
#include <ctime>
#include <fstream>
#include <iostream>   


BSplineSurface::BSplineSurface()
{
	SetRational(false);
	u_knots.clear();
	v_knots.clear();
	weights.clear();
	ctrlpoints.clear();
}

BSplineSurface::BSplineSurface(int udeg, int vdeg, std::vector<double> &uknots, std::vector<double> &vknots, std::vector<std::vector<double>> &w, std::vector<std::vector<Point>> &controlpoints)
{
	SetDegree(udeg, vdeg);
	SetKnotsU(uknots);
	SetKnotsV(vknots);
	SetRational(true);
	SetWeights(w);
	SetControlPoints(controlpoints);
	int m = GetNumOfCtrlptsU() - 1;
	int n = GetNumOfCtrlptsV() - 1;
	UV = Point4(uknots[u_degree], uknots[m + 1], vknots[v_degree], vknots[n + 1]); //2022.03.21
	Pw_ctrlpts.resize(m + 1, std::vector<Point4>(n + 1));
	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			Pw_ctrlpts[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
		}
	}

	PartialUCtrpts = Pw_ctrlpts;
	for (int i = 0; i <= m - 1; i++)
	{
		double diffu = u_knots[i + u_degree + 1] - u_knots[i + 1];
		if (diffu < DBL_EPSILON)
		{
			for (int j = 0; j <= n; j++)
			{
				setZero(PartialUCtrpts[i][j]);
			}	
			continue;
		}
		for (int j = 0; j <= n; j++)
		{
			PartialUCtrpts[i][j] = (PartialUCtrpts[i + 1][j] - PartialUCtrpts[i][j]) * (u_degree) / diffu;
		}
	}

	PartialVCtrpts = Pw_ctrlpts;
	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= n - 1; j++)
		{
			double diffv = v_knots[j + v_degree + 1] - v_knots[j + 1];
			if (diffv < DBL_EPSILON) setZero(PartialVCtrpts[i][j]);
			else PartialVCtrpts[i][j] = (PartialVCtrpts[i][j + 1] - PartialVCtrpts[i][j]) * (v_degree) / diffv;		
		}
	}

	PartialUUCtrpts = Pw_ctrlpts;
	for (int a = 1; a <= 2; a++)
	{
		for (int i = 0; i <= m - a; i++)
		{
			double diffu = u_knots[i + u_degree + 1] - u_knots[i + a];
			if (diffu < DBL_EPSILON)
			{
				for (int j = 0; j <= n; j++)
				{
					setZero(PartialUUCtrpts[i][j]);
				}
				continue;
			}
			for (int j = 0; j <= n; j++)
			{
				PartialUUCtrpts[i][j] = (PartialUUCtrpts[i + 1][j] - PartialUUCtrpts[i][j]) * (u_degree - a + 1) / diffu;
			}
		}
	}

	PartialUVCtrpts = Pw_ctrlpts;
	for (int i = 0; i <= m - 1; i++)
	{
		double diffu = u_knots[i + u_degree + 1] - u_knots[i + 1];
		if (diffu < DBL_EPSILON)
		{
			for (int j = 0; j <= n; j++)
			{
				setZero(PartialUVCtrpts[i][j]);
			}
			continue;
		}
		for (int j = 0; j <= n; j++)
		{
			PartialUVCtrpts[i][j] = (PartialUVCtrpts[i + 1][j] - PartialUVCtrpts[i][j]) * (u_degree) / diffu;
		}
	}
	for (int i = 0; i <= m - 1; i++)
	{
		for (int j = 0; j <= n - 1; j++)
		{
			double diffv = v_knots[j + v_degree + 1] - v_knots[j + 1];
			if (diffv < DBL_EPSILON) setZero(PartialUVCtrpts[i][j]);
			else PartialUVCtrpts[i][j] = (PartialUVCtrpts[i][j + 1] - PartialUVCtrpts[i][j]) * (v_degree) / diffv;			
		}
	}

	PartialVVCtrpts = Pw_ctrlpts;
	for (int b = 1; b <= 2; b++)
	{
		for (int i = 0; i <= m; i++)
		{
			for (int j = 0; j <= n - b; j++)
			{
				double diffv = v_knots[j + v_degree + 1] - v_knots[j + b];
				if (diffv < DBL_EPSILON) setZero(PartialVVCtrpts[i][j]);
				else PartialVVCtrpts[i][j] = (PartialVVCtrpts[i][j + 1] - PartialVVCtrpts[i][j]) * (v_degree - b + 1) / diffv;				
			}
		}
	}
}

BSplineSurface::BSplineSurface(int udeg, int vdeg, std::vector<double> &uknots, std::vector<double> &vknots, std::vector<std::vector<Point>> &controlpoints)
{
	SetDegree(udeg, vdeg);
	SetKnotsU(uknots);
	SetKnotsV(vknots);
	SetRational(false);
	SetControlPoints(controlpoints);
	weights.clear();
	int m = GetNumOfCtrlptsU() - 1;
	int n = GetNumOfCtrlptsV() - 1;
	UV = Point4(uknots[u_degree], uknots[m + 1], vknots[v_degree], vknots[n + 1]); //2022.03.21

	PartialUCtrpts_ = ctrlpoints;
	for (int i = 0; i <= m - 1; i++)
	{
		double diffu = u_knots[i + u_degree + 1] - u_knots[i + 1];
		if (diffu < DBL_EPSILON)
		{
			for (int j = 0; j <= n; j++)
			{
				setZero(PartialUCtrpts_[i][j]);
			}
			continue;
		}
		for (int j = 0; j <= n; j++)
		{
			PartialUCtrpts_[i][j] = (PartialUCtrpts_[i + 1][j] - PartialUCtrpts_[i][j]) * (u_degree) / diffu;
		}
	}

	PartialVCtrpts_ = ctrlpoints;
	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= n - 1; j++)
		{
			double diffv = v_knots[j + v_degree + 1] - v_knots[j + 1];
			if (diffv < DBL_EPSILON) setZero(PartialVCtrpts_[i][j]);
			else PartialVCtrpts_[i][j] = (PartialVCtrpts_[i][j + 1] - PartialVCtrpts_[i][j]) * (v_degree) / diffv;
		}
	}

	PartialUUCtrpts_ = ctrlpoints;
	for (int a = 1; a <= 2; a++)
	{
		for (int i = 0; i <= m - a; i++)
		{
			double diffu = u_knots[i + u_degree + 1] - u_knots[i + a];
			if (diffu < DBL_EPSILON)
			{
				for (int j = 0; j <= n; j++)
				{
					setZero(PartialUUCtrpts_[i][j]);
				}
				continue;
			}
			for (int j = 0; j <= n; j++)
			{
				PartialUUCtrpts_[i][j] = (PartialUUCtrpts_[i + 1][j] - PartialUUCtrpts_[i][j]) * (u_degree - a + 1) / diffu;
			}
		}
	}

	PartialUVCtrpts_ = ctrlpoints;
	for (int i = 0; i <= m - 1; i++)
	{
		double diffu = u_knots[i + u_degree + 1] - u_knots[i + 1];
		if (diffu < DBL_EPSILON)
		{
			for (int j = 0; j <= n; j++)
			{
				setZero(PartialUVCtrpts_[i][j]);
			}
			continue;
		}
		for (int j = 0; j <= n; j++)
		{
			PartialUVCtrpts_[i][j] = (PartialUVCtrpts_[i + 1][j] - PartialUVCtrpts_[i][j]) * (u_degree) / diffu;
		}
	}
	for (int i = 0; i <= m - 1; i++)
	{
		for (int j = 0; j <= n - 1; j++)
		{
			double diffv = v_knots[j + v_degree + 1] - v_knots[j + 1];
			if (diffv < DBL_EPSILON) setZero(PartialUVCtrpts_[i][j]);
			else PartialUVCtrpts_[i][j] = (PartialUVCtrpts_[i][j + 1] - PartialUVCtrpts_[i][j]) * (v_degree) / diffv;
		}
	}

	PartialVVCtrpts_ = ctrlpoints;
	for (int b = 1; b <= 2; b++)
	{
		for (int i = 0; i <= m; i++)
		{
			for (int j = 0; j <= n - b; j++)
			{
				double diffv = v_knots[j + v_degree + 1] - v_knots[j + b];
				if (diffv < DBL_EPSILON) setZero(PartialVVCtrpts_[i][j]);
				else PartialVVCtrpts_[i][j] = (PartialVVCtrpts_[i][j + 1] - PartialVVCtrpts_[i][j]) * (v_degree - b + 1) / diffv;
			}
		}
	}
}

BSplineSurface::~BSplineSurface()
{

}


const int BSplineSurface::FindSpan(int n, int p, double t, std::vector<double> K)const
{
	if (abs(t - K[n + 1]) < DBL_EPSILON)
		return n;

	int low = p, high = n + 1, mid = (low + high) / 2;
	while (t < K[mid] || t >= K[mid + 1])
	{
		if (t < K[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return(mid);
}

Point BSplineSurface::DeBoor(const double u, const double v) const
{
	if (u < u_knots.front() || u > u_knots.back()) return Point();
	if (v < v_knots.front() || v > v_knots.back()) return Point();

	if (isRational)
	{
		auto Pw = DeBoor(Pw_ctrlpts, u, v);
		return Pw.head(3) / Pw[3];
	}
	else
	{
		return DeBoor(ctrlpoints, u, v);
	}
}

Point BSplineSurface::PartialDerivativeU(const double u, const double v) const
{
	if (u_degree < 1) return Point(0.0, 0.0, 0.0);
	if (isRational)
	{
		auto Pw = DeBoor(Pw_ctrlpts, u, v);
		auto Aw_u = Derivative(PartialUCtrpts, u, v, 1, 0);
		double W_u = Aw_u(3);
		Point A_u = Aw_u.head(3);
		double W = Pw(3);
		Point S = Pw.head(3) / W;

		return (A_u - W_u * S) / W;
	}
	else
	{
		return Derivative(PartialUCtrpts_, u, v, 1, 0);
	}
}

Point BSplineSurface::PartialDerivativeV(const double u, const double v)const
{
	if (v_degree < 1) return Point(0.0, 0.0, 0.0);
	if (isRational)
	{
		auto Pw = DeBoor(Pw_ctrlpts, u, v);
		auto Aw_v = Derivative(PartialVCtrpts, u, v, 0, 1);
		double W_v = Aw_v(3);
		Point A_v = Aw_v.head(3);
		double W = Pw(3);
		Point S = Pw.head(3) / W;

		return (A_v - W_v * S) / W;
	}
	else
	{
		return Derivative(PartialVCtrpts_, u, v, 0, 1);
	}
}

Point BSplineSurface::PartialDerivativeUU(const double u, const double v) const
{
	if (u_degree < 2) return Point(0.0, 0.0, 0.0);
	if (isRational)
	{
		auto Pw = DeBoor(Pw_ctrlpts, u, v);
		auto Aw_u = Derivative(PartialUCtrpts, u, v, 1, 0);
		auto Aw_uu = Derivative(PartialUUCtrpts, u, v, 2, 0);

		double W_uu = Aw_uu(3);
		double W_u = Aw_u(3);
		double W = Pw(3);
		Point A_uu = Aw_uu.head(3);
		Point A_u = Aw_u.head(3);
		Point S = Pw.head(3) / W;
		Point S_u = (A_u - W_u * S) / W;

		return (A_uu - 2 * W_u * S_u - W_uu * S) / W;
	}
	else
	{
		return Derivative(PartialUUCtrpts_, u, v, 2, 0);
	}
}

Point BSplineSurface::PartialDerivativeUV(const double u, const double v)const
{
	if (u_degree < 1 || v_degree < 1) return Point(0.0, 0.0, 0.0);
	if (isRational)
	{
		auto Pw = DeBoor(Pw_ctrlpts, u, v);
		auto Aw_u = Derivative(PartialUCtrpts, u, v, 1, 0);
		auto Aw_v = Derivative(PartialVCtrpts, u, v, 0, 1);
		auto Aw_uv = Derivative(PartialUVCtrpts, u, v, 1, 1);

		double W_uv = Aw_uv(3);
		double W_u = Aw_u(3);
		double W_v = Aw_v(3);
		double W = Pw(3);
		Point A_uv = Aw_uv.head(3);
		Point A_u = Aw_u.head(3);
		Point A_v = Aw_v.head(3);
		Point S = Pw.head(3) / W;
		Point S_u = (A_u - W_u * S) / W;
		Point S_v = (A_v - W_v * S) / W;

		return (A_uv - W_u * S_v - W_v * S_u - W_uv * S) / W;
	}
	else
	{
		return Derivative(PartialUVCtrpts_, u, v, 1, 1);
	}
}

Point BSplineSurface::PartialDerivativeVV(const double u, const double v) const
{
	if (v_degree < 2) return Point(0.0, 0.0, 0.0);
	if (isRational)
	{
		auto Pw = DeBoor(Pw_ctrlpts, u, v);
		auto Aw_v = Derivative(PartialVCtrpts, u, v, 0, 1);
		auto Aw_vv = Derivative(PartialVVCtrpts, u, v, 0, 2);
		double W_vv = Aw_vv(3);
		double W_v = Aw_v(3);
		double W = Pw(3);
		Point A_vv = Aw_vv.head(3);
		Point A_v = Aw_v.head(3);
		Point S = Pw.head(3) / W;
		Point S_v = (A_v - W_v * S) / W;

		return (A_vv - 2 * W_v * S_v - W_vv * S) / W;
	}
	else
	{
		return Derivative(PartialVVCtrpts_, u, v, 0, 2);
	}
}

void BSplineSurface::PrincipalCurvature(const double u, const double v, double &k1, double &k2) const
{
	Eigen::Matrix2d Weingarten, value;
	Point ru = PartialDerivativeU(u, v);
	Point rv = PartialDerivativeV(u, v);
	double E = ru.dot(ru);
	double F = ru.dot(rv);
	double G = rv.dot(rv);
	Point n = (ru.cross(rv)).normalized();
	double L = PartialDerivativeUU(u, v).dot(n);
	double M = PartialDerivativeUV(u, v).dot(n);
	double N = PartialDerivativeVV(u, v).dot(n);
	Weingarten << L * G - M * F, M * E - L * F,
		M * G - N * F, N * E - M * F;
	Eigen::EigenSolver<Eigen::Matrix2d> es(Weingarten);
	value = (es.pseudoEigenvalueMatrix()) / (E*G - pow(F, 2));
	k1 = std::max(std::abs(value(0, 0)), std::abs(value(1, 1)));
	k2 = std::min(std::abs(value(0, 0)), std::abs(value(1, 1)));
}

void BSplineSurface::KnotInsertion(double uv, int k, DIRECTION dir, BSplineSurface & new_surface)
{
	int mp = GetNumOfKnotsU() - 1;
	int mq = GetNumOfKnotsV() - 1;
	int np = GetNumOfCtrlptsU() - 1;
	int nq = GetNumOfCtrlptsV() - 1;

	std::vector<double> new_knotsU;
	std::vector<double> new_knotsV;
	std::vector<std::vector<Point>> new_ctrlpts;

	if (dir == U_DIRECTION)
	{
		assert(uv >= u_knots[0] && uv <= u_knots[mp]);
		if (abs(uv - u_knots[0]) < DBL_EPSILON || abs(uv - u_knots[mp]) < DBL_EPSILON)	return;

		new_knotsU.resize(mp + 1 + k);
		new_knotsV.resize(mq + 1);
		new_ctrlpts.resize(np + 1 + k, std::vector<Point>(nq + 1));

		if (isRational)
		{
			std::vector<std::vector<double>> new_weights(np + 1 + k, std::vector<double>(nq + 1));
			std::vector<std::vector<Point4>> pts(np + 1 + k, std::vector<Point4>(nq + 1));
			std::vector<std::vector<Point4>> Rw(np + 1, std::vector<Point4>(nq + 1));
			for (int j = 0; j <= nq; j++)
			{
				for (int i = 0; i <= np; i++)
				{
					Rw[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
				}
			}

			KnotInsertion(uv, k, U_DIRECTION, u_knots, v_knots, Rw, new_knotsU, new_knotsV, pts);
			for (int j = 0; j <= nq; j++)
			{
				for (int i = 0; i <= np + k; i++)
				{
					new_weights[i][j] = pts[i][j](3);
					new_ctrlpts[i][j] = pts[i][j].head(3) / pts[i][j](3);
				}
			}

			new_surface.SetWeights(new_weights);
		}

		else
		{
			KnotInsertion(uv, k, U_DIRECTION, u_knots, v_knots, ctrlpoints, new_knotsU, new_knotsV, new_ctrlpts);
		}

	}

	else
	{
		assert(uv >= v_knots[0] && uv <= v_knots[mq]);
		if (abs(uv - v_knots[0]) < DBL_EPSILON || abs(uv - v_knots[mq]) < DBL_EPSILON)	return;

		new_knotsU.resize(mp + 1);
		new_knotsV.resize(mq + 1 + k);
		new_ctrlpts.resize(np + 1, std::vector<Point>(nq + 1 + k));

		if (isRational)
		{
			std::vector<std::vector<double>> new_weights(np + 1, std::vector<double>(nq + 1 + k));
			std::vector<std::vector<Point4>> pts(np + 1, std::vector<Point4>(nq + 1 + k));
			std::vector<std::vector<Point4>> Rw(np + 1, std::vector<Point4>(nq + 1));
			for (int j = 0; j <= nq; j++)
			{
				for (int i = 0; i <= np; i++)
				{
					Rw[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
				}
			}

			KnotInsertion(uv, k, V_DIRECTION, u_knots, v_knots, Rw, new_knotsU, new_knotsV, pts);
			for (int j = 0; j <= nq + k; j++)
			{
				for (int i = 0; i <= np; i++)
				{
					new_weights[i][j] = pts[i][j](3);
					new_ctrlpts[i][j] = pts[i][j].head(3) / pts[i][j](3);
				}
			}

			new_surface.SetWeights(new_weights);
		}

		else
		{
			KnotInsertion(uv, k, V_DIRECTION, u_knots, v_knots, ctrlpoints, new_knotsU, new_knotsV, new_ctrlpts);
		}
	}

	new_surface.SetDegree(u_degree, v_degree);
	new_surface.SetRational(isRational);
	new_surface.SetKnotsU(new_knotsU);
	new_surface.SetKnotsV(new_knotsV);
	new_surface.SetControlPoints(new_ctrlpts);
}

void BSplineSurface::Split(double uv, DIRECTION dir, BSplineSurface &left_s, BSplineSurface &right_s)
{
	int mp = GetNumOfKnotsU() - 1;
	int mq = GetNumOfKnotsV() - 1;
	int np = GetNumOfCtrlptsU() - 1;
	int nq = GetNumOfCtrlptsV() - 1;

	if (dir == U_DIRECTION)
	{
		assert(uv >= u_knots.front() && uv <= u_knots.back());
		if (abs(uv - u_knots.front()) < DBL_EPSILON)
		{
			left_s = *this;
			return;
		}
		else if (abs(uv - u_knots.back()) < DBL_EPSILON)
		{
			right_s = *this;
			return;
		}

		int mult = 0;
		int r = FindSpan(np, u_degree, uv, u_knots);
		for (int i = r; i > u_degree; i--)
		{
			if (abs(uv - u_knots[i]) < DBL_EPSILON)
				mult++;
		}
		int k = u_degree - mult;

		BSplineSurface tmp_surface;
		KnotInsertion(uv, k, U_DIRECTION, tmp_surface);
		auto tmp_pts = tmp_surface.GetControlPoints();
		int left_m = r - mult + (u_degree + 1);
		int right_m = mp - (r + 1) + (u_degree + 1);
		int left_n = left_m - (u_degree + 1);
		int right_n = right_m - (u_degree + 1);

		std::vector<std::vector<Point>> left_ctrlpts(left_n + 1, std::vector<Point>(nq + 1));
		std::vector<std::vector<Point>> right_ctrlpts(right_n + 1, std::vector<Point>(nq + 1));
		std::vector<double> left_uknots(left_m + 1);
		std::vector<double> right_uknots(right_m + 1);

		// control points;
		for (int j = 0; j <= nq; j++)
		{
			for (int i = 0; i <= left_n; i++)
			{
				left_ctrlpts[i][j] = tmp_pts[i][j];
			}

			for (int i = 0; i <= right_n; i++)
			{
				right_ctrlpts[i][j] = tmp_pts[i + left_n][j];
			}
		}

		// knots;
		for (int i = 0; i <= r - mult; i++)
		{
			left_uknots[i] = u_knots[i];
		}
		for (int i = r + 1; i <= mp; i++)
		{
			right_uknots[i - (r + 1) + (u_degree + 1)] = u_knots[i];
		}
		for (int i = 1; i <= u_degree + 1; i++)
		{
			left_uknots[i + r - mult] = uv;
			right_uknots[i - 1] = uv;
		}

		// weights;
		if (isRational)
		{
			auto tmp_weights = tmp_surface.GetWeights();
			std::vector<std::vector<double>> left_weights(left_n + 1, std::vector<double>(nq + 1));
			std::vector<std::vector<double>> right_weights(right_n + 1, std::vector<double>(nq + 1));

			for (int j = 0; j <= nq; j++)
			{
				for (int i = 0; i <= left_n; i++)
				{
					left_weights[i][j] = tmp_weights[i][j];
				}

				for (int i = 0; i <= right_n; i++)
				{
					right_weights[i][j] = tmp_weights[i + left_n][j];
				}
			}

			left_s.SetWeights(left_weights);
			right_s.SetWeights(right_weights);
		}

		left_s.SetKnotsU(left_uknots);
		right_s.SetKnotsU(right_uknots);
		left_s.SetKnotsV(v_knots);
		right_s.SetKnotsV(v_knots);
		left_s.SetControlPoints(left_ctrlpts);
		right_s.SetControlPoints(right_ctrlpts);

	}

	else
	{
		assert(uv >= v_knots.front() && uv <= v_knots.back());
		if (abs(uv - v_knots.front()) < DBL_EPSILON)
		{
			left_s = *this;
			return;
		}
		else if (abs(uv - v_knots.back()) < DBL_EPSILON)
		{
			right_s = *this;
			return;
		}

		int mult = 0;
		int r = FindSpan(nq, v_degree, uv, v_knots);
		for (int i = r; i > v_degree; i--)
		{
			if (abs(uv - v_knots[i]) < DBL_EPSILON)
				mult++;
		}
		int k = v_degree - mult;

		BSplineSurface tmp_surface;
		KnotInsertion(uv, k, V_DIRECTION, tmp_surface);
		auto tmp_pts = tmp_surface.GetControlPoints();
		int left_m = r - mult + (v_degree + 1);
		int right_m = mq - (r + 1) + (v_degree + 1);
		int left_n = left_m - (v_degree + 1);
		int right_n = right_m - (v_degree + 1);

		std::vector<std::vector<Point>> left_ctrlpts(np + 1, std::vector<Point>(left_n + 1));
		std::vector<std::vector<Point>> right_ctrlpts(np + 1, std::vector<Point>(right_n + 1));
		std::vector<double> left_vknots(left_m + 1);
		std::vector<double> right_vknots(right_m + 1);

		// control points;
		for (int i = 0; i <= np; i++)
		{
			for (int j = 0; j <= left_n; j++)
			{
				left_ctrlpts[i][j] = tmp_pts[i][j];
			}

			for (int j = 0; j <= right_n; j++)
			{
				right_ctrlpts[i][j] = tmp_pts[i][j + left_n];
			}
		}

		// knots;
		for (int i = 0; i <= r - mult; i++)
		{
			left_vknots[i] = v_knots[i];
		}
		for (int i = r + 1; i <= mq; i++)
		{
			right_vknots[i - (r + 1) + (v_degree + 1)] = v_knots[i];
		}
		for (int i = 1; i <= v_degree + 1; i++)
		{
			left_vknots[i + r - mult] = uv;
			right_vknots[i - 1] = uv;
		}

		// weights;
		if (isRational)
		{
			auto tmp_weights = tmp_surface.GetWeights();
			std::vector<std::vector<double>> left_weights(np + 1, std::vector<double>(left_n + 1));
			std::vector<std::vector<double>> right_weights(np + 1, std::vector<double>(right_n + 1));

			for (int i = 0; i <= np; i++)
			{
				for (int j = 0; j <= left_n; j++)
				{
					left_weights[i][j] = tmp_weights[i][j];
				}

				for (int j = 0; j <= right_n; j++)
				{
					right_weights[i][j] = tmp_weights[i][j + left_n];
				}
			}

			left_s.SetWeights(left_weights);
			right_s.SetWeights(right_weights);
		}

		left_s.SetKnotsU(u_knots);
		right_s.SetKnotsU(u_knots);
		left_s.SetKnotsV(left_vknots);
		right_s.SetKnotsV(right_vknots);
		left_s.SetControlPoints(left_ctrlpts);
		right_s.SetControlPoints(right_ctrlpts);

	}

	right_s.SetDegree(u_degree, v_degree);
	left_s.SetRational(isRational);
	right_s.SetRational(isRational);
}

void BSplineSurface::BSplineToBezier(std::vector<BezierSurface> &bezier_surfaces)
{
	bezier_surfaces.clear();
	int np = ctrlpoints.size() - 1;
	int nq = ctrlpoints[0].size() - 1;

	int num_ubezier; // number of bezier segments(-1) in u dierection
	int num_vbezier; // number of bezier segments(-1) in v direction
	std::vector<double> USpan;
	std::vector<double> VSpan;

	if (isRational)
	{
		std::vector<std::vector<Point4>> Pw(np + 1, std::vector<Point4>(nq + 1));
		std::vector<std::vector<std::vector<Point4>>> Qw;
		for (int i = 0; i <= np; i++)
		{
			for (int j = 0; j <= nq; j++)
			{
				Pw[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
			}
		}

		BSplineToBezier(U_DIRECTION, u_degree, u_knots, Pw, num_ubezier, USpan, Qw);

		std::vector<std::vector<std::vector<Point4>>> pts_vec;
		std::vector<std::vector<Point>> points(u_degree + 1, std::vector<Point>(v_degree + 1));
		std::vector<std::vector<double>> w(u_degree + 1, std::vector<double>(v_degree + 1));
		for (int i = 0; i <= num_ubezier; i++)
		{
			BSplineToBezier(V_DIRECTION, v_degree, v_knots, Qw[i], num_vbezier, VSpan, pts_vec);
			for (int j = 0; j <= num_vbezier; j++)
			{
				for (int k = 0; k <= u_degree; k++)
				{
					for (int l = 0; l <= v_degree; l++)
					{
						w[k][l] = pts_vec[j][k][l](3);
						points[k][l] = pts_vec[j][k][l].head(3) / w[k][l];
					}
				}
				BezierSurface surface_(u_degree, v_degree, USpan[i], USpan[i + 1], VSpan[j], VSpan[j + 1], w, points);
				bezier_surfaces.push_back(surface_);
			}
		}

	}

	else
	{
		std::vector<std::vector<std::vector<Point>>> Qw;
		BSplineToBezier(U_DIRECTION, u_degree, u_knots, ctrlpoints, num_ubezier, USpan, Qw);

		std::vector<std::vector<std::vector<Point>>> pts_vec;
		for (int i = 0; i <= num_ubezier; i++)
		{
			BSplineToBezier(V_DIRECTION, v_degree, v_knots, Qw[i], num_vbezier, VSpan, pts_vec);

			for (int j = 0; j <= num_vbezier; j++)
			{
				BezierSurface surface_(u_degree, v_degree, USpan[i], USpan[i + 1], VSpan[j], VSpan[j + 1], pts_vec[j]);
				bezier_surfaces.push_back(surface_);
			}
		}
	}
}


void BSplineSurface::Derivative(const double &u, const double &v, int &rx, int &ry, std::vector<double> &du, std::vector<double> &dv) const
{
	rx = 0;
	for (int i = u_degree; i + 1 <= ctrlpoints.size(); i++)
	{
		if (u >= u_knots[i])
		{
			rx = i;
		}
		if (u < u_knots[i + 1])
		{
			break;
		}
	}
	ry = 0;
	for (int i = v_degree; i + 1 <= ctrlpoints[0].size(); i++)
	{
		if (v >= v_knots[i])
		{
			ry = i;
		}
		if (v < v_knots[i + 1])
		{
			break;
		}
	}

	du.resize((u_degree + 1) * (v_degree + 1));
	dv.resize((u_degree + 1) * (v_degree + 1));
	std::vector<std::vector<double>> indicator(ctrlpoints.size(), std::vector<double>(ctrlpoints[0].size(), 0.0));
	for (size_t i = 0; i < u_degree + 1; i++)
	{
		for (size_t j = 0; j < v_degree + 1; j++)
		{
			indicator[i + rx - u_degree][j + ry - v_degree] = 1.0;
			BSplineSurface derivativeu;
			derivativeu.SetDegree(u_degree - 1, v_degree);
			derivativeu.SetNumberControlPoints(indicator.size() - 1, indicator[0].size());
			std::vector<std::vector<double>> ductrlpts(indicator.size() - 1, \
				std::vector<double>(indicator[0].size()));
			for (size_t i1 = 1; i1 < indicator.size(); i1++)
			{
				for (size_t j1 = 0; j1 < indicator[0].size(); j1++)
				{
					ductrlpts[i1 - 1][j1] = (indicator[i1][j1] - indicator[i1 - 1][j1]) / (u_knots[i1 + u_degree] - u_knots[i1]) * u_degree;
				}
			}
			BSplineSurface derivativev;
			derivativev.SetDegree(u_degree, v_degree - 1);
			derivativev.SetNumberControlPoints(indicator.size(), indicator[0].size() - 1);
			std::vector<std::vector<double>> dvctrlpts(indicator.size(), std::vector<double>(indicator[0].size() - 1));
			for (size_t i1 = 0; i1 < indicator.size(); i1++)
			{
				for (size_t j1 = 1; j1 < indicator[0].size(); j1++)
				{
					dvctrlpts[i1][j1 - 1] = (indicator[i1][j1] - indicator[i1][j1 - 1]) / (v_knots[j1 + v_degree] - v_knots[j1]) * v_degree;
				}
			}
			du[i * (v_degree + 1) + j] = derivativeu.DeBoor<double>(ductrlpts, u, v);
			dv[i * (v_degree + 1) + j] = derivativev.DeBoor<double>(dvctrlpts, u, v);
			indicator[i + rx - u_degree][j + ry - v_degree] = 0.0;
		}
	}
}

void BSplineSurface::PreprocessData(void)
{
	std::vector<Point> datanew;
	std::vector<double> lens;
	double totallen = 0.0;
	for (size_t i = 0; i + 1 < datapoints.size(); i++)
	{
		double len = (datapoints[i] - datapoints[i + 1]).norm();
		if (len > 1e-12)
		{
			datanew.push_back(datapoints[i]);
			lens.push_back(len);
			totallen += len;
		}
	}
	double len = (datapoints[datapoints.size() - 1] - datapoints[0]).norm();
	if (len > 1e-12)
	{
		datanew.push_back(datapoints[datapoints.size() - 1]);
		lens.push_back(len);
		totallen += len;
	}
	double curvelen = totallen / 4.0;
	size_t i1 = 0, i2 = 0, i3 = 0;
	double len1 = 0.0;
	for (size_t i = 0; i < datanew.size(); i++)
	{
		len1 += lens[i];
		if (len1 > curvelen)
		{

		}
	}
}

void BSplineSurface::Parameterization(std::vector<double>& t)
{
}

void BSplineSurface::Regression(void)
{
}

// 根据contorlpointnet把u,v归一化到[0,1]*[0,1] (uniform)
void BSplineSurface::UpdateKnots(void)
{
	size_t m = ctrlpoints.size();
	if (m == 0) return;
	size_t n = ctrlpoints[0].size();
	assert(m > u_degree && n > v_degree);
	u_knots.resize(m + u_degree + 1);
	v_knots.resize(n + v_degree + 1);
	size_t dtx = m - u_degree;
	size_t dty = n - v_degree;
	for (size_t i = 0; i < u_degree; i++)
	{
		u_knots[i] = 0.0;
	}
	for (size_t i = u_degree; i < m; i++)
	{
		u_knots[i] = static_cast<double>(i - u_degree) / dtx;
	}
	for (size_t i = m; i <= m + u_degree; i++)
	{
		u_knots[i] = 1.0;
	}
	for (size_t i = 0; i < v_degree; i++)
	{
		v_knots[i] = 0.0;
	}
	for (size_t i = v_degree; i < n; i++)
	{
		v_knots[i] = static_cast<double>(i - v_degree) / dty;
	}
	for (size_t i = n; i <= n + v_degree; i++)
	{
		v_knots[i] = 1.0;
	}
}

void BSplineSurface::UpdataControl_Auxiliary()
{
	int N_x_control = ctrlpoints.size();
	int N_y_control = ctrlpoints[0].size();

	control_list.resize(N_x_control*N_y_control);

	int mv = 0;
	for (int i = 0; i < N_x_control; i++)
	{
		for (int j = 0; j < N_y_control; j++)
		{
			control_list[mv] = OpenMesh::Vec3d(ctrlpoints[i][j][0], ctrlpoints[i][j][1], 0);
			mv++;
		}
	}

	Set_Face_control(Face_Control_index, N_x_control, N_y_control);

}

void BSplineSurface::Set_Face_control(std::vector<std::vector<int>>& face_control, int m_size, int n_size)
{
	face_control.resize((m_size - u_degree)*(n_size - v_degree));
	double N_face_x = m_size - u_degree; double N_face_y = n_size - v_degree;
	double u_c = 0.5 / N_face_x; double v_c = 0.5 / N_face_y;
	int mt = 0;
	int rx, ry;
	for (int i = 0; i < int(N_face_x); i++)
	{
		for (int j = 0; j < int(N_face_y); j++)
		{
			double u = i / N_face_x + u_c;
			double v = j / N_face_y + v_c;
			UV_to_Cell(u, v, rx, ry);
			std::vector<int>& FC = face_control[mt];
			FC.resize((u_degree + 1)*(v_degree + 1));

			int mf = 0;
			for (int fi = 0; fi < u_degree + 1; fi++)
			{
				for (int fj = 0; fj < v_degree + 1; fj++)
				{
					FC[mf] = (rx + fi)*m_size + (ry + fj);
					mf++;
				}
			}

			mt++;
		}
	}
}

void BSplineSurface::UV_to_Cell(const double & u, const double & v, int & rx, int & ry)
{
	rx = 0;
	for (size_t i = u_degree; i + 1 <= ctrlpoints.size(); i++)
	{
		if (u >= u_knots[i])
		{
			rx = static_cast<int>(i);
		}
		if (u < u_knots[i + 1])
		{
			break;
		}
	}
	ry = 0;
	for (size_t i = v_degree; i + 1 <= ctrlpoints[0].size(); i++)
	{
		if (v >= v_knots[i])
		{
			ry = static_cast<int>(i);
		}
		if (v < v_knots[i + 1])
		{
			break;
		}
	}

	rx = rx - u_degree;  // rx-kx+1：u涉及到的控制点最小index
	ry = ry - v_degree;
}

void BSplineSurface::UV_to_Cell(const double & u, const double & v, int & c_id)
{
	int rx = 0, ry = 0;
	UV_to_Cell(u, v, rx, ry);

	c_id = rx + (ctrlpoints.size() - u_degree)*ry;
}

void BSplineSurface::Compute_All_Basic(const double u, const double v, std::vector<double>& all_base)
{
	int m_size = ctrlpoints.size();
	int n_size = ctrlpoints[0].size();

	control.resize(m_size, std::vector<double>(n_size, 0));
	all_base.resize(m_size*n_size);

	int mb = 0;
	for (int i = 0; i < m_size; i++)
	{
		for (int j = 0; j < n_size; j++)
		{
			control[i][j] = 1;
			all_base[mb] = DeBoor(control, u, v);
			control[i][j] = 0;
			mb++;
		}
	}
}


bool BSplineSurface::SaveControlPoints(const std::string & filename) const
{
	std::ofstream ofs(filename);
	if (!ofs.is_open())
	{
		std::cerr << "Error: cannot save control points from file " << filename << std::endl;
		return false;
	}
	ofs.precision(16);
	ofs << ctrlpoints.size() << " " << ctrlpoints[0].size() << std::endl;
	for (const auto &pts : ctrlpoints)
	{
		for (const auto &pt : pts)
		{
			ofs << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
		}
	}
	ofs.close();
	return true;
}

bool BSplineSurface::LoadControlPoints(const std::string & filename)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open())
	{
		std::cerr << "Error: cannot load control points from file " << filename << std::endl;
		return false;
	}
	size_t m, n;
	ifs >> m >> n;
	SetDegree(3, 3);	// default degree(3,3)
	SetNumberControlPoints(m, n);
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			ifs >> ctrlpoints[i][j][0] >> ctrlpoints[i][j][1] >> ctrlpoints[i][j][2];
		}
	}
	ifs.close();
	return true;
}

void BSplineSurface::SetNumberControlPoints(const size_t & mp, const size_t & np)
{
	ctrlpoints.clear();
	ctrlpoints.resize(mp, std::vector<Point>(np, Point(0.0, 0.0, 0.0)));
	UpdateKnots();
}

void BSplineSurface::SetControlPoints(Eigen::VectorXd& rx, Eigen::VectorXd& ry)
{
	int m_size = ctrlpoints.size(); int n_size = ctrlpoints[0].size();

	int mt = 0;
	for (int i = 0; i < m_size; i++)
	{
		for (int j = 0; j < n_size; j++)
		{
			ctrlpoints[i][j] = Point(rx[mt], ry[mt], 0);

			mt++;
		}
	}
}

void BSplineSurface::SetControlPoints(Eigen::VectorXd& rx, Eigen::VectorXd& ry, Eigen::VectorXd& rz)
{
	int m_size = ctrlpoints.size(); int n_size = ctrlpoints[0].size();

	int mt = 0;
	for (int i = 0; i < m_size; i++)
	{
		for (int j = 0; j < n_size; j++)
		{
			ctrlpoints[i][j] = Point(rx[mt], ry[mt], rz[mt]);
			mt++;
		}
	}
}

void BSplineSurface::SetControlPoints(double ui, double vi, double x, double y, double z)
{
	ctrlpoints[ui][vi] = Point(x, y, z);
}


void BSplineSurface::Compute_All_Boundary_Basic(const double u, const double v, std::vector<double>& re)
{
	size_t N_1 = ctrlpoints.size(); size_t N_2 = ctrlpoints[0].size();
	control.clear();
	control.resize(N_1, std::vector<double>(N_2, 0));

	int i, j;
	re.clear();	re.reserve(2 * N_1 + 2 * N_2);
	//===========================================================
	//first
	i = 0;
	for (int j = 0; j < N_2; j++)
	{
		control[i][j] = 1;
		double r = DeBoor(control, u, v);
		re.push_back(r);
		control[i][j] = 0;
	}
	//===========================================================
	//second
	j = 0;
	for (int i = 1; i < N_1 - 1; i++)
	{
		control[i][j] = 1;
		double r = DeBoor(control, u, v);
		re.push_back(r);
		control[i][j] = 0;
	}
	//===========================================================
	//3
	j = N_2 - 1;
	for (int i = 1; i < N_1 - 1; i++)
	{
		control[i][j] = 1;
		double r = DeBoor(control, u, v);
		re.push_back(r);
		control[i][j] = 0;
	}
	//===========================================================
	//4
	i = N_1 - 1;
	for (int j = 0; j < N_2; j++)
	{
		control[i][j] = 1;
		double r = DeBoor(control, u, v);
		re.push_back(r);
		control[i][j] = 0;
	}
}

void BSplineSurface::Set_Boundary_ControlPoints(Eigen::VectorXd& rx, Eigen::VectorXd& ry)
{
	size_t N_1 = ctrlpoints.size(); size_t N_2 = ctrlpoints[0].size();

	for (int i = 0; i < N_1; i++)
	{
		for (int j = 0; j < N_2; j++)
		{
			ctrlpoints[i][j] = Point(0, 0, 0);
		}
	}

	int i, j;
	int time = 0;
	//===========================================================
	//first
	i = 0;
	for (int j = 0; j < N_2; j++)
	{
		ctrlpoints[i][j] = Point(rx[time], ry[time], 0);
		time++;
	}
	//===========================================================
	//second
	j = 0;
	for (int i = 1; i < N_1 - 1; i++)
	{
		ctrlpoints[i][j] = Point(rx[time], ry[time], 0);
		time++;
	}
	//===========================================================
	//3
	j = N_2 - 1;
	for (int i = 1; i < N_1 - 1; i++)
	{
		ctrlpoints[i][j] = Point(rx[time], ry[time], 0);
		time++;
	}
	//===========================================================
	//4
	i = N_1 - 1;
	for (int j = 0; j < N_2; j++)
	{
		ctrlpoints[i][j] = Point(rx[time], ry[time], 0);
		time++;
	}
}

void BSplineSurface::scale_control(double farc)
{
	int m = ctrlpoints.size(); int n = ctrlpoints[0].size();
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			ctrlpoints[i][j] = ctrlpoints[i][j] * farc;
		}
	}
}

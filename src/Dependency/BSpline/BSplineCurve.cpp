#include "BSplineCurve.h"
#include <ctime>

BSplineCurve::BSplineCurve(void)
	:degree(-1)
{
	SetRational(false);
	knots.clear();
	weights.clear();
	ctrlpoints.clear();
}


BSplineCurve::BSplineCurve(int deg, std::vector<double> knots, std::vector<double> weights, std::vector<Point> ctrlpoints)
{
	SetDegree(deg);
	SetRational(true);
	SetKnots(knots);
	SetWeights(weights);
	SetControlPoints(ctrlpoints);
}

BSplineCurve::BSplineCurve(int deg, std::vector<double> knots, std::vector<Point> ctrlpoints)
{
	SetDegree(deg);
	SetRational(false);
	bool re = isRational;
	SetKnots(knots);
	SetControlPoints(ctrlpoints);
}

BSplineCurve::~BSplineCurve()
{

}


/* Determine the knot span index */
const int BSplineCurve::FindSpan(int n, int p, double t, std::vector<double> K)const 
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


Point BSplineCurve::DeBoor(const double t) const
{
	if (isRational)
	{
		int n = ctrlpoints.size() - 1;
		std::vector<Point4> Pw_ctrlpts(n + 1);
		for (int i = 0; i <= n; i++)
		{
			Pw_ctrlpts[i] << ctrlpoints[i] * weights[i], weights[i];
		}
		auto Pw = DeBoor(Pw_ctrlpts, t);

		return Pw.head(3) / Pw[3];
	}
	else
	{
		return DeBoor(ctrlpoints, t);
	}
}

Point BSplineCurve::Derivative(const double t) const
{
	if (degree < 1) return Point(0.0, 0.0, 0.0);
	if (isRational)
	{
		int n = ctrlpoints.size() - 1;
		std::vector<Point4> Pw_ctrlpts(n + 1);
		for (int i = 0; i <= n; i++)
		{
			Pw_ctrlpts[i] << ctrlpoints[i] * weights[i], weights[i];
		}
		auto Pw = DeBoor(Pw_ctrlpts, t);
		auto Aw_t = Derivative(Pw_ctrlpts, t, 1);
		double W_t = Aw_t(3);
		Point A_t = Aw_t.head(3);
		double W = Pw(3);
		Point S = Pw.head(3) / W;

		return (A_t - W_t * S) / W;
	}
	else
	{
		return Derivative(ctrlpoints, t, 1);
	}
}

Point BSplineCurve::Derivative2(const double t) const
{
	if (degree < 2) return Point(0.0, 0.0, 0.0);
	if (isRational)
	{
		int n = ctrlpoints.size() - 1;
		std::vector<Point4> Pw_ctrlpts(n + 1);
		for (int i = 0; i <= n; i++)
		{
			Pw_ctrlpts[i] << ctrlpoints[i] * weights[i], weights[i];
		}
		auto Pw = DeBoor(Pw_ctrlpts, t);
		auto Aw_t = Derivative(Pw_ctrlpts, t, 1);
		auto Aw_tt = Derivative(Pw_ctrlpts, t, 2);

		double W_tt = Aw_tt(3);
		double W_t = Aw_t(3);
		double W = Pw(3);
		Point A_tt = Aw_tt.head(3);
		Point A_t = Aw_t.head(3);
		Point S = Pw.head(3) / W;
		Point S_t = (A_t - W_t * S) / W;

		return (A_tt - 2 * W_t * S_t - W_tt * S) / W;
	}
	else
	{
		return Derivative(ctrlpoints, t, 2);
	}
}


void BSplineCurve::KnotInsertion(double t, int k, BSplineCurve & new_curve)
{
	int m = GetNumOfKnots() - 1;
	int n = GetNumOfCtrlpts() - 1;
	assert(t > knots[0] && t < knots[m]);
	if (abs(t - knots[0]) < DBL_EPSILON || abs(t - knots[m]) < DBL_EPSILON)
		return;

	std::vector<double> new_knots(m + 1 + k);
	std::vector<Point> new_ctrlpts(n + 1 + k);

	if (isRational)
	{
		std::vector<double> new_weights(n + 1 + k);
		std::vector<Point4> pts(n + 1 + k);
		std::vector<Point4> Rw(n + 1);
		for (int i = 0; i <= n; i++)
		{
			Rw[i] << ctrlpoints[i] * weights[i], weights[i];
		}

		KnotInsertion<Point4>(t, k, knots, Rw, new_knots, pts);
		for (int i = 0; i <= n + k; i++)
		{
			new_weights[i] = pts[i](3);
			new_ctrlpts[i] = pts[i].head(3) / pts[i](3);
		}

		new_curve.SetWeights(new_weights);
	}

	else
	{
		KnotInsertion<Point>(t, k, knots, ctrlpoints, new_knots, new_ctrlpts);
	}

	new_curve.SetDegree(degree);
	new_curve.SetRational(isRational);
	new_curve.SetKnots(new_knots);
	new_curve.SetControlPoints(new_ctrlpts);
}

void BSplineCurve::Split(double t, BSplineCurve &left_c, BSplineCurve &right_c)
{
	int m = GetNumOfKnots() - 1;
	int n = GetNumOfCtrlpts() - 1;
	assert(t >= knots[0] && t <= knots[m]);
	if (abs(t - knots[0]) < DBL_EPSILON)		
	{
		left_c = *this;
		return;
	}
	else if (abs(t - knots[m]) < DBL_EPSILON)
	{
		right_c = *this;
		return;
	}

	int mult = 0;
	int r = FindSpan(n, degree, t, knots);
	for (int i = r; i > degree; i--)
	{
		if (abs(t - knots[i]) < DBL_EPSILON)
			mult++;
	}
	int k = degree - mult;

	BSplineCurve tmp_curve;
	KnotInsertion(t, k, tmp_curve);
	auto tmp_pts = tmp_curve.GetControlPoints();
	int left_m = r - mult + (degree + 1);
	int right_m = m - (r + 1) + (degree + 1);
	int left_n = left_m - (degree + 1);
	int right_n = right_m - (degree + 1);
	std::vector<Point> left_ctrlpts(tmp_pts.begin(), tmp_pts.begin() + left_n + 1);
	std::vector<Point> right_ctrlpts(tmp_pts.end() - right_n - 1, tmp_pts.end());

	std::vector<double> left_knots(left_m + 1);
	std::vector<double> right_knots(right_m + 1);
	for (int i = 0; i <= r - mult; i++)
	{
		left_knots[i] = knots[i];
	}
	for (int i = r + 1; i <= m; i++)
	{
		right_knots[i - (r + 1) + (degree + 1)] = knots[i];
	}
	for (int i = 1; i <= degree + 1; i++)
	{
		left_knots[i + r - mult] = t;
		right_knots[i - 1] = t;
	}
	
	if (isRational)
	{
		auto tmp_weights = tmp_curve.GetWeights();
		std::vector<double> left_weights(tmp_weights.begin(), tmp_weights.begin() + left_n + 1);
		std::vector<double> right_weights(tmp_weights.end() - right_n - 1, tmp_weights.end());
		left_c.SetWeights(left_weights);
		right_c.SetWeights(right_weights);
	}

	left_c.SetKnots(left_knots);
	right_c.SetKnots(right_knots);
	left_c.SetDegree(degree);
	right_c.SetDegree(degree);
	left_c.SetRational(isRational);
	right_c.SetRational(isRational);
	left_c.SetControlPoints(left_ctrlpts);
	right_c.SetControlPoints(right_ctrlpts);
}

void BSplineCurve::BSplineToBezier(std::vector<BezierCurve> &bezier_curves)
{
	int n = GetNumOfCtrlpts() - 1;
	std::vector<double> KnotsSpan;

	int num_bezier; // number of bezier curves

	if (isRational)
	{
		std::vector<Point4> Pw(n + 1);
		std::vector<std::vector<Point4>> pts_vec;
		for (int i = 0; i <= n; i++)
		{
			Pw[i] << ctrlpoints[i] * weights[i], weights[i];
		}
		BSplineToBezier(degree, knots, Pw, num_bezier, KnotsSpan, pts_vec);

		bezier_curves.resize(num_bezier + 1);
		for (int i = 0; i <= num_bezier; i++)
		{
			std::vector<double> w(degree + 1);
			std::vector<Point> pts(degree + 1);
			for (int j = 0; j <= degree; j++)
			{
				w[j] = pts_vec[i][j](3);
				pts[j] = pts_vec[i][j].head(3) / w[j];
			}
			bezier_curves[i].SetControlPoints(pts);
			bezier_curves[i].SetWeights(w);
		}

	}
	else
	{
		std::vector<std::vector<Point>> pts_vec;
		BSplineToBezier(degree, knots, ctrlpoints, num_bezier, KnotsSpan, pts_vec);

		bezier_curves.resize(num_bezier + 1);
		for (int i = 0; i <= num_bezier; i++)
		{
			std::vector<Point> pts(degree + 1);
			for (int j = 0; j <= degree; j++)
			{
				pts[j] = pts_vec[i][j];
			}
			bezier_curves[i].SetControlPoints(pts);
		}

	}

	for (int i = 0; i <= num_bezier; i++)
	{
		bezier_curves[i].SetDegree(degree);
		bezier_curves[i].SetParameterSpan(KnotsSpan[i], KnotsSpan[i+1]);
		bezier_curves[i].SetRational(isRational);
	}

}


// 参数t按datapoints的弦长归一化到[0,1]
void BSplineCurve::Parameterization(std::vector<double>& t)
{
	t.resize(datapoints.size());
	double totallen = 0.0;
	for (size_t i = 0; i + 1 < datapoints.size(); i++)
	{
		totallen += (datapoints[i] - datapoints[i + 1]).norm();
	}
	double currentlen = 0.0;
	t[0] = 0.0;
	for (size_t i = 0; i + 1 < datapoints.size(); i++)
	{
		currentlen += (datapoints[i] - datapoints[i + 1]).norm();
		t[i + 1] = currentlen / totallen;
	}
}

// t∈[0,1](uniform parameterization)
void BSplineCurve::UpdateKnots(void)
{
	size_t n = ctrlpoints.size();
	if (n == 0)	return;
	assert(n > degree);
	knots.resize(n + degree + 1);
	size_t dt = n - degree;
	for (size_t i = 0; i < degree; i++)
	{
		knots[i] = 0.0;
	}
	for (size_t i = degree; i < n; i++)
	{
		knots[i] = static_cast<double>(i - degree) / dt;
	}
	for (size_t i = n; i <= n + degree; i++)
	{
		knots[i] = 1.0;
	}
}

void BSplineCurve::Regression(void)
{
	if (datapoints.empty())
	{
		return;
	}
	std::vector<double> t;
	// Parameterization(t);
	Eigen::MatrixXd A(datapoints.size(), ctrlpoints.size());
	
	for (size_t i = 0; i < datapoints.size(); i++)
	{
		for (size_t j = 0; j < ctrlpoints.size(); j++)
		{
			size_t r = 0;
			for (size_t k = degree; k + 1 <= ctrlpoints.size(); k++)
			{
				if (t[i] >= knots[k])
				{
					r = k;
				}
				if (t[i] < knots[k + 1])
				{
					break;
				}
			}
			std::vector<double> d(degree + 1, 0.0);
			if (j + degree >= r && j <= r)
			{
				d[j + degree - r] = 1.0;
				for (size_t k = 1; k <= degree; k++)
				{
					for (size_t l = degree; l >= k; l--)
					{
						double alpha = (t[i] - knots[l + r - degree]) / (knots[l + 1 + r - k] - knots[l + r - degree]);
						d[l] = (1 - alpha) * d[l - 1] + alpha * d[l];
					}
				}
			}
			A(i, j) = d[degree];
		}
	}
	Eigen::VectorXd bx(datapoints.size()), by(datapoints.size());
	for (size_t i = 0; i < datapoints.size(); i++)
	{
		bx(i) = datapoints[i](0);
		by(i) = datapoints[i](1);
	}
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::VectorXd x = svd.solve(bx);
	Eigen::VectorXd y = svd.solve(by);
	for (size_t i = 0; i < ctrlpoints.size(); i++)
	{
		ctrlpoints[i][0] = x(i);
		ctrlpoints[i][1] = y(i);
		ctrlpoints[i][2] = 0.0;
	}
}

#include "BezierCurve.h"
#include <fstream>
#include <iostream>
#include <ctime>

BezierCurve::BezierCurve(): isRational(false)
{
	SetRational(false);
	weights.clear();
	ctrlpoints.clear();
}

BezierCurve::BezierCurve(int deg, double tmin, double tmax, std::vector<double>w, std::vector<Point> controlpoints)
{
	degree = deg;
	t_min = tmin;
	t_max = tmax;
	isRational = true;
	weights = w;
	ctrlpoints = controlpoints;
}

BezierCurve::BezierCurve(int deg, double tmin, double tmax, std::vector<Point> controlpoints)
{
	degree = deg;
	t_min = tmin;
	t_max = tmax;
	isRational = false;
	weights.clear();
	ctrlpoints = controlpoints;
}

BezierCurve::~BezierCurve()
{

}


Point BezierCurve::DeCasteljau(const double t)const
{
	if (isRational)
	{
		std::vector<Point4> Pw_ctrlpts(degree + 1);
		for (int i = 0; i <= degree; i++)
		{
			Pw_ctrlpts[i] << ctrlpoints[i] * weights[i], weights[i];
		}
		auto Pw = DeCasteljau(Pw_ctrlpts, t);

		return Pw.head(3) / Pw[3];
	}
	else
	{
		return DeCasteljau(ctrlpoints, t);
	}
}

Point BezierCurve::Derivative(const double t) const
{
	if (degree < 1) return Point(0.0, 0.0, 0.0);

	if (isRational) // C(t) = A(t) / W(t); C'(t) = (A'(t)-W'(t)*C(t)) / W(t) 
	{	
		std::vector<Point4> Pw_ctrlpts(degree + 1);
		for (int i = 0; i <= degree; i++)
		{
			Pw_ctrlpts[i] << ctrlpoints[i] * weights[i], weights[i];
		}
		
		auto Pw = DeCasteljau(Pw_ctrlpts, t);
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

Point BezierCurve::Derivative2(const double t) const
{	
	if (degree < 2) return Point(0.0, 0.0, 0.0);

	if (isRational) // C(t) = A(t) / W(t); C''(t) = (A''(t)-2W'(t)*C'(t) - W''(t)C(t)) / W(t) 
	{
		std::vector<Point4> Pw_ctrlpts(degree + 1);
		for (int i = 0; i <= degree; i++)
		{
			Pw_ctrlpts[i] << ctrlpoints[i] * weights[i], weights[i];
		}

		auto Pw = DeCasteljau(Pw_ctrlpts, t);
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


void BezierCurve::Split(double t, BezierCurve &left_c, BezierCurve &right_c)
{
	assert(t > t_min && t < t_max);
	left_c.SetDegree(degree);
	left_c.SetParameterSpan(t_min, t);
	left_c.SetRational(isRational);

	right_c.SetDegree(degree);
	right_c.SetParameterSpan(t, t_max);
	right_c.SetRational(isRational);

	double t_ = (t - t_min) / (t_max - t_min);
	std::vector<Point> left_pts(degree + 1), right_pts(degree + 1);
	left_pts[0] = ctrlpoints[0];
	right_pts[degree] = ctrlpoints[degree];

	if (isRational)
	{		
		std::vector<Point4> Rw(degree + 1);
		std::vector<double> left_weights(degree + 1), right_weights(degree + 1);
		for (int i = 0; i <= degree; i++)
		{
			Rw[i] << ctrlpoints[i] * weights[i], weights[i];
		}

		for (int j = 1; j <= degree; j++)
		{
			for (int i = degree; i >= j; i--)
			{
				Rw[i] = t_ * Rw[i] + (1.0 - t_) * Rw[i - 1];
			}

			left_weights[j] = Rw[j][3];
			right_weights[degree - j] = Rw[degree][3];
			left_pts[j] = Rw[j].head(3) / left_weights[j];
			right_pts[degree - j] = Rw[degree].head(3) / right_weights[degree - j];
		}

		left_weights[0] = weights[0];
		right_weights[degree] = weights[degree];
		left_c.SetWeights(left_weights);
		right_c.SetWeights(right_weights);
	}

	else
	{
		std::vector<Point> Rw(ctrlpoints);
		for (int j = 1; j <= degree; j++)
		{
			for (int i = degree; i >= j; i--)
			{
				Rw[i] = t_ * Rw[i] + (1.0 - t_) * Rw[i - 1];
			}

			left_pts[j] = Rw[j];
			right_pts[degree - j] = Rw[degree];
		}		
	}
	
	left_c.SetControlPoints(left_pts);
	right_c.SetControlPoints(right_pts);
}


bool BezierCurve::SaveControlPoints(const std::string & filename) const
{
	std::ofstream ofs(filename);
	if (!ofs.is_open())
	{
		std::cerr << "Error: cannot save control points from file " << filename << std::endl;
		return false;
	}
	ofs.precision(16);
	ofs << ctrlpoints.size() << std::endl;
	for (const auto &pt : ctrlpoints)
	{
		ofs << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
	}
	ofs.close();
	return true;
}

bool BezierCurve::LoadControlPoints(const std::string & filename)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open())
	{
		std::cerr << "Error: cannot load control points from file " << filename << std::endl;
		return false;
	}
	size_t n;
	ifs >> n;
	SetDegree(n);	
	for (size_t i = 0; i <= n; i++)
	{
		ifs >> ctrlpoints[i][0] >> ctrlpoints[i][1] >> ctrlpoints[i][2];
	}
	ifs.close();
	return true;
}

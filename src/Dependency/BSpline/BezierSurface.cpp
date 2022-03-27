#include "BezierSurface.h"
#include <fstream>
#include <iostream>
#include <ctime>

BezierSurface::BezierSurface() : isRational(false)
{
	SetRational(false);
	weights.clear();
	ctrlpoints.clear();
}

BezierSurface::BezierSurface(int u_deg, int v_deg, double umin, double umax, double vmin, double
	vmax, std::vector<std::vector<double>> w, std::vector<std::vector<Point>> controlpoints)
{
	SetDegree(u_deg, v_deg);
	SetParameterSpan(umin, umax, vmin, vmax);
	SetRational(true);
	SetWeights(w);
	SetControlPoints(controlpoints);
	UV = Point4(umin, umax, vmin, vmax);  //2022.03.21
}

BezierSurface::BezierSurface(int u_deg, int v_deg, double umin, double umax, double vmin, double vmax, std::vector<std::vector<Point>> controlpoints)
{
	SetDegree(u_deg, v_deg);
	SetParameterSpan(umin, umax, vmin, vmax);
	SetRational(false);
	SetControlPoints(controlpoints);
	UV = Point4(umin, umax, vmin, vmax);  //2022.03.21
}

BezierSurface::~BezierSurface()
{

}


Point BezierSurface::DeCasteljau(const double u, const double v) const
{
	if (isRational)
	{
		std::vector<std::vector<Point4>> Pw_ctrlpts(u_degree + 1, std::vector<Point4>(v_degree + 1));
		for (int i = 0; i <= u_degree; i++)	// Pw_(i,j)^0
		{
			for (int j = 0; j <= v_degree; j++)
			{
				Pw_ctrlpts[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
			}
		}
		auto Pw = DeCasteljau(Pw_ctrlpts, u, v);

		return Pw.head(3) / Pw[3];
	}

	else
	{
		return DeCasteljau(ctrlpoints, u, v);
	}
}

Point BezierSurface::PartialDerivativeU(const double u, const double v) const
{
	if (u_degree < 1)
		return Point(0.0, 0.0, 0.0);

	if (isRational) // S(u,v) = A(u,v)/W(u,v); S_u = (A_u-W_u*S)/W
	{
		std::vector<std::vector<Point4>> Pw_ctrlpts(u_degree + 1, std::vector<Point4>(v_degree + 1));
		for (int i = 0; i <= u_degree; i++)
		{
			for (int j = 0; j <= v_degree; j++)
			{
				Pw_ctrlpts[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
			}
		}
		auto Pw = DeCasteljau(Pw_ctrlpts, u, v);
		auto Aw_u = Derivative(Pw_ctrlpts, u, v, 1, 0);
		double W_u = Aw_u(3);
		Point A_u = Aw_u.head(3);
		double W = Pw(3);
		Point S = Pw.head(3) / W;

		return (A_u - W_u * S) / W;
	}

	else
	{		
		return Derivative(ctrlpoints, u, v, 1, 0);
	}
}

Point BezierSurface::PartialDerivativeV(const double u, const double v)const
{
	if (v_degree < 1)
		return Point(0.0, 0.0, 0.0);

	if (isRational) // S(u,v) = A(u,v)/W(u,v); S_v = (A_v-W_v*S)/W
	{
		std::vector<std::vector<Point4>> Pw_ctrlpts(u_degree + 1, std::vector<Point4>(v_degree + 1));
		for (int i = 0; i <= u_degree; i++)
		{
			for (int j = 0; j <= v_degree; j++)
			{
				Pw_ctrlpts[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
			}
		}
		auto Pw = DeCasteljau(Pw_ctrlpts, u, v);
		auto Aw_v = Derivative(Pw_ctrlpts, u, v, 0, 1);
		double W_v = Aw_v(3);
		Point A_v = Aw_v.head(3);
		double W = Pw(3);
		Point S = Pw.head(3) / W;

		return (A_v - W_v * S) / W;
	}

	else
	{
		return Derivative(ctrlpoints, u, v, 0, 1);
	}
}

Point BezierSurface::PartialDerivativeUU(const double u, const double v) const
{
	if (u_degree < 2)
		return Point(0.0, 0.0, 0.0);

	if (isRational)// S(u,v) = A(u,v)/W(u,v); S_uu = (A_uu - 2*W_u*S_u - W_uu*S)/W
	{
		std::vector<std::vector<Point4>> Pw_ctrlpts(u_degree + 1, std::vector<Point4>(v_degree + 1));
		for (int i = 0; i <= u_degree; i++)
		{
			for (int j = 0; j <= v_degree; j++)
			{
				Pw_ctrlpts[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
			}
		}

		auto Pw = DeCasteljau(Pw_ctrlpts, u, v);
		auto Aw_u = Derivative(Pw_ctrlpts, u, v, 1, 0);
		auto Aw_uu = Derivative(Pw_ctrlpts, u, v, 2, 0);

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
		return Derivative(ctrlpoints, u, v, 2, 0);
	}
}

Point BezierSurface::PartialDerivativeUV(const double u, const double v)const
{
	if (u_degree < 1 || v_degree < 1)
		return Point(0.0, 0.0, 0.0);

	if (isRational)// S(u,v) = A(u,v)/W(u,v); S_uv = (A_uv - W_u*S_v - W_v*S_u - W_uv*S)/W
	{
		std::vector<std::vector<Point4>> Pw_ctrlpts(u_degree + 1, std::vector<Point4>(v_degree + 1));
		for (int i = 0; i <= u_degree; i++)
		{
			for (int j = 0; j <= v_degree; j++)
			{
				Pw_ctrlpts[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
			}
		}

		auto Pw = DeCasteljau(Pw_ctrlpts, u, v);
		auto Aw_u = Derivative(Pw_ctrlpts, u, v, 1, 0);
		auto Aw_v = Derivative(Pw_ctrlpts, u, v, 0, 1);
		auto Aw_uv = Derivative(Pw_ctrlpts, u, v, 1, 1);

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
		return Derivative(ctrlpoints, u, v, 1, 1);
	}
}

Point BezierSurface::PartialDerivativeVV(const double u, const double v) const 
{
	if (v_degree < 2)
		return Point(0.0, 0.0, 0.0);

	if (isRational)// S(u,v) = A(u,v)/W(u,v); S_vv = (A_vv - 2*W_v*S_v - W_vv*S)/W
	{
		std::vector<std::vector<Point4>> Pw_ctrlpts(u_degree + 1, std::vector<Point4>(v_degree + 1));
		for (int i = 0; i <= u_degree; i++)
		{
			for (int j = 0; j <= v_degree; j++)
			{
				Pw_ctrlpts[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
			}
		}

		auto Pw = DeCasteljau(Pw_ctrlpts, u, v);
		auto Aw_v = Derivative(Pw_ctrlpts, u, v, 0, 1);
		auto Aw_vv = Derivative(Pw_ctrlpts, u, v, 0, 2);

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
		return Derivative(ctrlpoints, u, v, 0, 2);
	}
}


void BezierSurface::Split(double uv, DIRECTION direction, BezierSurface &left_s, BezierSurface &right_s)
{
	left_s.SetDegree(u_degree, v_degree);
	left_s.SetRational(isRational);

	right_s.SetDegree(u_degree, v_degree);
	right_s.SetRational(isRational);

	std::vector<std::vector<Point>> left_pts(u_degree + 1, std::vector<Point>(v_degree + 1));
	std::vector<std::vector<Point>> right_pts(u_degree + 1, std::vector<Point>(v_degree + 1));

	if (direction == U_DIRECTION)
	{
		assert(uv >= u_min && uv <= u_max);
		double alpha = (uv - u_min) / (u_max - u_min);
		left_s.SetParameterSpan(u_min, uv, v_min, v_max);
		right_s.SetParameterSpan(uv, u_max, v_min, v_max);

		if (isRational)
		{
			std::vector<std::vector<Point4>> Rw(u_degree + 1, std::vector<Point4>(v_degree + 1));
			std::vector<std::vector<double>> left_weights(u_degree + 1, std::vector<double>(v_degree + 1));
			std::vector<std::vector<double>> right_weights(u_degree + 1, std::vector<double>(v_degree + 1));

			for (int j = 0; j <= v_degree; j++)
			{
				for (int i = 0; i <= u_degree; i++)
				{
					Rw[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
				}
				left_weights[0][j] = weights[0][j];
				right_weights[u_degree][j] = weights[u_degree][j];

				left_pts[0][j] = ctrlpoints[0][j];
				right_pts[u_degree][j] = ctrlpoints[u_degree][j];				
			}

			for (int row = 0; row <= v_degree; row++)
			{
				for (int k = 1; k <= u_degree; k++)
				{
					for (int i = u_degree; i >= k; i--)
					{
						Rw[i][row] = alpha * Rw[i][row] + (1.0 - alpha) * Rw[i - 1][row];
					}

					left_weights[k][row] = Rw[k][row](3);
					right_weights[u_degree - k][row] = Rw[u_degree][row](3);

					left_pts[k][row] = Rw[k][row].head(3) / Rw[k][row](3);
					right_pts[u_degree - k][row] = Rw[u_degree][row].head(3) / Rw[u_degree][row](3);
				}
			}

			left_s.SetWeights(left_weights);
			right_s.SetWeights(right_weights);
		}

		else
		{
			std::vector<std::vector<Point>> Rw(ctrlpoints);
			for (int row = 0; row <= v_degree; row++)
			{
				left_pts[0][row] = ctrlpoints[0][row];
				right_pts[u_degree][row] = ctrlpoints[u_degree][row];

				for (int k = 1; k <= u_degree; k++)
				{
					for (int i = u_degree; i >= k; i--)
					{
						Rw[i][row] = alpha * Rw[i][row] + (1.0 - alpha) * Rw[i - 1][row];
					}

					left_pts[k][row] = Rw[k][row];
					right_pts[u_degree - k][row] = Rw[u_degree][row];
				}
			}
		}
	}

	else
	{
		assert(uv >= v_min && uv <= v_max);
		double alpha = (uv - v_min) / (v_max - v_min);
		left_s.SetParameterSpan(u_min, u_max, v_min, uv);
		right_s.SetParameterSpan(u_min, u_max, uv, v_max);

		if (isRational)
		{
			std::vector<std::vector<Point4>> Rw(u_degree + 1, std::vector<Point4>(v_degree + 1));
			std::vector<std::vector<double>> left_weights(u_degree + 1, std::vector<double>(v_degree + 1));
			std::vector<std::vector<double>> right_weights(u_degree + 1, std::vector<double>(v_degree + 1));

			for (int i = 0; i <= u_degree; i++)
			{
				for (int j = 0; j <= v_degree; j++)
				{
					Rw[i][j] << ctrlpoints[i][j] * weights[i][j], weights[i][j];
				}

				left_pts[i][0] = ctrlpoints[i][0];
				left_weights[i][0] = weights[i][0];

				right_pts[i][v_degree] = ctrlpoints[i][v_degree];
				right_weights[i][v_degree] = weights[i][v_degree];
			}
			
			for (int col = 0; col <= u_degree; col++)
			{
				for (int l = 1; l <= v_degree; l++)
				{
					for (int j = v_degree; j >= l; j--)
					{
						Rw[col][j] = alpha * Rw[col][j] + (1.0 - alpha) * Rw[col][j - 1];
					}
					left_weights[col][l] = Rw[col][l](3);
					right_weights[col][v_degree - l] = Rw[col][v_degree](3);

					left_pts[col][l] = Rw[col][l].head(3) / Rw[col][l](3);
					right_pts[col][v_degree - l] = Rw[col][v_degree].head(3) / Rw[col][v_degree](3);
				}
			}

			left_s.SetWeights(left_weights);
			right_s.SetWeights(right_weights);
		}

		else
		{
			std::vector<std::vector<Point>> Rw(ctrlpoints);
			for (int col = 0; col <= u_degree; col++)
			{
				left_pts[col][0] = Rw[col][0];
				right_pts[col][v_degree] = Rw[col][v_degree];

				for (int l = 1; l <= v_degree; l++)
				{
					for (int j = v_degree; j >= l; j--)
					{
						Rw[col][j] = alpha * Rw[col][j] + (1.0 - alpha) * Rw[col][j - 1];
					}

					left_pts[col][l] = Rw[col][l];
					right_pts[col][v_degree - l] = Rw[col][v_degree];
				}
			}
		}	
	}

	left_s.SetControlPoints(left_pts);
	right_s.SetControlPoints(right_pts);
}


bool BezierSurface::SaveControlPoints(const std::string & filename) const
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

bool BezierSurface::LoadControlPoints(const std::string & filename)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open())
	{
		std::cerr << "Error: cannot load control points from file " << filename << std::endl;
		return false;
	}
	size_t m, n;
	ifs >> m >> n;
	SetDegree(m, n);
	for (size_t i = 0; i <= n; i++)
	{
		for (size_t j = 0; j <= m; j++)
		{
			ifs >> ctrlpoints[i][j][0] >> ctrlpoints[i][j][1] >> ctrlpoints[i][j][2];
		}		
	}
	ifs.close();
	return true;
}

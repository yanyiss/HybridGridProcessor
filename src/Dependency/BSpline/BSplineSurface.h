#pragma once
#ifndef BSPLINE_S
#define BSPLINE_S
#include <Eigen\Eigen>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <vector>
#include "BezierSurface.h"
typedef Eigen::Vector3d Point;
typedef Eigen::Vector4d Point4;

class BSplineSurface
{
public:
	BSplineSurface();
	BSplineSurface(int udeg, int vdeg, std::vector<double> &uknots, std::vector<double> &vknots, std::vector<std::vector<double>> &w, std::vector<std::vector<Point>> &controlpoints);
	BSplineSurface(int udeg, int vdeg, std::vector<double> &uknots, std::vector<double> &vknots, std::vector<std::vector<Point>> &controlpoints);
	~BSplineSurface();

	void SetDegree(const int degu, const int degv) { u_degree = degu; v_degree = degv; }
	void SetRational(const bool rational) { isRational = rational; }
	void SetKnotsU(const std::vector<double> &uknots) { u_knots.clear(); u_knots = uknots; }
	void SetKnotsV(const std::vector<double> &vknots) { v_knots.clear(); v_knots = vknots; }
	void SetWeights(const std::vector<std::vector<double>> &w) 
	{ 
		weights.clear(); 
		weights = w; 
		isRational = true;
	}
	void SetControlPoints(const std::vector<std::vector<Point>> &points)
	{
		ctrlpoints.clear(); ctrlpoints = points;
	}
	void SetDataPoints(const std::vector<Point> &points) { datapoints.clear(); datapoints = points; }
	
	const int GetDegreeU(void) const { return u_degree; }
	const int GetDegreeV(void) const { return v_degree; }
	const bool GetRational(void) const { return isRational; }
	const int GetNumOfKnotsU(void) const { return u_knots.size(); }
	const int GetNumOfKnotsV(void) const { return v_knots.size(); }
	const int GetNumOfCtrlptsU(void) const { return ctrlpoints.size(); }
	const int GetNumOfCtrlptsV(void) const { return ctrlpoints[0].size(); }
	const std::vector<double> GetKnotsU(void) const { return u_knots; }
	const std::vector<double> GetKnotsV(void) const { return v_knots; }	
	const std::vector<std::vector<double>> GetWeights(void) const { return weights; }
	const std::vector<std::vector<Point>> GetControlPoints(void) const { return ctrlpoints; }
	const std::vector<Point> & DataPoints(void) const { return datapoints; }
	const int FindSpan(int n, int p, double u, std::vector<double> K) const;

	
	template<typename T>
	T DeBoor(std::vector<std::vector<T>> controlpoints, const double u, const double v) const
	{
		if (controlpoints.empty()) return T();
		if (u < u_knots.front() || u > u_knots.back()) return T();
		if (v < v_knots.front() || v > v_knots.back()) return T();
		//assert(u > u_knots.front() && u < u_knots.back());
		//assert(v > v_knots.front() && v < v_knots.back());

		int m = controlpoints.size() - 1;
		int n = controlpoints[0].size() - 1;
		int ru = FindSpan(m, u_degree, u, u_knots);
		int rv = FindSpan(n, v_degree, v, v_knots);

		std::vector<std::vector<T>> d(u_degree + 1, std::vector<T>(v_degree + 1));
		for (int i = 0; i <= u_degree; i++)	// P_(i,j)^0
		{
			for (int j = 0; j <= v_degree; j++)
			{
				d[i][j] = controlpoints[i + ru - u_degree][j + rv - v_degree];
			}
		}

		for (int j = 0; j <= v_degree; j++)	// DeBoor in u direction
		{
			for (int k = 1; k <= u_degree; k++)
			{
				for (int i = u_degree; i >= k; i--)
				{
					double alpha = (u - u_knots[i + ru - u_degree]) / (u_knots[i + 1 + ru - k] - u_knots[i + ru - u_degree]);
					d[i][j] = (1 - alpha) * d[i - 1][j] + alpha * d[i][j];
				}
			}
		}

		for (int l = 1; l <= v_degree; l++)	// DeBoor in v direction
		{
			for (int j = v_degree; j >= l; j--)
			{
				double alpha = (v - v_knots[j + rv - v_degree]) / (v_knots[j + 1 + rv - l] - v_knots[j + rv - v_degree]);
				d[u_degree][j] = (1 - alpha) * d[u_degree][j - 1] + alpha * d[u_degree][j];
			}
		}
		return d[u_degree][v_degree];
	}

	template<typename T>
	T Derivative(std::vector<std::vector<T>> controlpoints, const double u, const double v, const int k, const int l) const
	{
		if (controlpoints.empty()) return T();
		if (u_degree < k || v_degree < l)
			return T();

		if (u < u_knots.front() || u > u_knots.back()) return T();
		if (v < v_knots.front() || v > v_knots.back()) return T();
		int m = controlpoints.size() - 1;
		int n = controlpoints[0].size() - 1;

		std::vector<std::vector<T>> d(controlpoints);
		for (int a = 1; a <= k; a++)
		{
			for (int i = 0; i <= m - a; i++)
			{
				for (int j = 0; j <= n; j++)
				{
					d[i][j] = (d[i + 1][j] - d[i][j]) * (u_degree - a + 1) / (u_knots[i + u_degree + 1] - u_knots[i + a]);
				}
			}
		}

		for (int b = 1; b <= l; b++)
		{
			for (int i = 0; i <= m - k; i++)
			{
				for (int j = 0; j <= n - b; j++)
				{
					d[i][j] = (d[i][j + 1] - d[i][j]) * (v_degree - b + 1) / (v_knots[j + v_degree + 1] - v_knots[j + b]);
				}
			}
		}

		std::vector<std::vector<T>> dev_ctrlpts(m + 1 - k, std::vector<T>(n + 1 - l));
		for (int i = 0; i <= m - k; i++)
		{
			for (int j = 0; j <= n - l; j++)
			{
				dev_ctrlpts[i][j] = d[i][j];
			}
		}
		std::vector<double> du_knots(u_knots.begin() + k, u_knots.end() - k);
		std::vector<double> dv_knots(v_knots.begin() + l, v_knots.end() - l);

		BSplineSurface dev_surface;
		dev_surface.SetDegree(u_degree - k, v_degree - l);
		dev_surface.SetKnotsU(du_knots);
		dev_surface.SetKnotsV(dv_knots);

		return dev_surface.DeBoor(dev_ctrlpts, u, v);
	}

	Point operator()(const double u, const double v) const { return DeBoor(u, v); }
	Point DeBoor(const double u, const double v) const;
	Point PartialDerivativeU(const double u, const double v) const;
	Point PartialDerivativeV(const double u, const double v)const;
	Point PartialDerivativeUU(const double u, const double v) const;
	Point PartialDerivativeUV(const double u, const double v)const;
	Point PartialDerivativeVV(const double u, const double v) const;

	template<typename T>
	void KnotInsertion(double uv, int k, DIRECTION dir, const std::vector<double> &U, const std::vector<double> &V, const std::vector<std::vector<T>> &Pw, std::vector<double> &Unew, std::vector<double> &Vnew, std::vector<std::vector<T>> &Qw)
	{
		int mp = U.size() - 1;
		int mq = V.size() - 1;
		int np = Pw.size() - 1;
		int nq = Pw[0].size() - 1;

		if (dir == U_DIRECTION)
		{
			int mult = 0;
			int r = FindSpan(np, u_degree, uv, U);
			for (int i = r; i > u_degree; i--)
			{
				if (abs(uv - U[i]) < DBL_EPSILON)
					mult++;
			}

			Unew.resize(mp + 1 + k);
			Qw.resize(np + 1 + k, std::vector<T>(nq + 1));

			/* Load new knot vector */
			Vnew = V;
			for (int i = 0; i <= r; i++) Unew[i] = U[i];
			for (int i = 1; i <= k; i++) Unew[r + i] = uv;
			for (int i = r + 1; i <= mp; i++) Unew[i + k] = U[i];

			/* Save the alphas */
			std::vector<std::vector<double>> alpha(u_degree - mult, std::vector<double>(k));
			for (int j = 1; j <= k; j++)
			{
				int L = r - u_degree + j;
				for (int i = 0; i <= u_degree - j - mult; i++)
					alpha[i][j - 1] = (uv - U[i + L]) / (U[i + r + 1] - U[i + L]);
			}

			/* New Control Points */
			for (int row = 0; row <= nq; row++)
			{
				/* Save unaltered control points */
				for (int i = 0; i <= r - u_degree; i++) Qw[i][row] = Pw[i][row];
				for (int i = r - mult; i <= np; i++) Qw[i + k][row] = Pw[i][row];

				/* Load auxiliary control points. */
				std::vector<T> Rw(u_degree - mult + 1);
				for (int i = 0; i <= u_degree - mult; i++) Rw[i] = Pw[i + r - u_degree][row];
				for (int j = 1; j <= k; j++)
				{
					for (int i = 0; i <= u_degree - mult - j; i++)
					{
						Rw[i] = alpha[i][j - 1] * Rw[i + 1] + (1 - alpha[i][j - 1]) * Rw[i];
					}
					Qw[r - u_degree + j][row] = Rw[0];
					Qw[r - mult + k - j][row] = Rw[u_degree - mult - j];
				}

				/* Load the remaining control points */
				for (int i = 1; i < u_degree - mult - k; i++)
				{
					Qw[r - u_degree + k + i][row] = Rw[i];
				}
			}

		}

		else
		{
			int mult = 0;
			int r = FindSpan(nq, v_degree, uv, V);
			for (int i = r; i > v_degree; i--)
			{
				if (abs(uv - V[i]) < DBL_EPSILON)
					mult++;
			}

			Vnew.resize(mq + 1 + k);
			Qw.resize(np + 1, std::vector<T>(nq + 1 + k));

			/* Load new knot vector */
			Unew = U;
			for (int i = 0; i <= r; i++) Vnew[i] = V[i];
			for (int i = 1; i <= k; i++) Vnew[r + i] = uv;
			for (int i = r + 1; i <= mq; i++) Vnew[i + k] = V[i];

			/* Save the alphas */
			std::vector<std::vector<double>> alpha(v_degree - mult, std::vector<double>(k));
			for (int j = 1; j <= k; j++)
			{
				int L = r - v_degree + j;
				for (int i = 0; i <= v_degree - j - mult; i++)
					alpha[i][j - 1] = (uv - V[i + L]) / (V[i + r + 1] - V[i + L]);
			}

			/* New Control Points */
			for (int col = 0; col <= np; col++)
			{
				/* Save unaltered control points */
				for (int j = 0; j <= r - v_degree; j++) Qw[col][j] = Pw[col][j];
				for (int j = r - mult; j <= nq; j++) Qw[col][j + k] = Pw[col][j];

				/* Load auxiliary control points. */
				std::vector<T> Rw(v_degree + 1 - mult);
				for (int i = 0; i <= v_degree - mult; i++) Rw[i] = Pw[col][i + r - v_degree];
				for (int j = 1; j <= k; j++)
				{
					for (int i = 0; i <= v_degree - mult - j; i++)
					{
						Rw[i] = alpha[i][j - 1] * Rw[i + 1] + (1 - alpha[i][j - 1]) * Rw[i];
					}
					Qw[col][r - v_degree + j] = Rw[0];
					Qw[col][r - mult + k - j] = Rw[v_degree - mult - j];
				}

				/* Load the remaining control points */
				for (int i = 1; i < v_degree - mult - k; i++)
				{
					Qw[col][r - v_degree + k + i] = Rw[i];
				}
			}

		}

	}

	template<typename T>
	void BSplineToBezier(DIRECTION dir, int deg, std::vector<double> &Knot, std::vector<std::vector<T>> &Pw, int &nb, std::vector<double> &ParamSpan, std::vector<std::vector<std::vector<T>>> &Qw)
	{
		int m = Knot.size() - 1;
		int np = Pw.size() - 1;
		int nq = Pw[0].size() - 1;
		nb = 0;
		Qw.clear();
		ParamSpan.clear();

		int a = deg, b = deg + 1;
		ParamSpan.push_back(Knot[a]);

		if (dir == U_DIRECTION)
		{
			std::vector<std::vector<T>> tmp(deg + 1, std::vector<T>(nq + 1));
			for (int row = 0; row <= nq; row++)
			{
				for (int i = 0; i <= deg; i++)
					tmp[i][row] = Pw[i][row];
			}
			Qw.push_back(tmp);

			while (b < m)
			{
				int multi;
				int i = b;
				while (b < m && ((Knot[b + 1] - Knot[b]) < DBL_EPSILON))
					b++;
				multi = b - i + 1;


				std::vector<std::vector<T>> Qw_(deg + 1, std::vector<T>(nq + 1));
				if (multi < deg)
				{
					/* Compute alphas */
					std::vector<double> alphas(deg - multi);
					double numerator = Knot[b] - Knot[a];
					for (int j = 0; j < deg - multi; j++)
						alphas[j] = numerator / (Knot[a + multi + 1 + j] - Knot[a]);

					/* Insert Knot[b] r times */
					int r = deg - multi;
					for (int j = 1; j <= r; j++)
					{
						int save = r - j;
						int s = multi + j;
						for (int k = deg; k >= s; k--)
						{
							double alpha = alphas[k - s];
							for (int row = 0; row <= nq; row++)
								Qw[nb][k][row] = alpha * Qw[nb][k][row] + (1.0 - alpha) * Qw[nb][k - 1][row];
						}

						if (b < m)
							for (int row = 0; row <= nq; row++)
								Qw_[save][row] = Qw[nb][deg][row]; // control points of next segment
					}
				}
				ParamSpan.push_back(Knot[b]);

				if (b < m)
				{
					for (int i = deg - multi; i <= deg; i++)
						for (int row = 0; row <= nq; row++)
							Qw_[i][row] = Pw[b - deg + i][row];

					a = b;
					b = b + 1;
					Qw.push_back(Qw_);
					nb++;
				}
			}

		}

		else
		{
			std::vector<std::vector<T>> tmp(np + 1, std::vector<T>(deg + 1));
			for (int col = 0; col <= np; col++)
			{
				for (int i = 0; i <= deg; i++)
					tmp[col][i] = Pw[col][i];
				//Qw[nb][col][i] = Pw[col][i];
			}
			Qw.push_back(tmp);

			while (b < m)
			{
				int multi;
				int i = b;
				while (b < m && ((Knot[b + 1] - Knot[b]) < DBL_EPSILON))
					b++;
				multi = b - i + 1;

				std::vector<std::vector<T>> Qw_(np + 1, std::vector<T>(deg + 1));
				if (multi < deg)
				{
					/* Compute alphas */
					std::vector<double> alphas(deg - multi);
					double numerator = Knot[b] - Knot[a];
					for (int j = 0; j < deg - multi; j++)
						alphas[j] = numerator / (Knot[a + multi + 1 + j] - Knot[a]);

					/* Insert Knot[b] r times */
					int r = deg - multi;
					for (int j = 1; j <= r; j++)
					{
						int save = r - j;
						int s = multi + j;
						for (int k = deg; k >= s; k--)
						{
							double alpha = alphas[k - s];
							for (int col = 0; col <= np; col++)
								Qw[nb][col][k] = alpha * Qw[nb][col][k] + (1.0 - alpha) * Qw[nb][col][k - 1];
						}

						if (b < m)
							for (int col = 0; col <= np; col++)
								Qw_[col][save] = Qw[nb][col][deg]; // control points of next segment
					}
				}
				ParamSpan.push_back(Knot[b]);

				if (b < m)
				{
					for (int i = deg - multi; i <= deg; i++)
						for (int col = 0; col <= np; col++)
							Qw_[col][i] = Pw[col][b - deg + i];

					a = b;
					b = b + 1;
					Qw.push_back(Qw_);
					nb++;
				}
			}

		}
	}


	void KnotInsertion(double uv, int k, DIRECTION dir, BSplineSurface & new_surface);
	void Split(double uv, DIRECTION dir, BSplineSurface &left_s, BSplineSurface &right_s);
	void BSplineToBezier(std::vector<BezierSurface> &bezier_surfaces);

private:
	int u_degree;
	int v_degree;
	bool isRational;
	std::vector<double> u_knots;
	std::vector<double> v_knots;
	std::vector<std::vector<double>> weights;
	std::vector<std::vector<Point>> ctrlpoints;

	std::vector<Point> datapoints;
	std::vector<std::vector<double>> control; // 用于计算 基函数在特定点uv的值


public:
	void Derivative(const double &u, const double &v, int &rx, int &ry, std::vector<double> &du, std::vector<double> &dv) const;
	void PreprocessData(void);
	void Regression(void);
	void Parameterization(std::vector<double> &t);
	

	std::vector<std::vector<int>> Face_Control_index;
	std::vector<OpenMesh::Vec3d> control_list;

	void UpdataControl_Auxiliary();

	void Set_Face_control(std::vector<std::vector<int>>& face_control, int m_size, int n_size);
	void UV_to_Cell(const double & u, const double & v, int & rx, int & ry);
	void UV_to_Cell(const double & u, const double & v, int & c_id);

	void Compute_All_Basic(const double u, const double v, std::vector<double>& all_base);
	void Compute_All_Boundary_Basic(const double u, const double v, std::vector<double>& re);

	bool SaveControlPoints(const std::string & filename) const;
	bool LoadControlPoints(const std::string & filename);
	void SetNumberControlPoints(const size_t &mp, const size_t &np);
	void SetControlPoints(Eigen::VectorXd& rx, Eigen::VectorXd& ry);
	void SetControlPoints(Eigen::VectorXd& rx, Eigen::VectorXd& ry, Eigen::VectorXd& rz);
	void SetControlPoints(double ui, double vi, double x, double y, double z);
	void Set_Boundary_ControlPoints(Eigen::VectorXd& rx, Eigen::VectorXd& ry);

	void scale_control(double farc);
	void UpdateKnots(void);
};
#endif
/* 提取边界曲线 */
/* BSplineToBezier()还可以加速，U方向split完之后，得到的一系列Bezier-Bspline曲面，继续在V方向split，alpha和parameter span 重复计算 */
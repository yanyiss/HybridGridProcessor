#include "SquaredDistanceFunc_Surface.h"

SquaredDistanceFunc_Surface::SquaredDistanceFunc_Surface()
{
	isRational = false;
	weights.clear();
	coeffs.clear();
}

SquaredDistanceFunc_Surface::~SquaredDistanceFunc_Surface()
{

}

SquaredDistanceFunc_Surface::SquaredDistanceFunc_Surface(int u_deg, int v_deg, double umin, double umax, double vmin, double vmax,
	std::vector<std::vector<double>> w, std::vector<std::vector<double>> coefficients)
{
	SetDegree(u_deg, v_deg);
	SetParameterSpan(umin, umax, vmin, vmax);
	SetRational(true);
	SetWeights(w);
	SetCoeffs(coefficients);
}

SquaredDistanceFunc_Surface::SquaredDistanceFunc_Surface(int u_deg, int v_deg, double umin, double umax, double vmin, double vmax, std::vector<std::vector<double>> coefficients)
{
	SetDegree(u_deg, v_deg);
	SetParameterSpan(umin, umax, vmin, vmax);
	SetRational(false);
	weights.clear();
	SetCoeffs(coefficients);
}

SquaredDistanceFunc_Surface::SquaredDistanceFunc_Surface(BezierSurface &bezier_, Point P)
{
	// degree
	int p = bezier_.GetDegreeU();
	int q = bezier_.GetDegreeV();
	int deg_X = 2 * p;
	int deg_Y = 2 * q;
	SetDegree(deg_X, deg_Y);

	// parameter span
	bezier_.GetParameterSpan(u_min, u_max, v_min, v_max);

	// rational
	isRational = bezier_.GetRational();

	// weights and coefficients
	weights.clear();
	coeffs.clear();
	coeffs.resize(u_degree + 1, std::vector<double>(v_degree + 1));
	std::vector<std::vector<Point>> ctrlpts = bezier_.GetControlPoints();
	for (int i = 0; i <= p; i++)
	{
		for (int j = 0; j <= q; j++)
		{
			ctrlpts[i][j] -= P;
		}
	}
		
	std::vector<double> binominal_p = Binomial(p);
	std::vector<double> binominal_q = Binomial(q);
	std::vector<double> binominal_2p = Binomial(u_degree);
	std::vector<double> binominal_2q = Binomial(v_degree);

	if (isRational)
	{
		weights.resize(u_degree + 1, std::vector<double>(v_degree + 1));
		std::vector<std::vector<double>> w = bezier_.GetWeights();
		for (int k = 0; k <= u_degree; k++)
		{
			for (int l = 0; l <= v_degree; l++)
			{
				double temp_coeff = 0.0;
				double temp_weight = 0.0;
				for (int i = std::max(0, k - p); i <= std::min(k, p); i++)
				{
					for (int j = std::max(0, l - q); j <= std::min(l, q); j++)
					{
						double tmp = weights[i][j] * weights[k - i][l - j] * binominal_p[i] * binominal_p[k - i] * binominal_q[j] * binominal_q[l - j];
						temp_weight += tmp;
						temp_coeff += ctrlpts[i][j].dot(ctrlpts[k - i][l - j]) * tmp;
					}
				}
				weights[k][l] = temp_weight / binominal_2p[k] / binominal_2q[l];
				coeffs[k][l] = temp_coeff / binominal_2p[k] / binominal_2q[l];
			}
		}

	}

	else
	{
		for (int k = 0; k <= u_degree; k++)
		{
			for (int l = 0; l <= v_degree; l++)
			{
				double temp_coeff = 0.0;
				for (int i = std::max(0, k - p); i <= std::min(k, p); i++)
				{
					for (int j = std::max(0, l - q); j <= std::min(l, q); j++)
					{
						temp_coeff += ctrlpts[i][j].dot(ctrlpts[k - i][l - j]) * binominal_p[i] * binominal_p[k - i] * binominal_q[j] * binominal_q[l - j];
					}
				}
				coeffs[k][l] = temp_coeff / binominal_2p[k] / binominal_2q[l];
			}
		}

	}
}


double SquaredDistanceFunc_Surface::DeCasteljau(const std::vector<std::vector<double>> &coefficients, int degu, int degv, const double u, const double v) const
{
	if (u < u_min || u > u_max) return 0.0;
	if (v < v_min || v > v_max) return 0.0;

	double u_ = (u - u_min) / (u_max - u_min);
	double v_ = (v - v_min) / (v_max - v_min);

	std::vector<std::vector<double>> d(coefficients);
	for (int j = 0; j <= degv; j++) // P_(0,vj)^(u_degree)
	{
		for (int k = 1; k <= degu; k++)
		{
			for (int i = 0; i <= degu - k; i++)
			{
				d[i][j] = (1 - u_) * d[i][j] + u_ * d[i + 1][j];
			}
		}
	}

	for (int l = 1; l <= degv; l++)
	{
		for (int j = 0; j <= degv - l; j++)
		{
			d[0][j] = (1 - v_) * d[0][j] + v_ * d[0][j + 1];
		}
	}

	return d[0][0];
}

double SquaredDistanceFunc_Surface::Derivative(const std::vector<std::vector<double>> &coefficients, const double u, const double v, const int k, const int l) const
{
	// (k, l)th derivative
	if (u < u_min || u > u_max) return 0.0;
	if (v < v_min || v > v_max) return 0.0;

	std::vector<std::vector<double>> d(coefficients);
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

	std::vector<std::vector<double>> dev(u_degree + 1 - k, std::vector<double>(v_degree + 1 - l));
	for (int i = 0; i <= u_degree - k; i++)
	{
		for (int j = 0; j <= v_degree - l; j++)
		{
			dev[i][j] = d[i][j];
		}
	}

	return DeCasteljau(dev, u_degree - k, v_degree - l, u, v);
}

double SquaredDistanceFunc_Surface::DeCasteljau(const double u, const double v) const
{
	if (isRational)
	{
		return DeCasteljau(coeffs, u_degree, v_degree, u, v) / DeCasteljau(weights, u_degree, v_degree, u, v);
	}
	else
	{
		return DeCasteljau(coeffs, u_degree, v_degree, u, v);
	}
}

double SquaredDistanceFunc_Surface::PartialDerivativeU(const double u, const double v) const
{
	if (isRational) // F(u,v) = f(u,v)/W(u,v); F_u = (f_u*W-W_u*f)/W^2
	{
		double f = DeCasteljau(coeffs, u_degree, v_degree, u, v);
		double W = DeCasteljau(weights, u_degree, v_degree, u, v);
		double f_u = Derivative(coeffs, u, v, 1, 0);
		double W_u = Derivative(weights, u, v, 1, 0);

		return (f_u * W - W_u * f) / pow(W, 2);
	}

	else
	{
		return Derivative(coeffs, u, v, 1, 0);
	}
}

double SquaredDistanceFunc_Surface::PartialDerivativeV(const double u, const double v)const
{
	if (isRational) // F(u,v) = f(u,v)/W(u,v); F_u = (f_u*W-W_u*f)/W^2
	{
		double f = DeCasteljau(coeffs, u_degree, v_degree, u, v);
		double W = DeCasteljau(weights, u_degree, v_degree, u, v);
		double f_v = Derivative(coeffs, u, v, 0, 1);
		double W_v = Derivative(weights, u, v, 0, 1);

		return (f_v * W - W_v * f) / pow(W, 2);
	}

	else
	{
		return Derivative(coeffs, u, v, 0, 1);
	}
}

double SquaredDistanceFunc_Surface::PartialDerivativeUU(const double u, const double v) const
{
	if (isRational) // F(u,v) = f(u,v)/W(u,v); F_uu = [W*(f_uu*W - f*W_uu) - 2*W_u*(f_u*W- f*W_u)]/ W^3 
	{
		double f = DeCasteljau(coeffs, u_degree, v_degree, u, v);
		double f_u = Derivative(coeffs, u, v, 1, 0);
		double f_uu = Derivative(coeffs, u, v, 2, 0);

		double W = DeCasteljau(weights, u_degree, v_degree, u, v);	
		double W_u = Derivative(weights, u, v, 1, 0);		
		double W_uu = Derivative(weights, u, v, 2, 0);

		return (W * (f_uu * W - f * W_uu) - 2 * W_u * (f_u * W - f * W_u)) / pow(W, 3);
	}

	else
	{
		return Derivative(coeffs, u, v, 2, 0);
	}
}

double SquaredDistanceFunc_Surface::PartialDerivativeUV(const double u, const double v)const
{
	if (isRational) // F(u,v) = f(u,v)/W(u,v); F_uv = [W*(f_uu*W + f_u*W_v -f_v*W_u - fW_uu) - 2*W_v*(f_u*W- f*W_u)]/ W^3 
	{
		double f = DeCasteljau(coeffs, u_degree, v_degree, u, v);
		double f_u = Derivative(coeffs, u, v, 1, 0);
		double f_v = Derivative(coeffs, u, v, 0, 1);
		double f_uv = Derivative(coeffs, u, v, 1, 1);

		double W = DeCasteljau(weights, u_degree, v_degree, u, v);	
		double W_u = Derivative(weights, u, v, 1, 0);	
		double W_v = Derivative(weights, u, v, 0, 1);		
		double W_uv = Derivative(weights, u, v, 1, 1);

		return (W * (f_uv * W + f_u * W_v - f_v * W_u - f * W_uv) - 2 * W_v * (f_u * W - f * W_u)) / pow(W, 3);
	}

	else
	{
		return Derivative(coeffs, u, v, 1, 1);
	}
}

double SquaredDistanceFunc_Surface::PartialDerivativeVV(const double u, const double v) const
{
	if (isRational) // F(u,v) = f(u,v)/W(u,v); F_vv = [W*(f_vv*W - f*W_vv) - 2*W_v*(f_v*W- f*W_v)]/ W^3 
	{
		double f = DeCasteljau(coeffs, u_degree, v_degree, u, v);
		double f_v = Derivative(coeffs, u, v, 0, 1);
		double f_vv = Derivative(coeffs, u, v, 0, 2);

		double W = DeCasteljau(weights, u_degree, v_degree, u, v);		
		double W_v = Derivative(weights, u, v, 0, 1);		
		double W_vv = Derivative(weights, u, v, 0, 2);

		return (W * (f_vv * W - f * W_vv) - 2 * W_v * (f_v * W - f * W_v)) / pow(W, 3);
	}

	else
	{
		return Derivative(coeffs, u, v, 0, 2);
	}
}

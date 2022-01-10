#include "SquaredDistanceFunc_Curve.h"

SquaredDistanceFunc_Curve::SquaredDistanceFunc_Curve()
{
	isRational = false;
	weights.clear();
	coeffs.clear();	
}

SquaredDistanceFunc_Curve::~SquaredDistanceFunc_Curve()
{
}

SquaredDistanceFunc_Curve::SquaredDistanceFunc_Curve(int deg, double tmin, double tmax, std::vector<double> coefficients)
{
	degree = deg;
	t_min = tmin;
	t_max = tmax;
	isRational = false;
	weights.clear();
	coeffs = coefficients;
}

SquaredDistanceFunc_Curve::SquaredDistanceFunc_Curve(int deg, double tmin, double tmax, std::vector<double> w, std::vector<double> coefficients)
{
	degree = deg;
	t_min = tmin;
	t_max = tmax;
	isRational = true;
	weights = w;
	coeffs = coefficients;
}

SquaredDistanceFunc_Curve::SquaredDistanceFunc_Curve(BezierCurve &bezier_, Point P)
{
	// degree
	int p = bezier_.GetDegree();
	degree = 2 * p;

	// parameter span
	bezier_.GetParameterSpan(t_min, t_max);

	// rational
	isRational = bezier_.GetRational();

	// weights and coefficients
	weights.clear();
	coeffs.clear();
	coeffs.resize(degree + 1);
	std::vector<Point> ctrlpts = bezier_.GetControlPoints();
	for (int i = 0; i <= p; i++)
	{
		ctrlpts[i] -= P;	// control points of C(t) - P
	}

	std::vector<double> binominal_p = Binomial(p);
	std::vector<double> binominal_2p = Binomial(degree);
	if (isRational)
	{
		// 这里的weights仅为分母的系数，分子中不再包含weights
		// F(t) = sum(coeff[k] * B_k) / sum(weights[k] * B_k)
		weights.resize(degree + 1);
		std::vector<double> w = bezier_.GetWeights();
		for (int k = 0; k <= degree; k++)
		{
			double temp_coeff = 0.0;
			double temp_weight = 0.0;
			for (int i = std::max(0, k - p); i <= std::min(k, p); i++)
			{
				double tmp = weights[i] * weights[k - i] * binominal_p[i] * binominal_p[k - i];				
				temp_coeff += ctrlpts[i].dot(ctrlpts[k - i]) * tmp;
				temp_weight += tmp;
			}
			coeffs[k] = temp_coeff / binominal_2p[k];
			weights[k] = temp_weight / binominal_2p[k];
		}

	}
	else
	{
		for (int k = 0; k <= degree; k++)
		{
			double temp_coeff = 0.0;
			for (int i = std::max(0, k - p); i <= std::min(k, p); i++)
			{
				temp_coeff += ctrlpts[i].dot(ctrlpts[k - i]) * binominal_p[i] * binominal_p[k - i];
			}
			coeffs[k] = temp_coeff / binominal_2p[k];
		}
	}
}


double SquaredDistanceFunc_Curve::DeCasteljau(const std::vector<double> &coefficients, int deg, const double t)const
{
	if (t < t_min || t> t_max) return 0.0;
	double t_ = (t - t_min) / (t_max - t_min);

	std::vector<double> d(coefficients);
	for (int j = 1; j <= deg; j++)
	{
		for (int i = 0; i <= deg - j; i++)
		{
			d[i] = (1 - t_) * d[i] + t_ * d[i + 1];
		}
	}
	return d[0];
}

double SquaredDistanceFunc_Curve::Derivative(const std::vector<double> &coefficients, const double t, const int k) const
{
	// k: kth derivative
	if (t < t_min || t> t_max) return 0.0;
	
	std::vector<double> d(coefficients);
	for (int j = 1; j <= k; j++)
	{
		for (int i = 0; i <= degree - j; i++)
		{
			d[i] = (degree - j + 1) / (t_max - t_min) * (d[i + 1] - d[i]);
		}
	}
	std::vector<double> dev(d.begin(), d.end() - k);

	return DeCasteljau(dev, degree - k, t);
}


double SquaredDistanceFunc_Curve::DeCasteljau(const double t)const
{
	if (isRational)
	{
		return DeCasteljau(coeffs, degree, t) / DeCasteljau(weights, degree, t);
	}
	else
	{
		return DeCasteljau(coeffs, degree, t);
	}

}

double SquaredDistanceFunc_Curve::Derivative(const double t)const
{
	if (isRational) // F(t) = f(t) / W(t); F'(t) = f'/ W - f * W'/ W^2 
	{
		double f_t = Derivative(coeffs, t, 1);
		double W_t = Derivative(weights, t, 1);
		double f = DeCasteljau(coeffs, degree, t);
		double W = DeCasteljau(weights, degree, t);
		return (f_t / W - W_t * f / pow(W,2)); 
	}
	else
	{
		return Derivative(coeffs, t, 1);
	}
}

double SquaredDistanceFunc_Curve::Derivative2(const double t)const
{
	if (isRational) // F(t) = f(t) / W(t); F''(t) = [W*(f''*W - f*W'') - 2*W'*(f'*W- f*W')]/ W^3 
	{
		double f = DeCasteljau(coeffs, degree, t);
		double W = DeCasteljau(weights, degree, t);
		double f_t = Derivative(coeffs, t, 1);
		double W_t = Derivative(weights, t, 1);
		double f_tt = Derivative(coeffs, t, 2);
		double W_tt = Derivative(weights, t, 2);
		return (W * (f_tt * W - f * W_tt) - 2 * W_t * (f_t * W - f * W_t)) / pow(W, 3);
	}
	else
	{
		return Derivative(coeffs, t, 2);
	}
}


double Binomial(int n, int k)
{
	assert(n >= k && k >= 0);
	double C = 1.0;
	if (k > n / 2)	k = n - k;
	for (int i = 1; i <= k; i++)
	{
		C *= (n - i + 1) / i;
	}

	return C;
}

std::vector<double> Binomial(int n)
{
	assert(n >= 0);
	std::vector<double> Cn(n + 1);
	Cn[0] = 1.0;

	int mid = (n + 1) / 2;
	for (int i = 1; i <= mid; i++)
	{
		Cn[i] = Cn[i - 1] * (n - i + 1) / i;
	}

	for (int i = mid + 1; i <= n; i++)
	{
		Cn[i] = Cn[n - i];
	}

	return Cn;
}
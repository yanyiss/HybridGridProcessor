#include "NewtonRaphson.h"

NewtonRaphson::NewtonRaphson()
{
}

NewtonRaphson::~NewtonRaphson()
{
}


NewtonRaphson::NewtonRaphson(Point P, BezierCurve &curve, double &t)
{
	int iter_times = 50;
	double param_error = 100;
	double	tau = 0.25;
	double	gamma = 0.5;
	double t_min, t_max;
	double t_n = t;
	curve.GetParameterSpan(t_min, t_max);

	Point Ct, C_t, C_tt, Vec_PQ;
	double distance, dCt_norm;

	int i = 0;
	while (i < iter_times)
	{
		Ct = curve(t_n);
		C_t = curve.Derivative(t_n);
		C_tt = curve.Derivative2(t_n);
		Vec_PQ = Ct - P;
		distance = Vec_PQ.norm();
		dCt_norm = C_t.norm();
		//std::cout << std::endl << "t_n: " << t_n << "  distance: " << distance << "  dCt_norm: " << dCt_norm << "  angle: " << acos(C_t.dot(Vec_PQ) / dCt_norm / distance)/ M_PI << "  param_error:" << param_error << std::endl;

		//if ((distance < ERROR_BOUND) || (dCt_norm < ERROR_BOUND || (abs(C_t.dot(Vec_PQ) / dCt_norm / distance) < COSINE_BOUND)) || (param_error < ERROR_BOUND)) 
		if ((distance < ERROR_BOUND) || (abs(C_t.dot(Vec_PQ) / dCt_norm / distance) < ERROR_BOUND) || (param_error < ERROR_BOUND)) // 假设是正则曲线，一阶导不为0
		{
			t = t_n;
			return;
		}

		else
		{
			double dt;
			double t_pre = t_n;
			double df = Vec_PQ.dot(C_t);
			double ddf = C_t.squaredNorm() + Vec_PQ.dot(C_tt);
			if (ddf > eps)		// Newton's Method
			{
				dt = -df / ddf;
				UpdateParameter(t_n, dt, t_min, t_max);
			}
			else						// Steepest Descent Method: backtracking line search
			{
				dt = -df;
				UpdateParameter(t_n, dt, t_min, t_max);
				double alpha = (t_n - t_pre) / dt;

				double f_ = Vec_PQ.squaredNorm();
//				while ((curve(t_n) - P).squaredNorm() > (1 - tau * alpha) * f_)
//				while ((curve(t_n) - P).squaredNorm() > f_)
				double decrease = tau * df * dt;
				while ((curve(t_n) - P).squaredNorm() > f_ + alpha * decrease)
				{
					//std::cout << (curve(t_n) - P).squaredNorm() << "; " << f_ + alpha * decrease << std::endl;
					alpha *= gamma;
					t_n = t_pre + alpha * dt;
				}
			}

			//if (!((P - curve(t_n)).norm() < distance)) t_n = t_pre;  // 06.06改
			dt = t_n - t_pre;
			param_error = abs(dt * dCt_norm);
			//std::cout << "dt: " << dt << "  t_n: " << t_n << "  distance: " << distance << std::endl;
		}

		i++;
	}

	t = t_n;
}

NewtonRaphson::NewtonRaphson(Point P, BezierSurface &surface, double &u, double &v)
{
	int iter_times = 50;
	double param_error = 1;
	double	tau = 0.25;
	double	gamma = 0.5;
	double u_min, u_max, v_min, v_max;
	double u_n = u, v_n = v;
	surface.GetParameterSpan(u_min, u_max, v_min, v_max);

	Point Suv, S_u, S_v, S_uu, S_uv, S_vv, Vec_PQ;
	double distance, Su_norm, Sv_norm;

	int i = 0;
	while(i < iter_times)
	{
		Suv = surface(u_n, v_n);
		S_u = surface.PartialDerivativeU(u_n, v_n);
		S_v = surface.PartialDerivativeV(u_n, v_n);
		S_uu = surface.PartialDerivativeUU(u_n, v_n);
		S_uv = surface.PartialDerivativeUV(u_n, v_n);
		S_vv = surface.PartialDerivativeVV(u_n, v_n);


		Vec_PQ = Suv - P;
		distance = Vec_PQ.norm();
		Su_norm = S_u.norm();
		Sv_norm = S_v.norm();

	/*	if ((distance < ERROR_BOUND) || ((Su_norm < eps) || (abs(S_u.dot(Vec_PQ) / Su_norm / distance) <  COSINE_BOUND)) && ((Sv_norm < eps) || (abs(S_v.dot(Vec_PQ) / Sv_norm / distance) < COSINE_BOUND)) || (param_error < ERROR_BOUND))*/
		if ((distance < ERROR_BOUND) || ((abs(S_u.dot(Vec_PQ) / Su_norm / distance) <  COSINE_BOUND)) && ((abs(S_v.dot(Vec_PQ) / Sv_norm / distance) <  ERROR_BOUND)) || (param_error < ERROR_BOUND))
		{
			u = u_n;
			v = v_n;
			return;
		}
		else
		{
			double du, dv;
			double u_pre = u_n;
			double v_pre = v_n;
			double f_u = S_u.dot(Vec_PQ);
			double f_v = S_v.dot(Vec_PQ);
			double f_uu = S_u.squaredNorm() + S_uu.dot(Vec_PQ);
			double f_uv = S_u.dot(S_v) + S_uv.dot(Vec_PQ);
			double f_vv = S_v.squaredNorm() + S_vv.dot(Vec_PQ);

			Eigen::Vector2d gradient, duv;
			Eigen::Matrix2d hessian;
			gradient << f_u, f_v;
			hessian << f_uu, f_uv, f_uv, f_vv;

			double trace = hessian(0, 0) + hessian(1, 1);
			double determinant = hessian.determinant();
			if (std::min(trace, determinant) > eps)		// Newton's iteration(hessian is positive)
			{
				duv = -hessian.inverse() *  gradient;
				UpdateParameter(u_n, v_n, duv(0), duv(1), u_min, u_max, v_min, v_max);
			}

			else											// Steepest Descent Method
			{
				duv = -gradient;
				du = duv(0);
				dv = duv(1);
				UpdateParameter(u_n, v_n, du, dv, u_min, u_max, v_min, v_max);
				double alpha = (u_n - u_pre) / du;

				double f_ = Vec_PQ.squaredNorm();
				double decrease = tau * gradient.dot(duv);
				while ((surface(u_n, v_n) - P).squaredNorm() >  f_ + alpha * decrease)
				{
					alpha *= gamma;
					u_n = u_pre + alpha * du;
					v_n = v_pre + alpha * dv;
				}
			}

			//if (!((P - surface(u_n, v_n)).norm() < distance))	// 06.06改
			//{
			//	u_n = u_pre;
			//	v_n = v_pre;
			//}
			du = u_n - u_pre;
			dv = v_n - v_pre;
			param_error = (du * S_u + dv * S_v).norm();
		}

		i++;
	}

	u = u_n;
	v = v_n;
}

NewtonRaphson::NewtonRaphson(Point P, BSplineCurve &curve, double &t)
{
	int iter_times = 20;
	double param_error = 1;
	double	tau = 0.25;
	double	gamma = 0.5;
	double t_min, t_max;
	double t_n = t;
	auto knots = curve.GetKnots();
	t_min = knots.front();
	t_max = knots.back();

	Point Ct, C_t, C_tt, Vec_PQ;
	double distance, dCt_norm;

	int i = 0;
	while (i < iter_times)
	{
		Ct = curve(t_n);
		C_t = curve.Derivative(t_n);
		C_tt = curve.Derivative2(t_n);
		Vec_PQ = Ct - P;
		distance = Vec_PQ.norm();
		dCt_norm = C_t.norm();

		if ((distance < ERROR_BOUND) || (dCt_norm < ERROR_BOUND || (abs(C_t.dot(Vec_PQ) / dCt_norm / distance) <  ERROR_BOUND)) || (param_error < ERROR_BOUND)) // 假设是正则曲线，一阶导不为0
		{
			t = t_n;
			return;
		}

		else
		{
			double dt;
			double t_pre = t_n;
			double df = Vec_PQ.dot(C_t);
			double ddf = C_t.squaredNorm() + Vec_PQ.dot(C_tt);
			if (ddf > eps)		// Newton's Method
			{
				dt = -df / ddf;
				UpdateParameter(t_n, dt, t_min, t_max);
			}
			else						// Steepest Descent Method(解决除零的问题):backtracking line search
			{
				dt = -df;
				UpdateParameter(t_n, dt, t_min, t_max);
				double alpha = (t_n - t_pre) / dt;

				double f_ = Vec_PQ.squaredNorm();
				double decrease = tau * df * dt;
				while ((curve(t_n) - P).squaredNorm() > f_ + alpha * decrease)
				{
					alpha *= gamma;
					t_n = t_pre + alpha * dt;
				}
			}
			dt = t_n - t_pre;
			param_error = abs(dt * dCt_norm);
		}

		i++;
	}

	t = t_n;
}

NewtonRaphson::NewtonRaphson(Point P, BSplineSurface &surface, double &u, double &v)
{
	int iter_times = 20;
	double param_error = 1;
	double	tau = 0.25;
	double	gamma = 0.5;
	double u_min, u_max, v_min, v_max;
	double u_n = u, v_n = v;
	auto u_knots = surface.GetKnotsU();
	auto v_knots = surface.GetKnotsV();
	u_min = u_knots.front();
	u_max = u_knots.back();
	v_min = v_knots.front();
	v_max = v_knots.back();

	Point Suv, S_u, S_v, S_uu, S_uv, S_vv, Vec_PQ;
	double distance, Su_norm, Sv_norm;

	int i = 0;
	while (i < iter_times)
	{
		Suv = surface(u_n, v_n);
		S_u = surface.PartialDerivativeU(u_n, v_n);
		S_v = surface.PartialDerivativeV(u_n, v_n);
		S_uu = surface.PartialDerivativeUU(u_n, v_n);
		S_uv = surface.PartialDerivativeUV(u_n, v_n);;
		S_vv = surface.PartialDerivativeVV(u_n, v_n);
		Vec_PQ = Suv - P;
		distance = Vec_PQ.norm();
		Su_norm = S_u.norm();
		Sv_norm = S_v.norm();

		if ((distance < ERROR_BOUND) || ((Su_norm < eps) || (abs(S_u.dot(Vec_PQ) / Su_norm / distance) <   ERROR_BOUND)) && ((Sv_norm < eps) || (abs(S_v.dot(Vec_PQ) / Sv_norm / distance) <  ERROR_BOUND)) || (param_error < ERROR_BOUND))
		{
			u = u_n;
			v = v_n;
			return;
		}
		else
		{
			double du, dv;
			double u_pre = u_n;
			double v_pre = v_n;
			double f_u = S_u.dot(Vec_PQ);
			double f_v = S_v.dot(Vec_PQ);
			double f_uu = S_u.squaredNorm() + S_uu.dot(Vec_PQ);
			double f_uv = S_u.dot(S_v) + S_uv.dot(Vec_PQ);
			double f_vv = S_v.squaredNorm() + S_vv.dot(Vec_PQ);

			Eigen::Vector2d gradient, duv;
			Eigen::Matrix2d hessian;
			gradient << f_u, f_v;
			hessian << f_uu, f_uv, f_uv, f_vv;

			double trace = hessian(0, 0) + hessian(1, 1);
			double determinant = hessian.determinant();
			if (std::min(trace, determinant) > eps)		// Newton's iteration(hessian is positive)
			{
				duv = -hessian.inverse() *  gradient;
				UpdateParameter(u_n, v_n, duv(0), duv(1), u_min, u_max, v_min, v_max);
			}

			else											// Steepest Descent Method
			{
				duv = -gradient;
				du = duv(0);
				dv = duv(1);
				UpdateParameter(u_n, v_n, du, dv, u_min, u_max, v_min, v_max);
				double alpha = (u_n - u_pre) / du;

				double f_ = Vec_PQ.squaredNorm();
				double decrease = tau * gradient.dot(duv);
				while ((surface(u_n, v_n) - P).squaredNorm() >  f_ + alpha * decrease)
				{
					alpha *= gamma;
					u_n = u_pre + alpha * du;
					v_n = v_pre + alpha * dv;
				}
			}

			//if (!((P - surface(u_n, v_n)).norm() < distance))	// 06.06改
			//{
			//	u_n = u_pre;
			//	v_n = v_pre;
			//}
			du = u_n - u_pre;
			dv = v_n - v_pre;
			param_error = (du * S_u + dv * S_v).norm();
		}

		i++;
	}

	u = u_n;
	v = v_n;
}


void NewtonRaphson::UpdateParameter(double &t, double dt, double t_min, double t_max)
{
	double t_ = t + dt;
	while (t_ > t_max || t_ < t_min)
	{
		dt *= 0.5;
		t_ = t + dt;
	}
	t = t_;
}

void NewtonRaphson::UpdateParameter(double &u, double &v, double du, double dv, double u_min, double u_max, double v_min, double v_max)
{
	double u_ = u + du;
	double v_ = v + dv;
	while (u_ > u_max || u_ < u_min || v_>v_max || v_ < v_min)
	{
		du *= 0.5;
		dv *= 0.5;
		u_ = u + du;
		v_ = v + dv;
	}
	u = u_;
	v = v_;
}

// 迭代停止的判断条件到底是 && 还是 ||
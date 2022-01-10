#include "NewtonRaphsonIteration.h"
#include <iostream>


/*由于精度有限，alpha算出来再代回去乘迭代方向，得到的下一步的参数仍可能越界*/ 
NewtonRaphsonIteration::NewtonRaphsonIteration(Point P_, BSplineCurve &curve, double &t0)
{
	curve_ = curve;
	t = t0;
	P = P_;

	size_t ctrlpts_num = curve_.GetNumOfCtrlpts();
	std::vector<double> knots = curve_.GetKnots();
	size_t degree = curve_.GetDegree();
	t_max = knots[ctrlpts_num];
	t_min = knots[degree];

	iterations_C = 20;
	param_error_C = 1;

	int i = 0;
	for (; i < iterations_C; i++)
	{
		UpdateCurveValues();
//		std::cout << "distance = " << distance << std::endl << std::endl;

		Vec3d dCt = curve_Derivative;
		double dCt_norm = curve_Derivative.norm();
		//if ((distance < ERROR_BOUND) || ((dCt_norm > eps) && (abs(dCt.dot(Vec_PQ) / dCt_norm / distance) < COSINE_BOUND)) || (param_error_C < ERROR_BOUND)) 
		// 2021.4.27 
		if ((distance < ERROR_BOUND) && ((dCt_norm > eps) && (abs(dCt.dot(Vec_PQ) / dCt_norm / distance) < COSINE_BOUND)) || (param_error_C < ERROR_BOUND))
		{
			t0 = t;
			double dot_ = abs(dCt.dot(Vec_PQ));
			if (dot_ >= eps)
				cosine_error_C = dot_ / dCt_norm / distance;
			else
				cosine_error_C = dot_;
			std::cout << "The corresponding parameter of point projection is: " << t0 << ";  cosine = " << cosine_error_C << std::endl << std::endl;
			return;
		}
		else
		{
			double dt;
			double t_pre = t;
			double df = DistFuncDerivative_C();
			double ddf = DistFuncDerivative2_C();
			if (abs(ddf) > eps)		// Newton's Method
			{
				dt = -df / ddf;
				t += dt;
				CheckIfOutOfRange_C(t);
			}
			else					// Steepest Descent Method
			{
				dt = -df;			
				double  alpha = CheckIfOutOfRange_C(t, dt); // stepsize: backtracking line search
				t = t_pre + alpha * dt;
				CheckIfOutOfRange_C(t);

				double f_ = Vec_PQ.squaredNorm();
				while ((curve_(t) - P).squaredNorm() > (1 - tau * alpha) * f_)
				{
					alpha *= gamma;
					t = t_pre + alpha * dt;
				}
			}
			dt = t - t_pre;
			param_error_C = abs(dt * dCt_norm);
//			std::cout << "dt: " << dt << ";  " << std::endl << "t:" << t << std::endl << std::endl;
		}
	}

	if (i == iterations_C)
	{
		std::cout << "The initial parameter is not available." << std::endl << std::endl;
	}
}

NewtonRaphsonIteration::NewtonRaphsonIteration(Point P_, BSplineSurface &surface, Vec2d &uv_0)
{
	surface_ = surface;
	uv = uv_0;
	P = P_;
	surface_Derivative.resize(2, 3);
	surface_Derivative2.resize(2, 6);

	size_t u_ctrlpts_num = surface_.GetNumOfCtrlptsU();
	size_t v_ctrlpts_num = surface_.GetNumOfCtrlptsV();
	std::vector<double> knotsu = surface_.GetKnotsU();
	std::vector<double> knotsv = surface_.GetKnotsV();
	size_t u_degree = surface_.GetDegreeU();
	size_t v_degree = surface_.GetDegreeV();
	u_max = knotsu[u_ctrlpts_num];
	v_max = knotsv[v_ctrlpts_num];
	u_min = knotsu[u_degree];
	v_min = knotsv[v_degree];

	iterations_S = 20;
	param_error_S = 1;

	int i = 0;
	for (; i < iterations_S; i++)
	{
		UpdateSurfaceValues();
//		std::cout << "distance = " << distance << std::endl << std::endl;
		
		Vec3d Su = surface_Derivative.row(0);
		Vec3d Sv = surface_Derivative.row(1);
		double Su_norm = Su.norm();
		double Sv_norm = Sv.norm();

		/*if ((distance < ERROR_BOUND) || ((Su_norm > eps) && (Sv_norm > eps) && (abs(Su.dot(Vec_PQ) / Su_norm / distance) <  COSINE_BOUND) &&  (abs(Sv.dot(Vec_PQ) / Sv_norm / distance) < COSINE_BOUND)) || (param_error_S < ERROR_BOUND))*/
		// 2021.04.27
		if ((distance < ERROR_BOUND) && ((Su_norm > eps) && (Sv_norm > eps) && (abs(Su.dot(Vec_PQ) / Su_norm / distance) <  COSINE_BOUND) && (abs(Sv.dot(Vec_PQ) / Sv_norm / distance) < COSINE_BOUND)) || (param_error_S < ERROR_BOUND))
		{
			uv_0 = uv;
			Vec3d normal = Su.cross(Sv);
			double dot_ = abs(normal.dot(Vec_PQ));
			double sin_;
			if (dot_ >= eps)
			{
				double cos_ = dot_ / normal.norm() / distance;
				sin_ = sqrt(1 - cos_*cos_);
			}				
			else
				sin_ = dot_;
			std::cout << "The corresponding parameters of point orthogonal projection is: " << uv_0.transpose() << ";  sin = " << sin_ << std::endl << std::endl;
			return;
		}
		else
		{
			Vec2d duv;
			Vec2d uv_pre = uv;
			double alpha = 1;
			Vec2d gradient = DistFuncDerivative_S();
			Mat2d hessian = DistFuncDerivative2_S();

			double trace = hessian(0, 0) + hessian(1, 1);
			double determinant = hessian.determinant();

			if (std::min(trace, determinant) > eps)		// positive, Newton's iteration
			{
				duv = - hessian.inverse() *  gradient;
				alpha = CheckIfOutOfRange_S(uv, duv); // 不越界为1，否则压缩步长
				uv += alpha * duv;
				CheckIfOutOfRange_S(uv);
			}
			else							// Steepest Descent Method
			{
				duv = -gradient;
				alpha = CheckIfOutOfRange_S(uv, duv);
				uv = uv_pre + alpha * duv;
				CheckIfOutOfRange_S(uv);
				double f_ = Vec_PQ.squaredNorm();			

				while ((surface_(uv(0), uv(1)) - P).squaredNorm() > (1 - tau * alpha) * f_)
				{
					alpha *= gamma;
					uv = uv_pre + alpha * duv;
				}
			}
			duv = uv - uv_pre;
			param_error_S = (duv(0) * Su + duv(1) * Sv).norm();
//			std::cout << "duv: " << duv.transpose() << ";  " << std::endl << "uv:" << uv.transpose() << std::endl << std::endl;
			
		}
	}

	if (i == iterations_S)
	{
		std::cout << "The initial parameter is not available." << std::endl << std::endl;
	}
}

NewtonRaphsonIteration::~NewtonRaphsonIteration()
{
	
}


void NewtonRaphsonIteration::UpdateCurveValues()
{
	Q = curve_(t);
	Vec_PQ = Q - P;
	distance = Vec_PQ.norm();
	curve_Derivative = curve_.Derivative(t);
	curve_Derivative2 = curve_.Derivative2(t);
}

void NewtonRaphsonIteration::UpdateSurfaceValues()
{
	Q = surface_(uv(0), uv(1));
	Vec_PQ = Q - P;
	distance = Vec_PQ.norm();
	surface_Derivative.row(0) = surface_.PartialDerivativeU(uv(0), uv(1));
	surface_Derivative.row(1) = surface_.PartialDerivativeV(uv(0), uv(1));
	surface_Derivative2 << surface_.PartialDerivativeUU(uv(0), uv(1)), surface_.PartialDerivativeUV(uv(0), uv(1)), surface_.PartialDerivativeUV(uv(0), uv(1)), surface_.PartialDerivativeVV(uv(0), uv(1));

//	Vec3d tangent_u = surface_Derivative.row(0);
//	Vec3d tangent_v = surface_Derivative.row(1);
//	surface_normal = tangent_u.cross(tangent_v);
}


double NewtonRaphsonIteration::DistFuncDerivative_C()
{
	return Vec_PQ.dot(curve_Derivative);	
}

Vec2d NewtonRaphsonIteration::DistFuncDerivative_S()
{
	Vec2d gradient(0, 0);
	gradient(0) = (surface_Derivative.row(0)).dot(Vec_PQ);
	gradient(1) = (surface_Derivative.row(1)).dot(Vec_PQ);
	return gradient;
}


double NewtonRaphsonIteration::DistFuncDerivative2_C()
{
	return curve_Derivative.squaredNorm() + Vec_PQ.dot(curve_Derivative2);
}

Mat2d NewtonRaphsonIteration::DistFuncDerivative2_S(void)
{
	Mat2d hessian;
	Vec3d Suu = surface_Derivative2.topLeftCorner(1, 3).transpose();
	Vec3d Suv = surface_Derivative2.topRightCorner(1, 3).transpose();
	Vec3d Svv = surface_Derivative2.bottomRightCorner(1, 3).transpose();

	double fuu, fuv, fvv;
	fuu = Suu.dot(Vec_PQ) + surface_Derivative.row(0).squaredNorm();
	fuv = Suv.dot(Vec_PQ) + surface_Derivative.row(0).dot(surface_Derivative.row(1));
	fvv = Svv.dot(Vec_PQ) + surface_Derivative.row(1).squaredNorm();
	hessian << fuu, fuv,
			   fuv, fvv;
	return hessian;
}


double NewtonRaphsonIteration::CheckIfOutOfRange_C(double t, double dt)
{
//	std::cout << "t = " << t << ";  " << "dt = " << dt << std::endl;
	if (abs(dt) < eps)
		return 0;
	else
		return std::min((dt > 0 ? ((t_max - t) / dt) : ((t_min - t) / dt)), 1.0);
	//if (abs(dt) < eps)
	//	return 0;
	//else
	//	return std::min((dt > 0 ? ((t_max - t) / dt) : ((t_min - t) / dt)), 1.0);
}

double NewtonRaphsonIteration::CheckIfOutOfRange_S(Vec2d uv, Vec2d duv)
{
//	std::cout << "uv = " << uv << ";  " << "duv = " << duv << std::endl;
	if (abs(duv(0)) < eps)
	{
		if (abs(duv(1)) < eps)
			return 0;
		else
		{
			double rate_v = duv(1) > 0 ? ((v_max - uv(1)) / duv(1)) : ((v_min - uv(1)) / duv(1));
			return std::min(rate_v, 1.0);
		}
	}
	else
	{
		if (abs(duv(1)) < eps)
		{
			double rate_u = duv(0) > 0 ? ((u_max - uv(0)) / duv(0)) : (u_min - uv(0) / duv(0));
			return std::min(rate_u, 1.0);
		}			
		else
		{
			double rate_u = duv(0) > 0 ? ((u_max - uv(0)) / duv(0)) : (u_min - uv(0) / duv(0));
			double rate_v = duv(1) > 0 ? ((v_max - uv(1)) / duv(1)) : (v_min - uv(1) / duv(1));
			return std::min(std::min(rate_u, rate_v), 1.0);
		}
	}
}


void NewtonRaphsonIteration::CheckIfOutOfRange_C(double &t)
{
	if (t <= t_min)
	{
		t = t_min;
	}
	else if (t >= t_max)
	{
		t = t_max;
	}
}

void NewtonRaphsonIteration::CheckIfOutOfRange_S(Vec2d &uv)
{
	if (uv(0) <= u_min)
	{
		uv(0) = u_min;
	}
	else if (uv(0) >= u_max)
	{
		uv(0) = u_max;
	}
	if (uv(1) <= v_min)
	{
		uv(1) = v_min;
	}
	else if (uv(1) >= v_max)
	{
		uv(1) = v_max;
	}
}

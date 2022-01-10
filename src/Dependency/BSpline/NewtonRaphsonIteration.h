#pragma once
#include <Eigen/Eigen>
#include "Bsplinecurve.h"
#include "Bsplinesurface.h"

typedef Eigen::Vector3d Vec3d;
typedef Eigen::Vector3d Point;
typedef Eigen::Vector2d Vec2d;
typedef Eigen::Matrix2d Mat2d;

class NewtonRaphsonIteration
{
public:
	NewtonRaphsonIteration(Point P, BSplineCurve &curve, double &t0);
	NewtonRaphsonIteration(Point P, BSplineSurface &surface, Vec2d &uv_0);

	~NewtonRaphsonIteration();
	
	void	UpdateCurveValues(void);
	void	UpdateSurfaceValues(void);

	double	DistFuncDerivative_C(void);
	Vec2d	DistFuncDerivative_S(void);

	double	DistFuncDerivative2_C(void);
	Mat2d	DistFuncDerivative2_S(void);

	double	CheckIfOutOfRange_C(double t, double dt);
	double	CheckIfOutOfRange_S(Vec2d uv, Vec2d duv);

	void CheckIfOutOfRange_C(double &t);
	void CheckIfOutOfRange_S(Vec2d &uv);

private:
	Point			P;
	Point			Q;
	Vec3d			Vec_PQ;
	double			distance;

	// curve
	BSplineCurve	curve_;		
	double			t;						// 曲线参数值
	double			t_max;					// 曲线参数最大有效值
	double			t_min;					// 曲线参数最小有效值
	Vec3d			curve_Derivative;		// 曲线一阶导
	Vec3d			curve_Derivative2;		// 曲线二阶导
	Vec3d			curve_tangent;			// 曲线的切向
	int				iterations_C;			// 迭代次数
	double			cosine_error_C;			// 切线与距离向量夹角余弦
	double			param_error_C;			// 衡量迭代前后参数变化

	// surface
	BSplineSurface	surface_;
	Vec2d			uv;						// 曲面参数值(2D)
	double			u_max;					// 曲面参数u的最大值
	double			u_min;					// 曲面参数u的最小值
	double			v_max;					// 曲面参数v的最大值
	double			v_min;					// 曲面参数v的最小值
	Eigen::MatrixXd surface_Derivative;		// 曲面gradient
	Eigen::MatrixXd surface_Derivative2;	// 曲面Hessian阵
	int				iterations_S;			// 迭代次数
	double			cosine_error_Su;		// 切线与距离向量夹角余弦
	double			cosine_error_Sv;		// 切线与距离向量夹角余弦
	double			param_error_S;			// 衡量迭代前后参数变化

	double	ERROR_BOUND = 1e-10;
	double  COSINE_BOUND = cos(M_PI * 89 / 180);	// cos 89°= 0.01745
	double	tau = 0.25;
	double	gamma = 0.5;
	double	eps = 1e-14;
	// double	eta = 0.1;
};



// 输入是什么，输出是什么
// 本项目中确定唯一性后，不用最速下降，直接newton法迭代
// 注意越界的问题
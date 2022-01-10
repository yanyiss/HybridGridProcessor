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
	double			t;						// ���߲���ֵ
	double			t_max;					// ���߲��������Чֵ
	double			t_min;					// ���߲�����С��Чֵ
	Vec3d			curve_Derivative;		// ����һ�׵�
	Vec3d			curve_Derivative2;		// ���߶��׵�
	Vec3d			curve_tangent;			// ���ߵ�����
	int				iterations_C;			// ��������
	double			cosine_error_C;			// ��������������н�����
	double			param_error_C;			// ��������ǰ������仯

	// surface
	BSplineSurface	surface_;
	Vec2d			uv;						// �������ֵ(2D)
	double			u_max;					// �������u�����ֵ
	double			u_min;					// �������u����Сֵ
	double			v_max;					// �������v�����ֵ
	double			v_min;					// �������v����Сֵ
	Eigen::MatrixXd surface_Derivative;		// ����gradient
	Eigen::MatrixXd surface_Derivative2;	// ����Hessian��
	int				iterations_S;			// ��������
	double			cosine_error_Su;		// ��������������н�����
	double			cosine_error_Sv;		// ��������������н�����
	double			param_error_S;			// ��������ǰ������仯

	double	ERROR_BOUND = 1e-10;
	double  COSINE_BOUND = cos(M_PI * 89 / 180);	// cos 89��= 0.01745
	double	tau = 0.25;
	double	gamma = 0.5;
	double	eps = 1e-14;
	// double	eta = 0.1;
};



// ������ʲô�������ʲô
// ����Ŀ��ȷ��Ψһ�Ժ󣬲��������½���ֱ��newton������
// ע��Խ�������
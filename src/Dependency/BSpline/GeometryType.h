#pragma once
#include <Eigen\Eigen>
#include"BSplineCurve.h"

typedef Eigen::Vector3d Point;
typedef Eigen::Vector4d Point4;
#define PI 3.1415926535897932

class GeometryType
{
public:
	GeometryType() {};
	virtual ~GeometryType() {};
public:
	virtual Point PartialDerivativeU(const double u, const double v) const = 0;
	virtual Point PartialDerivativeV(const double u, const double v) const = 0;
	virtual Point PartialDerivativeUU(const double u, const double v) const = 0;
	virtual Point PartialDerivativeUV(const double u, const double v) const = 0;
	virtual Point PartialDerivativeVV(const double u, const double v) const = 0;
	virtual void PrincipalCurvature(const double u, const double v, double &k1, double &k2) const = 0;
	Point4 Getbounds() { return UV; }
protected:
	Point4 UV;
};

class PlaneType:public GeometryType
{
public:
	PlaneType(Point &origin, Point &xdir, Point &ydir) :
		Origin(origin), Xdir(xdir), Ydir(ydir) {UV = Point4(-DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX);};
	~PlaneType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return Xdir; }
	Point PartialDerivativeV(const double u, const double v) const { return Ydir; }
	Point PartialDerivativeUU(const double u, const double v) const { return Point(0, 0, 0); }
	Point PartialDerivativeUV(const double u, const double v) const { return Point(0, 0, 0); }
	Point PartialDerivativeVV(const double u, const double v) const { return Point(0, 0, 0); }
	void PrincipalCurvature(const double u, const double v, double &k1, double &k2) const { k1 = k2 = 0; }

private:
	Point Origin, Xdir, Ydir;	
};

class CylindricalType :public GeometryType
{
public:
	CylindricalType(Point &location, Point &xAxis, Point &yAxis, Point &zAxis, double &r) :
		Location(location), XAxis(xAxis), YAxis(yAxis), ZAxis(zAxis), R(r) {UV = Point4(-DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX);};
	~CylindricalType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return -std::sin(u)*R*XAxis + std::cos(u)*R*YAxis; }
	Point PartialDerivativeV(const double u, const double v) const { return ZAxis; }
	Point PartialDerivativeUU(const double u, const double v) const { return -std::cos(u)*R*XAxis - std::sin(u)*R*YAxis; }
	Point PartialDerivativeUV(const double u, const double v) const { return Point(0, 0, 0); }
	Point PartialDerivativeVV(const double u, const double v) const { return Point(0, 0, 0); }
	void PrincipalCurvature(const double u, const double v, double &k1, double &k2) const { k1 = 1/R; k2 = 0; }
private:
	Point Location, XAxis, YAxis, ZAxis;
	double R;
};

class ConicalType :public GeometryType
{
public:
	ConicalType(Point &origin, Point &xdir, Point &ydir, Point &zdir, double &r, double &ang) :
		Origin(origin), Xdir(xdir), Ydir(ydir), Zdir(zdir), R(r), Ang(ang) {UV = Point4(-DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX);}
	~ConicalType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return (R + v * std::sin(Ang)) * (-std::sin(u)*Xdir + std::cos(u)*Ydir);}
	Point PartialDerivativeV(const double u, const double v) const { return std::sin(Ang)*(std::cos(u)*Xdir + std::sin(u)*Ydir)+std::cos(Ang)*Zdir;}
	Point PartialDerivativeUU(const double u, const double v) const { return (R + v * std::sin(Ang)) * (-std::cos(u)*Xdir - std::sin(u)*Ydir); }
	Point PartialDerivativeUV(const double u, const double v) const { return std::sin(Ang) * (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeVV(const double u, const double v) const { return Point(0, 0, 0); }
	void PrincipalCurvature(const double u, const double v, double &k1, double &k2) const 
	{ 
		Eigen::Matrix2d Weingarten, value;
		Point ru = PartialDerivativeU(u, v);
		Point rv = PartialDerivativeV(u, v);
		double E = ru.dot(ru);
		double F = ru.dot(rv);
		double G = rv.dot(rv);
		Point n = (ru.cross(rv)).normalized();
		double L = PartialDerivativeUU(u, v).dot(n);
		double M = PartialDerivativeUV(u, v).dot(n);
		double N = 0;
		Weingarten << L * G - M * F, M * E - L * F,
			M * G - N * F, N * E - M * F;
		Eigen::EigenSolver<Eigen::Matrix2d> es(Weingarten);
		value = (es.pseudoEigenvalueMatrix()) / (E*G - pow(F, 2));
		//vec = es.pseudoEigenvectors();
		k1 = std::max(std::abs(value(0, 0)), std::abs(value(1, 1)));
		k2 = 0;
	}
private:
	Point Origin, Xdir, Ydir, Zdir;
	double R, Ang;
}; 

class SphericalType :public GeometryType
{
public:
	SphericalType(Point &origin, Point &xdir, Point &ydir, Point &zdir, double &r) :
		Origin(origin), Xdir(xdir), Ydir(ydir), Zdir(zdir), R(r) {UV = Point4(-DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX);};
	~SphericalType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return R * std::cos(v)* (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeV(const double u, const double v) const { return -R * std::sin(v)*(std::cos(u)*Xdir + std::sin(u)*Ydir) + R*std::cos(v)*Zdir; }
	Point PartialDerivativeUU(const double u, const double v) const { return R * std::cos(v)* (-std::cos(u)*Xdir - std::sin(u)*Ydir); }
	Point PartialDerivativeUV(const double u, const double v) const { return -R * std::sin(v) * (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeVV(const double u, const double v) const { return -R * std::cos(v)*(std::cos(u)*Xdir + std::sin(u)*Ydir)-R*std::sin(v)*Zdir; }
	void PrincipalCurvature(const double u, const double v, double &k1, double &k2) const { k1 = k2 = 1 / R; }
private:
	Point Origin, Xdir, Ydir, Zdir;
	double R;
};

class ToroidalType :public GeometryType
{
public:
	ToroidalType(Point &origin, Point &xdir, Point &ydir, Point &zdir, double &r1, double &r2):
		Origin(origin), Xdir(xdir), Ydir(ydir), Zdir(zdir), R(r1), r(r2) {UV = Point4(-DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX);};
	~ToroidalType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return (R+r * std::cos(v))* (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeV(const double u, const double v) const { return -r * std::sin(v)*(std::cos(u)*Xdir + std::sin(u)*Ydir) + r * std::cos(v)*Zdir; }
	Point PartialDerivativeUU(const double u, const double v) const { return (R + r * std::cos(v)) * (-std::cos(u)*Xdir - std::sin(u)*Ydir); }
	Point PartialDerivativeUV(const double u, const double v) const { return -r * std::sin(v) * (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeVV(const double u, const double v) const { return -r * std::cos(v)*(std::cos(u)*Xdir + std::sin(u)*Ydir) - r*std::sin(v)*Zdir; }
	void PrincipalCurvature(const double u, const double v, double &k1, double &k2) const 
	{
		Eigen::Matrix2d Weingarten, value;
		Point ru = PartialDerivativeU(u, v);
		Point rv = PartialDerivativeV(u, v);
		double E = ru.dot(ru);
		double F = ru.dot(rv);
		double G = rv.dot(rv);
		Point n = (ru.cross(rv)).normalized();
		double L = PartialDerivativeUU(u, v).dot(n);
		double M = PartialDerivativeUV(u, v).dot(n);
		double N = PartialDerivativeVV(u, v).dot(n);
		Weingarten << L * G - M * F, M * E - L * F,
			M * G - N * F, N * E - M * F;
		Eigen::EigenSolver<Eigen::Matrix2d> es(Weingarten);
		value = (es.pseudoEigenvalueMatrix()) / (E*G - pow(F, 2));
		k1 = std::max(std::abs(value(0, 0)), std::abs(value(1, 1)));
		k2 = std::min(std::abs(value(0, 0)), std::abs(value(1, 1)));
	}
private:
	Point Origin, Xdir, Ydir, Zdir;
	double R, r;
};

class SurfaceRevolutionType :public GeometryType
{
public:
	SurfaceRevolutionType(Point &Origin, Point &Xdir, BSplineCurve *Curve) :
		origin(Origin), xdir(Xdir)
	{
		curve = Curve;
		auto knots = curve->GetKnots();
		UV = Point4(-DBL_MAX, DBL_MAX, knots[curve->GetDegree()], knots[curve->GetNumOfCtrlpts()]);
		Point pos = curve->DeBoor((UV[0] + UV[1])*0.5);
		ydir = (pos - origin - (pos.dot(xdir) - origin.dot(xdir))*xdir).normalized();
		zdir = (xdir.cross(ydir)).normalized();
	};
	~SurfaceRevolutionType() 
	{
		if (curve)
		{
			delete curve;
			curve = nullptr;
		}
	};
public:
	Point PartialDerivativeU(const double u, const double v) const
	{
		Point pos = curve->DeBoor(v);
		double length = (pos - origin - (pos.dot(xdir) - origin.dot(xdir))*xdir).norm();
		return -length * std::sin(u)*ydir + length * std::cos(u)*zdir;
	}

	Point PartialDerivativeV(const double u, const double v) const
	{
		Point pos = curve->DeBoor(v);
		Point QP = pos - origin - (pos.dot(xdir) - origin.dot(xdir))*xdir;
		Point pos_v = curve->Derivative(v);
		double partial_len = QP.dot(pos_v - (pos_v.dot(xdir)*xdir)) / QP.norm();
		return partial_len * std::cos(u)*ydir + partial_len * std::sin(u)*zdir + pos_v.dot(xdir)*xdir;
	}

	Point PartialDerivativeUU(const double u, const double v) const
	{
		Point pos = curve->DeBoor(v);
		double length = (pos - origin - (pos.dot(xdir) - origin.dot(xdir))*xdir).norm();
		return -length * std::cos(u)*ydir - length * std::sin(u)*zdir;
	}

	Point PartialDerivativeUV(const double u, const double v) const
	{
		Point pos = curve->DeBoor(v);
		Point QP = pos - origin - (pos.dot(xdir) - origin.dot(xdir))*xdir;
		Point pos_v = curve->Derivative(v);
		double partial_len = QP.dot(pos_v - (pos_v.dot(xdir)*xdir)) / QP.norm();
		return -partial_len * std::sin(u)*ydir + partial_len * std::cos(u)*zdir;
	}

	Point PartialDerivativeVV(const double u, const double v) const
	{
		Point pos = curve->DeBoor(v);
		Point QP = pos - origin - (pos.dot(xdir) - origin.dot(xdir))*xdir;
		Point pos_v = curve->Derivative(v);
		Point pos_vv = curve->Derivative2(v);
		Point QP_v = pos_v - pos_v.dot(xdir)*xdir;
		double partial2_len = ((pos_v.dot(QP_v)+QP.dot(pos_vv-pos_vv.dot(xdir)*xdir))*QP.norm()-QP.dot(QP_v)*QP_v.dot(ydir)) / std::pow(QP.norm(),2);
		return partial2_len * std::cos(u)*ydir + partial2_len * std::sin(u)*zdir + pos_vv.dot(xdir)*xdir;
	}

	void PrincipalCurvature(const double u, const double v, double &k1, double &k2) const
	{
		Point pos = curve->DeBoor(v);   //旋转曲线对应的点的坐标
		Point pos_v = curve->Derivative(v);
		Point pos_vv = curve->Derivative2(v);
		Point QP = pos - origin - (pos.dot(xdir) - origin.dot(xdir))*xdir;  //v对应点到旋转轴的垂直向量
		Point QP_v = pos_v - pos_v.dot(xdir)*xdir;
		double length = QP.norm();
		double partial_len = QP.dot(pos_v - (pos_v.dot(xdir)*xdir)) / length;
		double partial2_len = ((pos_v.dot(QP_v) + QP.dot(pos_vv - pos_vv.dot(xdir)*xdir))*QP.norm() - QP.dot(QP_v)*QP_v.dot(ydir)) / std::pow(QP.norm(), 2);

		Point ru = -length * std::sin(u)*ydir + length * std::cos(u)*zdir;
		Point rv = partial_len * std::cos(u)*ydir + partial_len * std::sin(u)*zdir + pos_v.dot(xdir)*xdir;
		Point ruu = -length * std::cos(u)*ydir - length * std::sin(u)*zdir;
		Point ruv = -partial_len * std::sin(u)*ydir + partial_len * std::cos(u)*zdir;
		Point rvv = partial2_len * std::cos(u)*ydir + partial2_len * std::sin(u)*zdir + pos_vv.dot(xdir)*xdir;

		Eigen::Matrix2d Weingarten, value;
		double E = ru.dot(ru);
		double F = ru.dot(rv);
		double G = rv.dot(rv);
		Point n = (ru.cross(rv)).normalized();
		double L = ruu.dot(n);
		double M = ruv.dot(n);
		double N = rvv.dot(n);
		Weingarten << L * G - M * F, M * E - L * F,
			M * G - N * F, N * E - M * F;
		Eigen::EigenSolver<Eigen::Matrix2d> es(Weingarten);
		value = (es.pseudoEigenvalueMatrix()) / (E*G - pow(F, 2));
		k1 = std::max(std::abs(value(0, 0)), std::abs(value(1, 1)));
		k2 = std::min(std::abs(value(0, 0)), std::abs(value(1, 1)));
	}

private:
	Point origin, xdir, ydir, zdir;
	BSplineCurve *curve;
};

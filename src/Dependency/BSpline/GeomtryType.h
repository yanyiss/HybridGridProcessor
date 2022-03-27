#pragma once
#include <Eigen\Eigen>

typedef Eigen::Vector3d Point;
typedef Eigen::Vector4d Point4;
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
private:
	double umax, umin, vmax, vmin;
};

class PlaneType:public GeometryType
{
public:
	PlaneType(Point &origin, Point &xdir, Point &ydir):Origin(origin),Xdir(xdir),Ydir(ydir) {};
	~PlaneType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return Xdir; }
	Point PartialDerivativeV(const double u, const double v) const { return Ydir; }
	Point PartialDerivativeUU(const double u, const double v) const { return Point(0, 0, 0); }
	Point PartialDerivativeUV(const double u, const double v) const { return Point(0, 0, 0); }
	Point PartialDerivativeVV(const double u, const double v) const { return Point(0, 0, 0); }
private:
	Point Origin, Xdir, Ydir;
	
};

class CylindricalType :public GeometryType
{
public:
	CylindricalType() {};  //怎么赋给私有变量
	~CylindricalType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return -std::sin(u)*R*XAxis + std::cos(u)*R*YAxis; }
	Point PartialDerivativeV(const double u, const double v) const { return ZAxis; }
	Point PartialDerivativeUU(const double u, const double v) const { return -std::cos(u)*R*XAxis - std::sin(u)*R*YAxis; }
	Point PartialDerivativeUV(const double u, const double v) const { return Point(0, 0, 0); }
	Point PartialDerivativeVV(const double u, const double v) const { return Point(0, 0, 0); }
private:
	Point Location, XAxis, YAxis, ZAxis;
	double R;
};

class ConicalType :public GeometryType
{
public:
	ConicalType() {};  //怎么赋给私有变量
	~ConicalType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return (R + v * std::sin(Ang)) * (-std::sin(u)*Xdir + std::cos(u)*Ydir);}
	Point PartialDerivativeV(const double u, const double v) const { return std::sin(Ang)*(std::cos(u)*Xdir + std::sin(u)*Ydir)+std::cos(Ang)*Zdir;}
	Point PartialDerivativeUU(const double u, const double v) const { return (R + v * std::sin(Ang)) * (-std::cos(u)*Xdir - std::sin(u)*Ydir); }
	Point PartialDerivativeUV(const double u, const double v) const { return std::sin(Ang) * (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeVV(const double u, const double v) const { return Point(0, 0, 0); }
private:
	Point Origin, Xdir, Ydir, Zdir;
	double R, Ang;
}; 

class SphericalType :public GeometryType
{
public:
	SphericalType() {};  //怎么赋给私有变量
	~SphericalType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return R * std::cos(v)* (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeV(const double u, const double v) const { return -R * std::sin(v)*(std::cos(u)*Xdir + std::sin(u)*Ydir) + R*std::cos(v)*Zdir; }
	Point PartialDerivativeUU(const double u, const double v) const { return R * std::cos(v)* (-std::cos(u)*Xdir - std::sin(u)*Ydir); }
	Point PartialDerivativeUV(const double u, const double v) const { return -R * std::sin(v) * (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeVV(const double u, const double v) const { return -R * std::cos(v)*(std::cos(u)*Xdir + std::sin(u)*Ydir)-std::pow(R, 2)*cos(2*v)*Zdir; }
private:
	Point Origin, Xdir, Ydir, Zdir;
	double R;
};


class ToroidalType :public GeometryType
{
public:
	ToroidalType() {};  //怎么赋给私有变量
	~ToroidalType() {};
public:
	Point PartialDerivativeU(const double u, const double v) const { return (R+r * std::cos(v))* (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeV(const double u, const double v) const { return -r * std::sin(v)*(std::cos(u)*Xdir + std::sin(u)*Ydir) + r * std::cos(v)*Zdir; }
	Point PartialDerivativeUU(const double u, const double v) const { return (R + r * std::cos(v)) * (-std::cos(u)*Xdir - std::sin(u)*Ydir); }
	Point PartialDerivativeUV(const double u, const double v) const { return -r * std::sin(v) * (-std::sin(u)*Xdir + std::cos(u)*Ydir); }
	Point PartialDerivativeVV(const double u, const double v) const { return -r * std::cos(v)*(std::cos(u)*Xdir + std::sin(u)*Ydir) - std::pow(r, 2)*cos(2 * v)*Zdir; }
private:
	Point Origin, Xdir, Ydir, Zdir;
	double R, r;
};


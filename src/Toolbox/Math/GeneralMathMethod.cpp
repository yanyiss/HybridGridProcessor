#include "GeneralMathMethod.h"
#include <iostream>


namespace GeneralMathMethod {
	Eigen::Vector2d ComputePolygonInteriorPoint(Polygon &polygon)
	{
		int n = polygon.cols();
		double alpha = cos(PI * (n - 2) / n);

		int convex_point = 0;
		Vector2d frac_ray;
		while (convex_point < n)
		{
			Vector2d v0 = (polygon.col((convex_point + 1) % n) - polygon.col(convex_point)).normalized();
			Vector2d v1 = (polygon.col((convex_point + 2) % n) - polygon.col((convex_point + 1) % n)).normalized();
			if (v0(0)*v1(1) - v0(1)*v1(0) > 0)
			{
				if (-v0.dot(v1) >= alpha)
				{
					frac_ray = v1 - v0;
					break;
				}
			}
			convex_point++;
		}
		convex_point = (convex_point + 1) % n;
		Vector2d pc = polygon.col(convex_point);
		double c = frac_ray(0)*pc(1) - frac_ray(1)*pc(0);
		Vector2d intersect_point;
		double dis = DBL_MAX;
		for (int i = 0; i < n; i++) {
			if (i == convex_point || (i + 1) % n == convex_point) continue;
			auto &p0 = polygon.col(i);
			auto &p1 = polygon.col((i + 1) % n);
			if ((frac_ray(1)*p0(0) - frac_ray(0)*p0(1) + c)*(frac_ray(1)*p1(0) - frac_ray(0)*p1(1) + c) <= 0) {
				Matrix2d A(2, 2);
				A << frac_ray(1), -frac_ray(0), p0(1) - p1(1), p1(0) - p0(0);
				Vector2d b(2);
				b << -c, p1(0)*p0(1) - p0(0)*p1(1);
				b = A.inverse()*b;
				if (frac_ray(0)*(b(0) - pc(0)) > 0 || frac_ray(1)*(b(1) - pc(1)) > 0) {
					double d = (b - pc).norm();
					if (d < dis) {
						dis = d;
						intersect_point = b;
					}
				}
			}
		}
		return (pc + intersect_point) / 2;
	}

	bool IsInPolygon(Polygon &polygon, Eigen::Vector2d &p)
	{
		int n = polygon.cols();
		double sum = 0;
		for (int i = 0; i < n; i++) {
			//auto s = Vector2d(polygon(0, i), polygon(1, i)) - p;
			//auto t = Vector2d(polygon(0, (i + 1) % n), polygon(1, (i + 1) % n)) - p;
			auto s = polygon.col(i) - p;
			auto t = polygon.col((i + 1) % n) - p;
			sum += atan2(-s(1)*t(0) + s(0)*t(1), s(0)*t(0) + s(1)*t(1));
		}
		//cout << "the sum of angles: " << sum << endl;
		//if(sum>PI) cout << "the sum of angles: " << sum << endl;
		if (sum > PI)
			return true;
		else
			return false;
	}

	double ComputePolygonArea(Polygon &polygon)
	{
		double area = 0;
		int n = polygon.cols();
		for (int i = 0; i < n; i++) {
			auto &a = polygon.col(i);
			auto &b = polygon.col((i + 1) % n);
			//if (a(1)*b(0) > a(0)*b(1)) area -= Vector3d(a(0), a(1), 0).cross(Vector3d(b(0), b(1), 0)).norm();
			//else area += Vector3d(a(0), a(1), 0).cross(Vector3d(b(0), b(1), 0)).norm();
			area -= a(0)*b(1) - a(1)*b(0);
		}
		return std::fabs(area)*0.5;
	}


	double ComputePolygonPerimeter(Polygon &polygon)
	{
		double perimeter = (polygon.col(polygon.cols() - 1) - polygon.col(0)).norm();
		for (int i = 0; i < polygon.cols() - 1; i++) perimeter += (polygon.col(i) - polygon.col(i + 1)).norm();
		return perimeter;
	}

	Eigen::Vector4d ComputePolygonBound(Polygon &polygon)
	{
		Vector4d bound(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
		int n = polygon.cols();
		for (int i = 0; i < n; i++) {
			auto &pos = polygon.col(i);
			bound(0) = min(bound(0), pos(0));
			bound(1) = max(bound(1), pos(0));
			bound(2) = min(bound(2), pos(1));
			bound(3) = max(bound(3), pos(1));
		}
		return bound;
	}

	int ComputePolygonsOuterBoundIndex(Polygons &polygons)
	{
		if (polygons.size() == 1) return 0;
		for (int i = 0; i < polygons.size(); i++) {
			auto &pi = polygons[i];
			int n = pi.cols();
			double sum = 0;
			for (int j = 0; j < n; j++) {
				Vector2d p0 = pi.col((j + 2) % n) - pi.col((j + 1) % n);
				Vector2d p1 = pi.col((j + 1) % n) - pi.col(j);
				sum += atan2(p0(0)*p1(1) - p0(1)*p1(0), p0(0)*p1(0) + p0(1)*p1(1));
			}
			cout << "the sum should be 2*PI or -2*PI: " << sum << endl;
			if (sum > 0) return i;
		}
		return -1;
	}

	bool IsInPolygons(Polygons &polygons, Eigen::Vector2d &p)
	{
		for (int i = 0; i < polygons.size(); i++)
			if (!IsInPolygon(polygons[i], p))
				return false;
		return true;
	}

	double ComputePolygonsArea(Polygons &polygons)
	{
		double area = ComputePolygonArea(polygons[0]);
		for (int i = 1; i < polygons.size(); i++)
			area -= ComputePolygonArea(polygons[i]);
		return area;
	}

	void Find_Span(Matrix2Xd& dir, Matrix2Xd &para, bool direction, int &id, int num, double height)
	{
		id = -1;
		Vector2d p0, dir1, dir2;
		if (direction)
		{
			p0 = dir.col(dir.cols() - 1);
			dir1 = (dir.col(dir.cols() - 2) - p0).normalized();
			for (int i = 1; i < para.cols(); i++)
			{
				dir2 = para.col(i) - p0;
				if (abs(dir1[0] * dir2[1] - dir1[1] * dir2[0]) < height) continue;
				id = i - 1;
				break;
			}
			if (id < 0) return;
			for (int i = 0; i <= num; i++)
			{
				dir2 = ((num - i) * para.col(id) + i * para.col(id + 1))/num - p0;
				if (abs(dir1[0] * dir2[1] - dir1[1] * dir2[0]) < height) continue;
				id = i + id * num;
				break;
			}
		}
		else
		{
			p0 = dir.col(0);
			dir1 = (dir.col(1) - p0).normalized();
			for (int i = para.cols() - 2; i >= 0; i--)
			{
				dir2 = para.col(i) - p0;
				if (abs(dir1[0] * dir2[1] - dir1[1] * dir2[0]) < height) continue;
				id = i;
				break;
			}
			if (id < 0) return;
			for (int i = num; i >= 0; i--)
			{
				dir2 = ((num - i) * para.col(id) + i * para.col(id + 1)) / num - p0;
				if (abs(dir1[0] * dir2[1] - dir1[1] * dir2[0]) < height) continue;
				id = i + id * num;
				break;
			}
		}
	}

	/*void Find_Span(Matrix2Xd& dir, Matrix2Xd &para, bool direction, int &id, int num, double height, GeometryType *Surface)
	{
		id = -1;
		Vector2d p1, p2;
		if (direction)
		{
			p1 = dir.col(dir.cols() - 1);
			for (int i = 1; i < para.cols(); i++)
			{
				p2 = para.col(i);
				if (Riemanlen(Surface, p1, p2) < height) continue;
				id = i - 1;
				break;
			}
			if (id < 0) return;
			for (int i = 0; i <= num; i++)
			{
				p2 = ((num - i) * para.col(id) + i * para.col(id + 1)) / num;
				if (Riemanlen(Surface, p1, p2) < height) continue;
				id = i + id * num;
				break;
			}
		}
		else
		{
			p1 = dir.col(0);
			for (int i = para.cols() - 2; i >= 0; i--)
			{
				p2 = para.col(i);
				if (Riemanlen(Surface, p1, p2) < height) continue;
				id = i;
				break;
			}
			if (id < 0) return;
			for (int i = num; i >= 0; i--)
			{
				p2 = ((num - i) * para.col(id) + i * para.col(id + 1)) / num;
				if (Riemanlen(Surface, p1, p2) < height) continue;
				id = i + id * num;
				break;
			}
		}
	}*/

	/*double Riemanlen(GeometryType *Surface, Vector2d &p1, Vector2d &p2)
	{
		auto UV = Surface->Getbounds();
		double u1 = UV(0), u2 = UV(1), v1 = UV(2), v2 = UV(3);
		if (p1(0) < u1 || p1(0) > u2 || p1(1) < v1 || p1(1) > v2)
			return (p1 - p2).norm();
		if (p2(0) < u1 || p2(0) > u2 || p2(1) < v1 || p2(1) > v2)
			return (p1 - p2).norm();
		Eigen::MatrixXd xy = Eigen::MatrixXd::Zero(1, 2);
		xy << p1[0] - p2[0], p1[1] - p2[1];
		Eigen::Matrix2d M;
		Point ru = Surface->PartialDerivativeU(p1(0), p1(1));
		Point rv = Surface->PartialDerivativeV(p1(0), p1(1));
		double E = ru.dot(ru);
		double F = ru.dot(rv);
		double G = rv.dot(rv);
		ru = Surface->PartialDerivativeU(p2(0), p2(1));
		rv = Surface->PartialDerivativeV(p2(0), p2(1));
		E += ru.dot(ru);
		F += ru.dot(rv);
		G += rv.dot(rv);
		M << E, F,
			F, G;
		return pow(abs((xy * (M *0.5) * xy.transpose()).determinant()), 0.5);
	}*/

	double Binomial(int n, double p, int m)
	{
		if (p < 0 || p > 1 || m > n || m < 0) return 0;
		double a = 1, b = 1;
		for (int i = m + 1; i <= n; i++)
		{
			a *= i;
		}
		for (int i = 1; i <= n - m; i++)
		{
			b *= i;
		}
		return (a / b) * pow(p, m) * pow(1 - p, n - m);
	}

	void DataSet(Point4& UV, double eps, Matrix2Xd& pnts)
	{
		double u1 = UV(0), u2 = UV(1), v1 = UV(2), v2 = UV(3);
		for (int j = 0; j < pnts.cols(); j++)
		{
			double &u = pnts(0, j);
			double &v = pnts(1, j);
			if (u <= u1) u = u1 + (u2 - u1)*eps;
			else if (u >= u2) u = u2 - (u2 - u1)*eps;
			if (v <= v1) v = v1 + (v2 - v1)*eps;
			else if (v >= v2) v = v2 - (v2 - v1)*eps;
		}
	}

#ifdef OPENMESH_TRIMESH_ARRAY_KERNEL_HH
	double ComputeVoronoiArea(TriMesh* mesh, OpenMesh::VertexHandle v)//v不允许是边界点
	{
		double area = 0;
		for (auto tvoh : mesh->voh_range(v)) {
			double x = pow(mesh->calc_edge_length(tvoh), 2);
			double y = pow(mesh->calc_edge_length(mesh->next_halfedge_handle(tvoh)), 2);
			double z = pow(mesh->calc_edge_length(mesh->prev_halfedge_handle(tvoh)), 2);
			if (x + z <= y) area += 0.5*mesh->calc_face_area(mesh->face_handle(tvoh));
			else if (x + y <= z || y + z <= x) area += 0.25*mesh->calc_face_area(mesh->face_handle(tvoh));
			else {
				OpenMesh::Vec3d xv = mesh->calc_edge_vector(tvoh);
				OpenMesh::Vec3d yv = mesh->calc_edge_vector(mesh->next_halfedge_handle(tvoh));
				OpenMesh::Vec3d zv = mesh->calc_edge_vector(mesh->prev_halfedge_handle(tvoh));
				area -= (x*yv.dot(zv) / yv.cross(zv).norm() + z * xv.dot(yv) / xv.cross(yv).norm()) * 0.125;
			}
		}
		return area;
	}

	double ComputeGaussCurvature(TriMesh* mesh, OpenMesh::VertexHandle v)//v不允许是边界点
	{
		double curvature = 2 * PI;
		for (auto tvih : mesh->vih_range(v)) curvature -= mesh->calc_sector_angle(tvih);
		return curvature / ComputeVoronoiArea(mesh, v);
	}

	double ComputeMeanCurvature(TriMesh* mesh, OpenMesh::VertexHandle v)//v不允许是边界点
	{
		OpenMesh::Vec3d curvature(0.0, 0.0, 0.0);
		for (auto tvoh : mesh->voh_range(v))
			curvature += mesh->calc_edge_vector(tvoh) * (1.0 / tan(mesh->calc_sector_angle(mesh->next_halfedge_handle(tvoh))) +
				1.0 / tan(mesh->calc_sector_angle(mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(tvoh)))));
		return curvature.norm() / (4.0*ComputeVoronoiArea(mesh, v));
	}
#endif OPENMESH_TRIMESH_ARRAY_KERNEL_HH
}
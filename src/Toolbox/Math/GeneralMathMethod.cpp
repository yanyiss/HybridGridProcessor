#include "GeneralMathMethod.h"
#include <iostream>

namespace GeneralMathMethod {
	Eigen::Vector2d ComputePolygonInteriorPoint(Polygon &polygon)
	{
		int n = polygon.rows();
		double alpha = cos(PI * (n - 2) / n);

		int convex_point = 0;
		Vector2d frac_ray;
		while (convex_point < n)
		{
			Vector2d v0 = (polygon.row((convex_point + 1) % n) - polygon.row(convex_point)).normalized();
			Vector2d v1 = (polygon.row((convex_point + 2) % n) - polygon.row((convex_point + 1) % n)).normalized();
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
		Vector2d pc = polygon.row(convex_point);
		double c = frac_ray(0)*pc(1) - frac_ray(1)*pc(0);
		Vector2d intersect_point;
		double dis = DBL_MAX;
		for (int i = 0; i < n; i++) {
			if (i == convex_point || (i + 1) % n == convex_point) continue;
			auto &p0 = polygon.row(i);
			auto &p1 = polygon.row((i + 1) % n);
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
		int n = polygon.rows();
		double sum = 0;
		for (int i = 0; i < n; i++) {
			auto s = Vector2d(polygon(i, 0), polygon(i, 1)) - p;
			auto t = Vector2d(polygon((i + 1) % n, 0), polygon((i + 1) % n, 1)) - p;
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
		int n = polygon.rows();
		auto x = Vector2d::Zero();
		for (int i = 0; i < n; i++) {
			auto a = polygon.row(i);
			auto b = polygon.row((i + 1) % n);
			if (a(1)*b(0) > a(0)*b(1)) area -= Vector3d(a(0), a(1), 0).cross(Vector3d(b(0), b(1), 0)).norm();
			else area += Vector3d(a(0), a(1), 0).cross(Vector3d(b(0), b(1), 0)).norm();
		}
		return area / 2.0;
	}


	double ComputePolygonPerimeter(Polygon &polygon)
	{
		double perimeter = (polygon.row(polygon.rows() - 1) - polygon.row(0)).norm();
		for (int i = 0; i < polygon.rows() - 1; i++) perimeter += (polygon.row(i) - polygon.row(i + 1)).norm();
		return perimeter;
	}

	Eigen::Vector4d ComputePolygonBound(Polygon &polygon)
	{
		Vector4d bound(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
		int n = polygon.rows();
		for (int i = 0; i < n; i++) {
			auto pos = polygon.row(i);
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
			int n = pi.rows();
			double sum = 0;
			for (int j = 0; j < n; j++) {
				Vector2d p0 = pi.row((j + 2) % n) - pi.row((j + 1) % n);
				Vector2d p1 = pi.row((j + 1) % n) - pi.row(j);
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

#ifdef OPENMESH_TRIMESH_ARRAY_KERNEL_HH
	double ComputeVoronoiArea(TriMesh* mesh, OpenMesh::VertexHandle v)//v不允许是边界点
	{
		double area = 0;
		/*for (auto tvoh : mesh->voh_range(v)) {
			double x = mesh->calc_edge_length(tvoh);
			double y = mesh->calc_edge_length(mesh->next_halfedge_handle(tvoh));
			double z = mesh->calc_edge_length(mesh->prev_halfedge_handle(tvoh));
			double cosxz = (x*x + z * z - y * y) / (2 * x*z);
			if (cosxz < 0) { area += 0.25*x*z*sqrt(1.0 - cosxz * cosxz); cout << 0.25*x*z*sqrt(1.0 - cosxz * cosxz) << endl; continue; }
			else {
				double cosxy = (x*x + y * y - z * z) / (2 * x*y);
				if (cosxy < 0) { area += 0.125*x*z*sqrt(1.0 - cosxz * cosxz); continue; }
				double cosyz = (y*y + z * z - x * x) / (2 * y*z);
				if (cosyz < 0) { area += 0.125*x*z*sqrt(1.0 - cosxz * cosxz); continue; }
			}
			OpenMesh::Vec3d xv = mesh->calc_edge_vector(tvoh);
			OpenMesh::Vec3d yv = mesh->calc_edge_vector(mesh->next_halfedge_handle(tvoh));
			OpenMesh::Vec3d zv = mesh->calc_edge_vector(mesh->prev_halfedge_handle(tvoh));
			area -= (x*x*yv.dot(zv) / yv.cross(zv).norm() + z * z*xv.dot(yv) / xv.cross(yv).norm()) * 0.125;
		}*/
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
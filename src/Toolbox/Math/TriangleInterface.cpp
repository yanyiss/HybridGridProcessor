#include "TriangleInterface.h"
#include <iostream>
void triangulate(const Eigen::MatrixX2d &all_pts, const Eigen::MatrixX2i &E, const double area_threshold, TriMesh &mesh)
{
	int p_num = all_pts.rows();
	triangulateio in;

	in.numberofpoints = p_num;
	in.pointlist = (double *)calloc(all_pts.size(), sizeof(double));

	for (size_t i = 0; i < p_num; i++)
	{
		in.pointlist[2 * i] = all_pts(i, 0);
		in.pointlist[2 * i + 1] = all_pts(i, 1);
	}

	in.numberofpointattributes = 0;
	in.pointmarkerlist = (int *)calloc(p_num, sizeof(int));
	for (unsigned i = 0; i < p_num; ++i) in.pointmarkerlist[i] = 1;

	in.holelist = NULL;
	in.numberofholes = 0;
	in.trianglelist = NULL;
	in.numberoftriangles = 0;
	in.numberofcorners = 0;
	in.numberoftriangleattributes = 0;
	in.triangleattributelist = NULL;

	in.numberofsegments = E.size() ? E.rows() : 0;
	in.segmentlist = (int*)calloc(E.size(), sizeof(int));

	for (size_t i = 0; i < E.rows(); i++)
	{
		in.segmentlist[2 * i] = E(i, 0);
		in.segmentlist[2 * i + 1] = E(i, 1);
	}

	in.segmentmarkerlist = (int*)calloc(E.rows(), sizeof(int));
	for (unsigned i = 0; i < E.rows(); ++i) in.segmentmarkerlist[i] = 1;

	in.numberofholes = 0;
	in.numberofregions = 0;

	triangulateio out;
	out.pointlist = NULL;
	out.trianglelist = NULL;
	out.segmentlist = NULL;
	out.segmentmarkerlist = NULL;
	out.pointmarkerlist = NULL;
	std::string full_flags;
	
	if (area_threshold < 1.0e-15)
	{
		full_flags = "qYYQpz";
	}
	else
	{
		std::stringstream stream{};
		stream << std::fixed << std::setprecision(15) << area_threshold;
		full_flags = "qYYa" + stream.str() + "Qpz";
	}
	triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);


	mesh.clear();
	std::vector<TriMesh::VertexHandle> vhandle;
	for (size_t i = 0; i < out.numberofpoints; i++)
	{
		vhandle.push_back(mesh.add_vertex(TriMesh::Point(out.pointlist[2 * i], out.pointlist[2 * i + 1], 0)));
	}
	for (size_t i = 0; i < out.numberoftriangles; i++)
	{
		mesh.add_face(vhandle[out.trianglelist[3 * i]], vhandle[out.trianglelist[3 * i + 1]], vhandle[out.trianglelist[3 * i + 2]]);
	}

	// Cleanup in
	free(in.pointlist);
	free(in.pointmarkerlist);
	free(in.segmentlist);
	free(in.segmentmarkerlist);
	free(in.holelist);
	// Cleanup out
	free(out.pointlist);
	free(out.trianglelist);
	free(out.segmentlist);
	free(out.segmentmarkerlist);
	free(out.pointmarkerlist);
}

//若area_threshold为负值，则不允许在边界增加新点
//E[0]是外边界，其他是洞，洞内不生成网格
void triangulate(const Eigen::MatrixX2d &all_pts, const std::vector<Eigen::MatrixX2i> &E, const double area_threshold, TriMesh &mesh)
{
	if (E.size() == 1)
	{
		triangulate(all_pts, E[0], area_threshold, mesh);
		return;
	}
	int p_num = all_pts.rows();
	triangulateio in;

	in.numberofpoints = p_num;
	in.pointlist = (double *)calloc(all_pts.size(), sizeof(double));

	for (size_t i = 0; i < p_num; i++)
	{
		in.pointlist[2 * i] = all_pts(i, 0);
		in.pointlist[2 * i + 1] = all_pts(i, 1);
	}

	in.numberofpointattributes = 0;
	in.pointmarkerlist = (int *)calloc(p_num, sizeof(int));
	for (unsigned i = 0; i < p_num; ++i) in.pointmarkerlist[i] = 1;

	in.trianglelist = NULL;
	in.numberoftriangles = 0;
	in.numberofcorners = 0;
	in.numberoftriangleattributes = 0;
	in.triangleattributelist = NULL;

	int sum = 0;
	for (auto itr = E.begin(); itr != E.end(); itr++)
	{
		sum += (*itr).rows();
	}
	in.numberofsegments = sum;
	in.segmentlist = (int*)calloc(sum * 2, sizeof(int));

	sum = 0;
	for (auto itr = E.begin(); itr != E.end(); itr++)
	{
		for (size_t i = 0; i < (*itr).rows(); i++)
		{
			in.segmentlist[2 * sum] = (*itr)(i, 0);
			in.segmentlist[2 * sum + 1] = (*itr)(i, 1);
			sum++;
		}
	}
	in.segmentmarkerlist = (int*)calloc(sum, sizeof(int));
	for (unsigned i = 0; i < sum; ++i) in.segmentmarkerlist[i] = 1;

	sum = E.size() - 1;
	in.numberofholes = sum;
	in.holelist = (double*)calloc(sum * 2, sizeof(double));
	//由x坐标跨度判断外边界的序号，记为outer_flag
	std::vector<double> x_scale;
	for (int i = 0; i < E.size(); i++) {
		auto &e = E[i];
		double x_max = DBL_MIN;
		double x_min = DBL_MAX;
		for (int j = 0; j < e.rows(); j++) {
			double p = all_pts(e(j, 0), 0);
			x_max = std::max(x_max, p);
			x_min = std::min(x_min, p);
		}
		x_scale.push_back(x_max - x_min);
	}
	int outerflag = -1;
	double x_scale_max = DBL_MIN;
	for (int i = 0; i < x_scale.size(); i++) {
		if (x_scale[i] > x_scale_max) {
			x_scale_max = x_scale[i];
			outerflag = i;
		}
	}

	//根据逆向多边形外角和为360度判断外边界，记为outer_flag
	int itertimes = 0;
	for (int i = 0; i < E.size(); i++) {
		if (i == outerflag) continue;
		int r = E[i].rows() - 1;
		GeneralMathMethod::Polygon polygon(E[i].rows(), 2);
		for (int j = 0; j < E[i].rows(); j++) {
			polygon.row(j) = all_pts.row(E[i](r--, 0));
		}
		auto interior_point = GeneralMathMethod::ComputePolygonInteriorPoint(polygon);
		in.holelist[2 * itertimes] = interior_point(0);
		in.holelist[2 * itertimes + 1] = interior_point(1);
		itertimes++;
	}

	in.numberofregions = 0;

	triangulateio out;
	out.pointlist = NULL;
	out.trianglelist = NULL;
	out.segmentlist = NULL;
	out.segmentmarkerlist = NULL;
	out.pointmarkerlist = NULL;
	std::string full_flags;
	if (area_threshold < 1.0e-15) {
		full_flags = "qYYQpz";
	}
	else {
		std::stringstream stream{};
		stream << std::fixed << std::setprecision(15) << area_threshold;
		full_flags = "qYYa" + stream.str() + "Qpz";
	}
	triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);


	mesh.clear();
	std::vector<TriMesh::VertexHandle> vhandle;
	for (size_t i = 0; i < out.numberofpoints; i++)
	{
		vhandle.push_back(mesh.add_vertex(TriMesh::Point(out.pointlist[2 * i], out.pointlist[2 * i + 1], 0)));
	}
	for (size_t i = 0; i < out.numberoftriangles; i++)
	{
		mesh.add_face(vhandle[out.trianglelist[3 * i]], vhandle[out.trianglelist[3 * i + 1]], vhandle[out.trianglelist[3 * i + 2]]);
	}

	// Cleanup in
	free(in.pointlist);
	free(in.pointmarkerlist);
	free(in.segmentlist);
	free(in.segmentmarkerlist);
	free(in.holelist);
	// Cleanup out
	free(out.pointlist);
	free(out.trianglelist);
	free(out.segmentlist);
	free(out.segmentmarkerlist);
	free(out.pointmarkerlist);
}
#pragma once
#include <Eigen/Core>
#include "..\src\Dependency\triangle\triangle.h"
#include <sstream>
#include <iomanip>
#include "GeneralMathMethod.h"

void triangulate(const Eigen::MatrixX2d &all_pts, const Eigen::MatrixX2i &E, const double area_threshold, TriMesh &mesh);

void triangulate(const Eigen::MatrixX2d &all_pts, const std::vector<Eigen::MatrixX2i> &E, const double area_threshold, TriMesh &mesh);
#include"src\Dependency\BSpline\GeometryType.h"
#include"src\MeshViewer\MeshDefinition.h"
#include"src\Toolbox\dprinter\dprint.h"


class Riemannremesh
{
public:
	Riemannremesh( GeometryType* surface, TriMesh *paramesh)
	{
		B = surface;
		mesh = paramesh;
	}
	Eigen::Matrix2d Riemanndata(const TriMesh::Point &p);
	double Riemannlen(const TriMesh::VertexHandle &h1, const TriMesh::VertexHandle &h2);
	void split();
	void collapse();
	void flip();
	void updatepoint();
	void calulenth();
	void remesh();

	double highlenth;
	double lowlenth;
	GeometryType *B;
	TriMesh *mesh;
};




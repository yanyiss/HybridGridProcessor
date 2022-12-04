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
	bool if_reverse_face(TriMesh::FaceHandle& f);
	void split();
	void collapse();
	void equalize_valence();
	void updatepoint();
	void calclength();
	void remesh();

	double highlength;
	double lowlength;
	bool direction;
	GeometryType *B;
	TriMesh *mesh;
};




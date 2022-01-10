#ifndef OCCREADER_H
#define OCCREADER_H
#include "occheader.h"
#include "basic_def.h"
class OccReader {
public:
	OccReader(QString &fileName)
	{
		std::string filetype;
		if (fileName.endsWith(".stp") || fileName.endsWith(".STEP")) {
			reader = new STEPControl_Reader();
			dprint("CAD model from STEP file");
			filetype = "STEP";
		}
		else if (fileName.endsWith(".igs") || fileName.endsWith(".IGES")) {
			reader = new IGESControl_Reader();
			dprint("CAD model from IGES file");
			filetype = "IGES";
		}
		else
		{
			dprinterror("Is anything wrong? It can't be other file types");
			exit(0);
		}

		dprint(filetype + "file read beginning!\n");
		reader->ReadFile(fileName.toLatin1().data());
		//reader->TransferRoots();
		//Standard_Integer NbRoots = reader->NbRootsForTransfer();
		// gets the number of transferable roots
		Standard_Integer NbTrans = reader->TransferRoots();
		globalmodel.aShape = reader->OneShape();
		//aShape = reader->OneShape();
		dprint(filetype + "file read finished!\n");

		ComputeFaceAndEdge();
		Discrete_Bound();
	}
	~OccReader() {
		if (reader) { delete reader; reader = nullptr; }
	}

protected:
	XSControl_Reader *reader;
	

public:
	double expected_edge_length;
	double mu = 1.2;//三角形面积允许扩张系数
	vector<TriMesh> Surface_TriMeshes;
	vector<PolyMesh> Surface_PolyMeshes;

	void ComputeFaceAndEdge();
	void Discrete_Bound();
	void Set_TriMesh();
	//void Set_PolyMesh();

	//vector<ShapeFace> faceshape;
	//vector<ShapeEdge> edgeshape;
private:
};

#endif // !OCCREADER_H

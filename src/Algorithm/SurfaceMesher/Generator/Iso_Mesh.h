#pragma once
#include "OccReader.h"
class Iso_Mesh
{
public:
	Iso_Mesh(QString &fileName);
	~Iso_Mesh()
	{
		if (occ_reader) { delete occ_reader; occ_reader = nullptr; }
	};

private:
	//double expected_area;
	OccReader *occ_reader = nullptr;
	enum MergeType {
		Tri,Quad,Poly
	};
	void InitTree();


public:
	void MergeModel(MergeType mt);
	int EndId(vector<ShapeEdge> &edgeshape, int edge_id);
	void ResetFeature();

	std::string Write_Obj(TriMesh &aMesh);
	void Open_File(std::ofstream &file_writer);
};


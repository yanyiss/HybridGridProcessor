#pragma once
#include "OccReader.h"

namespace CADMesher
{
	class Iso_Mesh
	{
	public:
		Iso_Mesh(QString &fileName);
		~Iso_Mesh()
		{
			if (occ_reader) { delete occ_reader; occ_reader = nullptr; }
		};

	private:
		OccReader *occ_reader = nullptr;
		void InitTree();


	public:
		void MergeModel();
		int EndId(vector<ShapeEdge> &edgeshape, int edge_id);
		void ResetFeature();

		std::string Write_Obj(TriMesh &aMesh);
		void Open_File(std::ofstream &file_writer);
	};
}

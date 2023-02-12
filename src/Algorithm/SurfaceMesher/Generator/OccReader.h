#ifndef OCCREADER_H
#define OCCREADER_H
#include "occheader.h"
#include "basic_def.h"
#include "DomainRemesh.h"
#include "..\src\Toolbox\Math\GeneralMathMethod.h"

namespace CADMesher
{
	class OccReader {
	public:
		explicit OccReader(QString &file);
		OccReader(const OccReader& or) = delete;
		~OccReader() {
			if (reader) 
			{ 
				delete reader; reader = nullptr; 
			}
		}

	protected:
		XSControl_Reader *reader;


	public:
		OpenMesh::Vec3d bbmin, bbmax;
		double expected_edge_length;
		double mu = 1.2;        //三角形面积允许扩张系数
		double epsratio = 0.005;//网格和曲面的误差系数
		int offset_quad_num = 15;
		double offset_initial_ratio = 0.01;
		double offset_increase_ratio = 1.3;
		vector<TriMesh> Surface_TriMeshes;
		vector<PolyMesh> Surface_PolyMeshes;
		vector<int> quad_mun_info;
		vector<double> initial_ratio_info;
		vector<double> increase_ratio_info;

		void SetShape();
		void SetCADWireFrame();
		void ComputeFaceAndEdge();
		void Discrete_Edge();
		bool ProcessTangentialBoundary(int fid, int bid);
		void ClearBoundary(TriMesh &tm);
		bool Face_type();
		void C0_Feature();
		void Curvature_Feature();
		void Reset_feature();
		void Set_TriMesh();
		void Offset_lines(Matrix2Xd &parameters, vector<Matrix2Xd> &offset_pnts, int begin, int pntnum, int quadnum);
		void Set_PolyMesh();
		void Set_Offset_Grid();
		void Re_discrete(ShapeEdge &edge, int id, int discrete_num, int quad_num, bool direction);
		bool If_decline(vector<double>& curvature);
		void OccReader::MergeModel(TriMesh &model_mesh, std::vector<TriMesh> &surface_meshes);
		void OccReader::MergeModel(PolyMesh &model_mesh, std::vector<PolyMesh> &surface_meshes);
		int EndId(vector<ShapeEdge> &edgeshape, int edge_id);
		void SetTriFeature();
		void SetPolyFeature();

		template<typename T>
		void Set_Curvature(GeometryType * srf, T & aMesh)
		{
			double k1, k2;
			for (auto v : aMesh.vertices())
			{
				auto p = aMesh.point(v);
				if(srf->PrincipalCurvature(p[0], p[1], k1, k2))
					aMesh.data(v).GaussCurvature = std::max(std::fabs(k1), std::fabs(k2));		
				else aMesh.data(v).GaussCurvature = -1;
			}
			for (auto v : aMesh.vertices())
			{
				if (aMesh.data(v).GaussCurvature >= 0) continue;
				int count = 0;
				k1 = 0;
				for (auto vv : aMesh.vv_range(v))
				{
					if (aMesh.data(vv).GaussCurvature < 0) continue;
					k1 += aMesh.data(vv).GaussCurvature;
					count++;
				}
				if (!count)
				{
					dprint("Cannot compute the curvature of this surface");
					system("pause");
				}
				aMesh.data(v).GaussCurvature = k1 / count;
			}
		}


	//private:
		QString fileName;
		double initialRate = 0.004;
		//double degeneratedRate = 0.02;
	};

}
#endif // !OCCREADER_H


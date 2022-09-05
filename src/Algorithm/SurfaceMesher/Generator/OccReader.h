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
		explicit OccReader(QString &fileName)
		{
			std::string filetype;
			if (fileName.endsWith(".stp") || fileName.endsWith(".step") || fileName.endsWith(".STP") || fileName.endsWith(".STEP")) {
				reader = new STEPControl_Reader();
				dprint("CAD model from STEP file");
				filetype = "STEP";
			}
			else if (fileName.endsWith(".igs") || fileName.endsWith(".iges") || fileName.endsWith(".IGS") || fileName.endsWith(".IGES")) {
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
			Standard_Integer NbTrans = reader->TransferRoots();
			globalmodel.aShape = reader->OneShape();
			dprint(filetype + "file read finished\n");

			ComputeFaceAndEdge();
			narrow_surface();
			Discrete_Edge();
			Face_type();
			C0_Feature();
			curvature_feature();
		}
		OccReader(const OccReader& or) = delete;
		~OccReader() {
			if (reader) 
			{ 
				delete reader; reader = nullptr; 
				vector<ShapeFace> &faceshape = globalmodel.faceshape;
				for (int i = 0; i < faceshape.size(); i++)
				{
					if (!faceshape[i].if_exisited) continue;
					delete faceshape[i].Surface;
					faceshape[i].Surface = nullptr;
				}
			}
		}

	protected:
		XSControl_Reader *reader;


	public:
		double expected_edge_length;
		double mu = 1.2;        //三角形面积允许扩张系数
		double epsratio = 0.005;//网格和曲面的误差系数
		int offset_quad_num = 15;
		double offset_initial_ratio = 0.01;
		double offset_increase_ratio = 1.2;
		vector<TriMesh> Surface_TriMeshes;
		vector<PolyMesh> Surface_PolyMeshes;

		void ComputeFaceAndEdge();
		void narrow_surface();
		void Discrete_Edge();
		bool ProcessTangentialBoundary(int fid, int bid);
		void ClearBoundary(TriMesh &tm);
		void Face_type();
		void C0_Feature();
		void curvature_feature();
		void Set_TriMesh();
		void Offset_lines(Matrix2Xd &parameters, vector<Matrix2Xd> &offset_pnts, int begin, int pntnum, int quadnum);
		void Set_PolyMesh();
		void Set_Offset_Grid();
		void re_discrete(ShapeEdge &edge, int id, int discrete_num, int quad_num, bool direction);
		bool If_decline(vector<double>& curvature);

		template<typename T>
		void surface_curvature(GeometryType * srf, T & aMesh)
		{
			double k1, k2;
			for (auto v : aMesh.vertices())
			{
				auto p = aMesh.point(v);
				srf->PrincipalCurvature(p[0], p[1], k1, k2);
				//dprint("v.idx()", k1, k2);
				if (isnan(k1*k2))
				{
					dprint("bug", v.idx());
				}
				aMesh.data(v).GaussCurvature = std::max(std::fabs(k1), std::fabs(k2));
			}
		}

	private:
		double initialRate = 0.01;
		//double degeneratedRate = 0.02;
	};

}
#endif // !OCCREADER_H


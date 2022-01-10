#include "OccReader.h"
void OccReader::Set_TriMesh()
{
	vector<ShapeFace> &faceshape = globalmodel.faceshape;
    vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;

	Surface_TriMeshes.resize(faceshape.size());
	for (int i = 0; i < faceshape.size(); i++)
	//int i = 3;
	{
		//print(i);
		TriMesh &aMesh = Surface_TriMeshes[i];
		auto &wires = faceshape[i].wires;
		if (wires.empty())
		{
			continue;
		}

		auto &aface = faceshape[i].face;
		TopLoc_Location loc;
		//auto r = aface.TShape(); 
		Handle(Geom_Surface) asurface = BRep_Tool::Surface(aface, loc);
		Standard_Real x, y, z, w;
		asurface->Bounds(x, y, z, w);


		int pointsnumber = 0;
		vector<MatrixX2i> bnd;
		double x_step = 0;
		double y_step = 0;
		for (int m = 0; m < wires.size(); m++)
		{
			auto &edges = wires[m];
			int boundsum = 0;
			auto &end_paras = edgeshape[*edges.rbegin()].parameters;
			auto end_para = end_paras.row(end_paras.rows() - 1);
			for (int j = 0; j < edges.size(); j++)
			{
				auto &aedge = edgeshape[edges[j]];
				auto &boundpos = aedge.parameters;
				int rows = boundpos.rows() - 1;
				pointsnumber += rows;
				boundsum += rows;

				double temp;
				for (int r = 0; r < rows; r++)
				{
					temp = fabs(boundpos(r + 1, 0) - boundpos(r, 0));
					x_step = std::max(x_step, temp);
					temp = fabs(boundpos(r + 1, 1) - boundpos(r, 1));
					y_step = std::max(y_step, temp);
				}
			}
			MatrixX2i bound(boundsum, 2);
			for (int j = 0; j < boundsum - 1; j++)
			{
				bound.row(j) << j + pointsnumber - boundsum, j + pointsnumber - boundsum + 1;
			}
			bound.row(boundsum - 1) << pointsnumber - 1, pointsnumber - boundsum;
			bnd.push_back(bound);
		}
		

		MatrixX2d all_pnts(pointsnumber, 2);
		pointsnumber = 0;
		int s = 0;
		double if_reverse = aface.Orientation() ? -1.0 : 1.0;

		double ra = x_step / y_step * if_reverse;
		
		for (int j = 0; j < wires.size(); j++)
		{
			auto &edges = wires[j];
			for (int k = 0; k < edges.size(); k++)
			{
				auto &boundpos = edgeshape[edges[k]].parameters;
				int rows = boundpos.rows() - 1;
				/*all_pnts.block(s, 0, rows, 1) = boundpos.block(0, 0, rows, 1);
				all_pnts.block(s, 1, rows, 1) = ra * boundpos.block(0, 1, rows, 1);*/
				all_pnts.block(s, 0, rows, 2) = boundpos.block(0, 0, rows, 2);
				s += rows;
			}
			//cout << "\n\n" << endl;
			pointsnumber = s;
		}
		
		all_pnts.block(0, 1, all_pnts.rows(), 1) *= ra;

		triangulate(all_pnts, bnd, sqrt(3) * x_step * x_step / 4 * mu, aMesh);//sqrt(3) * x_step * x_step / 4 * mu
		//triangulate(all_pnts, bnd, 100, aMesh);
		for (auto tv : aMesh.vertices())
		{
			auto pos = aMesh.point(tv);
			auto v = asurface->Value(pos[0], pos[1] / ra);
			//aMesh.set_point(tv, Mesh::Point(pos[0], pos[1] / ra, 0));
			aMesh.set_point(tv, TriMesh::Point(v.X(), v.Y(), v.Z()));
		}
	}
	dprint("Piecewise TriMesh Done!");


	/*for (auto surface : faceshape)
	{
		TopLoc_Location loca;
		opencascade::handle<Geom_Surface> geom_surface = BRep_Tool::Surface(surface.face, loca);
		opencascade::handle<Standard_Type> type = geom_surface->DynamicType();

		if (type == STANDARD_TYPE(Geom_BSplineSurface)) {

			opencascade::handle<Geom_BSplineSurface> geom_bsplinesurface = Handle(Geom_BSplineSurface)::DownCast(geom_surface);
			TColStd_Array1OfReal uknotsequence = geom_bsplinesurface->UKnotSequence();
			TColStd_Array1OfReal vknotsequence = geom_bsplinesurface->VKnotSequence();
			vector<double> u;
			vector<double> v;
			for (auto itr = uknotsequence.begin(); itr != uknotsequence.end(); itr++)
				u.push_back(*itr);
			for (auto itr = vknotsequence.begin(); itr != vknotsequence.end(); itr++)
				v.push_back(*itr);

			TColgp_Array2OfPnt controlpoints = geom_bsplinesurface->Poles();


			const TColStd_Array2OfReal* weights = geom_bsplinesurface->Weights();
			if (weights) {
				dprint(weights->NbRows(), weights->NbColumns());
				vector<vector<double>> w(weights->NbRows());
				for (int r = 1; r <= weights->NbRows(); r++)
				{
					w[r - 1].reserve(weights->NbColumns());
					for (int c = 1; c <= weights->NbColumns(); c++)
						w[r - 1].push_back(weights->Value(r, c));
				}
				if (u.size() > 8 || v.size() > 8)
				{
					dprint(geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), "\n");
					dprint(u);
					dprint(v);
					dprint();
					for (int r = 1; r <= weights->NbRows(); r++)
						dprint(w[r - 1]);
					dprint();
					for (int r = 1; r <= controlpoints.NbRows(); r++)
						for (int c = 1; c <= controlpoints.NbColumns(); c++)
							dprint(controlpoints.Value(r, c));
				}
			}

		}
	}*/
}


//#include <kt84/util.h>
//#include <patchgen/generate_topology.h>
////#include "decl.h"
////#include "Patch.h"
//#include <patchgen_demo/decl.h>
//#include <patchgen_demo/Patch.h>
//#include <Eigen/Sparse>
//#include <Eigen/SparseLU>
//using namespace std;
//using namespace Eigen;
//using namespace kt84;
//void OccReader::Set_PolyMesh()
//{
//	vector<ShapeFace> &faceshape = globalmodel.faceshape;
//	vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;
//	int nf = faceshape.size();
//	Surface_PolyMeshes.clear();
//	//Surface_PolyMeshes.resize(nf);
//	for (int i = 0; i < nf; i++)
//	{
//		print(i);
//		//Mesh &aMesh = Surface_PolyMeshes[i];
//		auto &wires = faceshape[i].wires;
//		if (wires.empty())
//		{
//			continue;
//		}
//
//		auto &aface = faceshape[i].face;
//		TopLoc_Location loc;
//		//auto r = aface.TShape(); 
//		Handle(Geom_Surface) asurface = BRep_Tool::Surface(aface, loc);
//
//		auto &edges = faceshape[i].wires[0];
//		Eigen::VectorXi l(edges.size());
//		//int sum = 0;
//		for (int i = 0; i < edges.size(); i++)
//		{
//			l(i) = edgeshape[edges[i]].parameters.rows() - 1;
//			//sum += l(i);
//		}
//		//if (edges.size() > 6 || sum%2) continue;
//		//if (edges.size() > 6) continue;
//		if (edges.size() > 6)
//		{
//			int dif = edges.size() - 6;
//			l.resize(6);
//			int id = 0;
//			int i = 0;
//			for (; i < dif; i++)
//			{
//				l(i) = edgeshape[edges[id]].parameters.rows() + edgeshape[edges[id + 1]].parameters.rows() - 2;
//				id += 2;
//			}
//			for (; id < edges.size(); id++)
//			{
//				l(i++) = edgeshape[edges[id]].parameters.rows() - 1;
//			}
//		}
//		if (true)
//		{
//
//		}
//
//		//Surface_PolyMeshes.push_back(PolyMesh());
//		//PolyMesh* aMesh = &(Surface_PolyMeshes.back());
//
//
//		patchgen::PatchParam param;
//		demo::Patch patch;
//		patchgen::generate_topology(l, param, patch);
//
//		// fix boundary vertices
//		demo::Patch::HHandle h;
//		for (auto v : patch.vertices()) {
//			if (patch.data(v).patchgen.corner_index == 0) {
//				h = patch.halfedge_handle(v);
//				break;
//			}
//		}
//
//		int num_sides = l.size();
//		if (edges.size() > 6)
//		{
//
//		}
//		else 
//		{
//			for (int i = 0; i < num_sides; ++i) {
//				auto &aedge = edges[i];
//				for (int j = 0; j < l[i]; ++j) {
//					auto& vdata = patch.data(patch.from_vertex_handle(h)).laplaceDirect;
//					auto r = edgeshape[aedge].parameters.row(j);
//					vdata.value << r(0), r(1), 0;
//					h = patch.prev_halfedge_handle(h);
//				}
//			}
//		}
//		
//
//
//		//PolyMesh &pmesh = globalmodel.initial_polymesh;
//		PolyMesh pmesh;
//		//int nv = pmesh.n_vertices();
//		for (auto tv : patch.vertices())
//		{
//			auto p = patch.data(tv).laplaceDirect.value;
//			pmesh.add_vertex(OpenMesh::Vec3d(p[0],p[1],p[2]));
//		}
//		for (auto tf : patch.faces())
//		{
//			std::vector<OpenMesh::VertexHandle> vh;
//			//vector<int> vh;
//			for (auto tfv : patch.fv_range(tf))
//			{
//				vh.push_back(tfv);
//				//vh.push_back(pmesh.vertex_handle(tfv.idx() + nv));
//			}
//			pmesh.add_face(vh[0], vh[1], vh[2], vh[3]);
//		}
//		vector<Triplet<double>> triple;
//		Eigen::MatrixXd right = Eigen::MatrixXd::Zero(pmesh.n_vertices(), 2);
//
//		for (auto tv : pmesh.vertices())
//		{
//			int v_id = tv.idx();
//			if (pmesh.is_boundary(tv))
//			{
//				auto p = pmesh.point(tv);
//				right.row(v_id) << p[0], p[1];
//				triple.emplace_back(v_id, v_id, 1);
//			}
//			else
//			{
//				triple.emplace_back(v_id, v_id, pmesh.valence(tv));
//				for (auto tvv : pmesh.vv_range(tv))
//				{
//					triple.emplace_back(v_id, tvv.idx(), -1);
//				}
//			}
//		}
//
//		SparseMatrix<double> A(pmesh.n_vertices(), pmesh.n_vertices());
//		A.setFromTriplets(triple.begin(), triple.end());
//		SparseLU<SparseMatrix<double>> solver;
//		solver.compute(A);
//		/*if (!solver.info())
//		{
//			continue;
//		}*/
//		right = solver.solve(right);
//		for (auto tv : pmesh.vertices())
//		{
//			auto p = asurface->Value(right(tv.idx(), 0), right(tv.idx(), 1));
//			//pmesh.set_point(tv, OpenMesh::Vec3d(right(tv.idx(), 0), right(tv.idx(), 1), 0));
//			pmesh.set_point(tv,OpenMesh::Vec3d( p.X(),p.Y(),p.Z()));
//		}
//
//		PolyMesh &pm = globalmodel.initial_polymesh;
//		int nv = pm.n_vertices();
//		for (auto tv : pmesh.vertices())
//		{
//			pm.add_vertex(pmesh.point(tv));
//		}
//		for (auto tf : pmesh.faces())
//		{
//			vector<OV> vh;
//			for (auto tfv : pmesh.fv_range(tf))
//			{
//				vh.push_back(pm.vertex_handle(tfv.idx()+nv));
//			}
//			pm.add_face(vh[0], vh[1], vh[2], vh[3]);
//		}
//		//break;
//	}
//}

//{
//	//首先确定每个面内增加的点数
//	double MaxGaussCurvature = 0;
//	int fsize = faceshape.size();
//	vector<int> AddNum(fsize);
//	vector<vector<MatrixXd>> SurfacePolygons(fsize);
//	vector<vector<double>> CurvatureDistribution(fsize);
//	vector<Vector4d> Bound(fsize);
//	for (int i = 0; i < fsize; i++) {
//		//计算polygons
//		auto &aface = faceshape[i].face;
//		auto &polygons = SurfacePolygons[i];
//		auto &wires = faceshape[i].wires;
//		polygons.resize(wires.size());
//		for (int j = 0; j < wires.size(); j++) {
//			auto &polygon = polygons[j];
//			for (auto aedgeidx : wires[j]) {
//				auto &parameters = edgeshape[aedgeidx].parameters;
//				int pr = parameters.rows() - 1;
//				polygon.conservativeResize(polygon.rows() + pr, 2);
//				polygon.block(polygon.rows() - pr, 0, pr, 2) = parameters.block(0, 0, pr, 2);
//			}
//		}
//		//计算预测gauss curvature分布的网络
//		TopLoc_Location loc;
//		Handle(Geom_Surface) asurface = BRep_Tool::Surface(aface, loc);
//		Standard_Real umi, uma, vmi, vma;
//		asurface->Bounds(umi, uma, vmi, vma);
//		//cout << umi << "," << uma << "," << vmi << "," << vma << endl;
//		int x_num = 0;
//		int y_num = 0;
//		if (Precision::IsNegativeInfinite(umi) || Precision::IsPositiveInfinite(uma)
//			|| Precision::IsNegativeInfinite(vmi) || Precision::IsPositiveInfinite(vma)) {
//			x_num = 50;
//			y_num = 50;
//		}
//		else {
//			gp_Pnt corner0 = asurface->Value(umi, vmi);
//			gp_Pnt corner1 = asurface->Value(umi, vma);
//			gp_Pnt corner2 = asurface->Value(uma, vmi);
//			gp_Pnt corner3 = asurface->Value(uma, vma);
//			x_num = min(50, static_cast<int>((pow(corner0.X() - corner1.X(), 2) + pow(corner0.Y() - corner1.Y(), 2) + pow(corner0.Z() - corner1.Z(),
//				2) + pow(corner2.X() - corner3.X(), 2) + pow(corner2.Y() - corner3.Y(), 2) + pow(corner2.Z() - corner3.Z(), 2)) / expected_edge_length));
//			y_num = min(50, static_cast<int>((pow(corner0.X() - corner2.X(), 2) + pow(corner0.Y() - corner2.Y(), 2) + pow(corner0.Z() - corner2.Z(),
//				2) + pow(corner1.X() - corner3.X(), 2) + pow(corner1.Y() - corner3.Y(), 2) + pow(corner1.Z() - corner3.Z(), 2)) / expected_edge_length));
//		}
//		cout << i << "," << x_num << "," << y_num << endl;
//		//计算gauss curvature
//		int outerboundflag = ComputePolygonsOuterBoundIndex(polygons);
//		if (outerboundflag > 0) swap(polygons[0], polygons[outerboundflag]);
//		Vector4d &bound = Bound[i];
//		bound = ComputePolygonBound(polygons[0]);
//		double x_step = (bound(1) - bound(0)) / x_num;
//		double y_step = (bound(3) - bound(2)) / y_num;
//		GeomLProp_SLProps curvature_tool = GeomLProp_SLProps(asurface, 2, 1.0e-15);
//		auto &curvature_distribution = CurvatureDistribution[i];
//		for (int j = 0; j <= x_num; j++) {
//			for (int k = 0; k <= y_num; k++) {
//				if (IsInPolygons(polygons, Vector2d(j*x_step + bound(0), k*y_step + bound(1)))) {
//					curvature_tool.SetParameters(j*x_step + bound(0), k*y_step + bound(1));
//					double fabscurvature = fabs(curvature_tool.GaussianCurvature());
//					MaxGaussCurvature = max(MaxGaussCurvature, fabscurvature);
//					curvature_distribution.push_back(fabscurvature);
//				}
//			}
//		}
//	}
//	//在边界基础上增加内点，构建网格
//	Surface_Meshes.resize(fsize);
//	uniform_real_distribution<double> g_distrdouble(0, MaxGaussCurvature);
//	for (int i = 0; i < fsize; i++) {
//		double num = 0;
//		for (auto gc : CurvatureDistribution[i]) {
//			num += pow(gc * 6 / MaxGaussCurvature, 2);
//		}
//		int targetnum = static_cast<int>(num);
//		cout << i << "," << targetnum << endl;
//		TopLoc_Location loc;
//		Handle(Geom_Surface) asurface = BRep_Tool::Surface(faceshape[i].face, loc);
//		GeomLProp_SLProps curvature_tool = GeomLProp_SLProps(asurface, 2, 1.0e-15);
//		Vector4d &bound = Bound[i];
//
//		random_device rd;
//		mt19937_64 eng(rd());
//		uniform_real_distribution<double> x_distrdouble(bound(0), bound(1));
//		uniform_real_distribution<double> y_distrdouble(bound(2), bound(3));
//		int n = 0;
//		auto &polygons = SurfacePolygons[i];
//		/*for (auto polygon : polygons) {
//			for (int j = 0; j < polygon.rows(); j++) {
//				cout << j << "," << polygon.row(j) << endl;
//			}
//		}*/
//
//		int x_num = 13; int y_num = 50;
//		double x_step = (bound(1) - bound(0)) / x_num;
//		double y_step = (bound(3) - bound(2)) / y_num;
//		auto &curvature_distribution = CurvatureDistribution[i];
//		for (int j = 0; j <= x_num; j++) {
//			for (int k = 0; k <= y_num; k++) {
//				if (IsInPolygons(polygons, Vector2d(j*x_step + bound(0), k*y_step + bound(1)))) {
//					curvature_tool.SetParameters(j*x_step + bound(0), k*y_step + bound(1));
//					double fabscurvature = fabs(curvature_tool.GaussianCurvature());
//					MaxGaussCurvature = max(MaxGaussCurvature, fabscurvature);
//					curvature_distribution.push_back(fabscurvature);
//				}
//			}
//		}
//
//		polygons.push_back(MatrixXd::Zero(targetnum, 2));
//		MatrixXd &newpoints = *polygons.rbegin();
//		while (n < targetnum) {
//			double x = x_distrdouble(eng);
//			double y = y_distrdouble(eng);
//			double g = g_distrdouble(eng);
//			curvature_tool.SetParameters(x, y);
//			double gausscurvature = fabs(curvature_tool.GaussianCurvature());
//			if (gausscurvature < g) continue;
//			//cout << n << "," << x << "," << y << endl;
//			if (IsInPolygons(polygons, Vector2d(x, y))) {
//				newpoints.row(n++) << x, y;
//			}
//		}
//		MatrixXd all_pnts;
//		vector<MatrixXi> bnd(polygons.size() - 1);
//		for (int j = 0; j < polygons.size() - 1; j++) {
//			int ps = polygons[j].rows();
//			all_pnts.conservativeResize(all_pnts.rows() + ps, 2);
//			all_pnts.block(all_pnts.rows() - ps, 0, ps, 2) = polygons[j];
//			bnd[j].resize(ps, 2);
//			for (int k = 0; k < ps - 1; k++) {
//				bnd[j].row(k) << k, k + 1;
//			}
//			bnd[j].row(ps - 1) << ps - 1, 0;
//		}
//		
//		Mesh &aMesh = Surface_Meshes[i];
//		triangulate(all_pnts, bnd, 100, aMesh);
//		for (auto tv : aMesh.vertices()) {
//			auto pos = aMesh.point(tv);
//			gp_Pnt p = asurface->Value(pos[0], pos[1]);
//			aMesh.set_point(tv, Mesh::Point(p.X(), p.Y(), p.Z()));
//		}
//	}
//	cout << "Piecewise Mesh Done!" << endl;
//}


void OccReader::ComputeFaceAndEdge()
{
	TopoDS_Shape &aShape = globalmodel.aShape;

	vector<ShapeFace> &faceshape = globalmodel.faceshape;
	vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;
	faceshape.clear();
	edgeshape.clear();

	int faceSize = 0;
	int edgeSize = 0;
	for (TopExp_Explorer faceExp(aShape, TopAbs_FACE); faceExp.More(); faceExp.Next())
	{
		faceSize++;
		for (TopExp_Explorer wireExp(TopoDS::Face(faceExp.Current()), TopAbs_WIRE); wireExp.More(); wireExp.Next())
			edgeSize += wireExp.Current().NbChildren();
	}
	faceshape.reserve(faceSize);
	edgeshape.reserve(edgeSize);


	Vector3d ma(DBL_MIN, DBL_MIN, DBL_MIN);
	Vector3d mi(DBL_MAX, DBL_MAX, DBL_MAX);
	for (TopExp_Explorer vertexExp(aShape, TopAbs_VERTEX); vertexExp.More(); vertexExp.Next())
	{
		gp_Pnt v = BRep_Tool::Pnt(TopoDS::Vertex(vertexExp.Current()));
		Vector3d p(v.X(), v.Y(), v.Z());
		for (int i = 0; i < 3; i++)
		{
			ma(i) = std::max(ma(i), p(i));
			mi(i) = std::min(mi(i), p(i));
		}
	}
	expected_edge_length = 0.003 * sqrt(pow(ma(0) - mi(0), 2) + pow(ma(1) - mi(1), 2) + pow(ma(2) - mi(2), 2));
	//a new method to order the faces and edges
	int itertimes_face = 0;
	int itertimes_edge = 0;
	int numofholes = 0;
	for (TopExp_Explorer faceExp(aShape, TopAbs_FACE); faceExp.More(); faceExp.Next())
	{
		auto &aface = TopoDS::Face(faceExp.Current());
		/*if (!aface.Location().IsIdentity())
		{
			cout << itertimes_face << "th face location not identity\n\n\n\n\n Location Warning!!!" << endl;
			aface.Location().ShallowDump(cout);
			system("pause");
		}*/
		//faceshape.push_back(ShapeFace(itertimes_face, aface));
		faceshape.emplace_back(itertimes_face, aface);
		//if (aface.NbChildren() != 1)
		//{
		//	//cout << "The " << itertimes_face << "th face has " << aface.NbChildren() - 1 << " holes" << endl;
		//}
		numofholes += aface.NbChildren() - 1;
		for (TopExp_Explorer wireExp(aface, TopAbs_WIRE); wireExp.More(); wireExp.Next())
		{
			TopoDS_Wire awire = TopoDS::Wire(wireExp.Current());
			vector<int> edges;
			for (TopExp_Explorer edgeExp(awire, TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
			{
				auto &aedge = TopoDS::Edge(edgeExp.Current());
				if (BRep_Tool::Degenerated(aedge))
				{
					//cout << "a degenerated edge in the " << itertimes_face << "th face" << endl;
					continue;
				}
				if (!BRep_Tool::IsClosed(aedge))
				{
					gp_Pnt p0 = BRep_Tool::Pnt(TopExp::FirstVertex(aedge));
					gp_Pnt p1 = BRep_Tool::Pnt(TopExp::LastVertex(aedge));
					if (p0.IsEqual(p1, expected_edge_length * 1.0e-2)) continue;
				}

				//edgeshape.push_back(ShapeEdge(itertimes_edge, aedge));
				edgeshape.emplace_back(itertimes_edge, aedge);
				edgeshape[itertimes_edge].main_face = itertimes_face;
				edges.push_back(itertimes_edge);
				itertimes_edge++;
			}
			if (awire.Orientation() == TopAbs_REVERSED) edges = vector<int>(edges.rbegin(), edges.rend());
			if (edges.size()) faceshape[itertimes_face].wires.push_back(edges);
			//else cout << "the edges on a wire of the " << itertimes_face << "th face are all degenerated" << endl;
		}
		itertimes_face++;
	}
	for (int i = 0; i < edgeshape.size(); i++)
	{
		auto &aedge = edgeshape[i];
		if (aedge.secondary_face != -1) continue;
		for (int j = i + 1; j < edgeshape.size(); j++)
		{
			auto &redge = edgeshape[j];
			if (aedge.edge.IsSame(redge.edge))
			{
				aedge.reversed_edge = j;
				redge.reversed_edge = i;
				aedge.secondary_face = redge.main_face;
				redge.secondary_face = aedge.main_face;
				break;
			}
		}
	}
	dprint("CAD Info:\nModel Size: [", mi(0), ma(0), "],[", mi(1), ma(1), "],[", mi(2), ma(2),"]"
		"\nSurface Number: ", faceshape.size(), 
		"\nEdge Number: ", edgeshape.size(), 
		"\nHoles Number: ", numofholes, 
		"\n\nCompute Topology Done!");
}

void OccReader::Discrete_Bound()
{
	vector<ShapeFace> &faceshape = globalmodel.faceshape;
	vector<ShapeEdge> &edgeshape = globalmodel.edgeshape;

	const double GL_c[8] = { 0.1012285363,0.2223810345,0.3137066459,0.3626837834,0.3626837834 ,0.3137066459,0.2223810345,0.1012285363 };
	const double GL_x[8] = { -0.9602898565,-0.7966664774,-0.5255354099,-0.1834346425,0.1834346425,0.5255354099,0.7966664774,0.9602898565 };
	for (int i = 0; i < edgeshape.size(); i++)
	{
		//cout << i << endl;
		auto &edge = edgeshape[i];
		auto &aedge = edge.edge;
		auto &mainface = faceshape[edge.main_face].face;
		TopLoc_Location loc;
		Handle(Geom_Surface) asurface = BRep_Tool::Surface(mainface, loc);
		//loc.ShallowDump(cout);
		Standard_Real first = 0;
		Standard_Real last = 0;
		Handle(Geom_Curve) acurve = BRep_Tool::Curve(aedge, first, last);
		
		double curve_length = 0;
		gp_Pnt P;
		gp_Vec d1;
		//Gauss-Lengred积分公式
		double a = (last - first) / 2;
		double b = (last + first) / 2;
		for (int j = 0; j < 8; j++)
		{
			acurve->D1(a*GL_x[j] + b, P, d1);
			curve_length += a * GL_c[j] * d1.Magnitude();
		}
		//int seg = BRep_Tool::IsClosed(aedge) ? 4 : 3;
		int segment_number = std::max(4, int(curve_length / expected_edge_length));



		//*************************** for test
	//	if (segment_number % 2) segment_number++;
		//*************************** need to delete later




		edge.parameters.resize(segment_number + 1, 2);

		Handle_Geom2d_Curve thePCurve = BRep_Tool::CurveOnSurface(aedge, mainface, first, last);
		double step = (last - first) / segment_number;
		gp_Pnt2d uv;
		if(aedge.Orientation() == TopAbs_FORWARD)
		{
			for (int j = 0; j < segment_number; j++)
			{
				uv = thePCurve->Value(first + j * step);
				edge.parameters.row(j) << uv.X(), uv.Y();
			}
			uv = thePCurve->Value(last);
			edge.parameters.row(segment_number) << uv.X(), uv.Y();
		}
		else
		{
			for (int j = 0; j < segment_number; j++)
			{
				uv = thePCurve->Value(last - j * step);
				edge.parameters.row(j) << uv.X(), uv.Y();
			}
			uv = thePCurve->Value(first);
			edge.parameters.row(segment_number) << uv.X(), uv.Y();
		}
	}
	dprint("Discrete Boundary Done!");
}




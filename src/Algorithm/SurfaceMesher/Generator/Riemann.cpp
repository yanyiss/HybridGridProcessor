#include"Riemann.h"

Eigen::Matrix2d Riemannremesh::Riemanndata(const TriMesh::Point &p)
{
	Point ru = B->PartialDerivativeU(p[0], p[1]);
	Point rv = B->PartialDerivativeV(p[0], p[1]);
	double E = ru.dot(ru);
	double F = ru.dot(rv);
	double G = rv.dot(rv);
	Eigen::Matrix2d first_form;
	first_form << E, F,
		F, G;
	return first_form;
}

double Riemannremesh::Riemannlen(const TriMesh::VertexHandle &h1, const TriMesh::VertexHandle &h2)
{
	auto p1 = mesh->point(h1);
	auto p2 = mesh->point(h2);
	Eigen::MatrixXd xy = Eigen::MatrixXd::Zero(1, 2);
	xy << p1[0] - p2[0], p1[1] - p2[1];
	return pow((xy * ((mesh->data(h1).M + mesh->data(h2).M) / 2) * xy.transpose()).determinant(), 1.0 / 2);
}

void Riemannremesh::calulenth()
{
	mesh->request_edge_status();
	mesh->request_face_status();
	mesh->request_vertex_status();
	mesh->request_halfedge_status();
	TriMesh::Point p;
	for (auto v = mesh->vertices_begin(); v != mesh->vertices_end(); v++)
	{
		p = mesh->point(*v);
		if (p[0] < uu.front()) p[0] = uu.front();
		else if (p[0] > uu.back()) p[0] = uu.back();

		if (p[1] < vv.front()) p[1] = vv.front();
		else if (p[1] > vv.back()) p[1] = vv.back();

		mesh->set_point(*v, p);
		mesh->data(*v).M = Riemanndata(p);
	}
	double len = 0;
	TriMesh::HalfedgeHandle he;
	int count = 0;
	for (TriMesh::EdgeIter ee = mesh->edges_begin(); ee != mesh->edges_end(); ee++)
	{
		if (mesh->is_boundary(*ee))
		{
			count++;
			he = mesh->halfedge_handle(*ee, 0);
			len += Riemannlen(mesh->from_vertex_handle(he), mesh->to_vertex_handle(he));
		}
	}
	len /= count;
	highlenth = 4 * len / 3;
	lowlenth = 4 * len / 5;
}

void Riemannremesh::split()
{
	TriMesh::VertexHandle vh, vh1, vh2, vh3, vh4;
	TriMesh::EdgeHandle ee;
	TriMesh::HalfedgeHandle he;
	int edgeNum = mesh->n_edges();
	double len;
	std::vector<TriMesh::VertexHandle> facevhandle;
	for (int i = 0; i < edgeNum; i++)
	{
		ee = mesh->edge_handle(i);
		if (mesh->is_boundary(ee)) continue;
		he = mesh->halfedge_handle(ee, 0);
		len = Riemannlen(mesh->from_vertex_handle(he), mesh->to_vertex_handle(he));
		if (len > highlenth)
		{
			vh = mesh->add_vertex(mesh->calc_edge_midpoint(ee));
			mesh->split_edge(ee, vh);   //split²Ù×÷
			mesh->data(vh).M = Riemanndata(mesh->point(vh));
		}
	}
	mesh->garbage_collection();
}

void Riemannremesh::collapse()
{
	TriMesh::EdgeHandle e;
	TriMesh::HalfedgeHandle he;
	TriMesh::VertexHandle p1, p2;
	float len;
	bool is_collapse;
	for (int i = mesh->n_edges() - 1; i >= 0; i--)
	{
		if (i > mesh->n_edges() - 1) continue;
		e = mesh->edge_handle(i);
		he = mesh->halfedge_handle(e, 0);
		if (!mesh->is_collapse_ok(he)) continue;
		p1 = mesh->from_vertex_handle(he);
		p2 = mesh->to_vertex_handle(he);
		if (mesh->is_boundary(p1) || mesh->is_boundary(p2)) continue;
		len = Riemannlen(p1, p2);
		if (len < lowlenth)
		{
			is_collapse = 1;
			for (auto vv = mesh->vv_begin(p1); vv.is_valid(); vv++)
			{
				if (Riemannlen(*vv, p2) > highlenth)
				{
					is_collapse = 0;
					break;
				}
			}
			if (is_collapse)
			{
				mesh->collapse(he);
			}
		}
	}
	mesh->garbage_collection();
}

void Riemannremesh::flip()
{
	std::vector<int> targetv;
	TriMesh::HalfedgeHandle he, op;
	TriMesh::VertexHandle a, b, c, d;
	int pre, post;
	for (auto v = mesh->vertices_begin(); v != mesh->vertices_end(); v++)
	{
		if (mesh->is_boundary(*v))
		{
			targetv.push_back(4);
		}
		else
		{
			targetv.push_back(6);
		}
	}
	for (auto e = mesh->edges_begin(); e != mesh->edges_end(); e++)
	{
		if (mesh->is_boundary(*e) || !mesh->is_flip_ok(*e))
		{
			continue;
		}
		he = mesh->halfedge_handle(*e, 0);
		a = mesh->from_vertex_handle(he);
		b = mesh->to_vertex_handle(he);
		c = mesh->to_vertex_handle(mesh->next_halfedge_handle(he));
		op = mesh->opposite_halfedge_handle(he);
		d = mesh->to_vertex_handle(mesh->next_halfedge_handle(op));
		pre = abs(int(mesh->valence(a) - targetv[a.idx()])) +
			abs(int(mesh->valence(b) - targetv[b.idx()])) +
			abs(int(mesh->valence(c) - targetv[c.idx()])) +
			abs(int(mesh->valence(d) - targetv[d.idx()]));
		post = abs(int(mesh->valence(a) - 1 - targetv[a.idx()])) +
			abs(int(mesh->valence(b) - 1 - targetv[b.idx()])) +
			abs(int(mesh->valence(c) + 1 - targetv[c.idx()])) +
			abs(int(mesh->valence(d) + 1 - targetv[d.idx()]));
		if (pre > post)
		{
			mesh->flip(*e);
		}
	}
}

void Riemannremesh::updatepoint()
{
	float count;
	TriMesh::Point newpoint;
	for (auto v = mesh->vertices_begin(); v != mesh->vertices_end(); v++)
	{
		if (mesh->is_boundary(*v))
		{
			continue;
		}
		count = 0.0;
		newpoint = Mesh::Point(0, 0, 0);
		for (auto vv = mesh->vv_begin(*v); vv.is_valid(); vv++)
		{
			newpoint += mesh->point(*vv);
			count += 1;
		}
		mesh->set_point(*v, newpoint / count);
	}
}

void Riemannremesh::remesh()
{
	calulenth();
	for (int i = 0; i < 10; i++)
	{
		split();
		collapse();
		flip();
		updatepoint();
		if (i < 9)
		{
			for (auto v = mesh->vertices_begin(); v != mesh->vertices_end(); v++)
			{
				auto p = mesh->point(*v);
				mesh->set_point(*v, p);
				mesh->data(*v).M = Riemanndata(p);
			}
		}
	}
}

void Riemannremesh::curvature_feature(Eigen::Matrix2Xd &all_pnts, std::vector<bool> &curvature)
{
	Point ru, rv, n;
	Eigen::Matrix2d Weingarten, value;
	//std::cout << uu.front() << '\t' << uu.back() << '\t' << vv.front() << '\t' << vv.back() << '\n';
	for (int i = 0; i < all_pnts.cols(); i++)
	{
		double u =  all_pnts(0,i), v = all_pnts(1,i);
		//std::cout << i<<": "<<u << '\t' << v << '\t';
		if (u < uu.front()) 
		{
			u = uu.front();
			//std::cout << "1" << std::endl;
		}
		else if (u > uu.back()) 
		{
			u = uu.back();
			//std::cout << "2" << std::endl;
		}
		if (v < vv.front())
		{
			v = vv.front(); 
			//std::cout << "3" << std::endl;
		}
		else if (v > vv.back())
		{
			v = vv.back();
			//std::cout << "4" << std::endl;
		}
		//std::cout << u << '\t' << v << '\n';
		ru = B->PartialDerivativeU(u, v);
		rv = B->PartialDerivativeV(u, v);
		//std::cout << i << '\t' << rv[0] << '\t' << rv[1] << '\t' << rv[2] << '\n';
		double E = ru.dot(ru);
		double F = ru.dot(rv);
		double G = rv.dot(rv);
		n = (ru.cross(rv)).normalized();
		double L = (B->PartialDerivativeUU(u, v)).dot(n);
		double M = (B->PartialDerivativeUV(u, v)).dot(n);
		double N = (B->PartialDerivativeVV(u, v)).dot(n);
		//std::cout <<i<< ":"<<E << '\t' << F << '\t' << G << '\t' << L << '\t' << M << '\t' << N << '\t'<<'\n';
		Weingarten << L * G - M * F, M * E - L * F,
			M * G - N * F, N * E - M * F;
		//std::cout << Weingarten / (E*G - pow(F, 2)) << std::endl;
		Eigen::EigenSolver<Eigen::Matrix2d> es(Weingarten);
		value = (es.pseudoEigenvalueMatrix())/(E*G-pow(F,2));
		//std::cout << value << std::endl;
		double K1 = std::max(std::abs(value(0, 0)), std::abs(value(1, 1)));
		double K2 = std::min(std::abs(value(0, 0)), std::abs(value(1, 1)));
		//std::cout << i << '\t' << K1 << '\t' << K2 << std::endl<<std::endl;
		//if (K2 < 0.0008) continue;
		if (K2/K1 < 0.015)
		{
			curvature[i] = true;
		}
	}
}
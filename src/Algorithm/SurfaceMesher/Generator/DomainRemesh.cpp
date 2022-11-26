#include"DomainRemesh.h"
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
	return pow((xy * ((mesh->data(h1).M + mesh->data(h2).M) *0.5) * xy.transpose()).determinant(), 0.5);
}

void Riemannremesh::calclength()
{
	mesh->request_edge_status();
	mesh->request_face_status();
	mesh->request_vertex_status();
	mesh->request_halfedge_status();
	for (auto v : mesh->vertices())
	{
		mesh->data(v).M = Riemanndata(mesh->point(v));
	}
	double len = 0;
	TriMesh::HalfedgeHandle he;
	int count = 0;
	for (auto e : mesh->edges())
	{
		//if (!e.is_boundary()) continue;
		count++;
		he = mesh->halfedge_handle(e, 0);
		len += Riemannlen(mesh->from_vertex_handle(he), mesh->to_vertex_handle(he));
	}
	len /= count;
	highlength = 4 * len / 3;
	lowlength = 0.8 * len;
}

bool Riemannremesh::if_reverse_face(TriMesh::FaceHandle& f)
{
	auto itr = mesh->fh_begin(f);
	auto v1 = mesh->calc_edge_vector(*itr);
	itr++;
	auto v2 = mesh->calc_edge_vector(*itr);
	if ((v1[0] * v2[1] - v1[1] * v2[0] > 0) ^ direction) return true;
	else return false;
}

void Riemannremesh::split()
{
	TriMesh::VertexHandle vh, vh1, vh2;
	TriMesh::EdgeHandle ee;
	TriMesh::HalfedgeHandle he;
	int edgeNum = mesh->n_edges();
	double len;
	for (int i = 0; i < edgeNum; i++)
	{
		ee = mesh->edge_handle(i);
		if (mesh->is_boundary(ee)) continue;
		he = mesh->halfedge_handle(ee, 0);
		vh1 = mesh->from_vertex_handle(he);
		vh2 = mesh->to_vertex_handle(he);
		len = Riemannlen(vh1, vh2);
		if (len > highlength)
		{
			vh = mesh->add_vertex(mesh->calc_edge_midpoint(ee));
			mesh->split_edge(ee, vh);
			mesh->data(vh).M = (mesh->data(vh1).M + mesh->data(vh2).M) * 0.5;
		}
	}
	mesh->garbage_collection();
}

void Riemannremesh::collapse()
{
	TriMesh::EdgeHandle e;
	TriMesh::HalfedgeHandle he;
	TriMesh::VertexHandle p1, p2;
	for (int i = mesh->n_edges() - 1; i >= 0; i--)
	{
		if (i > mesh->n_edges() - 1) continue;
		e = mesh->edge_handle(i);
		he = mesh->halfedge_handle(e, 0);
		if (!mesh->is_collapse_ok(he)) continue;
		p1 = mesh->from_vertex_handle(he);
		p2 = mesh->to_vertex_handle(he);
		if (mesh->is_boundary(p1))
		{
			if (mesh->is_boundary(p2)) continue;
			he = mesh->opposite_halfedge_handle(he);
			p1 = mesh->from_vertex_handle(he);
			p2 = mesh->to_vertex_handle(he);
		}
		if (Riemannlen(p1, p2) >= lowlength) continue;
		bool is_collapse = true;
		for (auto vf : mesh->vf_range(p1))
		{
			int count = 0;
			for (auto fv : mesh->fv_range(vf))
			{
				if (fv.is_boundary() && fv.idx() != p2.idx()) count++;
			}
			if (count > 1)
			{
				is_collapse = false;
				break;
			}
		}
		if (!is_collapse) continue;
		for (auto vv = mesh->vv_begin(p1); vv.is_valid(); vv++)
		{
			if (Riemannlen(*vv, p2) > highlength)
			{
				is_collapse = false;
				break;
			}
		}
		if (!is_collapse) continue;

		auto pos = mesh->is_boundary(p2) ? mesh->point(p2) : mesh->calc_centroid(he);
		TriMesh::HalfedgeHandle flipH = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));
		auto& fromPos = mesh->point(mesh->from_vertex_handle(he));
		int endId = (mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(he))).idx();
		while (flipH.idx() != endId)
		{
			auto& p0 = mesh->point(mesh->to_vertex_handle(flipH));
			auto& p1 = mesh->point(mesh->to_vertex_handle(mesh->next_halfedge_handle(flipH)));
			if ((p0 - fromPos).cross(p1 - fromPos).dot((p0 - pos).cross(p1 - pos)) < 0)
			{
				is_collapse = false;
				break;
			}
			flipH = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(flipH));
		}
		if (!mesh->is_boundary(p2) && is_collapse)
		{
			auto& toPos = mesh->point(mesh->to_vertex_handle(he));
			flipH = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(he));
			endId = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(he)).idx();
			while (flipH.idx() != endId)
			{
				auto& p0 = mesh->point(mesh->from_vertex_handle(flipH));
				auto& p1 = mesh->point(mesh->from_vertex_handle(mesh->prev_halfedge_handle(flipH)));
				if ((p0 - toPos).cross(p1 - toPos).dot((p0 - pos).cross(p1 - pos)) < 0)
				{
					is_collapse = false;
					break;
				}
				flipH = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(flipH));
			}
		}

		if (is_collapse)
		{
			mesh->set_point(p2, pos);
			mesh->collapse(he);
		}
	}
	mesh->garbage_collection();
}
	
void Riemannremesh::equalize_valence()
{
	std::vector<int> targetv;
	TriMesh::VertexHandle a, b, c, d;
	int pre, post;
	for (auto v = mesh->vertices_begin(); v != mesh->vertices_end(); v++)
	{
		if (mesh->is_boundary(*v)) targetv.push_back(4);
		else targetv.push_back(6);
	}
	for (auto e : mesh->edges())
	{
		if (mesh->is_boundary(e) || !mesh->is_flip_ok(e)) continue;
		auto he = e.h0();
		a = mesh->from_vertex_handle(he);
		b = mesh->to_vertex_handle(he);
		c = mesh->to_vertex_handle(mesh->next_halfedge_handle(he));
		auto op = he.opp();
		d = mesh->to_vertex_handle(mesh->next_halfedge_handle(op));
		if (mesh->is_boundary(c) && mesh->is_boundary(d)) continue;
		pre = abs(int(mesh->valence(a) - targetv[a.idx()])) +
			abs(int(mesh->valence(b) - targetv[b.idx()])) +
			abs(int(mesh->valence(c) - targetv[c.idx()])) +
			abs(int(mesh->valence(d) - targetv[d.idx()]));
		post = abs(int(mesh->valence(a) - 1 - targetv[a.idx()])) +
			abs(int(mesh->valence(b) - 1 - targetv[b.idx()])) +
			abs(int(mesh->valence(c) + 1 - targetv[c.idx()])) +
			abs(int(mesh->valence(d) + 1 - targetv[d.idx()]));
		if (pre <= post) continue;
		auto edire = mesh->calc_edge_vector(he);
		double alpha = acos(-edire.dot(mesh->calc_edge_vector(he.prev()))) + acos(edire.dot(mesh->calc_edge_vector(op.next())));
		if (alpha > PI) continue;
		alpha = acos(-edire.dot(mesh->calc_edge_vector(he.next()))) + acos(edire.dot(mesh->calc_edge_vector(op.prev())));
		if (alpha > PI) continue;
		mesh->flip(e);
	}
	mesh->garbage_collection();
}

void Riemannremesh::updatepoint()
{
	int count;
	TriMesh::Point newpoint, oldpoint;
	for (auto v : mesh->vertices())
	{
		if (mesh->is_boundary(v)) continue;
		oldpoint = mesh->point(v);
		count = 0;
		newpoint = Mesh::Point(0, 0, 0);
		for (auto vv : mesh->vv_range(v))
		{
			newpoint += mesh->point(vv);
			count += 1;
		}
		mesh->set_point(v, newpoint / count);
		bool is_update = true;
		for (auto f : mesh->vf_range(v))
		{
			if (!if_reverse_face(f)) continue;
			is_update = false;
			break;
		}
		if (is_update)
		{
			auto& M = mesh->data(v).M;
			M.setZero();
			for (auto vv : mesh->vv_range(v)) M += mesh->data(vv).M;
			M /= count;
		}
		else mesh->set_point(v, oldpoint);
	}
}

void Riemannremesh::remesh()
{
	calclength();
	auto f = mesh->face_handle(0);
	auto itr = mesh->fh_begin(f);
	auto v1 = mesh->calc_edge_vector(*itr);
	itr++;
	auto v2 = mesh->calc_edge_vector(*itr);
	if (v1[0] * v2[1] - v1[1] * v2[0] > 0) direction = true;
	else direction = false;

	for (int i = 0; i < 5; i++)
	{
		split();
		collapse();
		equalize_valence();
		updatepoint();
	}
}


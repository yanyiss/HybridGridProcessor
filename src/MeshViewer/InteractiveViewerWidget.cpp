#include <QMouseEvent>
#include <QLineEdit>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtCore>
#include <QUrl>

#include "InteractiveViewerWidget.h"
#include "..\src\Toolbox\dprinter\dprint.h"

InteractiveViewerWidget::InteractiveViewerWidget(QWidget* parent /* = 0 */)
	:MeshViewerWidget(parent)
{
	draw_new_mesh = false;
	clearSelectedData();
	kdTree = NULL;
}

InteractiveViewerWidget::InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent)
:MeshViewerWidget(_fmt, _parent)
{
	draw_new_mesh = false;
	clearSelectedData();
	kdTree = NULL;
}

InteractiveViewerWidget::~InteractiveViewerWidget()
{
	if(kdTree) delete kdTree;
	if (occreader) { delete occreader; occreader = nullptr; }
	if (iso_mesh) { delete iso_mesh; iso_mesh = nullptr; }
	if (stripTree) { delete stripTree; stripTree = nullptr; }
}

void InteractiveViewerWidget::setMouseMode(int mm)
{
	if(mouse_mode_ != T2_MODE)
	{
		mouse_mode_ = mm;
		if( TRANS != mouse_mode_ )
		{ buildIndex(); }
		emit setMouseMode_signal(mm);
	}
}

void InteractiveViewerWidget::mousePressEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mousePressEvent(_event);
	}
	else
	{
		if(mouse_mode_ != T2_MODE)
		{
			pick_point( _event->x(), _event->y() );
			if(mouse_mode_ == VERTEXPICK)
			{
				pick_vertex( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == FACEPICK)
			{
				pick_face( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == EDGEPICK)
			{
				pick_edge( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == POINTPICK)
			{
			}
			else if (mouse_mode_ == CURVEPICK)
			{
				pick_curve(_event->x(), _event->y());
				dprint("iowjm");
			}
			else if( mouse_mode_ == MOVE )
			{
				pick_vertex( _event->x(), _event->y() );//set the selected handle
			}
			else if( mouse_mode_ == VECTOR_SET )
			{
				pick_point( _event->x(), _event->y() );//set the selected handle
			}
			else if(mouse_mode_ == EDGECOLLAPSE)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::HalfedgeHandle heh = mesh.halfedge_handle( mesh.edge_handle(desired_edge), 0 );
					OpenMesh::Vec3d from_p = mesh.point(mesh.from_vertex_handle(heh));
					OpenMesh::Vec3d to_p = mesh.point(mesh.to_vertex_handle(heh));
					OpenMesh::Vec3d sp(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
					bool collapse_ok = true;
					if( (sp-from_p).sqrnorm() > (to_p-sp).sqrnorm() )
					{
						if( mesh.is_collapse_ok(heh) )
						{
							mesh.collapse(heh);
						}
						else
						{
							collapse_ok = false;
							printf("[%d] Collapse Not OK!\n", desired_edge);
						}
					}
					else
					{
						heh = mesh.opposite_halfedge_handle(heh);
						if( mesh.is_collapse_ok(heh) )
						{
							mesh.collapse(heh);
						}
						else
						{
							collapse_ok = false;
							printf("[%d] Collapse Not OK!\n", desired_edge);
						}
					}
					if(collapse_ok)
					{
						mesh.garbage_collection();
						buildIndex();
						if( mesh_vector.size() - 1 > mesh_vector_index )
						{
							mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
						}
						mesh_vector.push_back( mesh ); mesh_vector_index += 1;
						emit set_edit_undo_enable_viewer_signal( true );
						emit set_edit_redo_enable_viewer_signal( false );
					}
					clearSelectedData();
				}
			}
			else if (mouse_mode_ == EDGEFLIP)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::EdgeHandle eh = mesh.edge_handle(desired_edge);
					if( is_flip_ok_openmesh(eh, mesh))
					{
						flip_openmesh(eh, mesh);
						if( mesh_vector.size() - 1 > mesh_vector_index )
						{
							mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
						}
						mesh_vector.push_back( mesh ); mesh_vector_index += 1;
						emit set_edit_undo_enable_viewer_signal( true );
						emit set_edit_redo_enable_viewer_signal( false );
					}
					else
					{
						printf("[%d] Flip Not OK!\n", desired_edge);
					}
					clearSelectedData();
				}
			}
			else if (mouse_mode_ == EDGESPLIT)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::EdgeHandle eh = mesh.edge_handle(desired_edge);
					Mesh::HalfedgeHandle heh = mesh.halfedge_handle( eh, 0 );
					Mesh::HalfedgeHandle heh_ = mesh.halfedge_handle( eh, 1 );
					Mesh::VertexHandle vh0 = mesh.to_vertex_handle(heh);
					Mesh::VertexHandle vh1 = mesh.to_vertex_handle(heh_);
					OpenMesh::Vec3d s = mesh.point( vh1 );
					OpenMesh::Vec3d e = mesh.point( vh0 );
					Mesh::VertexHandle vh = mesh.add_vertex( (s + e)*0.5 );
					std::vector<Mesh::VertexHandle> one_face(3);
					if(mesh.is_boundary(eh))
					{
						if(Mesh::InvalidFaceHandle != mesh.face_handle(heh))
						{
							Mesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
							mesh.delete_edge(eh, false); mesh.garbage_collection();
							one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh.add_face(one_face);
							one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh.add_face(one_face);
						}
						else
						{
							Mesh::VertexHandle vh3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh_));
							mesh.delete_edge(eh, false); mesh.garbage_collection();
							one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh.add_face(one_face);
							one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh.add_face(one_face);
						}
					}
					else
					{
						Mesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
						Mesh::VertexHandle vh3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh_));
						mesh.delete_edge(eh, false); mesh.garbage_collection();
						one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh.add_face(one_face);
					}

					mesh.update_normals();
					buildIndex();
					clearSelectedData();

					if( mesh_vector.size() - 1 > mesh_vector_index )
					{
						mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
					}
					mesh_vector.push_back( mesh ); mesh_vector_index += 1;
					emit set_edit_undo_enable_viewer_signal( true );
					emit set_edit_redo_enable_viewer_signal( false );
				}
			}
		}
	}
	updateGL();
}

void InteractiveViewerWidget::mouseMoveEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else
	{
		if( mouse_mode_ != T2_MODE)
		{
			if( mouse_mode_ == MOVE )
			{
				move_point_based_lastVertex( _event->x(), _event->y() );
				Mesh::Point P(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
				mesh.set_point( mesh.vertex_handle(lastestVertex), P );
				updateGL();
			}
			if (mouse_mode_ == VECTOR_SET)
			{
				move_point_based_lastVertex(_event->x(), _event->y());
				updateGL();
			}
		}
		else
		{
			
		}
		
	}
}

void InteractiveViewerWidget::mouseReleaseEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else
	{
		if(mouse_mode_ != T2_MODE )
		{
			if( mouse_mode_ == MOVE )
			{
				move_point_based_lastVertex( _event->x(), _event->y() );
				Mesh::Point P(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
				mesh.set_point( mesh.vertex_handle(lastestVertex), P );
				selectedVertex.clear();
				updateGL();
			}
		}
		else
		{
		}
	}
	
}

void InteractiveViewerWidget::wheelEvent(QWheelEvent* _event)
{
	if(mouse_mode_ != N_MODE && mouse_mode_ != T2_MODE)
	{
		MeshViewerWidget::wheelEvent(_event);
	}
}

void InteractiveViewerWidget::dragEnterEvent(QDragEnterEvent* event)
{
	if( event->mimeData()->hasFormat("text/uri-list") )
	{
		event->acceptProposedAction();
	}
}

void InteractiveViewerWidget::dropEvent(QDropEvent* event)
{
	QList<QUrl> urls = event->mimeData()->urls();
	if( urls.isEmpty() )
		return;
	QString fileName = urls.first().toLocalFile();
	if (fileName.isEmpty())
		return;

	if( fileName.endsWith(".off") || fileName.endsWith(".obj") || fileName.endsWith(".stl") || fileName.endsWith(".ply"))
	{
		if( openMesh(fileName.toLocal8Bit()))
		{
			emit(loadMeshOK(true,fileName));
			setDrawMode(FLAT_POINTS);
			setMouseMode(TRANS);
		}
		else
		{
			emit(loadMeshOK(false,"No Mesh"));
		}
	}
}

void InteractiveViewerWidget::pick_vertex(int x,int y)
{
	int r = find_vertex_using_selected_point();
	lastestVertex = r;
	//printf("Select Vertex : %d\n", r);
	dprint("Select Vertex:", r, "\tCurvature:", mesh.data(mesh.vertex_handle(r)).GaussCurvature,
		"\tLength:",mesh.data(mesh.vertex_handle(r)).get_targetlength(), "\tPosition:", mesh.point(mesh.vertex_handle(r)));


	static std::vector<double> K1, K2;
	static std::vector<OpenMesh::Vec3d> D1, D2;
	if (fei)
	{
		fei = false;
		TriMesh tm(mesh);
		K1.clear(); K2.clear();
		D1.clear(); D2.clear();
		compute_principal_curvature(&tm, K1, K2, D1, D2);
	}
	dprint("curvature:", K1[r], K2[r]);

	std::vector<int>::iterator it;
	if( (it = std::find(selectedVertex.begin(),selectedVertex.end(), r)) == selectedVertex.end() )
	{
		selectedVertex.push_back(r);
		metric_constraints.insert(std::make_pair<OpenMesh::VertexHandle, CADMesher::metric_info>(mesh.vertex_handle(r), CADMesher::metric_info()));
		double aa = 0;
		for (auto ve : mesh.ve_range(mesh.vertex_handle(r)))
		{
			aa += mesh.calc_edge_length(ve);
		}
		aa /= mesh.valence(mesh.vertex_handle(r));
		metric_constraints[mesh.vertex_handle(r)].avg = aa;
		double min_cur_ = 4.0 / (occreader->bbmax - occreader->bbmin).norm();
		double k1 = K1[r]; k1 = std::fabs(k1) < min_cur_ ? min_cur_ : k1;
		double k2 = K2[r]; k2 = std::fabs(k2) < min_cur_ ? min_cur_ : k2;
		metric_constraints[mesh.vertex_handle(r)].geo_cur[0] = k1;
		metric_constraints[mesh.vertex_handle(r)].geo_cur[1] = k2;
		metric_constraints[mesh.vertex_handle(r)].geo_dir[0] = D1[r];
		metric_constraints[mesh.vertex_handle(r)].geo_dir[1] = D2[r];

		OpenMesh::Vec3d n = D1[r].cross(D2[r]);
		Eigen::Matrix3d U;
		U(0, 0) = D1[r][0]; U(1, 0) = D1[r][1]; U(2, 0) = D1[r][2];
		U(0, 1) = D2[r][0]; U(1, 1) = D2[r][1]; U(2, 1) = D2[r][2];
		U(0, 2) = n[0]; U(1, 2) = n[1]; U(2, 2) = n[2];
		Eigen::Matrix3d D; D.setZero();
		D(0, 0) = std::fabs(k1); D(1, 1) = std::fabs(k2);
		auto H = U * D*U.transpose();
		Eigen::Vector3d d0(D1[r][0], D1[r][1], D1[r][2]);
		Eigen::Vector3d d1(D2[r][0], D2[r][1], D2[r][2]);
		metric_constraints[mesh.vertex_handle(r)].geo_dir[0] *= aa / sqrt(d0.transpose()*H*d0);
		metric_constraints[mesh.vertex_handle(r)].geo_dir[1] *= aa / sqrt(d1.transpose()*H*d1);

		metric_constraints[mesh.vertex_handle(r)].dir[0] = metric_constraints[mesh.vertex_handle(r)].geo_dir[0];
		metric_constraints[mesh.vertex_handle(r)].dir[1] = metric_constraints[mesh.vertex_handle(r)].geo_dir[1];

		MeshParam->set_red_length(metric_constraints[mesh.vertex_handle(r)].dir[0].norm());
		MeshParam->set_green_length(metric_constraints[mesh.vertex_handle(r)].dir[1].norm());
	}
	else
	{
		selectedVertex.erase(it);
		metric_constraints.erase(mesh.vertex_handle(r));
	}

	updateGL();
}
void InteractiveViewerWidget::pick_face(int x,int y)
{
	int desiredFace = find_face_using_selected_point();
	if(desiredFace < 0) return;
	lastestFace = desiredFace;
	printf("Select Face : %d\n", desiredFace);
	std::vector<int>::iterator it;
	if( (it = std::find(selectedFace.begin(),selectedFace.end(),desiredFace)) == selectedFace.end() )
	{
		selectedFace.push_back(desiredFace);
	}
	else
	{
		selectedFace.erase(it);
	}
	updateGL();
}
void InteractiveViewerWidget::pick_edge(int x,int y)
{
	int desiredEdge = find_edge_using_selected_point();
	if(desiredEdge < 0) return;
	lastestEdge = desiredEdge;
	printf("\nSelect Edge : %d\n", desiredEdge);
	dprint("Edge Length:", mesh.calc_edge_length(mesh.edge_handle(desiredEdge)));
	std::vector<int>::iterator it;
	if( (it = std::find(selectedEdge.begin(),selectedEdge.end(),desiredEdge)) == selectedEdge.end() )
	{
		selectedEdge.push_back(desiredEdge);
	}
	else
	{
		selectedEdge.erase(it);
	}
	updateGL();
}
void InteractiveViewerWidget::pick_curve(int x, int y)
{
	if (if_new_mesh)
		BuildCurveIndex();
	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	stripTree->annkSearch(tp, 1, nnIdx, dists, 4);

	int curveIndex = edgeshapeIndex[nnIdx[0]];
	std::vector<int>::iterator it;
	if ((it = std::find(selectedCurve.begin(), selectedCurve.end(), curveIndex)) == selectedCurve.end())
	{
		//selectedVertex.push_back(lastestVertex);
		selectedCurve.push_back(curveIndex);
	}
	else
	{
		//selectedVertex.erase(it);
		selectedCurve.erase(it);
	}

	updateGL();
}
void InteractiveViewerWidget::pick_point(int x,int y)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = double(x);
	GLdouble winY = double( height() - y );
	GLfloat winZ = 0.0;
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	gluUnProject(winX, winY, (GLdouble)winZ, &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

void InteractiveViewerWidget::move_point_based_lastVertex(int x,int y)
{
	if(lastestVertex<0 || lastestVertex>=mesh.n_vertices())
	{
		return;
	}
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = 0.0;
	GLdouble winY = 0.0;
	GLdouble winZ = 0.0;
	OpenMesh::Vec3d p = mesh.point(mesh.vertex_handle(lastestVertex));
	gluProject(p[0], p[1], p[2],  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &winX, &winY, &winZ);
	
	gluUnProject((GLdouble)(x), (GLdouble)( height() - y ), winZ,  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

int InteractiveViewerWidget::find_vertex_using_selected_point()
{
	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	kdTree->annkSearch(tp, 1, nnIdx, dists);
	return nnIdx[0];
}

int InteractiveViewerWidget::find_face_using_selected_point()
{
	int rv = find_vertex_using_selected_point();
	Mesh::VertexFaceIter vf_it = mesh.vf_iter( mesh.vertex_handle(rv) );
	int desiredFace = -1; //double minLen = 10*radius();
	std::vector<OpenMesh::Vec3d> tri_p(3); int tri_count = 0;
	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	for( vf_it; vf_it; ++vf_it )
	{
		tri_count = 0;
		for(Mesh::FaceVertexIter fv_it = mesh.fv_iter(vf_it.handle()); fv_it; ++fv_it)
		{
			tri_p[tri_count] = mesh.point(fv_it); ++tri_count;
		}
		if( check_in_triangle_face(tri_p, resultP) )
		{
			desiredFace = vf_it.handle().idx(); break;
		}
	}
	if(desiredFace < 0)
	{
		for(Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
		{
			tri_count = 0;
			for(Mesh::FaceVertexIter fv_it = mesh.fv_iter(f_it.handle()); fv_it; ++fv_it)
			{
				tri_p[tri_count] = mesh.point(fv_it); ++tri_count;
			}
			if( check_in_triangle_face(tri_p, resultP) )
			{
				desiredFace = f_it.handle().idx(); break;
			}
		}
	}

	return  desiredFace;
}

int InteractiveViewerWidget::find_edge_using_selected_point()
{
	int desiredFace = find_face_using_selected_point(); if(desiredFace < 0) return -1;
	Mesh::FaceHandle fh = mesh.face_handle(desiredFace);
	double min_len= 1e30; int desiredEdge = -1;
	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	for(Mesh::FaceHalfedgeIter fhe_it = mesh.fh_iter(fh); fhe_it; ++fhe_it)
	{
		OpenMesh::Vec3d s = mesh.point( mesh.from_vertex_handle(fhe_it) );
		OpenMesh::Vec3d e = mesh.point( mesh.to_vertex_handle(fhe_it) );
		double dis = OpenMesh::cross(resultP - s, resultP - e).norm() / (s - e).norm();
		if(dis < min_len){ min_len = dis; desiredEdge = mesh.edge_handle(fhe_it.handle()).idx(); }
	}
	
	return desiredEdge;
}

void InteractiveViewerWidget::buildIndex()
{
	if(mesh.n_vertices() == 0)
		return;

	Mesh::VertexIter v_it(mesh.vertices_begin());
	Mesh::VertexIter v_end(mesh.vertices_end());
	Mesh::Point p;
	unsigned nv = mesh.n_vertices();
	ANNpointArray dataPts = annAllocPts(nv, 3);
	int count = 0;
	for(; v_it != v_end; ++v_it)
	{
		p = mesh.point(v_it);
		dataPts[count][0] = p[0]; dataPts[count][1] = p[1]; dataPts[count][2] = p[2];
		++count;
	}

	if(kdTree) delete kdTree;
	kdTree = new ANNkd_tree(dataPts, nv, 3);
}

void InteractiveViewerWidget::BuildCurveIndex()
{
	const auto& edgeshape = CADMesher::globalmodel.edgeshape;
	int nv = 0;
	for (const auto& edge : edgeshape)
	{
		if (edge.id < edge.reversed_edge)
			continue;
		const auto& para = edge.parameters;
		nv += para.cols();
	}
	ANNpointArray dataPts = annAllocPts(nv, 3);
	edgeshapeIndex.clear();
	edgeshapeIndex.reserve(nv);
	int count = 0;
	for (const auto& edge : edgeshape)
	{
		if (edge.id < edge.reversed_edge)
			continue;
		const auto& para = edge.parameters;
		TopLoc_Location loc;
		Handle(Geom_Surface) asurface = BRep_Tool::Surface(CADMesher::globalmodel.faceshape[edge.main_face].face, loc);
		int col = para.cols();
		//Eigen::Matrix3Xd s(3, col);
		for (int i = 0; i < col; ++i)
		{
			auto pos = asurface->Value(para(0, i), para(1, i));
			//dataPts[count][0] = p[0]; dataPts[count][1] = p[1]; dataPts[count][2] = p[2];
			dataPts[count][0] = pos.X(); dataPts[count][1] = pos.Y(); dataPts[count][2] = pos.Z();
			edgeshapeIndex.push_back(edge.id);
			++count;
		}
	}
	if (stripTree) delete stripTree;
	stripTree = new ANNkd_tree(dataPts, nv, 3);
	if_new_mesh = false;
}

//with the first mesh
void InteractiveViewerWidget::draw_interactive_portion(int drawmode)
{
	glViewport ( 0,0, width(),height());
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( &ProjectionMatrix[0] );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( &ModelViewMatrix[0] );

	
	emit draw_from_out_signal();

	{
		//draw select vertex, face, edge.
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

		glPointSize(1);

		switch(mouse_mode_)
		{
		case POINTPICK:
			draw_selected_point();
			break;
		case VERTEXPICK:
			draw_selected_vertex();
			break;
		case FACEPICK:
			draw_selected_face();
			break;
		case EDGEPICK:
			draw_selected_edge();
			break;
		case CURVEPICK:
			draw_selected_curve();
			dprint("fio");
			break;
		case VECTOR_SET:
			draw_vector_set();
			break;
		default:
			draw_selected_vertex();
			draw_selected_face();
			draw_selected_edge();
			draw_selected_curve();
			break;
		}
	}
	
	//showFeature();
	drawMetricConstraint();
	if (draw_new_mesh)
	{
		draw_scene_mesh(drawmode);
	}
}

//with the second mesh
void InteractiveViewerWidget::draw_interactive_portion_mesh2()
{
	return;
}

void InteractiveViewerWidget::draw_selected_point()
{
	glColor3f(1.0, 0.5, 0.0);
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3d(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
	glEnd();
	glPointSize(1);
}

void InteractiveViewerWidget::draw_selected_vertex()
{
	if (selectedVertex.size() > 0)
	{
		Mesh::Point p;
		glColor3f(1.0, 0.5, 0.0);
		glPointSize(12);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < selectedVertex.size(); ++i)
		{
			p = mesh.point(mesh.vertex_handle(selectedVertex[i]));
			glVertex3dv(p.data());
		}
		glEnd();
		glPointSize(1);
	}
}

void InteractiveViewerWidget::draw_selected_face()
{
	if (selectedFace.size() > 0)
	{
		glColor3f(1.0, 0.5, 1.0);
		Mesh::Point p;
		Mesh::ConstFaceVertexIter fv_it;
		Mesh::FaceHandle f_handle;
		for (unsigned int i = 0; i < selectedFace.size(); ++i)
		{
			f_handle = mesh.face_handle(selectedFace[i]);
			fv_it = mesh.fv_iter(f_handle);
			glBegin(GL_POLYGON);
			for (fv_it; fv_it; ++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
	}
}

void InteractiveViewerWidget::draw_selected_edge()
{
	if (selectedEdge.size() > 0)
	{
		glColor3f(1.0, 0.5, 1.0);
		Mesh::Point p1; Mesh::Point p2;
		Mesh::EdgeHandle e_handle;
		Mesh::HalfedgeHandle he_handle;
		for (unsigned int i = 0; i < selectedEdge.size(); ++i)
		{
			e_handle = mesh.edge_handle(selectedEdge[i]);
			he_handle = mesh.halfedge_handle(e_handle, 0);
			p1 = mesh.point(mesh.from_vertex_handle(he_handle));
			p2 = mesh.point(mesh.to_vertex_handle(he_handle));
			glBegin(GL_LINES);
			glVertex3dv(p1.data());
			glVertex3dv(p2.data());
			glEnd();
		}
	}
}

void InteractiveViewerWidget::draw_selected_curve()
{
	if (selectedCurve.size() > 0 && drawCAD)
	{
		const auto& edgeshape = CADMesher::globalmodel.edgeshape;
		glColor3f(1.0, 0.5, 1.0);
		for (auto c : selectedCurve)
		{
			const auto& para = edgeshape[c].parameters;
			TopLoc_Location loc;
			Handle(Geom_Surface) asurface = BRep_Tool::Surface(CADMesher::globalmodel.faceshape[edgeshape[c].main_face].face, loc);
			int col = para.cols();
			//Eigen::Matrix3Xd s(3, col);
			glBegin(GL_LINE_STRIP);
			for (int i = 0; i < col; ++i)
			{
				auto pos = asurface->Value(para(0, i), para(1, i));
				glVertex3d(pos.X(), pos.Y(), pos.Z());
			}
			glEnd();
		}
	}
}

void InteractiveViewerWidget::drawMetricConstraint()
{
	glPointSize(5);
	glBegin(GL_POINTS);
	glColor3d(0.0, 0.1, 0.4);
	for (const auto& mec : metric_constraints)
	{
		glVertex3dv(mesh.point(mec.first).data());
	}
	glEnd();

	glLineWidth(5);
	glBegin(GL_LINES);
	glColor3d(0.9, 0.1, 0.0);
	for (const auto& mec : metric_constraints)
	{
		glVertex3dv(mesh.point(mec.first).data());
		glVertex3dv((mec.second.dir[0] + mesh.point(mec.first)).data());
	}
	glColor3d(0.1, 0.9, 0.0);
	for (const auto& mec : metric_constraints)
	{
		glVertex3dv(mesh.point(mec.first).data());
		glVertex3dv((mec.second.dir[1] + mesh.point(mec.first)).data());
	}

	glLineWidth(1);
	glBegin(GL_LINES);
	glColor3d(0.1, 0.1, 0.3);
	for (const auto& mec : metric_constraints)
	{
		auto xv = mec.second.geo_dir[0];// *mec.second.avg;
		auto yv = mec.second.geo_dir[1];// *mec.second.avg;
		auto pp = mesh.point(mec.first);
		glVertex3dv((xv + pp).data());
		for (double i = 0.0; i < 6.28; i += 0.0157)
		{
			auto fe = cos(i) * xv + sin(i) * yv + pp;
			glVertex3dv(fe.data());
			glVertex3dv(fe.data());
		}
		glVertex3dv((xv + pp).data());
	}
	glEnd();
}

void InteractiveViewerWidget::draw_vector_set()
{
	auto vh = mesh.vertex_handle(lastestVertex);
	if (metric_constraints.find(vh) != metric_constraints.end())
	{
		auto pos = mesh.point(vh);
		double avg_len = metric_constraints[vh].avg;
		OpenMesh::Vec3d zz(0, 0, 0);
		for (auto vf : mesh.vf_range(vh))
		{
			zz += mesh.calc_face_normal(vf);
		}
		zz.normalize();

		OpenMesh::Vec3d xx = metric_constraints[vh].dir[0];
		OpenMesh::Vec3d yy = metric_constraints[vh].dir[1];
		OpenMesh::Vec3d now_pos(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
		//dprint(now_pos);
		if (choose == 0)
		{
			xx = now_pos - pos;
			xx -= xx.dot(zz)*zz;
			yy = yy.norm()*zz.cross(xx.normalized());
		}
		else
		{
			yy = now_pos - pos;
			yy -= yy.dot(zz)*zz;
			xx = xx.norm()*yy.normalized().cross(zz);
		}

		metric_constraints[vh].dir[0] = xx;
		metric_constraints[vh].dir[1] = yy;
		MeshParam->set_red_length(xx.norm());
		MeshParam->set_green_length(yy.norm());
	}

	metric_constraints[vh].dir[0] *= MeshParam->get_red_length() / metric_constraints[vh].dir[0].norm();
	metric_constraints[vh].dir[1] *= MeshParam->get_green_length() / metric_constraints[vh].dir[1].norm();
}

void InteractiveViewerWidget::draw_scene(int drawmode)
{
	//if (!mesh.n_vertices()) { return; }
	draw_interactive_portion_mesh2();
	draw_interactive_portion(drawmode);

	if( !draw_new_mesh )
	{
		MeshViewerWidget::draw_scene(drawmode);
	}
}

void InteractiveViewerWidget::render_text_slot(OpenMesh::Vec3d pos, QString str)
{
	/*GLdouble  winX, winY, winZ;
	GLint     viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	gluProject(pos[0],pos[1],pos[2],&ModelViewMatrix[0][0],&ProjectionMatrix[0][0],viewport,&winX,&winY,&winZ);
	int x = (long)winX;
	int y = viewport[3]-(long)winY;
	render_text(x,y,str);*/
	render_text(pos[0],pos[1],pos[2],str);
}

void InteractiveViewerWidget::SetCADFileName(QString &fileName) {
	int id = fileName.lastIndexOf("/");
	assert(id != -1);
	CADFileName = fileName.right(fileName.length() - id - 1);
	id = CADFileName.lastIndexOf(".");
	assert(id != -1);
	CADFileName.truncate(id);
	dprint(fileName, "\nCAD name:", CADFileName);
	//BrepFileName = fileName;
};

void InteractiveViewerWidget::generateTriMesh(double ratio)
{
	if (!occreader) return;
	occreader->initialRate = ratio;
	occreader->Discrete_Edge();
	occreader->Face_type();
	occreader->C0_Feature();
	occreader->Curvature_Feature();
	occreader->Set_TriMesh();
	occreader->MergeModel(CADMesher::globalmodel.initial_trimesh, occreader->Surface_TriMeshes);
	CADMesher::globalmodel.init_trimesh_tree = new ClosestPointSearch::AABBTree(CADMesher::globalmodel.initial_trimesh);
	occreader->SetTriFeature();
	tri2poly(CADMesher::globalmodel.initial_trimesh, mesh, true);
	setMeshMode(N_MESH_MODES);
	initMeshStatusAndNormal(mesh);
	drawCAD = false;
	setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
	setMouseMode(InteractiveViewerWidget::TRANS);
	ifGenerateTriMesh = true;
	ifGeneratePolyMesh = false;
	dprint("Mesh Avg Length:", meshAverageLength(mesh));
	updateGL();
	//compute_principal_curvature(&CADMesher::globalmodel.initial_trimesh, K1, K2, D1, D2);
	fei = true;
}

void InteractiveViewerWidget::generatePolyMesh(double ratio, int quad_num, double initial_ratio, double increase_ratio)
{
	if (!occreader) return;
	occreader->initialRate = ratio;
	occreader->offset_quad_num = quad_num;
	occreader->offset_initial_ratio = initial_ratio;
	occreader->offset_increase_ratio = increase_ratio;
	occreader->Discrete_Edge();
	occreader->Face_type();
	occreader->C0_Feature();
	occreader->Curvature_Feature();
	//occreader->Set_Offset_Info(selectedCurve);
	occreader->Set_PolyMesh();
	occreader->MergeModel(CADMesher::globalmodel.initial_polymesh, occreader->Surface_PolyMeshes);
	TriMesh tm(CADMesher::globalmodel.initial_polymesh);
	CADMesher::globalmodel.init_trimesh_tree = new ClosestPointSearch::AABBTree(tm);
	occreader->SetPolyFeature();
	polycopy(CADMesher::globalmodel.initial_polymesh, mesh, true);
	initMeshStatusAndNormal(mesh);
	drawCAD = false;
	setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
	setMouseMode(InteractiveViewerWidget::TRANS);
	ifGenerateTriMesh = false;
	ifGeneratePolyMesh = true;
	dprint("Mesh Avg Length:", meshAverageLength(mesh));
	setMeshMode(2);
	updateGL();
	fei = true;
}

void InteractiveViewerWidget::showFeature()
{
	ifDrawFeature = !ifDrawFeature;
	setMouseMode(InteractiveViewerWidget::TRANS);
	updateGL();
}

void InteractiveViewerWidget::showIsotropicMesh(double tl)
{
	if (ifGenerateTriMesh)
	{
		dprint(tl);
		CADMesher::globalmodel.isotropic_trimesh = CADMesher::globalmodel.initial_trimesh;
		tmr = new CADMesher::TriangleMeshRemeshing(&(CADMesher::globalmodel.isotropic_trimesh), tl);
		tmr->run();
		tri2poly(CADMesher::globalmodel.isotropic_trimesh, mesh, true);
		initMeshStatusAndNormal(mesh);
	}
	else if (ifGeneratePolyMesh)
	{
		CADMesher::globalmodel.isotropic_polymesh = CADMesher::globalmodel.initial_polymesh;
		tmr = new CADMesher::TriangleMeshRemeshing(&(CADMesher::globalmodel.isotropic_polymesh));
		tmr->run();
		mesh = CADMesher::globalmodel.isotropic_polymesh;
		initMeshStatusAndNormal(mesh);
	}
	drawCAD = false;
	setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
	setMouseMode(InteractiveViewerWidget::TRANS);
	updateGL();
}

void InteractiveViewerWidget::offsetinfo(int quad_num, double initial_ratio, double increase_ratio)
{
	occreader->quad_mun_info.push_back(quad_num);
	occreader->initial_ratio_info.push_back(initial_ratio);
	occreader->increase_ratio_info.push_back(increase_ratio);
}

//#include "..\src\Algorithm\SurfaceMesher\Optimizer\AnisotropicMeshRemeshing.h"
void InteractiveViewerWidget::showAnisotropicMesh(double tl)
{
	if (ifGenerateTriMesh)
	{
		CADMesher::AnisotropicMeshRemeshing* amr = new CADMesher::AnisotropicMeshRemeshing();
		CADMesher::globalmodel.isotropic_trimesh = CADMesher::globalmodel.initial_trimesh;
		amr->SetMesh(&(CADMesher::globalmodel.isotropic_trimesh));
		amr->load_ref_mesh(&(CADMesher::globalmodel.initial_trimesh), tl, (occreader->bbmax - occreader->bbmin).norm());
		//double tl = amr->get_ref_mesh_ave_anisotropic_edge_length();
		dprint("anisotropic edge length:", tl);
		amr->set_metric(metric_constraints);
		amr->do_remeshing(amr->get_ref_mesh_ave_anisotropic_edge_length(), 1.5);

		tri2poly(CADMesher::globalmodel.isotropic_trimesh, mesh, true);
		initMeshStatusAndNormal(mesh);
	}
	drawCAD = false;
	setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
	setMouseMode(InteractiveViewerWidget::TRANS);
	metric_constraints.clear();
	updateGL();
}

#include "..\src\Toolbox\filesOperator.h"
#include "..\src\Algorithm\SurfaceMesher\Generator\Iso_Mesh.h"
#include "..\src\Algorithm\SurfaceMesher\Optimizer\TriangleMeshRemeshing.h"
//#include <fstream>
void InteractiveViewerWidget::showDebugTest()
{
#pragma region step files test
	{
		std::vector<std::string> allFileName;
		std::string path = "..\\model\\CAD";
		getFiles(path, allFileName);
		std::ofstream fileWriter;

#if 0   //导入各向同性数据
		int i = 0;
		fileWriter.open("C:\\Users\\1\\Desktop\\test\\test2\\New.csv", std::ios::app);
		for (; i < allFileName.size();)
		{
			auto fileName = allFileName[i];
			dprint("\n\n\nfile index:\t", i++, "\nfileName:\t", fileName);
			auto time_start = clock();
			CADMesher::OccReader occreader1(QString::fromStdString(fileName));
			occreader1.initialRate = 0.003;
			occreader1.Discrete_Edge();
			occreader1.Face_type();
			occreader1.C0_Feature();
			occreader1.Curvature_Feature();
			occreader1.Set_TriMesh();
			occreader1.MergeModel(CADMesher::globalmodel.initial_trimesh, occreader1.Surface_TriMeshes);
			CADMesher::globalmodel.init_trimesh_tree = new ClosestPointSearch::AABBTree(CADMesher::globalmodel.initial_trimesh);
			occreader1.SetTriFeature();
			CADMesher::globalmodel.isotropic_trimesh = CADMesher::globalmodel.initial_trimesh;
			CADMesher::TriangleMeshRemeshing trm1(&(CADMesher::globalmodel.isotropic_trimesh));
			trm1.run();
			double iso_time = (clock() - time_start) / 1000.0;
			truncateFilePath(fileName);
			truncateFileExtension(fileName);
			//if (!OpenMesh::IO::write_mesh(CADMesher::globalmodel.isotropic_trimesh, "C:\\Users\\1\\Desktop\\test\\test2\\ISO\\" + fileName + ".obj"));
			//{
			//	std::cerr << "fail";
			//}
			TriMeshQualityHelper tmqh(&CADMesher::globalmodel.isotropic_trimesh);
			fileWriter << fileName << ",";
			fileWriter << CADMesher::globalmodel.isotropic_trimesh.n_vertices() << ",";
			fileWriter << tmqh.getMinAngle() << "," << tmqh.getMaxAngle() << "," << tmqh.getAvgAngle() << ",";
			fileWriter << tmqh.getMinQuality() << "," << tmqh.getAvgQuality() << "," << iso_time;
			fileWriter << std::endl;
			CADMesher::globalmodel.clear();
		}
#else   //导入各向异性数据
		int i = 4;
		fileWriter.open("C:\\Users\\1\\Desktop\\test\\test3\\AnIsoRawData.csv", std::ios::app);
		for (; i < 20;)
		{
			auto fileName = allFileName[i];
			dprint("\n\n\nfile index:\t", i++, "\nfileName:\t", fileName);
			auto time_start = clock();
			CADMesher::OccReader occreader1(QString::fromStdString(fileName));
			occreader1.initialRate = 0.004;
			occreader1.Discrete_Edge();
			occreader1.Face_type();
			occreader1.C0_Feature();
			occreader1.Curvature_Feature();
			occreader1.Set_TriMesh();
			occreader1.MergeModel(CADMesher::globalmodel.initial_trimesh, occreader1.Surface_TriMeshes);
			CADMesher::globalmodel.init_trimesh_tree = new ClosestPointSearch::AABBTree(CADMesher::globalmodel.initial_trimesh);
			occreader1.SetTriFeature();
			compute_principal_curvature(&CADMesher::globalmodel.initial_trimesh, K1, K2, D1, D2);
			CADMesher::AnisotropicMeshRemeshing* amr = new CADMesher::AnisotropicMeshRemeshing();
			CADMesher::globalmodel.isotropic_trimesh = CADMesher::globalmodel.initial_trimesh;
			amr->SetMesh(&(CADMesher::globalmodel.isotropic_trimesh));
			double tl = meshAverageLength(CADMesher::globalmodel.isotropic_trimesh);
			amr->load_ref_mesh(&(CADMesher::globalmodel.initial_trimesh), tl, (occreader1.bbmax - occreader1.bbmin).norm());
			amr->set_metric(metric_constraints);
			amr->do_remeshing(amr->get_ref_mesh_ave_anisotropic_edge_length(), 1.5);
			double aniso_time = (clock() - time_start) / 1000.0;
			truncateFilePath(fileName);
			truncateFileExtension(fileName);
			if (!OpenMesh::IO::write_mesh(CADMesher::globalmodel.isotropic_trimesh, "C:\\Users\\1\\Desktop\\test\\test3\\ANISO\\" + fileName + ".obj"));
			{
				std::cerr << "fail";
			}
			TriMeshQualityHelper tmqh(&CADMesher::globalmodel.isotropic_trimesh);
			fileWriter << fileName << ",";
			fileWriter << CADMesher::globalmodel.isotropic_trimesh.n_vertices() << ",";
			fileWriter << amr->MQE << "," << amr->AQE << ",";
			fileWriter << (amr->MinAngle / 180.) * PI << "," << (amr->MaxAngle / 180.) * PI << "," << (amr->AveAngle / 180.) * PI << ",";
			fileWriter << amr->MinArea << "," << amr->MaxArea << ",";
			fileWriter << amr->MinRER << "," << amr->MaxRER << "," << aniso_time;
			fileWriter << std::endl;
			CADMesher::globalmodel.clear();
		}

#endif
		fileWriter.close();
	}
}



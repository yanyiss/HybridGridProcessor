#ifndef MESHPROCESSING_MAIN_VIEWWE_WIDGET_H
#define MESHPROCESSING_MAIN_VIEWWE_WIDGET_H

#include <QtGui>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
//main widget
#include "InteractiveViewerWidget.h"
#include "MeshParamDialog.h"

class MainViewerWidget : public QDialog
{
	Q_OBJECT
public:
	MainViewerWidget(QWidget* _parent=0);
	~MainViewerWidget();

	void setDrawMode(int dm)
	{
		MeshViewer->setDrawMode(dm);
	}
	void setMouseMode(int dm)
	{
		MeshViewer->setMouseMode(dm);
	}

	void set_show_BBox()
	{
		MeshViewer->set_draw_bbox_ok();
	}
	void set_show_mesh_boundary()
	{
		MeshViewer->set_draw_mesh_boundary_ok();
	}
	void openMesh_fromMain(char* filename)
	{
		QString str(filename);
		open_mesh_gui(filename);
	}

	void edit_undo()
	{
		MeshViewer->edit_undo_viewer();
	}

	void edit_redo()
	{
		MeshViewer->edit_redo_viewer();
	}

public slots:
	void open_query()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open mesh file"),
			tr("../model/CAD"),
			tr( "STEP Files (*.STEP;*.STP;*.step;*.stp);;"
				"IGES Files (*.IGES;*.IGS;*.iges;*.igs);;"
				"OBJ Files (*.obj);;""OFF Files (*.off);;""PLY Files (*.ply);;""STL Files (*.stl);;"
				"All Files (*)"));
		if (fileName.isEmpty())
			return;
		if (fileName.endsWith(".IGES") || fileName.endsWith(".STEP") ||
			fileName.endsWith(".IGS") || fileName.endsWith(".STP") ||
			fileName.endsWith(".iges") || fileName.endsWith(".step") ||
			fileName.endsWith(".igs") || fileName.endsWith(".stp"))
		{
			MeshViewer->SetCADFileName(fileName);
			CADMesher::globalmodel.clear();
			if (MeshViewer->occreader)
				delete MeshViewer->occreader;
			MeshViewer->occreader = new CADMesher::OccReader(fileName);
			auto occreader = MeshViewer->occreader;

			MeshViewer->calcCADStrip();
			MeshViewer->clearAllMesh();
			MeshViewer->set_scene_pos(0.5*(occreader->bbmax + occreader->bbmin), 0.5*(occreader->bbmax - occreader->bbmin).norm());

			MeshViewer->drawCAD = true;
			MeshViewer->setDrawMode(InteractiveViewerWidget::WIRE_FRAME);
			MeshViewer->setMouseMode(InteractiveViewerWidget::TRANS);

		}
		else
		{
			open_mesh_gui(fileName);
		}
	}
	void open_mesh_query()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open mesh file"),
			tr("../model/mesh"),
			tr("OBJ Files (*.obj);;"
			"OFF Files (*.off);;"
			"PLY Files (*.ply);;"
			"STL Files (*.stl);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			open_mesh_gui(fileName);
		}
	}
	void open_CAD_query()
	{
		open_query();
//		QString fileName = QFileDialog::getOpenFileName(this,
//			tr("Open mesh file"),
//			tr("../model/CAD"),
//			tr("IGES Files (*.IGES;*.IGS;*.iges;*.igs);;"
//				"STEP Files (*.STEP;*.STP;*.step;*.stp);;")
//		);
//		if (!fileName.isEmpty())
//		{
//			if (fileName.endsWith(".stp") || fileName.endsWith(".igs") ||
//				fileName.endsWith(".step") || fileName.endsWith(".iges") ||
//				fileName.endsWith(".IGS") || fileName.endsWith(".STP") ||
//				fileName.endsWith(".STEP") || fileName.endsWith(".IGES"))
//			{
//				MeshViewer->SetCADFileName(fileName);
//				CADMesher::globalmodel.clear();
//				timeRecorder tr;
//				CADMesher::Iso_Mesh iso_mesh(fileName);
//				tr.out("time of generating isotropic mesh:");
//#if 1
//				Mesh me;
//				tri2poly(CADMesher::globalmodel.initial_trimesh, me, true);
//#else
//				Mesh me = CADMesher::globalmodel.initial_polymesh;
//#endif
//				initMeshStatusAndNormal(me);
//				open_mesh_gui(me);
//				/*Mesh me;
//				bool read_OK = OpenMesh::IO::read_mesh(me, "step_to_obj.obj");
//				open_mesh_gui(me);*/
//			}
//		}
	}
	void save_mesh_query() 
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Save mesh file"),
			tr("../models/untitled.off"),
			tr("OBJ Files (*.obj);;"
			"OFF Files (*.off);;"
			"PLY Files (*.ply);;"
			"STL Files (*.stl);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			save_mesh_gui(fileName);
		}
	}
	void saveOpenGLScreen()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			("Save screen as image file"),
			("../Results/untitled.png"),
			("PNG Files (*.png);;BMP Files (*.bmp);;JPG Files (*.jpg);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			save_screen_gui(fileName);
		}
	}
	void save_opengl_screen(QString str)
	{
		MeshViewer->saveScreen(str.toLocal8Bit());
	}
	virtual void update_mesh()
	{
		if( MeshViewer->mesh_ref().n_vertices() != 0 )
		{
			MeshViewer->updateMesh();
		}
	}

	virtual void clear_all_mesh()
	{
		if(LoadMeshSuccess)
		{
			LoadMeshSuccess = false;
			MeshViewer->clearAllMesh();
		}
	}

	virtual void clear_all_selected()
	{
		if(LoadMeshSuccess)
		{
			MeshViewer->clearSelectedData();
			MeshViewer->updateGL();
		}
	}

	void LoadMeshFromInner(bool OK,QString fname)
	{
		LoadMeshSuccess = OK;
		if(LoadMeshSuccess)
		{
			SetMeshForALL();
		}
		emit( haveLoadMesh(fname) );
	};

private slots:
	void print_info();
	void submit_info();

signals:
	void haveLoadMesh(QString filePath);
	void setMouseMode_signal_main(int);
	void setDrawMode_signal_main(int);

	void set_edit_undo_enable_signal(bool);
	void set_edit_redo_enable_signal(bool);

protected:
	virtual void initViewerWindow();
	virtual void createParamIDialog();
	virtual void createViewerDialog();
	virtual void save_mesh_gui(QString fname);
	virtual void open_mesh_gui(QString fname);
	virtual void open_mesh_gui(Mesh &aMesh);
	virtual void save_screen_gui(QString fname);

public:
	void generateTriMesh()
	{
		MeshViewer->generateTriMesh(MeshParam->get_sample_ratio_AM());
		MeshParam->set_target_edge_length_AM(MeshViewer->getAngLen());
	}
	void generatePolyMesh()
	{
		MeshViewer->generatePolyMesh(MeshParam->get_sample_ratio_AM(), MeshParam->get_quad_num(), MeshParam->get_initial_ratio(), MeshParam->get_increase_ratio());
		MeshParam->set_target_edge_length_AM(MeshViewer->getAngLen());
	}
	void showFeature() {
		MeshViewer->showFeature();
	}
	void showIsotropicMesh() {
		MeshViewer->showIsotropicMesh(MeshParam->get_target_edge_length_AM());
	}
	void showAnisotropicMesh() {
		MeshViewer->showAnisotropicMesh();
	}
	void showDebugTest() {
		MeshViewer->showDebugTest();
	}
protected:
	bool LoadMeshSuccess;

private:
	InteractiveViewerWidget* MeshViewer;
	MeshParamDialog* MeshParam;
	
	void SetMeshForALL( )
	{
	}

#pragma region Auxiliary Function
public:
	void aux_inverse_mesh_connectivity()
	{
		MeshViewer->inverse_mesh_connectivity();
	}

	void aux_scale_mesh_using_BBox(int max_len)
	{
		MeshViewer->scale_mesh_using_BBox(max_len);
	}

	void aux_split_quad_mesh()
	{
		MeshViewer->split_quad_mesh();
	}

	void transform_mesh(const std::vector<double>& m)
	{
		MeshViewer->transform_mesh(m);
	}

	void aux_find_vertex_by_id(int id)
	{
		MeshViewer->find_vertex_by_id(id);
	}

	void aux_find_face_by_id(int id)
	{
		MeshViewer->find_face_by_id(id);
	}

	void aux_find_edge_by_id(int id)
	{
		MeshViewer->find_edge_by_id(id);
	}

	void aux_find_vertex_by_valance(int valance)
	{
		MeshViewer->find_vertex_by_valance( valance );
	}

	void aux_delete_vertex_valence_four()
	{
		MeshViewer->delete_vertex_valence_four( );
	}

	void aux_delete_vertex_valence_three()
	{
		MeshViewer->delete_vertex_valence_three( );
	}

	void aux_split_vertex_valence_eight()
	{
		MeshViewer->split_vertex_valence_eight( );
	}

#pragma endregion

};


#endif
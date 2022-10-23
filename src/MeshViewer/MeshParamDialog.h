#ifndef MESHPROCESSING_MESHPARAMDIALOG_H
#define MESHPROCESSING_MESHPARAMDIALOG_H

#include <QDialog>
#include <QtGui>
#include <QtWidgets>

class MeshParamDialog : public QDialog
{
	Q_OBJECT
public:
	MeshParamDialog(QWidget* parent=0);
	~MeshParamDialog();

	QSize sizeHint()
	{
		QRect rect = QApplication::desktop()->screenGeometry();
		return QSize( int( rect.width()*0.15), rect.height() );
	}

private:
	QTabWidget* tabWidget;

signals:
	void print_info_signal();

public:
	double get_target_edge_length_AM()
	{
		return target_edge_length_line_AM->text().toDouble();
	}
	double get_sample_ratio_AM()
	{
		return sample_ratio_line_AM->text().toDouble();
	}

signals:
	void load_ref_mesh_AM_signal();
	void do_remehsing_AM_signal();

private:
	QWidget* Basic_Operation_And_Information;
	QScrollArea *view_BOI;

	QLabel* leftLabel_BOI;
	QPushButton* print_info;

	QPushButton* load_ref_mesh_AM;
	QPushButton* do_remehsing_AM;
	QLabel* target_edge_length_AM; QLineEdit* target_edge_length_line_AM;
	QLabel* sample_ratio_AM; QLineEdit* sample_ratio_line_AM;

private:
	void create_Basic_Operation_Information_Widget();

private:
	void initDialog();
	void createWidget();
	void createLayout();

};
#ifndef MESHPROCESSING_MESHPARAMDIALOG_H
#define MESHPROCESSING_MESHPARAMDIALOG_H

#include <QDialog>
#include <QtGui>
#include <QtWidgets>

class MeshParamDialog : public QDialog
{
	Q_OBJECT
public:
	MeshParamDialog(QWidget* parent = 0);
	~MeshParamDialog();

	QSize sizeHint()
	{
		QRect rect = QApplication::desktop()->screenGeometry();
		return QSize(int(rect.width() * 0.15), rect.height());
	}

private:
	QTabWidget* tabWidget;

public:
	double get_target_edge_length_AM()
	{
		return target_edge_length_line_AM->text().toDouble();
	}
	double get_sample_ratio_AM()
	{
		return sample_ratio_line_AM->text().toDouble();
	}

signals:
	void load_ref_mesh_AM_signal();
	void do_remehsing_AM_signal();


private:
	QWidget* Basic_Operation_And_Information;
	QScrollArea* view_BOI;

	QLabel* leftLabel_BOI;
	QPushButton* print_info;

	QPushButton* load_ref_mesh_AM;
	QPushButton* do_remehsing_AM;
	QLabel* target_edge_length_AM; QLineEdit* target_edge_length_line_AM;
	QLabel* sample_ratio_AM; QLineEdit* sample_ratio_line_AM;

private:
	void create_Basic_Operation_Information_Widget();

private:
	void initDialog();
	void createWidget();
	void createLayout();

};

#endif
#endif
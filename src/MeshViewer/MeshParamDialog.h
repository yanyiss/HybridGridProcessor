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
		return QSize(int(rect.width()*0.15), rect.height());
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
	void set_target_edge_length_AM(double tel)
	{
		target_edge_length_line_AM->setText(QString::number(tel, 'g', 4));
	}
	void set_sample_ratio_AM(double tel)
	{
		sample_ratio_line_AM->setText(QString::number(tel, 'g', 4));
	}

	int get_quad_num()
	{
		int t = (int)(quad_num_line->text().toDouble());
		if (t < 3) return 3;
		if (t > 30) return 30;
		return t;
	}

	double get_initial_ratio()
	{
		double t = initial_ratio_line->text().toDouble();
		if (t < 1e-4 || t > 1) return 0.01;
		return t;
	}

	double get_increase_ratio()
	{
		double t = increase_ratio_line->text().toDouble();
		if (t < 1) return 1;
		if (t > 3) return 3;
		return t;
	}

	double get_red_length()
	{
		return red_length_line->text().toDouble();
	}
	void set_red_length(double r)
	{
		red_length_line->setText(QString::number(r, 'g', 4));
	}
	double get_green_length()
	{
		return green_length_line->text().toDouble();
	}
	void set_green_length(double g)
	{
		green_length_line->setText(QString::number(g, 'g', 4));
	}

signals:
	void load_ref_mesh_AM_signal();
	void do_remehsing_AM_signal();
	//void submit_info_signal();

private:
	QWidget* Basic_Operation_And_Information;
	QScrollArea *view_BOI;

	QLabel* leftLabel_BOI;
	QPushButton* print_info;

	QPushButton* load_ref_mesh_AM;
	QPushButton* do_remehsing_AM;
	QPushButton* submit_offset_info;
	QLabel* target_edge_length_AM; QLineEdit* target_edge_length_line_AM;
	QLabel* sample_ratio_AM; QLineEdit* sample_ratio_line_AM;
	QLabel* quad_num; QLineEdit* quad_num_line;
	QLabel* initial_ratio; QLineEdit* initial_ratio_line;
	QLabel* increase_ratio; QLineEdit* increase_ratio_line;
	QLabel* red_length; QLineEdit* red_length_line;
	QLabel* green_length; QLineEdit* green_length_line;

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
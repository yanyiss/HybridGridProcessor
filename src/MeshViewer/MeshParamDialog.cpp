#include "MeshParamDialog.h"
#include <QApplication>
#include <QDesktopWidget>

#include "..\src\Dependency\Common\CommonDefinitions.h"

MeshParamDialog::MeshParamDialog(QWidget* parent /* = 0 */)
	:QDialog(parent)
{
	initDialog();
}

MeshParamDialog::~MeshParamDialog()
{
}

void MeshParamDialog::initDialog()
{
	createWidget();
	createLayout();
}

void MeshParamDialog::createWidget()
{
	create_Basic_Operation_Information_Widget();
}

void MeshParamDialog::createLayout()
{
	tabWidget = new QTabWidget();
	tabWidget->addTab(view_BOI, "Parameter Setting");

	QGridLayout *layout = new QGridLayout();
	layout->addWidget(tabWidget, 0, 0, 1, 1);
	setLayout(layout);
}

void MeshParamDialog::create_Basic_Operation_Information_Widget()
{
	//print_info = new QPushButton("Print Mesh Information");
	load_ref_mesh_AM = new QPushButton("Load Ref Mesh");
	do_remehsing_AM = new QPushButton("Do Re meshing");
	//submit_offset_info = new QPushButton("Submit");
	leftLabel_BOI = new QLabel("");
	sample_ratio_AM = new QLabel("Model Ratio:");
	quad_num = new QLabel("Offset Num:");
	initial_ratio = new QLabel("Initial Ratio:");
	increase_ratio = new QLabel("Increase Ratio");
	sample_ratio_line_AM = new QLineEdit("0.004");
	sample_ratio_line_AM->setValidator(new QDoubleValidator(0.0, 1000.0, 10));
	target_edge_length_AM = new QLabel("Target Edge Length:");
	target_edge_length_line_AM = new QLineEdit("1.0");
	target_edge_length_line_AM->setValidator(new QDoubleValidator(0.0, 1000000.0, 10));
	quad_num_line = new QLineEdit("15");
	initial_ratio_line = new QLineEdit("0.01");
	increase_ratio_line = new QLineEdit("1.3");

	//QGridLayout* mainLayout = new QGridLayout(); int main_index = 0;
	//mainLayout->addWidget(print_info, main_index, 0, 1, 2); main_index += 1;
	//mainLayout->addWidget(leftLabel_BOI, main_index, 0, 1, 40);




	QGridLayout* LCOT_layout = new QGridLayout(); int LCOT_layout_index = 0;
	LCOT_layout->setSizeConstraint(QLayout::SetFixedSize);
	//LCOT_layout->addWidget(print_info, LCOT_layout_index, 0, 1, 20); LCOT_layout_index += 1;
	//LCOT_layout->addWidget(leftLabel_BOI, LCOT_layout_index, 0, 1, 40);
	//LCOT_layout->addWidget(load_ref_mesh_AM, LCOT_layout_index, 0, 1, 2); LCOT_layout_index += 1;
	LCOT_layout->addWidget(sample_ratio_AM, LCOT_layout_index, 0, 1, 1);
	LCOT_layout->addWidget(sample_ratio_line_AM, LCOT_layout_index, 1, 1, 1); LCOT_layout_index += 1;
	LCOT_layout->addWidget(target_edge_length_AM, LCOT_layout_index, 0, 1, 1);
	LCOT_layout->addWidget(target_edge_length_line_AM, LCOT_layout_index, 1, 1, 1); LCOT_layout_index += 1;
	LCOT_layout->addWidget(quad_num, LCOT_layout_index, 0, 1, 1);
	LCOT_layout->addWidget(quad_num_line, LCOT_layout_index, 1, 1, 1); LCOT_layout_index += 1;
	LCOT_layout->addWidget(initial_ratio, LCOT_layout_index, 0, 1, 1);
	LCOT_layout->addWidget(initial_ratio_line, LCOT_layout_index, 1, 1, 1); LCOT_layout_index += 1;
	LCOT_layout->addWidget(increase_ratio, LCOT_layout_index, 0, 1, 1);
	LCOT_layout->addWidget(increase_ratio_line, LCOT_layout_index, 1, 1, 1); LCOT_layout_index += 1;
	//LCOT_layout->addWidget(submit_offset_info, LCOT_layout_index, 0, 1, 2); LCOT_layout_index += 1;
	//LCOT_layout->addWidget(do_remehsing_AM, LCOT_layout_index, 0, 1, 2); LCOT_layout_index += 1;





	Basic_Operation_And_Information = new QWidget();
	Basic_Operation_And_Information->setLayout(LCOT_layout);
	//Basic_Operation_And_Information->setAutoFillBackground(true);
	//Basic_Operation_And_Information->setBackgroundRole(QPalette::Background); 
	//Basic_Operation_And_Information->setPalette(QPalette(QColor(255, 255, 0)));

	view_BOI = new QScrollArea;
	view_BOI->setFocusPolicy(Qt::NoFocus);
	view_BOI->setFrameStyle(QFrame::NoFrame);
	view_BOI->setWidget(Basic_Operation_And_Information);
	view_BOI->setWidgetResizable(true);

	//connect(print_info, SIGNAL(clicked()), SIGNAL(print_info_signal()));
	connect(load_ref_mesh_AM, SIGNAL(clicked()), SIGNAL(load_ref_mesh_AM_signal()));
	connect(do_remehsing_AM, SIGNAL(clicked()), SIGNAL(do_remehsing_AM_signal()));
	//connect(submit_offset_info, SIGNAL(clicked()), SIGNAL(submit_info_signal()));
}

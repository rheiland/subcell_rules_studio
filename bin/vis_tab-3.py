import sys
import os
import time
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
from pathlib import Path
# from ipywidgets import Layout, Label, Text, Checkbox, Button, BoundedIntText, HBox, VBox, Box, \
    # FloatText, Dropdown, SelectMultiple, RadioButtons, interactive
# import matplotlib.pyplot as plt
# from matplotlib.colors import BoundaryNorm
# from matplotlib.ticker import MaxNLocator
# from matplotlib.collections import LineCollection
# from matplotlib.patches import Circle, Ellipse, Rectangle
# from matplotlib.collections import PatchCollection
# import matplotlib.colors as mplc
# from matplotlib import gridspec
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
# from collections import deque
import glob

from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QFrame,QApplication,QWidget,QTabWidget,QFormLayout,QLineEdit, QHBoxLayout,QVBoxLayout,QRadioButton,QLabel,QCheckBox,QComboBox,QScrollArea,  QMainWindow,QGridLayout, QPushButton, QFileDialog, QMessageBox

import numpy as np
import scipy.io
# import matplotlib
# matplotlib.use('Qt5Agg')
# import matplotlib.pyplot as plt
#from QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
# from mpl_toolkits.axes_grid1 import make_axes_locatable

# from PyQt5 import QtCore, QtWidgets

# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
# from matplotlib.figure import Figure

class Vis(QWidget):

    def __init__(self, nanohub_flag):
        super().__init__()
        # global self.config_params

        self.nanohub_flag = nanohub_flag

        self.animating_flag = False

        self.xml_root = None
        self.current_svg_frame = 0
        self.timer = QtCore.QTimer()
        # self.t.timeout.connect(self.task)
        self.timer.timeout.connect(self.play_plot_cb)

        # self.tab = QWidget()
        # self.tabs.resize(200,5)
        
        self.num_contours = 15
        self.num_contours = 25
        self.num_contours = 50
        self.fontsize = 5

        self.plot_svg_flag = True
        # self.plot_svg_flag = False
        self.field_index = 4  # substrate (0th -> 4 in the .mat)

        self.use_defaults = True
        self.title_str = ""

        # self.config_file = "mymodel.xml"
        self.reset_model_flag = True
        self.xmin = -80
        self.xmax = 80
        self.x_range = self.xmax - self.xmin

        self.ymin = -50
        self.ymax = 100
        self.y_range = self.ymax - self.ymin

        self.aspect_ratio = 0.7

        self.show_nucleus = False
        self.show_edge = False
        self.alpha = 0.7

        basic_length = 12.0
        self.figsize_width_substrate = 15.0  # allow extra for colormap
        self.figsize_height_substrate = basic_length

        self.figsize_width_2Dplot = basic_length
        self.figsize_height_2Dplot = basic_length

        # self.width_substrate = basic_length  # allow extra for colormap
        # self.height_substrate = basic_length

        self.figsize_width_svg = basic_length
        self.figsize_height_svg = basic_length

        # self.output_dir = "/Users/heiland/dev/PhysiCell_V.1.8.0_release/output"
        # self.output_dir = "output"
        # self.output_dir = "../tmpdir"   # for nanoHUB
        # self.output_dir = "tmpdir"   # for nanoHUB
        self.output_dir = "."   # for nanoHUB


        # do in create_figure()?
        # xlist = np.linspace(-3.0, 3.0, 50)
        # print("len(xlist)=",len(xlist))
        # ylist = np.linspace(-3.0, 3.0, 50)
        # X, Y = np.meshgrid(xlist, ylist)
        # Z = np.sqrt(X**2 + Y**2) + 10*np.random.rand()
        # self.cmap = plt.cm.get_cmap("viridis")
        # self.cs = self.ax0.contourf(X, Y, Z, cmap=self.cmap)
        # self.cbar = self.figure.colorbar(self.cs, ax=self.ax)


        #-------------------------------------------
        label_width = 110
        domain_value_width = 100
        value_width = 60
        label_height = 20
        units_width = 70

        self.substrates_cbox = QComboBox(self)

        self.myscroll = QScrollArea()  # might contain centralWidget
        self.create_figure()

        self.config_params = QWidget()

        self.main_layout = QVBoxLayout()

        self.vbox = QVBoxLayout()
        self.vbox.addStretch(0)

        #------------------
        controls_hbox = QHBoxLayout()
        w = QPushButton("Directory")
        w.clicked.connect(self.open_directory_cb)
        # if self.nanohub_flag:
        if self.nanohub_flag:
            w.setEnabled(False)  # for nanoHUB
        controls_hbox.addWidget(w)

        # self.output_dir = "/Users/heiland/dev/PhysiCell_V.1.8.0_release/output"
        self.output_dir_w = QLineEdit()
        self.output_dir_w.setFixedWidth(domain_value_width)
        # w.setText("/Users/heiland/dev/PhysiCell_V.1.8.0_release/output")
        self.output_dir_w.setText(self.output_dir)
        if self.nanohub_flag:
            self.output_dir_w.setEnabled(False)  # for nanoHUB
        # w.textChanged[str].connect(self.output_dir_changed)
        # w.textChanged.connect(self.output_dir_changed)
        controls_hbox.addWidget(self.output_dir_w)

        self.first_button = QPushButton("<<")
        self.first_button.clicked.connect(self.first_plot_cb)
        controls_hbox.addWidget(self.first_button)

        self.back_button = QPushButton("<")
        self.back_button.clicked.connect(self.back_plot_cb)
        controls_hbox.addWidget(self.back_button)

        self.forward_button = QPushButton(">")
        self.forward_button.clicked.connect(self.forward_plot_cb)
        controls_hbox.addWidget(self.forward_button)

        self.last_button = QPushButton(">>")
        self.last_button.clicked.connect(self.last_plot_cb)
        controls_hbox.addWidget(self.last_button)

        self.play_button = QPushButton("Play")
        self.play_button.setStyleSheet("background-color : lightgreen")
        # self.play_button.clicked.connect(self.play_plot_cb)
        self.play_button.clicked.connect(self.animate)
        controls_hbox.addWidget(self.play_button)

        # self.prepare_button = QPushButton("Prepare")
        # self.prepare_button.clicked.connect(self.prepare_plot_cb)
        # controls_hbox.addWidget(self.prepare_button)

        self.cells_checkbox = QCheckBox('Cells')
        self.cells_checkbox.setChecked(True)
        self.cells_checkbox.clicked.connect(self.cells_toggle_cb)
        self.cells_checked_flag = True

        self.substrates_checkbox = QCheckBox('Substrates')
        self.substrates_checkbox.setChecked(False)
        self.substrates_checkbox.clicked.connect(self.substrates_toggle_cb)
        self.substrates_checked_flag = False
        

        hbox = QHBoxLayout()
        hbox.addWidget(self.cells_checkbox)
        hbox.addWidget(self.substrates_checkbox)
        controls_hbox.addLayout(hbox)

        #-------------------
        controls_hbox2 = QHBoxLayout()
        visible_flag = True

        label = QLabel("xmin")
        label.setFixedWidth(label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        controls_hbox2.addWidget(label)
        self.my_xmin = QLineEdit()
        self.my_xmin.textChanged.connect(self.change_plot_range)
        self.my_xmin.setFixedWidth(domain_value_width)
        self.my_xmin.setValidator(QtGui.QDoubleValidator())
        controls_hbox2.addWidget(self.my_xmin)
        self.my_xmin.setVisible(visible_flag)

        label = QLabel("xmax")
        label.setFixedWidth(label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        controls_hbox2.addWidget(label)
        self.my_xmax = QLineEdit()
        self.my_xmax.textChanged.connect(self.change_plot_range)
        self.my_xmax.setFixedWidth(domain_value_width)
        self.my_xmax.setValidator(QtGui.QDoubleValidator())
        controls_hbox2.addWidget(self.my_xmax)
        self.my_xmax.setVisible(visible_flag)

        label = QLabel("ymin")
        label.setFixedWidth(label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        controls_hbox2.addWidget(label)
        self.my_ymin = QLineEdit()
        self.my_ymin.textChanged.connect(self.change_plot_range)
        self.my_ymin.setFixedWidth(domain_value_width)
        self.my_ymin.setValidator(QtGui.QDoubleValidator())
        controls_hbox2.addWidget(self.my_ymin)
        self.my_ymin.setVisible(visible_flag)

        label = QLabel("ymax")
        label.setFixedWidth(label_width)
        label.setAlignment(QtCore.Qt.AlignRight)
        controls_hbox2.addWidget(label)
        self.my_ymax = QLineEdit()
        self.my_ymax.textChanged.connect(self.change_plot_range)
        self.my_ymax.setFixedWidth(domain_value_width)
        self.my_ymax.setValidator(QtGui.QDoubleValidator())
        controls_hbox2.addWidget(self.my_ymax)
        self.my_ymax.setVisible(visible_flag)

        w = QPushButton("Reset")
        w.clicked.connect(self.reset_plot_range)
        controls_hbox2.addWidget(w)

        self.my_xmin.setText(str(self.xmin))
        self.my_xmax.setText(str(self.xmax))
        self.my_ymin.setText(str(self.ymin))
        self.my_ymax.setText(str(self.ymax))

        #-------------------
        # self.substrates_cbox = QComboBox(self)
        # self.substrates_cbox.setGeometry(200, 150, 120, 40)
  
        # self.substrates_cbox.addItem("substrate1")
        # self.substrates_cbox.addItem("substrate2")
        self.substrates_cbox.currentIndexChanged.connect(self.substrates_cbox_changed_cb)
        controls_hbox.addWidget(self.substrates_cbox)

        controls_vbox = QVBoxLayout()
        controls_vbox.addLayout(controls_hbox)
        controls_vbox.addLayout(controls_hbox2)

        #==================================================================
        self.config_params.setLayout(self.vbox)

        self.myscroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.myscroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.myscroll.setWidgetResizable(True)

        # self.myscroll.setWidget(self.config_params) # self.config_params = QWidget()
        self.myscroll.setWidget(self.canvas) # self.config_params = QWidget()
        self.layout = QVBoxLayout(self)
        # self.layout.addLayout(controls_hbox)
        # self.layout.addLayout(controls_hbox2)
        self.layout.addLayout(controls_vbox)
        self.layout.addWidget(self.myscroll)

        # self.create_figure()


    def reset_plot_range(self):
        try:  # due to the initial callback
            self.my_xmin.setText(str(self.xmin))
            self.my_xmax.setText(str(self.xmax))
            self.my_ymin.setText(str(self.ymin))
            self.my_ymax.setText(str(self.ymax))

            self.plot_xmin = float(self.xmin)
            self.plot_xmax = float(self.xmax)
            self.plot_ymin = float(self.ymin)
            self.plot_ymax = float(self.ymax)
        except:
            pass

        self.update_plots()


    def init_plot_range(self, config_tab):
        print("----- init_plot_range:")
        try:
            # beware of widget callback 
            self.my_xmin.setText(config_tab.xmin.text())
            self.my_xmax.setText(config_tab.xmax.text())
            self.my_ymin.setText(config_tab.ymin.text())
            self.my_ymax.setText(config_tab.ymax.text())
        except:
            pass

    def change_plot_range(self):
        print("----- change_plot_range:")
        # print("----- my_xmin= ",self.my_xmin.text())
        # print("----- my_xmax= ",self.my_xmax.text())
        try:  # due to the initial callback
            self.plot_xmin = float(self.my_xmin.text())
            self.plot_xmax = float(self.my_xmax.text())
            self.plot_ymin = float(self.my_ymin.text())
            self.plot_ymax = float(self.my_ymax.text())
        except:
            pass

        self.update_plots()

    def update_plots(self):
        # self.ax0.cla()
        # if self.substrates_checked_flag:
        #     self.plot_substrate(self.current_svg_frame)
        # if self.cells_checked_flag:
        #     self.plot_svg(self.current_svg_frame)

        # self.canvas.update()
        # self.canvas.draw()
        return

    def fill_substrates_combobox(self, substrate_list):
        print("vis_tab.py: ------- fill_substrates_combobox")
        print("substrate_list = ",substrate_list )
        self.substrates_cbox.clear()
        for s in substrate_list:
            # print(" --> ",s)
            self.substrates_cbox.addItem(s)
        # self.substrates_cbox.setCurrentIndex(2)  # not working; gets reset to oxygen somehow after a Run

    def substrates_cbox_changed_cb(self,idx):
        print("----- vis_tab.py: substrates_cbox_changed_cb: idx = ",idx)
        self.field_index = 4 + idx # substrate (0th -> 4 in the .mat)
        self.update_plots()


    def open_directory_cb(self):
        dialog = QFileDialog()
        # self.output_dir = dialog.getExistingDirectory(self, 'Select an output directory')
        tmp_dir = dialog.getExistingDirectory(self, 'Select an output directory')
        print("open_directory_cb:  tmp_dir=",tmp_dir)
        if tmp_dir == "":
            return

        self.output_dir = tmp_dir
        self.output_dir_w.setText(self.output_dir)
        self.reset_model()

    def reset_model(self):
        print("\n--------- vis_tab: reset_model ----------")
        # Verify initial.xml and at least one .svg file exist. Obtain bounds from initial.xml
        # tree = ET.parse(self.output_dir + "/" + "initial.xml")
        xml_file = Path(self.output_dir, "initial.xml")
        if not os.path.isfile(xml_file):
            print("vis_tab: Warning: Expecting initial.xml, but does not exist.")
            # msgBox = QMessageBox()
            # msgBox.setIcon(QMessageBox.Information)
            # msgBox.setText("Did not find 'initial.xml' in the output directory. Will plot a dummy substrate until you run a simulation.")
            # msgBox.setStandardButtons(QMessageBox.Ok)
            # msgBox.exec()
            return

        tree = ET.parse(Path(self.output_dir, "initial.xml"))
        xml_root = tree.getroot()

        bds_str = xml_root.find(".//microenvironment//domain//mesh//bounding_box").text
        bds = bds_str.split()
        print('bds=',bds)
        self.xmin = float(bds[0])
        self.xmax = float(bds[3])
        print('reset_model(): self.xmin, xmax=',self.xmin, self.xmax)
        self.x_range = self.xmax - self.xmin
        self.plot_xmin = self.xmin
        self.plot_xmax = self.xmax

        try:
            self.my_xmin.setText(str(self.plot_xmin))
            self.my_xmax.setText(str(self.plot_xmax))
            self.my_ymin.setText(str(self.plot_ymin))
            self.my_ymax.setText(str(self.plot_ymax))
        except:
            pass

        self.ymin = float(bds[1])
        self.ymax = float(bds[4])
        self.y_range = self.ymax - self.ymin
        self.plot_ymin = self.ymin
        self.plot_ymax = self.ymax

        xcoords_str = xml_root.find(".//microenvironment//domain//mesh//x_coordinates").text
        xcoords = xcoords_str.split()
        print('reset_model(): xcoords=',xcoords)
        print('reset_model(): len(xcoords)=',len(xcoords))
        self.numx =  len(xcoords)
        self.numy =  len(xcoords)
        print("reset_model(): self.numx, numy = ",self.numx,self.numy)

        #-------------------
        vars_uep = xml_root.find(".//microenvironment//domain//variables")
        if vars_uep:
            sub_names = []
            for var in vars_uep:
            # self.substrate.clear()
            # self.param[substrate_name] = {}  # a dict of dicts

            # self.tree.clear()
                idx = 0
            # <microenvironment_setup>
		    #   <variable name="food" units="dimensionless" ID="0">
                # print(cell_def.attrib['name'])
                if var.tag == 'variable':
                    substrate_name = var.attrib['name']
                    print("substrate: ",substrate_name )
                    sub_names.append(substrate_name)
                self.substrates_cbox.clear()
                print("sub_names = ",sub_names)
                self.substrates_cbox.addItems(sub_names)


        # and plot 1st frame (.svg)
        self.current_svg_frame = 0
        # self.forward_plot_cb("")  

    def reset_axes(self):
        print("--------- vis_tab: reset_axes ----------")
        # Verify initial.xml and at least one .svg file exist. Obtain bounds from initial.xml
        # tree = ET.parse(self.output_dir + "/" + "initial.xml")
        xml_file = Path(self.output_dir, "initial.xml")
        if not os.path.isfile(xml_file):
            print("Expecting initial.xml, but does not exist.")
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText("Did not find 'initial.xml' in this directory.")
            msgBox.setStandardButtons(QMessageBox.Ok)
            msgBox.exec()
            return

        tree = ET.parse(Path(self.output_dir, "initial.xml"))
        xml_root = tree.getroot()

        bds_str = xml_root.find(".//microenvironment//domain//mesh//bounding_box").text
        bds = bds_str.split()
        print('bds=',bds)
        self.xmin = float(bds[0])
        self.xmax = float(bds[3])
        self.x_range = self.xmax - self.xmin

        self.ymin = float(bds[1])
        self.ymax = float(bds[4])
        self.y_range = self.ymax - self.ymin

        # and plot 1st frame (.svg)
        self.current_svg_frame = 0
        # self.forward_plot_cb("")  


    # def output_dir_changed(self, text):
    #     self.output_dir = text
    #     print(self.output_dir)

    def first_plot_cb(self, text):
        if self.reset_model_flag:
            self.reset_model()
            self.reset_model_flag = False

        self.current_svg_frame = 0
        self.update_plots()

    def last_plot_cb(self, text):
        if self.reset_model_flag:
            self.reset_model()
            self.reset_model_flag = False

        print('cwd = ',os.getcwd())
        print('self.output_dir = ',self.output_dir)
        # xml_file = Path(self.output_dir, "initial.xml")
        # xml_files = glob.glob('tmpdir/output*.xml')
        xml_files = glob.glob('output*.xml')
        if len(xml_files) == 0:
            return
        xml_files.sort()
        svg_files = glob.glob('snapshot*.svg')
        svg_files.sort()
        print('xml_files = ',xml_files)
        num_xml = len(xml_files)
        print('svg_files = ',svg_files)
        num_svg = len(svg_files)
        print('num_xml, num_svg = ',num_xml, num_svg)
        last_xml = int(xml_files[-1][-12:-4])
        last_svg = int(svg_files[-1][-12:-4])
        print('last_xml, _svg = ',last_xml,last_svg)
        self.current_svg_frame = last_xml
        if last_svg < last_xml:
            self.current_svg_frame = last_svg
        self.update_plots()

    def back_plot_cb(self, text):
        if self.reset_model_flag:
            self.reset_model()
            self.reset_model_flag = False

        self.current_svg_frame -= 1
        if self.current_svg_frame < 0:
            self.current_svg_frame = 0
        # print('svg # ',self.current_svg_frame)

        self.update_plots()


    def forward_plot_cb(self, text):
        if self.reset_model_flag:
            self.reset_model()
            self.reset_model_flag = False

        self.current_svg_frame += 1
        # print('svg # ',self.current_svg_frame)

        self.update_plots()


    # def task(self):
            # self.dc.update_figure()

    # used by animate
    def play_plot_cb(self):
        for idx in range(1):
            self.current_svg_frame += 1
            # print('svg # ',self.current_svg_frame)

            fname = "snapshot%08d.svg" % self.current_svg_frame
            full_fname = os.path.join(self.output_dir, fname)
            # print("full_fname = ",full_fname)
            # with debug_view:
                # print("plot_svg:", full_fname) 
            # print("-- plot_svg:", full_fname) 
            if not os.path.isfile(full_fname):
                # print("Once output files are generated, click the slider.")   
                print("play_plot_cb():  Reached the end (or no output files found).")
                # self.timer.stop()
                # self.current_svg_frame -= 1
                self.animating_flag = True
                self.current_svg_frame = 0
                self.animate()
                return

            self.update_plots()


    def cells_toggle_cb(self,bval):
        self.cells_checked_flag = bval

        self.update_plots()


    def substrates_toggle_cb(self,bval):
        self.substrates_checked_flag = bval

        self.update_plots()


    def animate(self):
        if not self.animating_flag:
            self.animating_flag = True
            self.play_button.setText("Halt")
            self.play_button.setStyleSheet("background-color : red")

            if self.reset_model_flag:
                self.reset_model()
                self.reset_model_flag = False

            # self.current_svg_frame = 0
            self.timer.start(1)

        else:
            self.animating_flag = False
            self.play_button.setText("Play")
            self.play_button.setStyleSheet("background-color : lightgreen")
            self.timer.stop()


    # def play_plot_cb0(self, text):
    #     for idx in range(10):
    #         self.current_svg_frame += 1
    #         print('svg # ',self.current_svg_frame)
    #         self.plot_svg(self.current_svg_frame)
    #         self.canvas.update()
    #         self.canvas.draw()
    #         # time.sleep(1)
    #         # self.ax0.clear()
    #         # self.canvas.pause(0.05)

    def prepare_plot_cb(self, text):
        self.current_svg_frame += 1
        print('\n\n   ====>     prepare_plot_cb(): svg # ',self.current_svg_frame)

        self.update_plots()


    def create_figure(self):
        print("\n--------- create_figure(): ------- creating figure, canvas, ax0")
        self.canvas = QWidget()
        self.vl = QVBoxLayout(self.canvas)
        # self.setCentralWidget(self.canvas)
        # self.resize(640, 480)

        self.ren = vtk.vtkRenderer()
        # vtk_widget = QVTKRenderWindowInteractor(rw=render_window)
        # self.vtkWidget = QVTKRenderWindowInteractor(self.ren)
        self.vtkWidget = QVTKRenderWindowInteractor(self.canvas)
        self.vl.addWidget(self.vtkWidget)
        # self.vtkWidget.Initialize()
        # self.vtkWidget.Start()

        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

        # Create source
        source = vtk.vtkSphereSource()
        source.SetCenter(0, 0, 0)
        source.SetRadius(5.0)

        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        self.ren.AddActor(actor)
        self.ren.ResetCamera()
        # self.frame.setLayout(self.vl)
        # self.setCentralWidget(self.frame)
        self.show()
        self.iren.Initialize()
        self.iren.Start()

        # self.figure = plt.figure()
        # self.canvas = FigureCanvasQTAgg(self.figure)
        # self.canvas.setStyleSheet("background-color:transparent;")

        # Adding one subplot for image
        # self.ax0 = self.figure.add_subplot(111)
        # self.ax0 = self.figure.add_subplot(111, adjustable='box', aspect=1.2)
        # self.ax0 = self.figure.add_subplot(111, adjustable='box', aspect=self.aspect_ratio)
        # self.ax0 = self.figure.add_subplot(111, adjustable='box')
        
        # self.ax0.get_xaxis().set_visible(False)
        # self.ax0.get_yaxis().set_visible(False)
        # plt.tight_layout()

        self.reset_model()

        # np.random.seed(19680801)  # for reproducibility
        # N = 50
        # x = np.random.rand(N) * 2000
        # y = np.random.rand(N) * 2000
        # colors = np.random.rand(N)
        # area = (30 * np.random.rand(N))**2  # 0 to 15 point radii
        # self.ax0.scatter(x, y, s=area, c=colors, alpha=0.5)

        # if self.plot_svg_flag:
        # if False:
        #     self.plot_svg(self.current_svg_frame)
        # else:
        #     self.plot_substrate(self.current_svg_frame)

        # print("create_figure(): ------- creating dummy contourf")
        # xlist = np.linspace(-3.0, 3.0, 50)
        # print("len(xlist)=",len(xlist))
        # ylist = np.linspace(-3.0, 3.0, 50)
        # X, Y = np.meshgrid(xlist, ylist)
        # Z = np.sqrt(X**2 + Y**2) + 10*np.random.rand()

        # self.cmap = plt.cm.get_cmap("viridis")
        # self.mysubstrate = self.ax0.contourf(X, Y, Z, cmap=self.cmap)
        # # if self.field_index > 4:
        # #     # plt.contour(xgrid, ygrid, M[self.field_index, :].reshape(self.numy,self.numx), [0.0])
        # #     plt.contour(X, Y, Z, [0.0])

        # self.cbar = self.figure.colorbar(self.mysubstrate, ax=self.ax0)
        # self.cbar.ax.tick_params(labelsize=self.fontsize)

        # # substrate_plot = self.ax0.contourf(xgrid, ygrid, M[self.field_index, :].reshape(self.numy,self.numx), num_contours, cmap='viridis')  # self.colormap_dd.value)

        # print("------------create_figure():  # axes = ",len(self.figure.axes))

        # # self.imageInit = [[255] * 320 for i in range(240)]
        # # self.imageInit[0][0] = 0

        # # Init image and add colorbar
        # # self.image = self.ax0.imshow(self.imageInit, interpolation='none')
        # # divider = make_axes_locatable(self.ax0)
        # # cax = divider.new_vertical(size="5%", pad=0.05, pack_start=True)
        # # self.colorbar = self.figure.add_axes(cax)
        # # self.figure.colorbar(self.image, cax=cax, orientation='horizontal')

        # # plt.subplots_adjust(left=0, bottom=0.05, right=1, top=1, wspace=0, hspace=0)

        # # self.plot_substrate(self.current_svg_frame)
        # # self.plot_svg(self.current_svg_frame)
        self.plot_cells3D(self.current_svg_frame)
        # # self.canvas.draw()

    #------------------------------------------------------------
    # def plot_svg(self, frame, rdel=''):
    def plot_cells3D(self, frame):
        return

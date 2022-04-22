import sys
import os
import time
import random
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
# import vtk
from vtk import *
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
# from collections import deque
import glob

from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QFrame,QApplication,QWidget,QTabWidget,QFormLayout,QLineEdit, QHBoxLayout,QVBoxLayout,QRadioButton,QLabel,QCheckBox,QComboBox,QScrollArea,  QMainWindow,QGridLayout, QPushButton, QFileDialog, QMessageBox

import numpy as np
import scipy.io
from pyMCDS_cells import pyMCDS_cells
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
        self.output_dir = "../tmpdir"   # for nanoHUB
        self.output_dir = "tmpdir"   # for nanoHUB
        # self.output_dir = "."   # for nanoHUB


        # do in create_figure()?
        # xlist = np.linspace(-3.0, 3.0, 50)
        # print("len(xlist)=",len(xlist))
        # ylist = np.linspace(-3.0, 3.0, 50)
        # X, Y = np.meshgrid(xlist, ylist)
        # Z = np.sqrt(X**2 + Y**2) + 10*np.random.rand()
        # self.cmap = plt.cm.get_cmap("viridis")
        # self.cs = self.ax0.contourf(X, Y, Z, cmap=self.cmap)
        # self.cbar = self.figure.colorbar(self.cs, ax=self.ax)

        # rendering objects
        self.points = vtkPoints()
        self.cellID = vtkFloatArray()
        self.cellVolume = vtkFloatArray()
        self.polydata = vtkPolyData()
        self.colors = vtkUnsignedCharArray()
        self.sphereSource = vtkSphereSource()
        nres = 20
        self.sphereSource.SetPhiResolution(nres)
        self.sphereSource.SetThetaResolution(nres)
        # self.sphereSource.SetRadius(0.02)
        # self.sphereSource.SetRadius(2.6)  # rwh ??
        self.sphereSource.SetRadius(8.412)  # rwh ??
        self.glyph = vtkGlyph3D()
        self.cells_mapper = vtkPolyDataMapper()
        self.cells_mapper.SetInputConnection(self.glyph.GetOutputPort())

        self.cells_actor = vtkActor()
        self.cells_actor.SetMapper(self.cells_mapper)

        xmax = 400
        ymax = 400
        zmax = 150
        self.cyl = vtkCylinderSource()
        self.cyl.CappingOff()
        self.cyl.SetCenter(0.0, 0.0, 300.0)
        self.cyl.SetRadius(xmax)
        self.cyl.SetHeight(800.0)
        self.cyl.SetResolution(50)
        self.cyl.Update()

        # define x-z plane
        self.plane_source = vtkPlaneSource()
        xorig = -100
        zBM = -20
        self.plane_source.SetOrigin(xorig,0,zBM)
        plength = 300.0
        self.plane_source.SetPoint1(xorig+plength,0,zBM)
        self.plane_source.SetPoint2(xorig,plength,zBM)
        self.plane_source.SetResolution(20,20)
        self.plane_source.Update()
        self.plane_mapper = vtkDataSetMapper()
        self.plane_mapper.SetInputData(self.plane_source.GetOutput())
        self.plane_actor = vtkActor()
        self.plane_actor.SetMapper(self.plane_mapper)
        self.plane_actor.GetProperty().SetColor(1.0, 0.5, 0.5)
        self.plane_actor.GetProperty().SetDiffuse(0.5)
        self.plane_actor.GetProperty().SetAmbient(0.4)
        self.plane_actor.GetProperty().SetRepresentationToWireframe()

        # clip the cylinder
        # https://kitware.github.io/vtk-examples/site/Python/Meshes/CapClip/
        self.plane = vtkPlane()
        print("cyl center= ",self.cyl.GetOutput().GetCenter())
        self.plane.SetOrigin(self.cyl.GetOutput().GetCenter())
        self.plane.SetOrigin(0,0,300)
        self.plane.SetNormal(0.0, 0.0, -1.0)

        self.clipper = vtkClipPolyData()
        self.clipper.SetInputData(self.cyl.GetOutput())
        self.clipper.SetClipFunction(self.plane)
        self.clipper.SetValue(100)
        self.clipper.Update()

        self.polyData = self.clipper.GetOutput()

        self.clipMapper = vtkDataSetMapper()
        self.clipMapper.SetInputData(self.polyData)

        # cyl_mapper = vtkPolyDataMapper()
        # cyl_mapper.SetInputConnection(cyl.GetOutputPort())

        self.cyl_actor = vtkActor()
        # self.cyl_actor.SetMapper(cyl_mapper)
        self.cyl_actor.SetMapper(self.clipMapper)
        self.cyl_actor.GetProperty().SetColor(1.0, 0.5, 0.5)
        self.cyl_actor.GetProperty().SetDiffuse(0.5)
        self.cyl_actor.GetProperty().SetAmbient(0.4)
        self.cyl_actor.GetProperty().SetRepresentationToWireframe()

        #--------
        self.box_outline = vtkOutlineSource()
        bds = [-xmax,xmax, -ymax,ymax, -zmax,zmax]    # {xmin,xmax,ymin,ymax,zmin,zmax} via SetBounds()
        self.box_outline.SetBounds(bds)

        self.box_mapper = vtkPolyDataMapper()
        self.box_mapper.SetInputConnection(self.box_outline.GetOutputPort())

        self.box_actor = vtkActor()
        self.box_actor.SetMapper(self.box_mapper)
        self.box_actor.GetProperty().SetColor(1.0, 1.0, 1.0)



        #-------------------------------------------
        label_width = 110
        domain_value_width = 100
        value_width = 60
        label_height = 20
        units_width = 70

        self.substrates_cbox = QComboBox(self)
        self.substrates_cbox.setEnabled(False)

        self.myscroll = QScrollArea()  # might contain centralWidget
        self.create_figure(False)

        # self.ren.AddActor(self.cyl_actor)
        # self.ren.AddActor(self.plane_actor)
        self.show_domain_box = True
        self.ren.AddActor(self.box_actor)
        self.ren.AddActor(self.cells_actor)


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
        self.substrates_checkbox.setEnabled(False)
        self.substrates_checkbox.setChecked(False)
        self.substrates_checkbox.clicked.connect(self.substrates_toggle_cb)
        self.substrates_checked_flag = False
        

        hbox = QHBoxLayout()
        hbox.addWidget(self.cells_checkbox)
        hbox.addWidget(self.substrates_checkbox)
        controls_hbox.addLayout(hbox)

        #-------------------
        # controls_hbox2 = QHBoxLayout()
        # visible_flag = True

        # label = QLabel("xmin")
        # label.setFixedWidth(label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # controls_hbox2.addWidget(label)
        # self.my_xmin = QLineEdit()
        # self.my_xmin.textChanged.connect(self.change_plot_range)
        # self.my_xmin.setFixedWidth(domain_value_width)
        # self.my_xmin.setValidator(QtGui.QDoubleValidator())
        # controls_hbox2.addWidget(self.my_xmin)
        # self.my_xmin.setVisible(visible_flag)

        # label = QLabel("xmax")
        # label.setFixedWidth(label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # controls_hbox2.addWidget(label)
        # self.my_xmax = QLineEdit()
        # self.my_xmax.textChanged.connect(self.change_plot_range)
        # self.my_xmax.setFixedWidth(domain_value_width)
        # self.my_xmax.setValidator(QtGui.QDoubleValidator())
        # controls_hbox2.addWidget(self.my_xmax)
        # self.my_xmax.setVisible(visible_flag)

        # label = QLabel("ymin")
        # label.setFixedWidth(label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # controls_hbox2.addWidget(label)
        # self.my_ymin = QLineEdit()
        # self.my_ymin.textChanged.connect(self.change_plot_range)
        # self.my_ymin.setFixedWidth(domain_value_width)
        # self.my_ymin.setValidator(QtGui.QDoubleValidator())
        # controls_hbox2.addWidget(self.my_ymin)
        # self.my_ymin.setVisible(visible_flag)

        # label = QLabel("ymax")
        # label.setFixedWidth(label_width)
        # label.setAlignment(QtCore.Qt.AlignRight)
        # controls_hbox2.addWidget(label)
        # self.my_ymax = QLineEdit()
        # self.my_ymax.textChanged.connect(self.change_plot_range)
        # self.my_ymax.setFixedWidth(domain_value_width)
        # self.my_ymax.setValidator(QtGui.QDoubleValidator())
        # controls_hbox2.addWidget(self.my_ymax)
        # self.my_ymax.setVisible(visible_flag)

        # w = QPushButton("Reset")
        # w.clicked.connect(self.reset_plot_range)
        # controls_hbox2.addWidget(w)

        # self.my_xmin.setText(str(self.xmin))
        # self.my_xmax.setText(str(self.xmax))
        # self.my_ymin.setText(str(self.ymin))
        # self.my_ymax.setText(str(self.ymax))

        #-------------------
        # self.substrates_cbox = QComboBox(self)
        # self.substrates_cbox.setGeometry(200, 150, 120, 40)
  
        # self.substrates_cbox.addItem("substrate1")
        # self.substrates_cbox.addItem("substrate2")
        self.substrates_cbox.currentIndexChanged.connect(self.substrates_cbox_changed_cb)
        controls_hbox.addWidget(self.substrates_cbox)

        controls_vbox = QVBoxLayout()
        controls_vbox.addLayout(controls_hbox)
        # controls_vbox.addLayout(controls_hbox2)

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


    def toggle_domain_box(self):
        self.show_domain_box = not self.show_domain_box
        if self.show_domain_box:
            # self.ren.AddActor(self.box_actor)
            self.box_actor.VisibilityOn()
        else:
            # self.ren.RemoveActor(self.box_actor)
            self.box_actor.VisibilityOff()

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
        self.plot_cells3D(self.current_svg_frame)
        # self.plot_cells3D(0)

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

    def reset_model2(self):
        return
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

        # svg_files = glob.glob('snapshot*.svg')
        # svg_files.sort()
        mat_files = glob.glob('output*_cells_physicell.mat')
        mat_files.sort()
        # print("mat_files: ",mat_files)

        # print('xml_files = ',xml_files)
        num_xml = len(xml_files)

        # print('svg_files = ',svg_files)
        # num_svg = len(svg_files)
        # print('num_xml, num_svg = ',num_xml, num_svg)

        last_xml = int(xml_files[-1][-12:-4])
        # last_svg = int(svg_files[-1][-12:-4])
        last_mat = int(mat_files[-1][6:14])  # output00000012_cells_physicell.mat

        # print('last_xml, _svg = ',last_xml,last_svg)
        print('last_xml, _mat = ',last_xml,last_mat)

        self.current_svg_frame = last_xml
        # if last_svg < last_xml:
            # self.current_svg_frame = last_svg
        if last_mat < last_xml:
            self.current_svg_frame = last_mat

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

            # fname = "snapshot%08d.svg" % self.current_svg_frame
            fname = "output%08d.xml" % self.current_svg_frame
            # full_fname = os.path.join(self.output_dir, fname)
            full_fname = fname
            print("play_plot_cb(): full_fname = ",full_fname)
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


    def create_figure(self, plot_flag):
        print("\n--------- create_figure(): ------- creating figure, canvas, ax0")
        self.canvas = QWidget()
        self.vl = QVBoxLayout(self.canvas)
        # self.setCentralWidget(self.canvas)
        # self.resize(640, 480)

        self.ren = vtkRenderer()
        # vtk_widget = QVTKRenderWindowInteractor(rw=render_window)
        # self.vtkWidget = QVTKRenderWindowInteractor(self.ren)
        self.vtkWidget = QVTKRenderWindowInteractor(self.canvas)
        self.vl.addWidget(self.vtkWidget)
        # self.vtkWidget.Initialize()
        # self.vtkWidget.Start()

        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

        # Create source
        # # # source = vtkSphereSource()
        # # # source.SetCenter(0, 0, 0)
        # # # source.SetRadius(5.0)

        # # # # Create a mapper
        # # # mapper = vtkPolyDataMapper()
        # # # mapper.SetInputConnection(source.GetOutputPort())

        # # # # Create an actor
        # # # actor = vtkActor()
        # # # actor.SetMapper(mapper)
        # # self.ren.AddActor(actor)
        # self.ren.ResetCamera()
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
        # self.plot_cells3D(self.current_svg_frame)
        if plot_flag:
            # self.plot_cells3D(0)
            self.plot_cells3D(self.current_svg_frame)
        # # self.canvas.draw()

    #------------------------------------------------------------
    # def plot_svg(self, frame, rdel=''):
    def plot_cells3D(self, frame):
        print("plot_cells3D:  self.output_dir= ",self.output_dir)
        print("plot_cells3D:  frame= ",frame)
        # xml_file = Path(self.output_dir, "output00000000.xml")
        # xml_file = "output00000000.xml"
        xml_file = "output%08d.xml" % frame
        print("plot_cells3D: xml_file = ",xml_file)

        if not os.path.isfile(xml_file):
            print("plot_cells3D(): file not found, return. ", xml_file)
            return

        # self.iren.ReInitialize()
        # self.iren.GetRenderWindow().Render()

        mcds = pyMCDS_cells(xml_file, '.')  
        # mcds = pyMCDS_cells(xml_file, 'tmpdir')  
        print('time=', mcds.get_time())

        print(mcds.data['discrete_cells'].keys())

        ncells = len(mcds.data['discrete_cells']['ID'])
        print('total_volume= ',mcds.data['discrete_cells']['total_volume'])
        print('ncells=', ncells)

        # global xyz
        xyz = np.zeros((ncells, 3))
        xyz[:, 0] = mcds.data['discrete_cells']['position_x']
        xyz[:, 1] = mcds.data['discrete_cells']['position_y']
        xyz[:, 2] = mcds.data['discrete_cells']['position_z']
        #xyz = xyz[:1000]
        # print("position_x = ",xyz[:,0])
        xmin = min(xyz[:,0])
        xmax = max(xyz[:,0])
        print("xmin = ",xmin)
        print("xmax = ",xmax)

        ymin = min(xyz[:,1])
        ymax = max(xyz[:,1])
        print("ymin = ",ymin)
        print("ymax = ",ymax)

        zmin = min(xyz[:,2])
        zmax = max(xyz[:,2])
        print("zmin = ",zmin)
        print("zmax = ",zmax)

        # cell_type = mcds.data['discrete_cells']['cell_type']
        # cell_custom_ID = mcds.data['discrete_cells']['cell_ID']
        # # print(type(cell_type))
        # # print(cell_type)
        # unique_cell_type = np.unique(cell_type)
        # unique_cell_custom_ID = np.unique(cell_custom_ID)
        # print("\nunique_cell_type = ",unique_cell_type )
        # print("\nunique_cell_custom_ID = ",unique_cell_custom_ID )

        #------------
        # colors = vtkNamedColors()

        self.points.Reset()
        self.cellID.Reset()
        self.cellVolume.Reset()
        for idx in range(ncells):
            x= mcds.data['discrete_cells']['position_x'][idx]
            y= mcds.data['discrete_cells']['position_y'][idx]
            z= mcds.data['discrete_cells']['position_z'][idx]
            # id = mcds.data['discrete_cells']['cell_type'][idx]
            # id = mcds.data['discrete_cells']['cell_ID'][idx]
            self.points.InsertNextPoint(x, y, z)
            self.cellVolume.InsertNextValue(2.2)  # actually want radius?
            # self.cellID.InsertNextValue(id)

        # print("min (parent)cell_ID = ",min(mcds.data['discrete_cells']['cell_ID']))
        # print("max (parent)cell_ID = ",max(mcds.data['discrete_cells']['cell_ID']))
        # cell_flavor = mcds.data['discrete_cells']['cell_ID']
        # unique_cell_flavor = np.unique(cell_flavor)
        # print("\nunique_cell_flavor = ",unique_cell_flavor )

        # self.polydata = vtkPolyData()
        self.polydata.SetPoints(self.points)
        self.polydata.GetPointData().SetScalars(self.cellVolume)
        # self.polydata.GetPointData().SetScalars(self.cellID)

        cellID_color_dict = {}
        # for idx in range(ncells):
        random.seed(42)
        # for utype in unique_cell_type:
        # for utype in unique_cell_custom_ID:
        #     # colors.InsertTuple3(0, randint(0,255), randint(0,255), randint(0,255)) # reddish
        #     cellID_color_dict[utype] = [random.randint(0,255), random.randint(0,255), random.randint(0,255)]
        # cellID_color_dict[0.0]=[255,255,0]  # yellow basement membrane
        # cellID_color_dict[1.]=[255,255,0]  # yellow basement membrane
        cellID_color_dict[0.]=[255,255,0]  # yellow basement membrane
        print("color dict=",cellID_color_dict)

        # self.colors = vtkUnsignedCharArray()
        # self.colors.Reset()
        # self.colors.SetNumberOfComponents(3)
        # self.colors.SetNumberOfTuples(self.polydata.GetNumberOfPoints())  # ncells

        # for idx in range(ncells):
        # for idx in range(len(unique_cell_type)):
            # colors.InsertTuple3(idx, randint(0,255), randint(0,255), randint(0,255)) 
            # if idx < 5:
                # print(idx,cellID_color_dict[cell_type[idx]])
            # self.colors.InsertTuple3(idx, cellID_color_dict[cell_type[idx]][0], cellID_color_dict[cell_type[idx]][1], cellID_color_dict[cell_type[idx]][2])

            # self.colors.InsertTuple3(idx, cellID_color_dict[cell_custom_ID[idx]][0], cellID_color_dict[cell_custom_ID[idx]][1], cellID_color_dict[cell_custom_ID[idx]][2])

        # self.polydata.GetPointData().SetScalars(self.colors)

        # self.sphereSource = vtkSphereSource()
        # nres = 20
        # self.sphereSource.SetPhiResolution(nres)
        # self.sphereSource.SetThetaResolution(nres)
        # self.sphereSource.SetRadius(20.0)

        # self.glyph = vtkGlyph3D()
        self.glyph.SetSourceConnection(self.sphereSource.GetOutputPort())
        self.glyph.SetInputData(self.polydata)
        # self.glyph.SetColorModeToColorByScalar()
        self.glyph.SetScaleModeToScaleByScalar()

        # using these 2 results in fixed size spheres
        # self.glyph.SetScaleModeToDataScalingOff()  # results in super tiny spheres without 'ScaleFactor'
        # self.glyph.SetScaleFactor(170)  # overall (multiplicative) scaling factor

        # glyph.SetScaleModeToDataScalingOn()
        # glyph.ScalingOn()
        self.glyph.Update()

        # Visualize
        # self.cells_mapper = vtkPolyDataMapper()
        # self.cells_mapper.SetInputConnection(self.glyph.GetOutputPort())
        self.cells_mapper.Update()

        # self.cells_actor = vtkActor()
        # self.cells_actor.SetMapper(self.cells_mapper)

        # actor.GetProperty().SetInterpolationToPBR()
        # actor.GetProperty().SetColor(colors.GetColor3d('Salmon'))
        # print("-- actor defaults:")
        # print("-- diffuse:",self.actor.GetProperty().GetDiffuse())  # 1.0
        # print("-- specular:",self.actor.GetProperty().GetSpecular())  # 0.0
        # print("-- roughness:",self.actor.GetProperty().GetCoatRoughness ())  # 0.0
        # self.ren.AddActor(self.cyl_actor)
        # self.cells_actor.GetProperty().SetAmbient(0.3)
        # self.cells_actor.GetProperty().SetDiffuse(0.5)
        # self.cells_actor.GetProperty().SetSpecular(0.2)
        # actor.GetProperty().SetCoatRoughness (0.5)
        # actor.GetProperty().SetCoatRoughness (0.2)
        # actor.GetProperty().SetCoatRoughness (1.0)

        # renderer = vtkRenderer()
        # amval = 1.0  # default
        # renderer.SetAmbient(amval, amval, amval)

        # renderWindow = vtkRenderWindow()
        # renderWindow.SetPosition(100,100)
        # renderWindow.SetSize(1400,1200)
        # renderWindow.AddRenderer(renderer)
        # renderWindowInteractor = vtkRenderWindowInteractor()
        # renderWindowInteractor.SetRenderWindow(renderWindow)

        # renderer.AddActor(actor)
        #------------------------------------------
        # renderer.SetBackground(colors.GetColor3d('SlateGray'))  # Background Slate Gray

        # self.iren.ReInitialize()
        self.iren.GetRenderWindow().Render()

    # renderWindow.SetWindowName('PhysiCell model')
    # renderWindow.Render()
    # renderWindowInteractor.Start()
        self.ren.GetActiveCamera().ParallelProjectionOn()
        # self.ren.ResetCamera()

        return

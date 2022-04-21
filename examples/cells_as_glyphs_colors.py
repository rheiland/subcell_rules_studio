# https://kitware.github.io/vtk-examples/site/Python/Filtering/Glyph3D/
#!/usr/bin/env python

# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2

from vtk import *
# from vtkmodules.vtkCommonColor import vtkNamedColors
# from vtkmodules.vtkCommonCore import vtkPoints, vtkFloatArray
# from vtkmodules.vtkCommonDataModel import vtkPolyData,vtkUnsignedCharArray
# from vtkmodules.vtkFiltersCore import vtkGlyph3D
# from vtkmodules.vtkFiltersSources import vtkSphereSource
# from vtkmodules.vtkRenderingCore import (
#     vtkActor,
#     vtkPolyDataMapper,
#     vtkRenderWindow,
#     vtkRenderWindowInteractor,
#     vtkRenderer
# )

#from pyMCDS_cells import pyMCDS_cells
#from vtk.util import numpy_support
import numpy as np
from numpy import genfromtxt
import random


def main():
    csv_file = sys.argv[1]
    pts = genfromtxt(csv_file, delimiter=',')
    print(pts.shape)
    ncells = len(pts[:,0])
    print("# pts = ",ncells)

    #------------
    colors = vtkNamedColors()

    points = vtkPoints()
    # points.InsertNextPoint(0, 0, 0)
    # points.InsertNextPoint(1, 1, 1)
    # points.InsertNextPoint(2, 2, 2)
    cellID = vtkFloatArray()
    cellVolume = vtkFloatArray()
    for idx in range(ncells):
        x= pts[idx,0]
        y= pts[idx,1]
        z= pts[idx,2]
        id = pts[idx,3]
        points.InsertNextPoint(x, y, z)
        # cellVolume.InsertNextValue(30.0)
        cellID.InsertNextValue(id)

    polydata = vtkPolyData()
    polydata.SetPoints(points)
    # polydata.GetPointData().SetScalars(cellVolume)
    polydata.GetPointData().SetScalars(cellID)

    unique_cell_type = np.unique(pts[:,3])
    print("\nunique_cell_type = ",unique_cell_type )

    cellID_color_dict = {}
    # for idx in range(ncells):
    random.seed(42)
    for utype in unique_cell_type:
        # colors.InsertTuple3(0, randint(0,255), randint(0,255), randint(0,255)) # reddish
        cellID_color_dict[utype] = [random.randint(0,255), random.randint(0,255), random.randint(0,255)]
    cellID_color_dict[0.0]=[255,255,0]  # yellow basement membrane
    print("color dict=",cellID_color_dict)

    colors = vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.SetNumberOfTuples(polydata.GetNumberOfPoints())  # ncells
    for idx in range(ncells):
    # for idx in range(len(unique_cell_type)):
        # colors.InsertTuple3(idx, randint(0,255), randint(0,255), randint(0,255)) 
        # if idx < 5:
            # print(idx,cellID_color_dict[cell_type[idx]])
        # colors.InsertTuple3(idx, cellID_color_dict[cell_type[idx]][0], cellID_color_dict[cell_type[idx]][1], cellID_color_dict[cell_type[idx]][2])
        colors.InsertTuple3(idx, cellID_color_dict[pts[idx,3]][0], cellID_color_dict[pts[idx,3]][1], cellID_color_dict[pts[idx,3]][2]) 
        # colors.InsertTuple3(idx, 255, 0, 0)

    polydata.GetPointData().SetScalars(colors)

    sphereSource = vtkSphereSource()
    nres = 20
    sphereSource.SetPhiResolution(nres)
    sphereSource.SetThetaResolution(nres)
    sphereSource.SetRadius(0.1)

    glyph = vtkGlyph3D()
    glyph.SetSourceConnection(sphereSource.GetOutputPort())
    glyph.SetInputData(polydata)
    glyph.SetColorModeToColorByScalar()
    # glyph.SetScaleModeToScaleByScalar()

    # using these 2 results in fixed size spheres
    glyph.SetScaleModeToDataScalingOff()  # results in super tiny spheres without 'ScaleFactor'
    glyph.SetScaleFactor(170)  # overall (multiplicative) scaling factor

    # glyph.SetScaleModeToDataScalingOn()
    # glyph.ScalingOn()
    glyph.Update()


    # Visualize
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())

    actor = vtkActor()
    actor.SetMapper(mapper)
#    actor.GetProperty().SetInterpolationToPBR()
    # actor.GetProperty().SetColor(colors.GetColor3d('Salmon'))
    print("-- actor defaults:")
    print("-- ambient:",actor.GetProperty().GetAmbient())  # 1.0
    print("-- diffuse:",actor.GetProperty().GetDiffuse())  # 1.0
    print("-- specular:",actor.GetProperty().GetSpecular())  # 0.0
    print("-- roughness:",actor.GetProperty().GetCoatRoughness ())  # 0.0
    # actor.GetProperty().SetSpecular(0.2)
#    actor.GetProperty().SetCoatRoughness (0.5)
#    actor.GetProperty().SetCoatRoughness (0.2)
#    actor.GetProperty().SetCoatRoughness (1.0)

    actor.GetProperty().SetDiffuse(0.8)
    actor.GetProperty().SetSpecular(0.05)
#    actor.GetProperty().SetAmbient(0.2)  # 1.0

    renderer = vtkRenderer()
    amval = 0.9  # default
#    renderer.SetAmbient(amval, amval, amval)

    renderWindow = vtkRenderWindow()
    renderWindow.SetPosition(100,100)
    renderWindow.SetSize(1400,1200)
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(actor)
    # renderer.SetBackground(colors.GetColor3d('SlateGray'))  # Background Slate Gray

    renderWindow.SetWindowName('PhysiCell model')
    renderWindow.Render()
    renderWindowInteractor.Start()


if __name__ == '__main__':
    main()

#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'FractalMammaryLobule\\tests'))
	print(os.getcwd())
except:
	pass

#%%
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
 
# create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
 
# create source
source = vtk.vtkCylinderSource()
source.SetCenter(0,0,0)
source.SetRadius(5.0)
source.SetHeight(7.0)
source.SetResolution(100.0)

apd=vtk.vtk

source2 = vtk.vtkCylinderSource()
source2.SetCenter(20,0,0)
source2.SetRadius(2.0)
source2.SetHeight(10.0)
source2.SetResolution(100.0)

# mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInput(source.GetOutput())

mapper2 = vtk.vtkPolyDataMapper()
mapper2.SetInput(source2.GetOutput())
 
# actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)

actor2= vtk.vtkActor()
actor2.SetMapper(mapper2)
 
    
# assign actor to the renderer
ren.AddActor(actor)
ren.AddActor(actor2)

plane=vtk.vtkPlane()
plane.SetOrigin(0,0,5)
plane.SetNormal(0,0,1)
 
#create cutter
cutter=vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(source.GetOutputPort())
cutter.Update()
cutterMapper=vtk.vtkPolyDataMapper()
cutterMapper.SetInputConnection( cutter.GetOutputPort())
 

# enable user interface interactor
iren.Initialize()
renWin.Render()
iren.Start()


#%%
#!/usr/bin/env python

# In this example vtkClipPolyData is used to cut a polygonal model
# of a cow in half. In addition, the open clip is closed by triangulating
# the resulting complex polygons.

import vtk
from vtk.util.misc import vtkGetDataRoot
from vtk.util.colors import peacock, tomato
VTK_DATA_ROOT = vtkGetDataRoot()

# First start by reading a cow model. We also generate surface normals for
# prettier rendering.
source = vtk.vtkCylinderSource()
source.SetCenter(0,0,0)
source.SetRadius(5.0)
source.SetHeight(7.0)
source.SetResolution(100.0)

source2 = vtk.vtkCylinderSource()
source2.SetCenter(20,0,0)
source2.SetRadius(2.0)
source2.SetHeight(10.0)
source2.SetResolution(100.0)
    
cowNormals = vtk.vtkPolyDataNormals()
cowNormals.SetInputConnection(source.GetOutputPort())


cowNormals2 = vtk.vtkPolyDataNormals()
cowNormals2.SetInputConnection(source2.GetOutputPort())

# We clip with an implicit function. Here we use a plane positioned near
# the center of the cow model and oriented at an arbitrary angle.
plane = vtk.vtkPlane()
plane.SetOrigin(0, 3, 0)
plane.SetNormal(0, 1, 0)

# Here we are cutting the cow. Cutting creates lines where the cut
# function intersects the model. (Clipping removes a portion of the
# model but the dimension of the data does not change.)
#
# The reason we are cutting is to generate a closed polygon at the
# boundary of the clipping process. The cutter generates line
# segments, the stripper then puts them together into polylines. We
# then pull a trick and define polygons using the closed line
# segements that the stripper created.
cutEdges = vtk.vtkCutter()
cutEdges.SetInputConnection(cowNormals.GetOutputPort())
cutEdges.SetCutFunction(plane)
cutEdges.GenerateCutScalarsOn()
cutEdges.SetValue(0, 0.5)
cutStrips = vtk.vtkStripper()
cutStrips.SetInputConnection(cutEdges.GetOutputPort())
cutStrips.Update()
cutPoly = vtk.vtkPolyData()
cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
cutPoly.SetPolys(cutStrips.GetOutput().GetLines())

cutEdges2 = vtk.vtkCutter()
cutEdges2.SetInputConnection(cowNormals2.GetOutputPort())
cutEdges2.SetCutFunction(plane)
cutEdges2.GenerateCutScalarsOn()
cutEdges2.SetValue(0, 0.5)
cutStrips2 = vtk.vtkStripper()
cutStrips2.SetInputConnection(cutEdges2.GetOutputPort())
cutStrips2.Update()
cutPoly2 = vtk.vtkPolyData()
cutPoly2.SetPoints(cutStrips2.GetOutput().GetPoints())
cutPoly2.SetPolys(cutStrips2.GetOutput().GetLines())

# Triangle filter is robust enough to ignore the duplicate point at
# the beginning and end of the polygons and triangulate them.
cutTriangles = vtk.vtkTriangleFilter()
cutTriangles.SetInput(cutPoly)
cutMapper = vtk.vtkPolyDataMapper()
cutMapper.SetInput(cutPoly)
cutMapper.SetInputConnection(cutTriangles.GetOutputPort())
cutActor = vtk.vtkActor()
cutActor.SetMapper(cutMapper)
cutActor.GetProperty().SetColor(peacock)

cutTriangles2 = vtk.vtkTriangleFilter()
cutTriangles2.SetInput(cutPoly2)
cutMapper2 = vtk.vtkPolyDataMapper()
cutMapper2.SetInput(cutPoly2)
cutMapper2.SetInputConnection(cutTriangles2.GetOutputPort())
cutActor2 = vtk.vtkActor()
cutActor2.SetMapper(cutMapper2)
cutActor2.GetProperty().SetColor(tomato)

# Create graphics stuff
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Add the actors to the renderer, set the background and size
ren.AddActor(cutActor)
ren.AddActor(cutActor2)
ren.SetBackground(1, 1, 1)
ren.ResetCamera()
ren.GetActiveCamera().Azimuth(30)
ren.GetActiveCamera().Elevation(30)
ren.GetActiveCamera().Dolly(1.5)
ren.ResetCameraClippingRange()

renWin.SetSize(300, 300)
iren.Initialize()

# Lets you move the cut plane back and forth by invoking the function
# Cut with the appropriate plane value (essentially a distance from
# the original plane).  This is not used in this code but should give
# you an idea of how to define a function to do this.
def Cut(v):
    cutEdges.SetValue(0, v)
    cutStrips.Update()
    cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
    cutPoly.SetPolys(cutStrips.GetOutput().GetLines())
    cutMapper.Update()
    cutEdges2.SetValue(0, v)
    cutStrips2.Update()
    cutPoly2.SetPoints(cutStrips2.GetOutput().GetPoints())
    cutPoly2.SetPolys(cutStrips2.GetOutput().GetLines())
    cutMapper2.Update()
    renWin.Render()

renWin.Render()
iren.Start()


#%%



#%%




#!/usr/bin/env python
# coding: utf-8

# In[4]:


import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import pi
import vtk
from vtk.util.colors import peacock, tomato


# In[37]:


#functions for creating vtk representation

def radToGrad(angle):
    return angle*180/pi
            
def createSourceFromPoint(node):
    if node.parent is None:
        px = node.x
        py = node.y
        pz = node.z-node.length
    else:
        px = node.parent.x
        py = node.parent.y
        pz = node.parent.z
    source = vtk.vtkLineSource()
    source.SetPoint1(px,py,pz)
    source.SetPoint2(node.x, node.y, node.z)
    tube=vtk.vtkTubeFilter()
    tube.SetRadius(2*node.model.length_ratio**node.level)
    tube.SetNumberOfSides(100)
    tube.CappingOn()
    tube.SetInputConnection(source.GetOutputPort())
    return tube

def createActorFromPoint(node, color):
    source = createSourceFromPoint(node)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(source.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color[0], color[1], color[2])
    return actor

def isCutting(node, z):
    if node.parent is None:
        px = node.x
        py = node.y
        pz = node.z-node.length
    else:
        px = node.parent.x
        py = node.parent.y
        pz = node.parent.z
    minz = min(pz, node.z)
    maxz = max(pz, node.z)
    return minz < z < maxz;

def createCutFromPoint(node, color, plane):
    #z=plane.GetOrigin()[2]
    #if not isCutting(node, z):
    #   return None
    
    segment_source = createSourceFromPoint(node)
    segment_normal = vtk.vtkPolyDataNormals()
    segment_normal.SetInputConnection(segment_source.GetOutputPort())
    cutEdges = vtk.vtkCutter()
    cutEdges.SetInputConnection(segment_normal.GetOutputPort())
    cutEdges.SetCutFunction(plane)
    cutEdges.GenerateCutScalarsOn()
    cutEdges.SetValue(0, 0.5)
    
    cutStrips = vtk.vtkStripper()
    cutStrips.SetInputConnection(cutEdges.GetOutputPort())
    cutStrips.Update()
    
    cutPoly = vtk.vtkPolyData()
    cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
    cutPoly.SetPolys(cutStrips.GetOutput().GetLines())
    
    cutTriangles = vtk.vtkTriangleFilter()
    cutTriangles.SetInput(cutPoly)
    cutMapper = vtk.vtkPolyDataMapper()
    cutMapper.SetInput(cutPoly)
    cutMapper.SetInputConnection(cutTriangles.GetOutputPort())
    
    cutActor = vtk.vtkActor()
    cutActor.SetMapper(cutMapper)
    cutActor.GetProperty().SetColor(color[0], color[1], color[2])

    return cutActor

def createNiceTreeSource(model, radius, radius_ratio):
    
    points = vtk.vtkPoints()
    scalars = vtk.vtkDoubleArray()
    lines = vtk.vtkCellArray()
    root = model.nodes[0]
    start_point = np.array(root.coords)
    start_point[2] -= root.length
    start_point_visual = points.InsertNextPoint(start_point)
    for node in model.nodes:
        node.visual = points.InsertNextPoint(node.coords)
        segment_radius = radius*(radius_ratio**node.level)
        scalars.InsertNextValue(segment_radius)
        prev_visual = node.parent.visual if node.parent else start_point_visual
        lines.InsertNextCell(2)
        lines.InsertCellPoint(prev_visual)
        lines.InsertCellPoint(node.visual)
    
    skeleton = vtk.vtkPolyData()
    skeleton.SetPoints(points)
    skeleton.SetLines(lines)
    skeleton.GetPointData().SetScalars(scalars)
    vtkTubeFilter = vtk.vtkTubeFilter()
    vtkTubeFilter.SetInput(skeleton)
    vtkTubeFilter.SetVaryRadiusToVaryRadiusByScalar() 
    vtkTubeFilter.SetRadius(0.15)
    vtkTubeFilter.SetRadiusFactor(2)
    vtkTubeFilter.SetNumberOfSides(20)
    vtkTubeFilter.CappingOn()
    return vtkTubeFilter
    


# In[20]:


class duct_model:
    def __init__(self, alpha, beta, length, length_ratio, n_levels):
        self.alpha = alpha
        self.beta = beta
        self.length = length
        self.length_ratio = length_ratio
        self.nodes = []
        self.n_levels = n_levels
    def add_node(self, node):
        node.idx = len(self.nodes)
        self.nodes.append(node)
        if node.level < self.n_levels:
            self.branch(node)
    def branch(self, parent):
        self.add_node(node(parent, -1))
        self.add_node(node(parent, 1))
    def create_ducts(self):
        self.add_node(node(model=self))
    

def rotation_matrix(axis, theta):
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
rotation_matrix

def rotate_point(p, axis, theta):
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)
    rotmat = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

    return np.dot(rotmat,p)
      
class node:
    def __init__(self, parent=None, sign=0, model=None):
        self.parent=parent
        if (parent is None):
            self.model = model
            self.length = model.length
            self.level = 0
            self.beta = pi/2
            self.alpha = 0
            bc = np.array([0,0,-model.length])
            self.axis = np.array([0, 0, 1])
            self.coords = self.length*self.axis
        else:
            self.model = parent.model
            self.level = parent.level + 1
            self.beta = parent.beta + sign * self.model.beta
            self.alpha = parent.alpha + alpha
            self.length = parent.length * self.model.length_ratio
            self.axis = np.array([0, cos(self.beta), sin(self.beta)])
            p = self.parent
            ancestors = []
            while p:
                ancestors.append(p)
                p = p.parent
            while len(ancestors):
                a = ancestors.pop()
                self.axis = rotate_point(self.axis, a.axis, self.model.alpha)
            self.coords = parent.coords + self.length*self.axis
   
    x = property(lambda self: self.coords[0])
    y = property(lambda self: self.coords[1])
    z = property(lambda self: self.coords[2])
    
class node_rotmat:
    def __init__(self, parent=None, sign=0, model=None):
        self.parent=parent
        if (parent is None):
            self.model = model
            self.length = model.length
            self.level = 0
            self.beta = pi/2
            self.alpha = 0
            bc = np.array([0,0,-model.length])
            self.rotation = np.identity(3)
            self.axis = np.array([0, 0, 1])
            self.coords = self.length*self.axis
        else:
            self.model = parent.model
            self.level = parent.level + 1
            self.beta = parent.beta + sign * self.model.beta
            self.alpha = parent.alpha + alpha
            self.length = parent.length * self.model.length_ratio
            self.rotation = np.dot(rotation_matrix(parent.axis, alpha), parent.rotation)
            self.axis = np.array([0, cos(self.beta), sin(self.beta)])
            self.axis = np.dot(self.rotation, self.axis)
            self.coords = parent.coords + self.length*self.axis
        #p = self.parent
        #ancestors = []
        #while p:
        #    ancestors.append(p)
        #    p = p.parent
        #while len(ancestors):
        #    a = ancestors.pop()
        #    self.axis = rotate_point(self.axis, a.axis, self.model.alpha)
        #self.coords = bc + self.length*self.axis
   
    x = property(lambda self: self.coords[0])
    y = property(lambda self: self.coords[1])
    z = property(lambda self: self.coords[2])


# In[21]:


#model parameters
l = 4
len_ratio = 0.88
alpha = pi/3
beta = pi/6
n_levels = 3


# In[22]:


model = duct_model(alpha, beta, l, len_ratio, n_levels)
model.create_ducts()


# In[23]:


# Create an object per segment (
ren = vtk.vtkRenderer()

lut = vtk.vtkLookupTable()
lutNum = model.n_levels
lut.SetNumberOfTableValues(lutNum)
ctf = vtk.vtkColorTransferFunction()
ctf.SetColorSpaceToDiverging()
ctf.AddRGBPoint(0.0, 0, 0, 1.0)
ctf.AddRGBPoint(1.0, 1.0, 0, 0)
for ii,ss in enumerate([float(xx)/float(lutNum) for xx in range(lutNum)]):
	cc = ctf.GetColor(ss)
	lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)

for n in model.nodes:
    color=lut.GetTableValue(n.level)
    ren.AddActor(createActorFromPoint(n, color))

vtkRenderWindow = vtk.vtkRenderWindow()
vtkRenderWindow.AddRenderer(ren)

vtkRenderWindowInteractor = vtk.vtkRenderWindowInteractor()
vtkRenderWindowInteractor.SetRenderWindow(vtkRenderWindow)
vtkRenderWindowInteractor.Initialize()
vtkRenderWindow.Render()
vtkRenderWindowInteractor.Start()


# In[10]:


# Create Cuts
ren = vtk.vtkRenderer()

lut = vtk.vtkLookupTable()
lutNum = model.n_levels
lut.SetNumberOfTableValues(lutNum)
ctf = vtk.vtkColorTransferFunction()
ctf.SetColorSpaceToDiverging()
ctf.AddRGBPoint(0.0, 0, 0, 1.0)
ctf.AddRGBPoint(1.0, 1.0, 0, 0)
for ii,ss in enumerate([float(xx)/float(lutNum) for xx in range(lutNum)]):
	cc = ctf.GetColor(ss)
	lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)

plane = vtk.vtkPlane()
plane.SetOrigin(0, 0, 15)
plane.SetNormal(0, 0, 1)
    
for n in model.nodes:
    color=lut.GetTableValue(n.level)
    actor = createCutFromPoint(n, color, plane)
    if actor:
        ren.AddActor(actor)

vtkRenderWindow = vtk.vtkRenderWindow()
vtkRenderWindow.AddRenderer(ren)

vtkRenderWindowInteractor = vtk.vtkRenderWindowInteractor()
vtkRenderWindowInteractor.SetRenderWindow(vtkRenderWindow)
vtkRenderWindowInteractor.Initialize()
vtkRenderWindow.Render()
vtkRenderWindowInteractor.Start()


# In[ ]:


# Tree at the Ric's way

source = createNiceTreeSource(model, 3, 0.88)

vtkPolyDataMapper = vtk.vtkPolyDataMapper()
vtkPolyDataMapper.SetInputConnection(source.GetOutputPort())

vtkActor = vtk.vtkActor()
vtkActor.SetMapper(vtkPolyDataMapper)

vtkRenderer = vtk.vtkRenderer()
vtkRenderer.AddActor(vtkActor)

vtkRenderWindow = vtk.vtkRenderWindow()
vtkRenderWindow.AddRenderer(vtkRenderer)

vtkRenderWindowInteractor = vtk.vtkRenderWindowInteractor()
vtkRenderWindowInteractor.SetRenderWindow(vtkRenderWindow)
vtkRenderWindowInteractor.Initialize()
vtkRenderWindow.Render()
vtkRenderWindowInteractor.Start()


# In[19]:





# In[ ]:





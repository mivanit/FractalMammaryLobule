#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import vtk
from vtk.util.colors import peacock, tomato


# In[5]:


#model parameters
l = 4
len_ratio = 0.88
alpha = 0
beta = pi/6
n_levels = 5


# In[2]:


def doScatterPlot(x, y, z, colors):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax = fig.gca()

    plt.hold(True)
    cm = matplotlib.cm.get_cmap('RdYlBu')

    ax.scatter(x, y, z, c = colors, cmap=cm)


# In[3]:


#original way of creating points

def successive_points(previous_p, l, theta_alpha, theta_beta):
    next_point_0 = [x for x in previous_p]
    alpha_angle = previous_p[3] + theta_alpha
    beta_angle = previous_p[4] + theta_beta
    
    next_point_0[0] += l * np.cos(beta_angle)*np.sin(alpha_angle)  # x   
    next_point_0[1] += l * np.sin(beta_angle)                # y
    next_point_0[2] += l * np.cos(beta_angle)*np.cos(alpha_angle)  # z
    next_point_0[3] = alpha_angle
    next_point_0[4] = beta_angle
    
    next_point_1 = [x for x in previous_p]
    alpha_angle = previous_p[3] - theta_alpha
    beta_angle = previous_p[4] - theta_beta
    next_point_1[0] += l * np.cos(beta_angle)*np.sin(alpha_angle) # x
    next_point_1[1] += l * np.sin(beta_angle)                # y
    next_point_1[2] += l * np.cos(beta_angle)*np.cos(alpha_angle) # z
    next_point_1[3] = alpha_angle
    next_point_1[4] = beta_angle
    
    return [next_point_0, next_point_1]


# In[4]:


p0 = [0, 0, 0, 0, 0]
points = [[p0]]

for I in range(1, n_levels+1):
    next_points = [successive_points(point, l, alpha, beta) for point in points[I-1]]
    points.append([point for tup_points in next_points for point in tup_points])
    l *= len_ratio # Exponentially decaying Length


# In[7]:


x = []
y = []
z = []
plot_colors = []

for point_tup, color in zip(points, range(len(points)) ):
    for point in point_tup:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
        plot_colors.append(color)

doScatterPlot(x, y, z, plot_colors)


# In[187]:


class duct_model:
    LEFT=0
    RIGHT=1
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
        node1 = node(parent, duct_model.LEFT)
        self.add_node(node1)
        node2 = node(parent, duct_model.RIGHT)
        self.add_node(node(parent, duct_model.RIGHT
    def create_ducts(self):
        self.add_node(node(model=self))

def rotate(p, axis, theta):
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)
    rotmat = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

    return np.dot(rotmat,p)
      
class node:
    def __init__(self, parent=None, node_side=None, model=None):
        self.parent=parent
        if (parent is None):
            self.model = model
            self.length = model.length
            self.level = 0
            self.beta = pi/2
            self.axis = np.array([0,0,1])
            bc = np.array([0,0,-model.length])
        else:
            self.model = parent.model
            self.level = parent.level + 1
            sign = 1 if node_side==duct_model.LEFT else -1
            self.beta = parent.beta + sign * self.model.beta
            self.length = parent.length * self.model.length_ratio
            bc = parent.coords
        self.axis = np.array([0, cos(self.beta), sin(self.beta)])
        if parent is not None:
            self.axis = rotate(self.axis, parent.axis, alpha*self.level)
        self.coords = bc + self.length*self.axis

    x = property(lambda self: self.coords[0])
    y = property(lambda self: self.coords[1])
    z = property(lambda self: self.coords[2])


# In[188]:


model = duct_model(alpha, beta, l, len_ratio, n_levels)
model.create_ducts()


# In[189]:


#scatter plot 
x = [p.x for p in model.nodes]
y = [p.y for p in model.nodes]
z = [p.z for p in model.nodes]
plot_colors = [p.level for p in model.nodes]

doScatterPlot(x, y, z, plot_colors)


# In[12]:


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
    tube.SetRadius(0.5*node.model.length_ratio**node.level)
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


# In[194]:


# Create a tree as single object (single color)
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
apd=vtk.vtkAppendPolyData() 
for n in model.nodes:
    apd.AddInput(createSourceFromPoint(n).GetOutput())

mapper = vtk.vtkPolyDataMapper() 
mapper.SetInput(apd.GetOutput())             
actor = vtk.vtkActor() 
actor.SetMapper(mapper) 
actor.GetProperty().SetColor(peacock)
ren.AddActor(actor)

ren.GetActiveCamera().Elevation(30)
ren.ResetCameraClippingRange()
renWin.SetSize(300, 300)
iren.Initialize()

renWin.Render()
iren.Start()


# In[166]:


# Create an object per segment (
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
apd=vtk.vtkAppendPolyData() 

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

ren.GetActiveCamera().Elevation(30)
ren.ResetCameraClippingRange()
renWin.SetSize(300, 300)
iren.Initialize()

renWin.Render()
iren.Start()


# In[200]:





# In[201]:





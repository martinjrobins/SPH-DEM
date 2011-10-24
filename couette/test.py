from mayavi import ivtk
from vtk import *
c = vtkConeSource()
m = vtkPolyDataMapper()
m.SetInput(c.GetOutput())
a = vtkActor()
a.SetMapper(m)
v = ivtk.create_viewer() # or ivtk.viewer()
# this creates the easy to use render window that can be used from
# the interpreter.  It has several useful menus.

v.AddActors(a)    # add actor(s) to viewer
v.config(c)       # pops up a GUI configuration for object.
v.doc(c)          # pops up class documentation for object.
v.help_browser()  # pops up a help browser where you can search!
v.RemoveActors(a) # remove actor(s) from viewer.

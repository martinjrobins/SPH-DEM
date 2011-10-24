import mayavi

theMayavi = mayavi.mayavi()
theMayavi.open_vtk("resultCouetteGrid0000002.vtk",0)
theMayavi.load_module("GridPlane",0)
theMayavi.open_vtk("resultCouette0000002.vtk",0)
theMayavi.load_module("Glyph",0)


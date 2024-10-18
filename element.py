import gmsh
import numpy as np

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("Serendipity_16x16_Mesh")

# Define geometry: Create a rectangle
length = 1.0  # Length of the rectangle (x-direction)
width = 1.0   # Width of the rectangle (y-direction)
Nx = 16       # Number of divisions in x-direction
Ny = 16       # Number of divisions in y-direction

# Create points at the corners
p1 = gmsh.model.geo.addPoint(0, 0, 0)
p2 = gmsh.model.geo.addPoint(length, 0, 0)
p3 = gmsh.model.geo.addPoint(length, width, 0)
p4 = gmsh.model.geo.addPoint(0, width, 0)

# Create lines between points to form the boundary
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# Define a loop and create a plane surface
loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
surface = gmsh.model.geo.addPlaneSurface([loop])

# Synchronize CAD with Gmsh
gmsh.model.geo.synchronize()

# Set meshing parameters for a 16x16 mesh
gmsh.option.setNumber("Mesh.MeshSizeMax", length / Nx)
gmsh.option.setNumber("Mesh.MeshSizeMin", length / Nx)

# Set element order to 2 for quadratic serendipity elements
gmsh.model.mesh.setOrder(2)  # Generates 8-node serendipity elements
gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)  # Optimize for high-order serendipity elements

# Generate the 2D mesh
gmsh.model.mesh.generate(2)

# Save the mesh to a file
gmsh.write("serendipity_16x16_mesh.msh")

# Optional: Visualize the mesh using Gmsh's GUI
gmsh.fltk.run()

# Finalize the Gmsh API session
gmsh.finalize()

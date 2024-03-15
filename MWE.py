from fenics import *

inlet_id = 8
outlet_id = 7
fluid_id = 6

mesh = Mesh()
XDMFFile("mesh_cells.xdmf").read(mesh)

print("Succesfully load mesh with " + str(len(mesh.cells())) + " cells")

# create dx from the volume markers
volume_markers = MeshFunction("size_t", mesh, mesh.topology().dim())
XDMFFile("mesh_cells.xdmf").read(volume_markers)
dx = Measure("dx", domain=mesh, subdomain_data=volume_markers)

print("Succesfully load mesh with " + str(len(volume_markers)) + " cells")
import numpy as np

print(np.unique(volume_markers.array()))
surface_markers = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
XDMFFile("mesh_facets.xdmf").read(surface_markers, "f")
surface_markers = MeshFunction("size_t", mesh, surface_markers)

V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
v = TestFunction(V)

inlet = DirichletBC(V, Constant(1), surface_markers, inlet_id)
outlet = DirichletBC(V, Constant(0), surface_markers, outlet_id)
bcs = [inlet, outlet]

F = dot(2 * grad(u), grad(v)) * dx(fluid_id)  # doesn't work
# F = 2 * dot(grad(u), grad(v)) * dx(fluid_id)  # doesn't works
# F = dot(grad(u), grad(v)) * dx(fluid_id)  # works
# F = dot(2 * grad(u), grad(v)) * dx  # works

solve(F == 0, u, bcs)

XDMFFile("solution.xdmf").write(u)

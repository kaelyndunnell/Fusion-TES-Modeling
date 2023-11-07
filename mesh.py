from fenics import plot, Point, CompiledSubDomain, MeshFunction, XDMFFile
from mshr import Rectangle, generate_mesh


height_fluid = 3e-2
pipe_thickness = 4e-3

length = float(input('length ='))

refinement = 200

p1 = Point(0, 0)
p2 = Point(length, height_fluid)
fluid_rectangle = Rectangle(p1, p2)

p1 = Point(0, height_fluid)
p2 = Point(length, height_fluid + pipe_thickness)
pipe_rectangle = Rectangle(p1, p2)

domain = fluid_rectangle + pipe_rectangle

domain.set_subdomain(1, fluid_rectangle)
domain.set_subdomain(2, pipe_rectangle)
mesh = generate_mesh(domain, refinement)

# marking physical groups (volumes and surfaces)
volume_markers = MeshFunction("size_t", mesh, mesh.topology().dim())
volume_markers.set_all(1)

# mark subdomains
from fenics import SubDomain

tol = 1e-14

id_fluid = 5
id_pipe_walls = 6


class Fluid(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] <= height_fluid + tol


class Pipe(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] >= height_fluid - tol


fluid = Fluid()
pipe = Pipe()
# marking volumes
fluid.mark(volume_markers, id_fluid)
pipe.mark(volume_markers, id_pipe_walls)

tol = 1e-14

inlet_surface = CompiledSubDomain(
    "on_boundary && near(x[0], 0, tol) && x[1] < height_fluid + tol",
    tol=tol,
    height_fluid=height_fluid,
)
outlet_surface = CompiledSubDomain(
    "on_boundary && near(x[0], outlet_x, tol) && x[1] < height_fluid + tol",
    tol=tol,
    height_fluid=height_fluid,
    outlet_x=length,
)
bottom_surface = CompiledSubDomain("on_boundary && near(x[1], 0, tol)", tol=tol)
top_surface = CompiledSubDomain(
    "on_boundary && near(x[1], top_y, tol)",
    tol=tol,
    top_y=height_fluid + pipe_thickness,
)

surface_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
surface_markers.set_all(0)

inlet_id = 1
top_id = 2
outlet_id = 3
bottom_id = 4
inlet_surface.mark(surface_markers, inlet_id)
outlet_surface.mark(surface_markers, outlet_id)
top_surface.mark(surface_markers, top_id)
bottom_surface.mark(surface_markers, bottom_id)

output_file = XDMFFile("surface_markers.xdmf")
output_file.write(surface_markers)

output_file2 = XDMFFile("volume_markers.xdmf")
output_file2.write(volume_markers)
# modeling extractor efficiency

import festim as F
import matplotlib.pyplot as plt
import fenics
from dolfin import *

my_model = F.Simulation()

inlet_id = 1
top_id = 2
outlet_id = 3
bottom_id = 4

id_fluid = 5
id_pipe_walls = 6

temperature = 1200
# temperature = input('temperature = ')

# Read mesh file
my_model.mesh = F.MeshFromXDMF(
    volume_file="volume_markers.xdmf",
    boundary_file="surface_markers.xdmf",
    type="cylindrical",
)

# verify these numbers and insert with your own materials

eurofer = F.Material(  # if have time, compare efficnecy between PAV and without membrane
        id=id_pipe_walls,
        D_0=4.57e-7,
        E_D=0.23,
        S_0=1.35e22,
        E_S=0.16,
    )
lipb = F.Material(
    id=id_fluid,
    D_0=4.68e-08,
    E_D=0.2021,
    S_0=1.62e21,
    E_S=0.17,
)

my_model.materials = F.Materials([eurofer, lipb])

my_model.T = F.Temperature(temperature)

my_model.boundary_conditions = [
    F.DirichletBC(field="solute", surfaces=top_id, value=0),
    # F.DirichletBC(field="solute", surfaces=bottom_id, value=1e18),
    F.DirichletBC(field="solute", surfaces=inlet_id, value=1e18),
]

c_out = F.AverageSurface("solute", outlet_id)
c_in = F.AverageSurface("solute", inlet_id)

derived_quantities = F.DerivedQuantities([c_in, c_out])

my_model.exports = F.Exports(
    [F.XDMFExport("solute", folder="task3/", mode=1), derived_quantities]
)

my_model.settings = F.Settings(
    absolute_tolerance=1e10,
    relative_tolerance=1e-10,
    transient=False,
    chemical_pot=True,
)

my_model.initialise()

# adding advection term:
from fenics import interpolate, Expression, VectorFunctionSpace, SubMesh, XDMFFile, Function, TestFunction, solve, inner

mesh_sub = SubMesh(my_model.mesh.mesh, my_model.mesh.volume_markers, id_fluid)

functionspace = VectorFunctionSpace(mesh_sub, "CG", 1)

velocity = Expression(("6e-05*(x[1] - 0.5)*(x[1] + 0.5)", "0"), degree=2)

velocity = interpolate(velocity, functionspace)

V = VectorFunctionSpace(my_model.mesh.mesh, "CG", 1)
u = Function(V)
v = TestFunction(V)

form = inner(u, v) * my_model.mesh.dx
form += inner(velocity, v) * my_model.mesh.dx(id_fluid)
solve(form == 0, u, bcs=[])

velocity = u

XDMFFile("velocity_field.xdmf").write(u)

# velocity profile
import matplotlib.pyplot as plt
from fenics import plot

# CS = plot(velocity)
# plt.colorbar(CS)
# plt.show()

print(my_model.h_transport_problem.F)

# add advection term to FESTIM model

from fenics import inner, dot, grad

my_model.initialise()  # reinitialisation is needed

(
    hydrogen_concentration,
    _,
) = my_model.h_transport_problem.mobile.get_concentration_for_a_given_material(
    lipb, my_model.T
)
test_function_mobile = my_model.h_transport_problem.mobile.test_function

advection_term = inner(
    dot(grad(hydrogen_concentration), u), test_function_mobile
) * my_model.mesh.dx(id_fluid)

my_model.h_transport_problem.F += advection_term
# print(my_model.h_transport_problem.F)

# eff_lst = []

my_model.run()
# for i in range(0,2):
#     my_model.run()
#     c_out_value = c_out.data[0]
#     c_in_value = c_in.data[0]   
#     n = 1 - c_out_value / c_in_value
#     eff_lst.append(n)
#     print(n)
#     # temperature += 100
#     # my_model.T = F.Temperature(temperature)


post_processing_mobile_conc = (
    my_model.h_transport_problem.mobile.post_processing_solution
)

plot(post_processing_mobile_conc)
XDMFFile("mobile_conc_with_advection.xdmf").write(post_processing_mobile_conc)

c_out_value = c_out.data[0]
c_in_value = c_in.data[0] 

n = 1 - c_out_value / c_in_value  # efficiency

print(f'efficiency is {n}')
# print(f'temperature is {temperature}')
# print(f'length is {length}')

# plt.show()

# print(eff_lst)


import festim as F
import matplotlib.pyplot as plt
import fenics
import numpy as np
from dolfin import *
from fenics import (
    interpolate,
    Expression,
    VectorFunctionSpace,
    SubMesh,
    XDMFFile,
    Function,
    TestFunction,
    solve,
    inner,
    set_log_level,
)
from fenics import inner, dot, grad

from mesh_function import (
    create_mesh,
    id_fluid,
    id_pipe_walls,
    inlet_id,
    outlet_id,
    top_id,
)


def run_model(temperature, length, height_fluid, pipe_thickness):
    my_model = F.Simulation()

    # create mesh
    create_mesh(length, height_fluid, pipe_thickness)

    my_model.mesh = F.MeshFromXDMF(
        volume_file="volume_markers.xdmf",
        boundary_file="surface_markers.xdmf",
        # type="cylindrical",
    )
    # materials
    eurofer = (
        F.Material(  # if have time, compare efficnecy between PAV and without membrane
            id=id_pipe_walls,
            D_0=4.57e-7,
            E_D=0.23,
            S_0=1.35e22,
            E_S=0.16,
        )
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
        F.DirichletBC(field="solute", surfaces=inlet_id, value=1e18),
    ]

    c_in = F.AverageSurface("solute", inlet_id)
    c_out = F.AverageSurface("solute", outlet_id)
    

    derived_quantities = F.DerivedQuantities([c_in, c_out])

    my_model.exports = F.Exports(
        [F.XDMFExport("solute", folder="results/", mode=1), derived_quantities]
    )

    my_model.settings = F.Settings(
        absolute_tolerance=1e10,
        relative_tolerance=1e-10,
        transient=False,
        chemical_pot=True,
    )
    my_model.initialise()
    # adding advection term:
    mesh_sub = SubMesh(my_model.mesh.mesh, my_model.mesh.volume_markers, id_fluid)

    functionspace = VectorFunctionSpace(mesh_sub, "CG", 1)

    velocity = Expression(
        ("6e2*(x[1] - height_fluid)*(x[1] + height_fluid)", "0"),
        height_fluid=height_fluid,
        degree=2,
    )

    velocity = interpolate(velocity, functionspace)

    V = VectorFunctionSpace(my_model.mesh.mesh, "CG", 1)
    u = Function(V)
    v = TestFunction(V)

    # set_log_level(20)

    form = inner(u, v) * my_model.mesh.dx
    form += inner(velocity, v) * my_model.mesh.dx(id_fluid)
    solve(form == 0, u, bcs=[])

    velocity = u

    XDMFFile("velocity_field.xdmf").write(u)
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

    my_model.run()

    return c_out, c_in


def compute_efficiency(temperature, length, height_fluid, pipe_thickness):
    c_out, c_in = run_model(temperature, length, height_fluid, pipe_thickness)
    c_out_value = c_out.data[0]
    c_in_value = c_in.data[0]
    efficiency = 1 - c_out_value / c_in_value
    return efficiency


if __name__ == "__main__":
    # run_model(temperature=700, length=0.3, height_fluid=1e-3, pipe_thickness=4e-3)
    print(compute_efficiency(temperature=800, length=0.3, height_fluid=1e-2, pipe_thickness=4e-3))

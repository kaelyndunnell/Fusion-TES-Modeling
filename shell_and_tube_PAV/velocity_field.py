import festim as F
from dolfinx import *
import dolfinx as fe
import tqdm.autonotebook
import numpy as np
import h_transport_materials as htm

# writes velocity field to an xdmf so that doesn't have to solve for it each time in the festim simulation

fluid_id = 6
inlet_id = 7
outlet_id = 8
vacuum_id = 9
walls_id = 10


def rho_lipb(T):  # units in kg/(m**3)
    return 10520.35 - 1.19051 * T


def visc_lipb(T):  # units (Pa s)
    return (
        0.01555147189 - 4.827051855e-05 * T + 5.641475215e-08 * T**2 - 2.2887e-11 * T**3
    )


def fluid_dynamics_sim_chorin(
    mesh, volume_markers, surface_markers, id_inlet, id_outlet, id_walls
):
    """
    Solves the Navier-Stokes equations using the Chorin's projection method
    See https://fenicsproject.org/pub/tutorial/html/._ftut1009.html#ftut1:NS for more details

    Args:
        mesh (fenics.Mesh): the mesh
        volume_markers (fenics.MeshFunction): the volume markers (subdomains)
        surface_markers (fenics.MeshFunction): the surface markers (boundaries)
        id_inlet (int): the id of the inlet boundary
        id_outlet (int): the id of the outlet boundary
        id_walls (int or list): the id(s) of the walls

    Returns:
        fenics.Function: the velocity field
        fenics.Function: the pressure field
    """
    V = fe.VectorFunctionSpace(mesh, "CG", 2)
    Q = fe.FunctionSpace(mesh, "CG", 1)

    u = fe.TrialFunction(V)
    p = fe.TrialFunction(Q)
    v = fe.TestFunction(V)
    q = fe.TestFunction(Q)

    u_n = fe.Function(V)
    u_ = fe.Function(V)
    p_n = fe.Function(Q)
    p_ = fe.Function(Q)

    # ##### Boundary conditions ##### #

    inlet_velocity = 2e-03  # units: m s-1
    outlet_pressure = 0  # units: Pa

    inflow = fe.DirichletBC(
        V, fe.Constant((0.0, 0.0, -inlet_velocity)), surface_markers, id_inlet
    )

    # make sure id_walls is a list
    if isinstance(id_walls, int):
        id_walls = [id_walls]

    # iterate through the walls
    walls = []
    for id_wall in id_walls:
        walls.append(
            fe.DirichletBC(V, fe.Constant((0.0, 0.0, 0.0)), surface_markers, id_wall)
        )

    pressure_outlet = fe.DirichletBC(
        Q, fe.Constant(outlet_pressure), surface_markers, id_outlet
    )
    bcu = [inflow] + walls
    bcp = [pressure_outlet]

    # ##### Solver ##### #
    dx = fe.Measure("dx", subdomain_data=volume_markers)
    ds = fe.Measure("ds", subdomain_data=surface_markers)

    t = 0
    total_time = 20
    dt = 0.1  # Time step size
    num_steps = int(total_time / dt)

    k = fe.Constant(dt)
    n = fe.FacetNormal(mesh)
    U = 0.5 * (u_n + u)

    # LiPb
    T = 700  # TODO make this more generic
    mu = visc_lipb(T)
    rho = rho_lipb(T)

    def epsilon(u):
        return fe.sym(fe.nabla_grad(u))

    def sigma(u, p):
        return 2 * mu * epsilon(u) - p * fe.Identity(len(u))

    # ##### Solver ##### #

    # Tentative velocity step
    F1 = rho * fe.dot((u - u_n) / k, v) * dx
    F1 += rho * fe.dot(fe.dot(u_n, fe.nabla_grad(u_n)), v) * dx
    F1 += fe.inner(sigma(U, p_n), epsilon(v)) * dx
    F1 += fe.dot(p_n * n, v) * ds
    F1 -= fe.dot(mu * fe.nabla_grad(U) * n, v) * ds

    a1 = fe.lhs(F1)
    L1 = fe.rhs(F1)

    solver1 = fe.KrylovSolver("gmres", "jacobi")
    solver1.parameters["absolute_tolerance"] = 1e-08
    solver1.parameters["relative_tolerance"] = 1e-08
    solver1.parameters["maximum_iterations"] = 1000
    solver1.parameters["report"] = False
    solver1.parameters["monitor_convergence"] = False

    # Pressure update
    a2 = fe.dot(fe.nabla_grad(p), fe.nabla_grad(q)) * dx
    L2 = fe.dot(fe.nabla_grad(p_n), fe.nabla_grad(q)) * dx
    L2 -= (1 / k) * fe.div(u_) * q * dx

    solver2 = fe.KrylovSolver("bicgstab", "hypre_amg")
    solver2.parameters["absolute_tolerance"] = 1e-08
    solver2.parameters["relative_tolerance"] = 1e-08
    solver2.parameters["maximum_iterations"] = 1000
    solver2.parameters["report"] = False
    solver2.parameters["monitor_convergence"] = False

    # Velocity update
    a3 = fe.dot(u, v) * dx
    L3 = fe.dot(u_, v) * dx
    L3 -= k * fe.dot(fe.nabla_grad(p_ - p_n), v) * dx

    solver3 = fe.KrylovSolver("cg", "sor")
    solver3.parameters["absolute_tolerance"] = 1e-08
    solver3.parameters["relative_tolerance"] = 1e-08
    solver3.parameters["maximum_iterations"] = 1000
    solver3.parameters["report"] = False
    solver3.parameters["monitor_convergence"] = False

    # Assemble matrices
    A1 = fe.assemble(a1)
    A2 = fe.assemble(a2)
    A3 = fe.assemble(a3)

    velocity_file = fe.XDMFFile("velocity_field.xdmf")
    pressure_file = fe.XDMFFile("pressure_field.xdmf")

    max_u = []

    # Time-stepping
    progress = tqdm.autonotebook.tqdm(desc="Solving Navier-Stokes", total=num_steps)
    append = False
    for i in range(num_steps):
        progress.update(1)
        # Update current time step
        t += dt

        # Compute tentative velocity step
        b1 = fe.assemble(L1)
        [bc.apply(A1, b1) for bc in bcu]
        solver1.solve(A1, u_.vector(), b1)

        # Pressure correction
        b2 = fe.assemble(L2)
        [bc.apply(A2, b2) for bc in bcp]
        solver2.solve(A2, p_.vector(), b2)

        # Velocity correction
        b3 = fe.assemble(L3)
        [bc.apply(A3, b3) for bc in bcu]
        solver3.solve(A3, u_.vector(), b3)

        # Move to next time step
        u_n.assign(u_)
        p_n.assign(p_)

        max_u.append(u_.vector().max())
        np.savetxt("3D_case_max_u.txt", np.array(max_u))

        velocity_file.write_checkpoint(
            u_,
            "velocity",
            t,
            XDMFFile.Encoding.HDF5,
            append=append,
        )
        pressure_file.write_checkpoint(
            p_,
            "pressure",
            t,
            XDMFFile.Encoding.HDF5,
            append=append,
        )
        append = True
    return u_, p_


if __name__ == "__main__":
    my_sim = F.Simulation()

    my_sim.mesh = F.MeshFromXDMF(
        volume_file="mesh_cells.xdmf", boundary_file="mesh_facets.xdmf"
    )

    fluid_dynamics_sim_chorin(
        mesh=my_sim.mesh.mesh,
        volume_markers=my_sim.mesh.volume_markers,
        surface_markers=my_sim.mesh.surface_markers,
        id_inlet=inlet_id,
        id_outlet=outlet_id,
        id_walls=[walls_id, vacuum_id],
    )

import festim as F
from fenics import *
import fenics as fe
import tqdm.autonotebook
import numpy as np
import h_transport_materials as htm


# TODO find original refs for these but for now see James' paper
# TODO is T in K or degC?
# def rho_lipb(T):  # units in kg/(m**3)
#     return 10520.35 - 1.19051 * T


# def visc_lipb(T):  # units (Pa s)
#     return (
#         0.01555147189 - 4.827051855e-05 * T + 5.641475215e-08 * T**2 - 2.2887e-11 * T**3
#     )


# def fluid_dynamics_sim(
#     mesh, volume_markers, surface_markers, id_inlet, id_outlet, id_walls
# ):
#     """
#     Solves the Navier-Stokes equations in steady state.
#     Doesn't work for large meshes cause too memory intensive

#     Args:
#         mesh (fenics.Mesh): the mesh
#         volume_markers (fenics.MeshFunction): the volume markers (subdomains)
#         surface_markers (fenics.MeshFunction): the surface markers (boundaries)
#         id_inlet (int): the id of the inlet boundary
#         id_outlet (int): the id of the outlet boundary
#         id_walls (int): the id of the walls
#     """

#     V_ele = fe.VectorElement("CG", mesh.ufl_cell(), 2)
#     Q_ele = fe.FiniteElement("CG", mesh.ufl_cell(), 1)
#     W = fe.FunctionSpace(mesh, fe.MixedElement([V_ele, Q_ele]))

#     # ##### Boundary conditions ##### #

#     inlet_velocity = 2e-04  # units: m s-1
#     outlet_pressure = 0  # units: Pa

#     # Simulation boundary conditions
#     inflow = fe.DirichletBC(
#         W.sub(0),
#         fe.Constant((inlet_velocity, 0.0, 0.0)),
#         surface_markers,
#         id_inlet,
#     )
#     walls = fe.DirichletBC(
#         W.sub(0), fe.Constant((0.0, 0.0, 0.0)), surface_markers, id_walls
#     )
#     pressure_outlet = fe.DirichletBC(
#         W.sub(1), fe.Constant(outlet_pressure), surface_markers, id_outlet
#     )
#     bcs = [inflow, pressure_outlet, walls]

#     # ##### CFD --> Define Variational Parameters ##### #

#     v, q = fe.TestFunctions(W)
#     up = fe.Function(W)
#     u, p = fe.split(up)

#     # ##### CFD --> Fluid Materials properties ##### #

#     # Fluid properties
#     rho = 1  # rho_0_lipb
#     mu = 1  # visc_lipb

#     # ##### Solver ##### #
#     dx = fe.Measure("dx", subdomain_data=volume_markers)
#     ds = fe.Measure("ds", subdomain_data=surface_markers)

#     F = (
#         #           momentum
#         rho * fe.inner(fe.grad(u), fe.grad(v)) * dx
#         - fe.inner(p, fe.div(v)) * dx
#         + mu * fe.inner(fe.dot(fe.grad(u), u), v) * dx
#         #           continuity
#         + fe.inner(fe.div(u), q) * dx
#     )
#     print("Solving Navier-Stokes")
#     fe.solve(
#         F == 0,
#         up,
#         bcs,
#         solver_parameters={
#             "newton_solver": {"linear_solver": "mumps"},
#         },
#     )

#     u, p = up.split()

#     fe.XDMFFile("Results/velocity_field.xdmf").write_checkpoint(
#         u, "u", 0, fe.XDMFFile.Encoding.HDF5, append=False
#     )
#     return u, p


# def fluid_dynamics_sim_chorin(
#     mesh, volume_markers, surface_markers, id_inlet, id_outlet, id_walls
# ):
#     """
#     Solves the Navier-Stokes equations using the Chorin's projection method
#     See https://fenicsproject.org/pub/tutorial/html/._ftut1009.html#ftut1:NS for more details

#     Args:
#         mesh (fenics.Mesh): the mesh
#         volume_markers (fenics.MeshFunction): the volume markers (subdomains)
#         surface_markers (fenics.MeshFunction): the surface markers (boundaries)
#         id_inlet (int): the id of the inlet boundary
#         id_outlet (int): the id of the outlet boundary
#         id_walls (int or list): the id(s) of the walls

#     Returns:
#         fenics.Function: the velocity field
#         fenics.Function: the pressure field
#     """
#     V = fe.VectorFunctionSpace(mesh, "CG", 2)
#     Q = fe.FunctionSpace(mesh, "CG", 1)

#     u = fe.TrialFunction(V)
#     p = fe.TrialFunction(Q)
#     v = fe.TestFunction(V)
#     q = fe.TestFunction(Q)

#     u_n = fe.Function(V)
#     u_ = fe.Function(V)
#     p_n = fe.Function(Q)
#     p_ = fe.Function(Q)

#     # ##### Boundary conditions ##### #

#     inlet_velocity = 2e-03  # units: m s-1
#     outlet_pressure = 0  # units: Pa

#     inflow = fe.DirichletBC(
#         V, fe.Constant((inlet_velocity, 0.0, 0.0)), surface_markers, id_inlet
#     )

#     # make sure id_walls is a list
#     if isinstance(id_walls, int):
#         id_walls = [id_walls]

#     # iterate through the walls
#     walls = []
#     for id_wall in id_walls:
#         walls.append(
#             fe.DirichletBC(V, fe.Constant((0.0, 0.0, 0.0)), surface_markers, id_wall)
#         )

#     pressure_outlet = fe.DirichletBC(
#         Q, fe.Constant(outlet_pressure), surface_markers, id_outlet
#     )
#     bcu = [inflow] + walls
#     bcp = [pressure_outlet]

#     # ##### Solver ##### #
#     dx = fe.Measure("dx", subdomain_data=volume_markers)
#     ds = fe.Measure("ds", subdomain_data=surface_markers)

#     t = 0
#     total_time = 20
#     dt = 0.1  # Time step size
#     num_steps = int(total_time / dt)
#     num_steps = 10

#     k = fe.Constant(dt)
#     n = fe.FacetNormal(mesh)
#     U = 0.5 * (u_n + u)

#     # LiPb
#     T = 700  # TODO make this more generic
#     mu = visc_lipb(T)
#     rho = rho_lipb(T)

#     def epsilon(u):
#         return fe.sym(fe.nabla_grad(u))

#     def sigma(u, p):
#         return 2 * mu * epsilon(u) - p * fe.Identity(len(u))

#     # ##### Solver ##### #

#     # Tentative velocity step
#     F1 = rho * fe.dot((u - u_n) / k, v) * dx
#     F1 += rho * fe.dot(fe.dot(u_n, fe.nabla_grad(u_n)), v) * dx
#     F1 += fe.inner(sigma(U, p_n), epsilon(v)) * dx
#     F1 += fe.dot(p_n * n, v) * ds
#     F1 -= fe.dot(mu * fe.nabla_grad(U) * n, v) * ds

#     a1 = fe.lhs(F1)
#     L1 = fe.rhs(F1)

#     solver1 = fe.KrylovSolver("bicgstab", "hypre_amg")
#     solver1.parameters["absolute_tolerance"] = 1e-08
#     solver1.parameters["relative_tolerance"] = 1e-08
#     solver1.parameters["maximum_iterations"] = 1000
#     solver1.parameters["report"] = False
#     solver1.parameters["monitor_convergence"] = False

#     # Pressure update
#     a2 = fe.dot(fe.nabla_grad(p), fe.nabla_grad(q)) * dx
#     L2 = fe.dot(fe.nabla_grad(p_n), fe.nabla_grad(q)) * dx
#     L2 -= (1 / k) * fe.div(u_) * q * dx

#     solver2 = fe.KrylovSolver("bicgstab", "hypre_amg")
#     solver2.parameters["absolute_tolerance"] = 1e-08
#     solver2.parameters["relative_tolerance"] = 1e-08
#     solver2.parameters["maximum_iterations"] = 1000
#     solver2.parameters["report"] = False
#     solver2.parameters["monitor_convergence"] = False

#     # Velocity update
#     a3 = fe.dot(u, v) * dx
#     L3 = fe.dot(u_, v) * dx
#     L3 -= k * fe.dot(fe.nabla_grad(p_ - p_n), v) * dx

#     solver3 = fe.KrylovSolver("cg", "sor")
#     solver3.parameters["absolute_tolerance"] = 1e-08
#     solver3.parameters["relative_tolerance"] = 1e-08
#     solver3.parameters["maximum_iterations"] = 1000
#     solver3.parameters["report"] = False
#     solver3.parameters["monitor_convergence"] = False

#     # Assemble matrices
#     A1 = fe.assemble(a1)
#     A2 = fe.assemble(a2)
#     A3 = fe.assemble(a3)

#     [bc.apply(A1) for bc in bcu]
#     [bc.apply(A2) for bc in bcp]

#     results_folder = "Results/velocity_fields/"
#     velocity_file = fe.XDMFFile(results_folder + "u_3D_sub.xdmf")
#     pressure_file = fe.XDMFFile("Results/3D_results/p_3D_sub.xdmf")

#     max_u = []

#     # Time-stepping
#     progress = tqdm.autonotebook.tqdm(desc="Solving Navier-Stokes", total=num_steps)
#     for i in range(num_steps):
#         progress.update(1)
#         # Update current time step
#         t += dt

#         # Compute tentative velocity step
#         b1 = fe.assemble(L1)
#         [bc.apply(A1, b1) for bc in bcu]
#         solver1.solve(A1, u_.vector(), b1)

#         # Pressure correction
#         b2 = fe.assemble(L2)
#         [bc.apply(A2, b2) for bc in bcp]
#         solver2.solve(A2, p_.vector(), b2)

#         # Velocity correction
#         b3 = fe.assemble(L3)
#         [bc.apply(A3, b3) for bc in bcu]
#         solver3.solve(A3, u_.vector(), b3)

#         # Move to next time step
#         u_n.assign(u_)
#         p_n.assign(p_)

#         velocity_file.write(u_, t)
#         # pressure_file.write(p_, t)

#         max_u.append(u_.vector().max())
#         np.savetxt(results_folder + "3D_case_max_u.txt", np.array(max_u))

#     return u_, p_


def run_simple_sim():
    """
    Runs FESTIM coupled to the fluid dynamics solver
    on a mesh with a single fluid domain

    Returns:
        festim.Simulation: the simulation object
    """

    # TODO find a way to make this more generic
    fluid_id = 6
    inlet_id = 7
    outlet_id = 8
    vacuum_id = 9
    walls_id = 10
    my_sim = F.Simulation()

    my_sim.mesh = F.MeshFromXDMF(
        volume_file="mesh_cells.xdmf", boundary_file="mesh_facets.xdmf"
    )

    u = "velocity_field.xdmf"

    # u, p = fluid_dynamics_sim_chorin(
    #     mesh=my_sim.mesh.mesh,
    #     volume_markers=my_sim.mesh.volume_markers,
    #     surface_markers=my_sim.mesh.surface_markers,
    #     id_inlet=inlet_id,
    #     id_outlet=outlet_id,
    #     id_walls=[walls_id, vacuum_id],
    # )

    D_LiPb = htm.diffusivities.filter(material=htm.LIPB).mean()
    print(D_LiPb)
    # factor = 100  # TODO remove this by improving mesh
    my_sim.materials = F.Material(
        id=fluid_id,
        D_0=D_LiPb.pre_exp.magnitude,
        E_D=D_LiPb.act_energy.magnitude,
    )

    inlet_BC = F.DirichletBC(surfaces=inlet_id, value=1e15, field="solute")
    recombination_on_vacuum = F.DirichletBC(surfaces=vacuum_id, value=0, field="solute")
    my_sim.boundary_conditions = [
        inlet_BC,
        recombination_on_vacuum,
    ]

    my_sim.T = F.Temperature(700)

    my_sim.settings = F.Settings(1e-10, 1e-10, transient=False)

    my_sim.exports = [F.XDMFExport(field="solute")]
    my_sim.log_level = 20

    my_sim.initialise()

    hydrogen_concentration = my_sim.h_transport_problem.mobile.mobile_concentration()
    test_function_solute = my_sim.h_transport_problem.mobile.test_function

    advection_term = (
        fe.inner(fe.dot(fe.grad(hydrogen_concentration), u), test_function_solute)
        * my_sim.mesh.dx
    )

    my_sim.h_transport_problem.F += advection_term
    print(my_sim.h_transport_problem.F)
    my_sim.run()
    return my_sim


if __name__ == "__main__":
    run_simple_sim()

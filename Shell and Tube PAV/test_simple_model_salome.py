import festim as F
from fenics import *
import fenics as fe
import h_transport_materials as htm


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

    fenics_mesh = my_sim.mesh.mesh
    V = fe.VectorFunctionSpace(fenics_mesh, "CG", 2)
    u = fe.Function(V)
    fe.XDMFFile("velocity_field.xdmf").read_checkpoint(u, "velocity", -1)

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

    average_inlet = F.AverageSurface("solute", inlet_id)
    average_outlet = F.AverageSurface("solute", outlet_id)
    extracted_flux = F.SurfaceFlux("solute", vacuum_id)
    derived_quantities = F.DerivedQuantities(
        [average_inlet, average_outlet, extracted_flux]
    )

    my_sim.exports = [F.XDMFExport(field="solute"), derived_quantities]
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

    extraction_efficiency = (
        average_inlet.data[-1] - average_outlet.data[-1]
    ) / average_inlet.data[-1]
    print(f"Extraction efficiency: {extraction_efficiency:.2%}")
    print(f"Extracted flux: {extracted_flux.data[-1]:.2e} T/s")

    return my_sim


if __name__ == "__main__":
    run_simple_sim()

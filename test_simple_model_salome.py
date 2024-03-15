import festim as F


def run_simple_sim():

    # TODO find a way to make this more generic
    fluid_id = 6
    inlet_id = 7
    outlet_id = 8
    my_sim = F.Simulation()

    my_sim.mesh = F.MeshFromXDMF(
        volume_file="mesh_cells.xdmf", boundary_file="mesh_facets.xdmf"
    )

    my_sim.materials = F.Material(id=fluid_id, D_0=2, E_D=0)

    inlet_BC = F.DirichletBC(surfaces=inlet_id, value=1, field="solute")
    outlet_BC = F.DirichletBC(surfaces=outlet_id, value=0, field="solute")
    my_sim.boundary_conditions = [
        inlet_BC,
        outlet_BC,
    ]

    my_sim.T = F.Temperature(500)

    my_sim.settings = F.Settings(1e-10, 1e-10, transient=False)

    my_sim.exports = [F.XDMFExport(field="solute")]

    my_sim.initialise()
    my_sim.run()


if __name__ == "__main__":
    run_simple_sim()

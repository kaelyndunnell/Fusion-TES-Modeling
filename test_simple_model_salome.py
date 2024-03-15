import festim as F
from fenics import *

from convert_mesh import convert_mesh


def run_simple_sim():
    correspondance_dict = convert_mesh("mesh.med")

    fluid_id = correspondance_dict["fluid"]
    inlet_id = correspondance_dict["inlet"]
    outlet_id = correspondance_dict["outlet"]
    walls_id = correspondance_dict["walls"]

    my_sim = F.Simulation()

    my_sim.mesh = F.MeshFromXDMF(
        volume_file="mesh_cells.xdmf", boundary_file="mesh_facets.xdmf"
    )

    my_sim.materials = F.Material(id=fluid_id, D_0=1, E_D=0)

    inlet_BC = F.DirichletBC(surfaces=inlet_id, value=1, field="solute")
    outlet_BC = F.DirichletBC(surfaces=outlet_id, value=0, field="solute")
    my_sim.boundary_conditions = [inlet_BC, outlet_BC]

    my_sim.T = F.Temperature(500)

    my_sim.settings = F.Settings(1e-10, 1e-10, transient=False)

    my_sim.exports = [F.XDMFExport(field="solute")]

    my_sim.initialise()
    my_sim.run()


if __name__ == "__main__":
    run_simple_sim()

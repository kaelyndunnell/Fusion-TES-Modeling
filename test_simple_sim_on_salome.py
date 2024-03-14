import festim as F
from fenics import *
from convert_mesh import convert_med_to_xdmf
import numpy as np


def run_simple_sim():
    correspondance_dict = convert_med_to_xdmf("Mesh_1.med")

    print(correspondance_dict)


    fluid_id = 6#np.abs(correspondance_dict["fluid"])
    inlet_id = 8#np.abs(correspondance_dict["inlet"])
    outlet_id = 7#np.abs(correspondance_dict["outlet"])
    walls_id = 9#np.abs(correspondance_dict["walls"])

    my_sim = F.Simulation()

    my_sim.mesh = F.MeshFromXDMF(volume_file='mesh_cells.xdmf', boundary_file='mesh_facets.xdmf')

    my_sim.materials = F.Material(id=fluid_id, D_0=1, E_D=0)

    inlet_BC = F.DirichletBC(surfaces=walls_id, value=1, field="solute")
    outlet_BC = F.DirichletBC(surfaces=outlet_id, value=1, field="solute")
    my_sim.boundary_conditions = [
        inlet_BC,
        outlet_BC,
    ]

    my_sim.T = F.Temperature(500)

    my_sim.settings = F.Settings(1e-10, 1e-10, transient=False)

    my_sim.exports = [F.XDMFExport(field="solute")]
    my_sim.log_level = 20
    my_sim.initialise()
    my_sim.run()

if __name__=="__main__":
    run_simple_sim()
import meshio
import numpy as np
import argparse
from dolfinx.io import XDMFFile
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI


### SCRIPT TO CONVERT VOLUME SINGLE ELEMENT TET MESHES PRODUCED WITH CUBIT TO .gmsh MESHES FOR OPENFOAM OR FESTIM ###


def convert_mesh(openfoam: bool):
    """converts an exodus mesh produced from CUBIT to a gmsh .msh file for use
    in openfoam or FESTIM."""

    if openfoam:
        # openfoam one vol path
        mesh_file_path = "meshing/CAD/cubit_files/inlet_tank.e"
        write_path = "meshing/inlet_tank.msh"
    else:
        # festim two vols paths
        mesh_file_path = "meshing/CAD/cubit_files/test_two_vols.e"
        write_path = "meshing/test_festim.msh"

    mesh = meshio.read(mesh_file_path)
    print("cell types:", [b.type for b in mesh.cells])

    #### MESH CONVERSION ####

    points = mesh.points

    tet_keys = {"tetra", "TETRA3"}
    tri_keys = {"triangle", "TRIANGLE3"}

    # assigning names to each block idk from cubit
    if openfoam:
        physical_groups = {  # ids start from 0 here but 1 in cubit
            0: (1, "fluid", 3),
            1: (2, "tank_inlet", 2),
            2: (3, "tank_outlet", 2),
            3: (4, "tank_walls", 2),
        }
    else:
        physical_groups = {  # ids start from 0 here but 1 in cubit
            0: (1, "fluid", 3),
            1: (2, "membrane", 3),
            2: (3, "inlet", 2),
            3: (4, "outlet", 2),
            4: (5, "walls", 2),
            5: (6, "vacuum", 2),
            6: (7, "interface", 2),
        }

    # SURFACES
    all_cells = []  # output cell data
    all_cells_phys = []  # tagging for gmsh physical groups (also geometry)
    all_points = []

    for block_idx, cell_block in enumerate(
        mesh.cells
    ):  # block_idx is equal to block id in cubit

        ctype = cell_block.type.lower()
        if ctype not in tri_keys:
            continue

        tag, name, _ = physical_groups[block_idx]
        data = cell_block.data
        n = len(data)

        # assign tag based on physical groups
        phys = np.full(n, tag, dtype=np.int32)
        unique_phys = np.unique(phys)

        all_cells.append(("triangle", data))
        all_cells_phys.append(phys)

    # VOLUMES
    all_tets = []  # same logic as above

    for block_idx, cell_block in enumerate(mesh.cells):
        ctype = cell_block.type.lower()
        data = cell_block.data
        n = len(data)

        if ctype in tet_keys:
            phys = np.full(n, block_idx + 1, dtype=np.int32)

            all_tets.append(data)
            all_cells.append(("tetra", data))
            all_cells_phys.append(phys)

    combined_tets = np.vstack(all_tets)

    print(f"total tets in mesh is {len(combined_tets)}.")

    field_data = {
        name: np.array([tag, dim]) for _, (tag, name, dim) in physical_groups.items()
    }

    point_data = {"gmsh:dim_tags": np.array(all_points)}

    ### WRITE MESH ###
    out_mesh = meshio.Mesh(
        points=points,
        cells=all_cells,
        cell_data={
            "gmsh:physical": all_cells_phys,
            "gmsh:geometrical": all_cells_phys,
        },
        field_data=field_data,
    )

    print(f"writing mesh to {write_path}.")

    meshio.write(
        write_path, out_mesh, file_format="gmsh22", binary=False
    )  # has to be gmsh22 for both (gmsh v4 has errors with meshio)
    print("finished.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Determines which program to convert mesh for (OpenFOAM or FESTIM)."
    )
    parser.add_argument(
        "-program",
        type=str,
        help="Program to convert mesh for, openfoam for OpenFOAM or festim for FESTIM.",
    )

    args = parser.parse_args()
    program = args.program

    if program == "openfoam":
        openfoam = True
    else:
        openfoam = False

    convert_mesh(openfoam=openfoam)

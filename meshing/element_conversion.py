import sys
import meshio
import numpy as np
import gmsh

mesh_file_path = "meshing/CAD/mixed_mesh.e"
write_path = "meshing/tes_openfoam.msh"

mesh = meshio.read(mesh_file_path)
print("cell types:", [b.type for b in mesh.cells])

QUAD_TO_TRIS = np.array(  # conversion array to make quads --> tris (2 tris per quad)
    [
        [0, 1, 2],
        [0, 2, 3],
    ],
    dtype=np.int64,
)

WEDGE_TO_TETS = np.array(  # conversion array to make wedges --> tets (3 tets per wedge)
    [
        [0, 1, 2, 5],
        [0, 1, 5, 4],
        [0, 4, 5, 3],
    ],
    dtype=np.int64,
)


def quad_to_tris(quad_conn: np.ndarray) -> np.ndarray:
    """converts (n,4) quads to (n*2, 3) tris."""
    n = quad_conn.shape[0]
    tris = np.empty((n * 2, 3), dtype=np.int64)
    for i, local in enumerate(QUAD_TO_TRIS):
        tris[i::2] = quad_conn[:, local]
    return tris


def wedge_cells_to_tets(wedge_conn: np.ndarray) -> np.ndarray:
    """converts a (n,6) wedge to an (n*3, 4) tet."""
    n = wedge_conn.shape[0]
    tets = np.empty((n * 3, 4), dtype=np.int64)
    for i, local in enumerate(WEDGE_TO_TETS):
        tets[i::3] = wedge_conn[:, local]
    return tets


#### MESH CONVERSION ####

points = mesh.points

wedge_keys = {"wedge", "wedge6", "penta6"}
tet_keys = {"tetra", "tetra4", "tet4"}
tri_keys = {"triangle", "triangle3", "tri3"}
surf_keys = {"triangle", "triangle3", "tri3", "quad"}

# assigning names to each block idk from cubit
physical_groups = {
    1: (1, "fluid", 3),
    2: (2, "inlet", 2),
    3: (2, "inlet", 2),
    4: (3, "outlet", 2),
    5: (3, "outlet", 2),
    6: (4, "walls", 2),
}

# SURFACES
all_tris = []  # output cell data
all_tris_phys = []  # tagging for gmsh physical groups (also geometry)

for block_idx, cell_block in enumerate(
    mesh.cells
):  # block_idx is equal to block id in cubit

    ctype = cell_block.type.lower()
    if ctype not in surf_keys:
        continue

    print(ctype, block_idx)

    tag, name, _ = physical_groups[block_idx]
    data = cell_block.data
    n = len(data)

    # assign tag based on physical groups
    phys = np.full(n, tag, dtype=np.int32)
    unique_phys = np.unique(phys)

    if ctype in tri_keys:
        print(f"found {n} {ctype} in block {block_idx}, {name}.")
        all_tris.append(("triangle", data))
        all_tris_phys.append(phys)

    elif "quad" in ctype:
        converted = quad_to_tris(data)
        print(
            f"found {n} quads --> converted to {len(converted)} triangles in block {block_idx}."
        )
        phys_rep = np.repeat(phys, 2)
        all_tris.append(("triangle", converted))
        all_tris_phys.append(phys_rep)

# VOLUMES
all_tets = []  # same logic as above
all_tet_phys = []

for block_idx, cell_block in enumerate(mesh.cells):
    ctype = cell_block.type.lower()
    data = cell_block.data
    n = len(data)

    def get_vol_tag():  # wedges are given id 0 and tets are given id 1, manually setting all vol ids as 1 as in cubit
        if block_idx == 0:
            return np.full(len(data), block_idx + 1, dtype=np.int32)
        else:
            return np.full(len(data), block_idx, dtype=np.int32)

    if ctype in tet_keys:
        phys = get_vol_tag()
        print(f"found {n} tets in block {block_idx}.")
        all_tets.append(data)
        all_tet_phys.append(phys)

    elif ctype in wedge_keys:
        phys = get_vol_tag()
        converted = wedge_cells_to_tets(data)
        phys_rep = np.repeat(
            phys, 3
        )  # convert each wedge to 3 tets, have to repeat this in these objects
        print(
            f"found {len(data)} wedges --> converted to {len(converted)} tets (in block {phys[0]})."
        )
        all_tets.append(converted)
        all_tet_phys.append(phys_rep)

combined_tets = np.vstack(all_tets)
combined_tet_phys = np.concatenate(all_tet_phys)

print(f"total tets in new mesh is {len(combined_tets)}.")

all_tris.append(("tetra", combined_tets))
all_tris_phys.append(combined_tet_phys)

field_data = {
    name: np.array([tag, dim]) for _, (tag, name, dim) in physical_groups.items()
}

### WRITE MESH ###
out_mesh = meshio.Mesh(
    points=points,
    cells=all_tris,
    cell_data={
        "gmsh:physical": all_tris_phys,
        "gmsh:geometrical": all_tris_phys,
    },
    field_data=field_data,
)

print("new mesh cell types:", [b.type for b in out_mesh.cells])

print(f"writing mesh to {write_path}.")
meshio.write(write_path, out_mesh, file_format="gmsh22", binary=False)
print("finished.")

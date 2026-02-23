import festim as F
from foam2dolfinx import OpenFOAMReader
from dolfinx import fem
from festim.helpers import nmm_interpolate
from dolfinx.io import VTXWriter
from mpi4py import MPI
import numpy as np
from dolfinx.mesh import meshtags, exterior_facet_indices
from scipy.spatial import cKDTree


def read_openfoam_data(file_name, final_time):
    """
    Read OpenFOAM data from a file and return the pressure, velocity, and viscosity (if it exists) fields.
    """
    print("Reading OpenFOAM data...")
    openfoam_reader = OpenFOAMReader(filename=file_name, cell_type=10)
    p = openfoam_reader.create_dolfinx_function(t=final_time, name="p")
    u = openfoam_reader.create_dolfinx_function(t=final_time, name="U")
    mesh = openfoam_reader.dolfinx_meshes_dict["default"]

    try:
        # read turbulent viscosity if it exists
        nut = openfoam_reader.create_dolfinx_function(t=final_time, name="nut")
        print("'nut' field found and read.")
    except Exception:
        nut = None
        print("No 'nut' field found.")

    facet_meshtags, volume_meshtags = define_meshtags(openfoam_reader)

    return p, u, mesh, nut, facet_meshtags, volume_meshtags


def export_openfoam_data(p, u):
    """
    Export OpenFOAM data to VTX files.
    """

    print("Exporting OpenFOAM data")
    writer_p = VTXWriter(
        MPI.COMM_WORLD,
        "OpenFOAM/pressure.bp",
        p,
        "BP5",
    )
    writer_u = VTXWriter(
        MPI.COMM_WORLD,
        "OpenFOAM/velocity.bp",
        u,
        "BP5",
    )

    writer_p.write(t=0)
    writer_u.write(t=0)


def tag_boundary_patch(dolfinx_mesh, patch_dataset, patch_id, tol=1e-6):
    fdim = dolfinx_mesh.topology.dim - 1
    dolfinx_mesh.topology.create_connectivity(fdim, 0)
    dolfinx_mesh.topology.create_connectivity(0, fdim)
    dolfinx_mesh.topology.create_connectivity(fdim, dolfinx_mesh.topology.dim)

    facet_indices = exterior_facet_indices(dolfinx_mesh.topology)
    x = dolfinx_mesh.geometry.x
    patch_points = patch_dataset.points
    tree = cKDTree(x)
    matched_vertex_indices = tree.query_ball_point(patch_points, tol)
    matched_vertex_indices = list(set(i for sub in matched_vertex_indices for i in sub))

    matched_facets = []
    for facet in facet_indices:
        vertices = dolfinx_mesh.topology.connectivity(fdim, 0).links(facet)
        if all(v in matched_vertex_indices for v in vertices):
            matched_facets.append(facet)

    # print(f"Tagging {len(matched_facets)} facets for patch ID {patch_id}")
    return np.array(matched_facets, dtype=np.int32), np.full(
        len(matched_facets), patch_id, dtype=np.int32
    )


def define_meshtags(cfd_reader):
    OF_multiblock = cfd_reader.reader.read()
    boundary = OF_multiblock["boundary"]
    mesh = cfd_reader.dolfinx_meshes_dict["default"]

    all_facets = np.array([], dtype=np.int32)
    all_tags = np.array([], dtype=np.int32)

    for i, name in enumerate(boundary.keys()):
        if name not in ["inlet", "outlet", "wall", "probe"]:
            continue
        facets, tags = tag_boundary_patch(mesh, boundary[name], i + 1)
        all_facets = np.concatenate([all_facets, facets])
        all_tags = np.concatenate([all_tags, tags])

        # print(f"Tagging {len(facets)} facets for patch {name} with ID {i + 1}")

    facet_tags = meshtags(
        mesh,
        mesh.topology.dim - 1,
        all_facets,
        all_tags,
    )

    num_cells = mesh.topology.index_map(mesh.topology.dim).size_local
    mesh_cell_indices = np.arange(num_cells, dtype=np.int32)
    tags_volumes = np.full(num_cells, 1, dtype=np.int32)
    volume_meshtags = meshtags(mesh, mesh.topology.dim, mesh_cell_indices, tags_volumes)

    num_cells = (
        mesh.topology.index_map(mesh.topology.dim).size_local
        + mesh.topology.index_map(mesh.topology.dim).num_ghosts
    )

    print(f"number of cells in mesh: {num_cells}")
    return (
        facet_tags,
        volume_meshtags,
    )


if __name__ == "__main__":
    # read openfoam data
    p, u = read_openfoam_data("OpenFOAM/probe-case/case.foam", final_time=100)

    export_openfoam_data(p, u)

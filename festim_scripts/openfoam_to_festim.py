from foam2dolfinx import OpenFOAMReader
from dolfinx.io import VTXWriter
from mpi4py import MPI
import io4dolfinx


def read_openfoam_data(file_name):
    """
    Read OpenFOAM data from a file and return the pressure, velocity, and viscosity (if it exists) fields.
    """
    print("Reading OpenFOAM data...")
    openfoam_reader = OpenFOAMReader(filename=file_name, cell_type=12)

    final_time = max(openfoam_reader.times)

    # p = openfoam_reader.create_dolfinx_function_with_cell_data(t=final_time, name="p")
    u = openfoam_reader.create_dolfinx_function_with_cell_data(t=final_time, name="U")
    mesh = openfoam_reader.dolfinx_meshes_dict["default"]

    try:
        # read turbulent viscosity if it exists
        nut = openfoam_reader.create_dolfinx_function_with_cell_data(
            t=final_time, name="nut"
        )
        print("'nut' field found and read.")
    except Exception:
        nut = None
        print("No 'nut' field found.")

    facet_meshtags = openfoam_reader.create_facet_meshtags()
    volume_meshtags = openfoam_reader.create_cell_meshtags()

    checkpoint_file = "test_cp.bp"
    io4dolfinx.write_mesh(checkpoint_file, mesh)

    io4dolfinx.write_function(checkpoint_file, u, time=0.0, name="U")
    io4dolfinx.write_function(checkpoint_file, nut, time=0.0, name="nut")

    io4dolfinx.write_meshtags(
        checkpoint_file, mesh, facet_meshtags, meshtag_name="facet_tags"
    )
    io4dolfinx.write_meshtags(
        checkpoint_file, mesh, volume_meshtags, meshtag_name="cell_tags"
    )

    return u, mesh, nut, facet_meshtags, volume_meshtags


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


if __name__ == "__main__":
    # read openfoam data
    p, u, mesh, nut, facet_meshtags, volume_meshtags = read_openfoam_data(
        "OpenFOAM/probe-case/case.foam"
    )

    export_openfoam_data(p, u)

import festim as F  # using festim2
import numpy as np
from scifem import assemble_scalar
from dolfinx import fem
from dolfinx.io import VTXWriter, XDMFFile
import ufl
from dolfinx import cpp as _cpp
from openfoam_to_festim import read_openfoam_data
from dolfinx.log import set_log_level, LogLevel
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI
from basix.ufl import element
import h_transport_materials as htm
from dolfinx.io import XDMFFile
from mpi4py import MPI


def evaluate_stabalisation_term(mesh, u, delta):
    """See more at https://www.comsol.com/blogs/understanding-stabilization-methods"""

    # evaluate Cell size
    tdim = mesh.topology.dim
    num_cells = mesh.topology.index_map(tdim).size_local
    cells = np.arange(num_cells, dtype=np.int32)
    mesh_ = _cpp.mesh.Mesh_float64(
        mesh.comm, mesh.topology._cpp_object, mesh.geometry._cpp_object
    )
    h = _cpp.mesh.h(mesh_, tdim, cells)
    V0 = fem.functionspace(mesh, ("DG", 0))
    h_as_function = fem.Function(V0)
    h_as_function.x.array[:] = h

    # Compute magnitude of velocity
    v_mag = ufl.sqrt(ufl.dot(u, u))

    D_art = delta * v_mag * h_as_function

    return D_art


def build_festim_model(
    openfoam_data_file,
    openfoam_final_time,
    festim_mesh_file,
    breeder_temperature,
    results_folder,
    visualize_fields=True,
):

    # READ OPENFOAM MESH
    p, u, openfoam_mesh, nut, facet_meshtags, volume_meshtags = read_openfoam_data(
        openfoam_data_file, final_time=openfoam_final_time
    )

    # READ GMSH MESH
    model_rank = 0
    festim_mesh_data = gmshio.read_from_msh(
        festim_mesh_file, MPI.COMM_WORLD, model_rank, gdim=3
    )
    festim_mesh = festim_mesh_data.mesh

    # DEFINE & INITIALIZE MODEL

    print("Building FESTIM model...")

    my_model = F.HydrogenTransportProblemDiscontinuous()

    my_model.mesh = F.Mesh(festim_mesh)
    my_model.facet_meshtags = festim_mesh_data.facet_tags
    my_model.volume_meshtags = festim_mesh_data.cell_tags

    # interpolate OpenFOAM velocity field onto FESTIM mesh
    el = element(
        "Lagrange",
        festim_mesh.topology.cell_name(),
        1,
        shape=(festim_mesh.geometry.dim,),
    )
    V_openfoam = fem.functionspace(openfoam_mesh, el)
    V_festim = fem.functionspace(festim_mesh, el)

    u_openfoam = fem.Function(V_openfoam)
    u_openfoam.interpolate(
        u
    )  # u is a fem.function.Function! paraview visualization of u_openfoam is accurate with this formulation
    festim_velocity = fem.Function(V_festim)

    festim_cells = festim_mesh_data.cell_tags.find(1)  # breeder cells to interpolate to

    interpolation_data = fem.create_interpolation_data(
        V_to=V_festim, V_from=V_openfoam, cells=festim_cells
    )

    festim_velocity.interpolate_nonmatching(
        u_openfoam, cells=festim_cells, interpolation_data=interpolation_data
    )

    # interpolate OpenFOAM nut field onto FESTIM mesh
    N_openfoam = fem.functionspace(openfoam_mesh, ("CG", 1))
    N_festim = fem.functionspace(festim_mesh, ("CG", 1))

    nut_openfoam = fem.Function(N_openfoam)
    nut_openfoam.interpolate(nut)
    festim_nut = fem.Function(N_festim)

    interpolation_data = fem.create_interpolation_data(
        V_to=N_festim, V_from=N_openfoam, cells=festim_cells
    )

    festim_nut.interpolate_nonmatching(
        nut_openfoam, cells=festim_cells, interpolation_data=interpolation_data
    )

    nut_field_array = festim_nut.x.array
    nut_field_array[nut_field_array < 0.0] = 0.0  # ensure no negative eddy viscosity
    festim_nut.x.array[:] = nut_field_array

    # LIPB

    lipb_diffusivity = (
        htm.diffusivities.filter(material=htm.LIPB)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )

    D_0_PbLi = lipb_diffusivity.pre_exp.magnitude  # m2/s,
    E_D_PbLi = lipb_diffusivity.act_energy.magnitude  # eV

    D_fick = D_0_PbLi * ufl.exp(-E_D_PbLi / (F.k_B * breeder_temperature))

    # add stabilization term for diffusion
    D_art = evaluate_stabalisation_term(mesh=festim_mesh, u=festim_velocity, delta=1)

    # add turbulent diffusion term
    Sc = 0.7
    D_turb = festim_nut / Sc

    D_expr = D_fick + D_turb + D_art
    V = fem.functionspace(festim_mesh, ("CG", 1))
    D_pbli = fem.Function(V)
    D_pbli.interpolate(fem.Expression(D_expr, V.element.interpolation_points))

    if visualize_fields:
        my_writer_v_festim = VTXWriter(
            MPI.COMM_WORLD,
            results_folder + "/velocity_field.bp",
            festim_velocity,
            "BP5",
        )
        my_writer_v_festim.write(t=0)
        my_writer_D_field = VTXWriter(
            MPI.COMM_WORLD, results_folder + "/D_field.bp", D_pbli, "BP5"
        )
        my_writer_D_field.write(t=0)
        my_writer_nut = VTXWriter(
            MPI.COMM_WORLD, results_folder + "/festim_nut.bp", festim_nut, "BP5"
        )
        my_writer_nut.write(t=0)
        my_writer_v_openfoam = VTXWriter(
            MPI.COMM_WORLD, results_folder + "/openfoam_velocity.bp", u_openfoam, "BP5"
        )
        my_writer_v_openfoam.write(t=0)

    lipb_solubility = (
        htm.solubilities.filter(material=htm.LIPB)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )

    # STEEL

    steel_diffusivity = (
        htm.diffusivities.filter(material=htm.STEEL_316L)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )

    steel_solubility = (
        htm.solubilities.filter(material=htm.STEEL_316L)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )

    # MATERIALS

    breeder_material = F.Material(
        D=D_pbli,
        K_S_0=lipb_solubility.pre_exp.magnitude,
        E_K_S=lipb_solubility.act_energy.magnitude,
    )

    membrane_material = F.Material(
        D_0=steel_diffusivity.pre_exp.magnitude,
        E_D=steel_diffusivity.act_energy.magnitude,
        K_S_0=steel_solubility.pre_exp.magnitude,
        E_K_S=steel_solubility.act_energy.magnitude,
    )

    # SET DOMAINS

    # use same tags as gmsh markers
    breeder_marker = 1
    membrane_marker = 2

    breeder = F.VolumeSubdomain(id=breeder_marker, material=breeder_material)
    membrane = F.VolumeSubdomain(id=membrane_marker, material=membrane_material)

    inlet_marker = 3
    outlet_marker = 4
    walls_marker = 5
    vacuum_marker = 6

    inlet = F.SurfaceSubdomain(id=inlet_marker)
    outlet = F.SurfaceSubdomain(id=outlet_marker)
    walls = F.SurfaceSubdomain(id=walls_marker)
    vacuum = F.SurfaceSubdomain(id=vacuum_marker)

    my_model.subdomains = [
        breeder,
        membrane,
        inlet,
        outlet,
        walls,
        vacuum,
    ]

    my_model.surface_to_volume = (
        {  # anything defined for BC needs to be here (and exports)
            inlet: breeder,
            outlet: breeder,
            vacuum: membrane,
        }
    )

    T = F.Species("T", subdomains=my_model.volume_subdomains)
    my_model.species = [T]

    interface_marker = 7

    my_model.method_interface = F.InterfaceMethod.penalty
    my_model.interfaces = [
        F.Interface(
            id=interface_marker,
            subdomains=[breeder, membrane],
            penalty_term=1e25,
        ),
    ]

    # SET TEMP AND BOUNDARY CONDITIONS

    advection = F.AdvectionTerm(velocity=festim_velocity, subdomain=breeder, species=T)
    my_model.advection_terms = [advection]

    my_model.temperature = breeder_temperature  # K

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=inlet, value=1e25, species=T),
        F.FixedConcentrationBC(subdomain=vacuum, value=0, species=T),
    ]

    # SETTINGS

    my_model.settings = F.Settings(
        atol=1e-10,
        rtol=1e-10,
        transient=False,
    )

    # EXPORTS

    concentration_field_breeder = F.VTXSpeciesExport(
        filename=f"{results_folder}/T_breeder.bp", field=T, subdomain=breeder
    )
    concentration_field_probe = F.VTXSpeciesExport(
        filename=f"{results_folder}/T_membrane.bp", field=T, subdomain=membrane
    )
    c_out = F.TotalSurface(
        field=T, surface=outlet, filename=f"{results_folder}/c_out.csv"
    )
    c_in = F.TotalSurface(field=T, surface=inlet, filename=f"{results_folder}/c_in.csv")

    permeation_flux = F.SurfaceFlux(
        field=T, surface=vacuum, filename=f"{results_folder}/permeation_flux.csv"
    )

    my_model.exports = [
        c_out,
        c_in,
        concentration_field_breeder,
        concentration_field_probe,
        permeation_flux,
    ]

    return my_model


if __name__ == "__main__":

    my_model = build_festim_model(
        openfoam_data_file="openfoam/lipb_simple/tes.foam",
        openfoam_final_time=5600,
        festim_mesh_file="meshing/test_festim.msh",
        breeder_temperature=603.15,
        results_folder="lipb_festim_results/benchmark",
    )

    # INITIALISE AND RUN
    my_model.initialise()
    set_log_level(LogLevel.INFO)
    my_model.run()

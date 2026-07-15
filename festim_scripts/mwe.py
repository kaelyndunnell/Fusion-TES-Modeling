import festim as F  # using festim2
import numpy as np
from dolfinx import fem
from scifem import assemble_scalar
from dolfinx.io import VTXWriter
import ufl
from dolfinx import cpp as _cpp
from openfoam_to_festim import read_openfoam_data
from dolfinx.log import set_log_level, LogLevel
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI
from basix.ufl import element
import h_transport_materials as htm
import os
import re
import meshio
from Nb_recombination import nb_recomb
from festim.helpers import nmm_interpolate

N_A = 6.0221e23  # atms/mol


class SurfaceAdvectionFlux(F.SurfaceFlux):
    """Computes the advection flux of a field on a given surface

    Args:
        field (festim.Species): species for which the surface flux is computed
        surface (festim.SurfaceSubdomain1D): surface subdomain
        filename (str, optional): name of the file to which the surface flux is exported

    Attributes:
        see `festim.SurfaceFlux`
    """

    def __init__(self, field, surface, filename, velocity_field):

        super().__init__(field=field, surface=surface, filename=filename)
        self.velocity_field = velocity_field

    @property
    def title(self):
        return f"{self.field.name} advection flux surface {self.surface.id}"

    def compute(self, u, ds: ufl.Measure, entity_maps=None):
        if isinstance(u, ufl.indexed.Indexed):
            mesh = self.field.sub_function_space.mesh
        else:
            mesh = u.function_space.mesh

        n = ufl.FacetNormal(mesh)

        # Diffusive flux
        surface_flux = assemble_scalar(
            fem.form(
                -self.D * ufl.dot(ufl.grad(u), n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )

        vel_local = self.velocity_field

        advective_flux = assemble_scalar(
            fem.form(
                u * ufl.inner(vel_local, n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )

        self.value = surface_flux + advective_flux
        self.data.append(self.value)


def convert_med_to_xdmf(
    med_file,
    cell_file="mesh_domains.xdmf",
    facet_file="mesh_boundaries.xdmf",
    cell_type="hexahedron",
    facet_type="quad",
):
    """Converts a MED mesh to XDMF
    Args:
        med_file (str): the name of the MED file
        cell_file (str, optional): the name of the file containing the
            volume markers. Defaults to "mesh_domains.xdmf".
        facet_file (str, optional): the name of the file containing the
            surface markers.. Defaults to "mesh_boundaries.xdmf".
        cell_type (str, optional): The topology of the cells. Defaults to "tetra".
        facet_type (str, optional): The topology of the facets. Defaults to "triangle".
    Returns:
        dict, dict: the correspondance dict, the cell types
    """
    msh = meshio.read(med_file)

    correspondance_dict = {-k: v for k, v in msh.cell_tags.items()}

    cell_data_types = msh.cell_data_dict["cell_tags"].keys()

    for mesh_block in msh.cells:
        if mesh_block.type == cell_type:
            meshio.write_points_cells(
                cell_file,
                msh.points,
                [mesh_block],
                cell_data={"f": [-1 * msh.cell_data_dict["cell_tags"][cell_type]]},
            )
        elif mesh_block.type == facet_type:
            meshio.write_points_cells(
                facet_file,
                msh.points,
                [mesh_block],
                cell_data={"f": [-1 * msh.cell_data_dict["cell_tags"][facet_type]]},
            )

    return correspondance_dict, cell_data_types


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
    breeder,
    openfoam_data_folder,
    festim_mesh_file,
    results_folder,
    penalty_term,
    visualize_fields=True,
):

    if os.path.isdir(results_folder):
        print(f"Results folder: {results_folder}")
    else:
        os.mkdir(results_folder)
        print(f"Results folder: {results_folder}")

    # READ OPENFOAM MESH
    p, openfoam_velocity, openfoam_mesh, nut, facet_meshtags, volume_meshtags = (
        read_openfoam_data(openfoam_data_folder + "/tes.foam")
    )

    nut.x.array[nut.x.array < 0.0] = 0.0  # ensure no negative eddy viscosity

    # READ FESTIM MESH
    correspondance_dict, cell_data_types = convert_med_to_xdmf(
        festim_mesh_file,
        cell_file=f"{results_folder}/mesh_domains.xdmf",
        facet_file=f"{results_folder}/mesh_boundaries.xdmf",
    )

    for index, label in correspondance_dict.items():
        print(f"{index}: {label[0]}")

    festim_mesh = F.MeshFromXDMF(
        volume_file=f"{results_folder}/mesh_domains.xdmf",
        facet_file=f"{results_folder}/mesh_boundaries.xdmf",
    )

    # DEFINE & INITIALIZE MODEL

    print("Building FESTIM model...")

    my_model = F.HydrogenTransportProblemDiscontinuous()

    my_model.mesh = festim_mesh
    mesh = my_model.mesh.mesh
    my_model.facet_meshtags = festim_mesh.define_surface_meshtags()
    my_model.volume_meshtags = festim_mesh.define_volume_meshtags()

    # interpolate OpenFOAM velocity field onto FESTIM mesh
    el = element(
        "Lagrange",
        mesh.topology.cell_name(),
        1,
        shape=(mesh.geometry.dim,),
    )
    V_openfoam = fem.functionspace(openfoam_mesh, el)
    V_festim = fem.functionspace(mesh, el)

    u_openfoam = fem.Function(V_openfoam)
    u_openfoam.interpolate(openfoam_velocity)
    festim_velocity = fem.Function(V_festim)

    festim_cells = my_model.volume_meshtags.find(6)  # breeder cells to interpolate to

    interpolation_data = fem.create_interpolation_data(
        V_to=V_festim, V_from=V_openfoam, cells=festim_cells
    )

    festim_velocity.interpolate_nonmatching(
        u_openfoam, cells=festim_cells, interpolation_data=interpolation_data
    )

    # interpolate OpenFOAM nut field onto FESTIM mesh
    N_openfoam = fem.functionspace(openfoam_mesh, ("CG", 1))
    N_festim = fem.functionspace(mesh, ("CG", 1))

    nut_openfoam = fem.Function(N_openfoam)
    nut_openfoam.interpolate(nut)
    festim_nut = fem.Function(N_festim)

    interpolation_data = fem.create_interpolation_data(
        V_to=N_festim, V_from=N_openfoam, cells=festim_cells
    )

    festim_nut.interpolate_nonmatching(
        nut_openfoam, cells=festim_cells, interpolation_data=interpolation_data
    )

    # BREEDER MATERIAL

    breeder_diffusivity = (
        htm.diffusivities.filter(material=htm.LIPB)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )
    breeder_solubility = (
        htm.solubilities.filter(material=htm.LIPB)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )
    breeder_solubility_law = "SIEVERT"
    breeder_temperature = 723.15  # K from https://doi.org/10.3390/en16073022

    # membrane material (Nb for LiPb)

    membrane_diffusivity = (
        htm.diffusivities.filter(material=htm.NIOBIUM)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")[0]
    )

    membrane_solubility = htm.solubilities.filter(material=htm.NIOBIUM)[0]
    membrane_solubility_law = "SIEVERT"

    u = htm.ureg
    membrane_recombo = nb_recomb

    my_model.method_interface = F.InterfaceMethod.penalty

    D_0_breeder = breeder_diffusivity.pre_exp.magnitude  # m2/s,
    E_D_breeder = breeder_diffusivity.act_energy.magnitude  # eV

    D_fick = D_0_breeder * ufl.exp(-E_D_breeder / (F.k_B * breeder_temperature))

    # add turbulent diffusion term
    Sc = 0.7
    D_turb = festim_nut / Sc

    # add stabilization term for diffusion
    D_art = evaluate_stabalisation_term(mesh=mesh, u=festim_velocity, delta=0.1)

    D_expr = D_fick + D_turb  # + D_art
    V = fem.functionspace(mesh, ("CG", 1))
    D_tot = fem.Function(V)
    D_tot.interpolate(fem.Expression(D_expr, V.element.interpolation_points))

    if visualize_fields:
        my_writer_v_festim = VTXWriter(
            MPI.COMM_WORLD,
            results_folder + "/velocity_field.bp",
            openfoam_velocity,
            "BP5",
        )
        my_writer_v_festim.write(t=0)
        my_writer_v_festim = VTXWriter(
            MPI.COMM_WORLD,
            results_folder + "/festim_velocity_field.bp",
            festim_velocity,
            "BP5",
        )
        my_writer_v_festim.write(t=0)
        my_writer_D_field = VTXWriter(
            MPI.COMM_WORLD, results_folder + "/D_field.bp", D_tot, "BP5"
        )
        my_writer_D_field.write(t=0)
        my_writer_nut = VTXWriter(
            MPI.COMM_WORLD, results_folder + "/festim_nut.bp", festim_nut, "BP5"
        )
        my_writer_nut.write(t=0)

    # MATERIALS

    breeder_material = F.Material(
        D=D_tot,
        K_S_0=breeder_solubility.pre_exp.to("mol/m**3/(Pa**0.5)").magnitude,
        E_K_S=breeder_solubility.act_energy.magnitude,
        solubility_law=breeder_solubility_law,
    )

    membrane_material = F.Material(
        D_0=membrane_diffusivity.pre_exp.magnitude,
        E_D=membrane_diffusivity.act_energy.magnitude,
        K_S_0=membrane_solubility.pre_exp.to("mol/m**3/(Pa**0.5)").magnitude,
        E_K_S=membrane_solubility.act_energy.magnitude,
        solubility_law=membrane_solubility_law,
    )

    # my_model.method_interface = F.InterfaceMethod.penalty

    # SET DOMAINS

    # use same tags as xdmf markers
    breeder_marker = 6
    membrane_marker = 7

    breeder = F.VolumeSubdomain(id=breeder_marker, material=breeder_material)
    membrane = F.VolumeSubdomain(id=membrane_marker, material=membrane_material)

    inlet_marker = 8
    outlet_marker = 9
    walls_marker = 10
    vacuum_marker = 11

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

    my_model.surface_to_volume = {  # anything defined for BC needs to be here (and exports)
        inlet: breeder,
        outlet: breeder,
        vacuum: membrane,
        walls: membrane,
    }

    T = F.Species("T", subdomains=my_model.volume_subdomains)
    my_model.species = [T]

    interface_marker = 12

    my_model.interfaces = [
        F.Interface(
            id=interface_marker,
            subdomains=[breeder, membrane],
            penalty_term=penalty_term,
        ),
    ]

    # SET TEMP AND BOUNDARY CONDITIONS

    advection = F.AdvectionTerm(velocity=festim_velocity, subdomain=breeder, species=T)
    my_model.advection_terms = [advection]

    my_model.temperature = 723  # K

    vacuum_TT_rxn = F.SurfaceReactionBC(
        reactant=[T, T],
        gas_pressure=0.0,  # residual pressure from surrounding pipes
        k_r0=(membrane_recombo.pre_exp.to("m**4/(mol*s)")).magnitude,
        E_kr=membrane_recombo.act_energy.magnitude,
        k_d0=0,
        E_kd=0,
        subdomain=vacuum,
    )

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=inlet, value=1, species=T),
        # F.FixedConcentrationBC(subdomain=vacuum, value=0, species=T),
        vacuum_TT_rxn,
    ]

    # SETTINGS
    my_model.settings = F.Settings(
        atol=1e-20, rtol=1e-9, transient=False, max_iterations=100
    )

    # EXPORTS

    inlet_flux = SurfaceAdvectionFlux(
        field=T,
        surface=inlet,
        velocity_field=festim_velocity,
        filename=f"{results_folder}/inlet_flux.csv",
    )
    outlet_flux = SurfaceAdvectionFlux(
        field=T,
        surface=outlet,
        velocity_field=festim_velocity,
        filename=f"{results_folder}/outlet_flux.csv",
    )

    concentration_field_breeder = F.VTXSpeciesExport(
        filename=f"{results_folder}/T_breeder.bp", field=T, subdomain=breeder
    )
    concentration_field_membrane = F.VTXSpeciesExport(
        filename=f"{results_folder}/T_membrane.bp", field=T, subdomain=membrane
    )

    permeation_flux = F.SurfaceFlux(
        field=T, surface=vacuum, filename=f"{results_folder}/permeation_flux.csv"
    )
    c_in = F.AverageSurface(
        field=T, surface=inlet, filename=f"{results_folder}/c_in.csv"
    )

    c_out = F.AverageSurface(
        field=T, surface=outlet, filename=f"{results_folder}/c_out.csv"
    )

    my_model.exports = [
        permeation_flux,
        inlet_flux,
        outlet_flux,
        concentration_field_breeder,
        concentration_field_membrane,
        c_in,
        c_out,
    ]

    return my_model


if __name__ == "__main__":
    to_run = {
        1.00e0: [1.93, 0.14, 1.71],
    }

    for c_in, lst in to_run.items():
        for penalty_term in [1e7]:
            print(penalty_term)

            my_model = build_festim_model(
                breeder="lipb",
                openfoam_data_folder=f"openfoam/velocity_parametrization/lipb/pipe_{lst[0]:.2f}m_s_r{lst[1]:.2f}_l{lst[2]:.2f}",
                festim_mesh_file=f"meshing/parametric_meshes/two_vol_0.0046_{lst[1]:.2f}_{lst[2]:.2f}.med",
                results_folder=f"lipb_festim_results/test/in_{c_in}_vel_{lst[0]:.2f}_r{lst[1]:.2f}_l{lst[2]:.2f}",
                penalty_term=penalty_term,
            )

            # INITIALISE AND RUN
            my_model.initialise()
            set_log_level(LogLevel.INFO)
            my_model.run()
            print(f"inlet flux is {my_model.exports[1].value}")
            print(f"outlet flux is {my_model.exports[2].value}")
            print(f"permeation flux is {my_model.exports[0].value}")
            # if abs(my_model.exports[1].value) >= my_model.exports[2].value + my_model.exports[0].value:
            #     print('mass conserved')
            #     if my_model.exports[0].value < 0:
            #         print('flux negative')
            #     else:
            #         break
            # else:
            #     print('mass not conserved')

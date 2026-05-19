import festim as F  # using festim2
import numpy as np
from dolfinx import fem
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

N_A = 6.0221e23  # atms/mol


def findDir(basePath):
    names = []
    for root, dirs, files in os.walk(basePath):
        for name in dirs:
            try:
                if re.match(r"\d", name):
                    names.append(int(name))
            except:
                continue
    return np.max(names)


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


def build_festim_model(
    c_inlet,
    residual_pressure,
    breeder,
    openfoam_data_folder,
    festim_mesh_file,
    results_folder,
    visualize_fields=True,
    penalty_term=1e20,
):

    # READ OPENFOAM MESH
    print("Reading OpenFOAM data...")
    openfoam_final_time = findDir(openfoam_data_folder)

    p, openfoam_velocity, openfoam_mesh, nut, facet_meshtags, volume_meshtags = (
        read_openfoam_data(
            openfoam_data_folder + "/tes.foam", final_time=openfoam_final_time
        )
    )
    print("OpenFOAM mesh read successfully.")

    nut.x.array[nut.x.array < 0.0] = 0.0  # ensure no negative eddy viscosity

    if os.path.isdir(results_folder):
        print(f"Results folder: {results_folder}")
    else:
        os.mkdir(results_folder)
        print(f"Results folder: {results_folder}")

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

    # interpolate OpenFOAM nut field onto FESTIM mesh
    festim_cells = my_model.volume_meshtags.find(6)  # breeder cells to interpolate to
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
    breeder_temperature = 603.15  # K

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

    # penalty_term = 1e43

    if c_inlet < 1e20 / N_A:
        penalty_term = penalty_term / 100

    D_0_breeder = breeder_diffusivity.pre_exp.magnitude  # m2/s,
    E_D_breeder = breeder_diffusivity.act_energy.magnitude  # eV

    D_fick = D_0_breeder * ufl.exp(-E_D_breeder / (F.k_B * breeder_temperature))

    # add turbulent diffusion term
    Sc = 0.7
    D_turb = festim_nut / Sc

    D_expr = D_fick + D_turb
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
        K_S_0=breeder_solubility.pre_exp.magnitude / N_A,
        E_K_S=breeder_solubility.act_energy.magnitude,
        solubility_law=breeder_solubility_law,
    )

    membrane_material = F.Material(
        D_0=membrane_diffusivity.pre_exp.magnitude,
        E_D=membrane_diffusivity.act_energy.magnitude,
        K_S_0=membrane_solubility.pre_exp.magnitude / N_A,
        E_K_S=membrane_solubility.act_energy.magnitude,
        solubility_law=membrane_solubility_law,
    )

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

    advection = F.AdvectionTerm(
        velocity=openfoam_velocity, subdomain=breeder, species=T
    )
    my_model.advection_terms = [advection]

    my_model.temperature = breeder_temperature  # K

    vacuum_TT_rxn = F.SurfaceReactionBC(
        reactant=[T, T],
        gas_pressure=residual_pressure,  # residual pressure from surrounding pipes
        k_r0=membrane_recombo.pre_exp.magnitude * N_A,
        E_kr=membrane_recombo.act_energy.magnitude,
        k_d0=0,
        E_kd=0,
        subdomain=vacuum,
    )

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=inlet, value=c_inlet, species=T),
        # F.FixedConcentrationBC(subdomain=vacuum, value=0, species=T),
        vacuum_TT_rxn,
    ]

    # SETTINGS

    my_model.settings = F.Settings(
        atol=1e-10, rtol=1e-10, transient=False, max_iterations=100
    )

    # EXPORTS

    concentration_field_breeder = F.VTXSpeciesExport(
        filename=f"{results_folder}/T_breeder.bp", field=T, subdomain=breeder
    )
    concentration_field_membrane = F.VTXSpeciesExport(
        filename=f"{results_folder}/T_membrane.bp", field=T, subdomain=membrane
    )
    c_in = F.AverageSurface(
        field=T, surface=inlet, filename=f"{results_folder}/c_in.csv"
    )

    c_out = F.AverageSurface(
        field=T, surface=outlet, filename=f"{results_folder}/c_out.csv"
    )
    permeation_flux = F.SurfaceFlux(
        field=T, surface=vacuum, filename=f"{results_folder}/permeation_flux.csv"
    )

    my_model.exports = [
        c_out,
        permeation_flux,
        c_in,
        concentration_field_breeder,
        concentration_field_membrane,
    ]

    return my_model


if __name__ == "__main__":
    for c_in in [1e20 / N_A]:
        for v_in in [0.05]:
            my_model = build_festim_model(
                c_inlet=c_in,
                residual_pressure=0,
                breeder="lipb",
                openfoam_data_folder=f"openfoam/velocity_parametrization/lipb_new/pipe_{v_in:.2f}m_s_r0.10_l1.50",  # test one
                festim_mesh_file="meshing/parametric_meshes/two_vol_0.0046_0.10_1.50.med",
                results_folder=f"lipb_festim_results/with_meshes/in_{c_in}_vel_{v_in}_r0.10_l1.50_fixedBC",
                penalty_term=1e-3,
            )

            # INITIALISE AND RUN
            my_model.initialise()
            set_log_level(LogLevel.INFO)
            my_model.run()

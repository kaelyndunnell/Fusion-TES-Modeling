import festim as F
from dolfinx.mesh import create_unit_square, CellType
from mpi4py import MPI
import numpy as np
import h_transport_materials as htm
from dolfinx.io import gmsh as gmshio
from dolfinx.log import set_log_level, LogLevel

N_A = 6.022e23  # atms/mol


# def build_festim_model(c_inlet, penalty_term):
#     fenics_mesh = create_unit_square(
#         MPI.COMM_WORLD, 10, 10, cell_type=CellType.quadrilateral
#     )
#     festim_mesh = F.Mesh(fenics_mesh)

#     my_model = F.HydrogenTransportProblemDiscontinuous()

#     breeder_diffusivity = (
#         htm.diffusivities.filter(material=htm.FLIBE)
#         .filter(exclude=True, isotope="H")
#         .filter(exclude=True, isotope="D")
#         .mean()
#     )
#     breeder_solubility = (
#         htm.solubilities.filter(material=htm.FLIBE)
#         .filter(exclude=True, isotope="H")
#         .filter(exclude=True, isotope="D")
#         .mean()
#     )
#     breeder_solubility_law = "HENRY"
#     breeder_temperature = 900  # K

#     membrane_diffusivity = (
#         htm.diffusivities.filter(material=htm.STEEL_316L)
#         .filter(exclude=True, isotope="H")
#         .filter(exclude=True, isotope="D")
#         .mean()
#     )
#     membrane_solubility = (
#         htm.solubilities.filter(material=htm.STEEL_316L)
#         .filter(exclude=True, isotope="H")
#         .filter(exclude=True, isotope="D")
#         .mean()
#     )
#     membrane_solubility_law = "SIEVERT"

#     membrane_recombo = htm.recombination_coeffs.filter(material=htm.STEEL_316L).mean()
#     membrane_diss = htm.dissociation_coeffs.filter(material=htm.STEEL_316L).mean()

#     material_top = F.Material(
#         D_0=breeder_diffusivity.pre_exp.magnitude / N_A,
#         E_D=breeder_diffusivity.act_energy.magnitude,
#         K_S_0=breeder_solubility.pre_exp.magnitude / N_A,
#         E_K_S=breeder_solubility.act_energy.magnitude,
#         solubility_law=breeder_solubility_law,
#     )
#     material_bottom = F.Material(
#         D_0=membrane_diffusivity.pre_exp.magnitude / N_A,
#         E_D=membrane_diffusivity.act_energy.magnitude,
#         K_S_0=membrane_solubility.pre_exp.magnitude / N_A,
#         E_K_S=membrane_solubility.act_energy.magnitude,
#         solubility_law=membrane_solubility_law,
#     )

#     top_volume = F.VolumeSubdomain(
#         id=3, material=material_top, locator=lambda x: x[1] >= 0.5
#     )
#     bottom_volume = F.VolumeSubdomain(
#         id=4, material=material_bottom, locator=lambda x: x[1] <= 0.5
#     )

#     top_surface = F.SurfaceSubdomain(id=1, locator=lambda x: np.isclose(x[1], 1.0))
#     bottom_surface = F.SurfaceSubdomain(id=2, locator=lambda x: np.isclose(x[1], 0.0))

#     my_model.mesh = festim_mesh
#     my_model.subdomains = [top_surface, bottom_surface, top_volume, bottom_volume]

#     my_model.surface_to_volume = {
#         top_surface: top_volume,
#         bottom_surface: bottom_volume,
#     }

#     my_model.interfaces = [
#         F.Interface(5, (bottom_volume, top_volume), penalty_term=penalty_term)
#     ]

#     H = F.Species("H")
#     my_model.species = [H]
#     H.subdomains = [top_volume, bottom_volume]

#     my_model.temperature = breeder_temperature

#     vacuum_TT_rxn = F.SurfaceReactionBC(
#         reactant=[H, H],
#         gas_pressure=0.0,  # residual pressure from surrounding pipes
#         k_r0=membrane_recombo.pre_exp.magnitude / N_A,
#         E_kr=membrane_recombo.act_energy.magnitude,
#         k_d0=membrane_diss.pre_exp.magnitude / N_A,
#         E_kd=membrane_diss.act_energy.magnitude,
#         subdomain=bottom_surface,
#     )

#     my_model.boundary_conditions = [
#         F.FixedConcentrationBC(subdomain=top_surface, value=c_inlet, species=H),
#         # vacuum_TT_rxn,
#         F.FixedConcentrationBC(subdomain=bottom_surface, value=0.0, species=H),
#     ]

#     my_model.settings = F.Settings(
#         atol=1e-10, rtol=1e-10, transient=False, max_iterations=100
#     )

#     flux_in = F.SurfaceFlux(surface=top_surface, field=H)
#     flux_out = F.SurfaceFlux(surface=bottom_surface, field=H)
#     concentration_field_top = F.VTXSpeciesExport(
#         filename=f"test/T_breeder.bp", field=H, subdomain=top_volume
#     )
#     concentration_field_bottom = F.VTXSpeciesExport(
#         filename=f"test/T_membrane.bp", field=H, subdomain=bottom_volume
#     )

#     my_model.exports = [
#         flux_in,
#         flux_out,
#         concentration_field_bottom,
#         concentration_field_top,
#     ]

#     return my_model


def build_festim_model(c_inlet, penalty_term):

    # READ GMSH MESH
    model_rank = 0
    festim_mesh_data = gmshio.read_from_msh(
        "meshing/test_festim.msh", MPI.COMM_WORLD, model_rank, gdim=3
    )
    festim_mesh = festim_mesh_data.mesh

    # DEFINE & INITIALIZE MODEL

    print("Building FESTIM model...")

    my_model = F.HydrogenTransportProblemDiscontinuous()

    my_model.mesh = F.Mesh(festim_mesh)
    my_model.facet_meshtags = festim_mesh_data.facet_tags
    my_model.volume_meshtags = festim_mesh_data.cell_tags

    breeder_diffusivity = (
        htm.diffusivities.filter(material=htm.FLIBE)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )
    breeder_solubility = (
        htm.solubilities.filter(material=htm.FLIBE)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )
    breeder_solubility_law = "HENRY"
    breeder_temperature = 900  # K

    membrane_diffusivity = (
        htm.diffusivities.filter(material=htm.STEEL_316L)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )
    membrane_solubility = (
        htm.solubilities.filter(material=htm.STEEL_316L)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )
    membrane_solubility_law = "SIEVERT"

    membrane_recombo = htm.recombination_coeffs.filter(material=htm.STEEL_316L).mean()
    membrane_diss = htm.dissociation_coeffs.filter(material=htm.STEEL_316L).mean()

    breeder_material = F.Material(
        D_0=breeder_diffusivity.pre_exp.magnitude,
        E_D=breeder_diffusivity.act_energy.magnitude,
        K_S_0=breeder_solubility.pre_exp.magnitude,
        E_K_S=breeder_solubility.act_energy.magnitude,
        solubility_law=breeder_solubility_law,
    )
    membrane_material = F.Material(
        D_0=membrane_diffusivity.pre_exp.magnitude,
        E_D=membrane_diffusivity.act_energy.magnitude,
        K_S_0=membrane_solubility.pre_exp.magnitude,
        E_K_S=membrane_solubility.act_energy.magnitude,
        solubility_law=membrane_solubility_law,
    )

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

    H = F.Species("H", subdomains=my_model.volume_subdomains)
    my_model.species = [H]

    interface_marker = 7

    my_model.interfaces = [
        F.Interface(
            id=interface_marker,
            subdomains=[breeder, membrane],
            penalty_term=penalty_term,
        ),
    ]

    my_model.temperature = breeder_temperature

    vacuum_TT_rxn = F.SurfaceReactionBC(
        reactant=[H, H],
        gas_pressure=0.0,  # residual pressure from surrounding pipes
        k_r0=membrane_recombo.pre_exp.magnitude,
        E_kr=membrane_recombo.act_energy.magnitude,
        k_d0=membrane_diss.pre_exp.magnitude,
        E_kd=membrane_diss.act_energy.magnitude,
        subdomain=vacuum,
    )

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=inlet, value=c_inlet, species=H),
        vacuum_TT_rxn,
        # F.FixedConcentrationBC(subdomain=bottom_surface, value=0.0, species=H),
    ]

    my_model.settings = F.Settings(
        atol=1e10, rtol=1e-10, transient=False, max_iterations=100
    )

    flux_in = F.SurfaceFlux(surface=inlet, field=H)
    flux_out = F.SurfaceFlux(surface=vacuum, field=H)
    concentration_field_top = F.VTXSpeciesExport(
        filename=f"test/T_breeder.bp", field=H, subdomain=breeder
    )
    concentration_field_bottom = F.VTXSpeciesExport(
        filename=f"test/T_membrane.bp", field=H, subdomain=membrane
    )

    my_model.exports = [
        flux_in,
        flux_out,
        concentration_field_bottom,
        concentration_field_top,
    ]

    return my_model


def check_no_negative_concentrations(model):
    overall_min = np.inf
    worst_subdomain = None

    for species in model.species:
        for subdomain, u in species.subdomain_to_post_processing_solution.items():
            if not hasattr(u, "x"):
                continue
            local_min = float(u.x.array.min())
            if local_min < overall_min:
                overall_min = local_min
                worst_subdomain = getattr(subdomain, "id", str(subdomain))

    passed = overall_min >= 0.0
    return passed, overall_min, worst_subdomain


def evaluate_penalty(penalty_term, c_inlet=1e20):
    model = build_festim_model(c_inlet=c_inlet, penalty_term=penalty_term)
    model.initialise()
    model.run()

    flux_in = model.exports[0].data[-1] if model.exports[0].data else None
    flux_out = model.exports[1].data[-1] if model.exports[1].data else None

    conc_passed, min_conc, worst_subdomain = check_no_negative_concentrations(model)

    return flux_in, flux_out, conc_passed, min_conc, worst_subdomain


def signs_are_correct(flux_in, flux_out):
    return flux_in is not None and flux_out is not None and flux_in < 0 and flux_out > 0


if __name__ == "__main__":
    set_log_level(LogLevel.WARNING)

    penalty_candidates = np.logspace(15, 35, 20)

    found = False
    for penalty in penalty_candidates:
        print(f"Running penalty_term = {penalty:.1e}")
        try:
            flux_in, flux_out, conc_passed, min_conc, worst_subdomain = (
                evaluate_penalty(penalty)
            )
            print(f"flux_in  = {flux_in}")
            print(f"flux_out = {flux_out}")
            print(
                f"concentration check passed: {conc_passed}  (min value = {min_conc:.3e}.')"
            )

            if signs_are_correct(flux_in, flux_out):
                print(f"passed!!! penalty_term = {penalty:.6e} solves:")
                print(f"flux_in  = {flux_in}  (negative)")
                print(f"flux_out = {flux_out}  (positive)")
                print(f"all concentrations non-negative (min = {min_conc:.3e})")
                found = True
                break
        except Exception as e:
            print(f"failed")

    if not found:
        print("no passing penalty terms found.")

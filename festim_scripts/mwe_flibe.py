import festim as F
from dolfinx.mesh import create_unit_square, CellType
from mpi4py import MPI
import numpy as np
import h_transport_materials as htm
from dolfinx.io import gmsh as gmshio
from dolfinx.log import set_log_level, LogLevel

N_A = 6.022e23  # atms/mol


def build_festim_model(c_inlet, penalty_term):
    fenics_mesh = create_unit_square(
        MPI.COMM_WORLD, 10, 10, cell_type=CellType.quadrilateral
    )
    festim_mesh = F.Mesh(fenics_mesh)

    my_model = F.HydrogenTransportProblemDiscontinuous()

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

    material_top = F.Material(
        D_0=breeder_diffusivity.pre_exp.magnitude / N_A,
        E_D=breeder_diffusivity.act_energy.magnitude,
        K_S_0=breeder_solubility.pre_exp.magnitude / N_A,
        E_K_S=breeder_solubility.act_energy.magnitude,
        solubility_law=breeder_solubility_law,
    )
    material_bottom = F.Material(
        D_0=membrane_diffusivity.pre_exp.magnitude / N_A,
        E_D=membrane_diffusivity.act_energy.magnitude,
        K_S_0=membrane_solubility.pre_exp.magnitude / N_A,
        E_K_S=membrane_solubility.act_energy.magnitude,
        solubility_law=membrane_solubility_law,
    )

    top_volume = F.VolumeSubdomain(
        id=3, material=material_top, locator=lambda x: x[1] >= 0.5
    )
    bottom_volume = F.VolumeSubdomain(
        id=4, material=material_bottom, locator=lambda x: x[1] <= 0.5
    )

    top_surface = F.SurfaceSubdomain(id=1, locator=lambda x: np.isclose(x[1], 1.0))
    bottom_surface = F.SurfaceSubdomain(id=2, locator=lambda x: np.isclose(x[1], 0.0))

    my_model.mesh = festim_mesh
    my_model.subdomains = [top_surface, bottom_surface, top_volume, bottom_volume]

    my_model.surface_to_volume = {
        top_surface: top_volume,
        bottom_surface: bottom_volume,
    }

    my_model.interfaces = [
        F.Interface(5, (bottom_volume, top_volume), penalty_term=penalty_term)
    ]

    H = F.Species("H")
    my_model.species = [H]
    H.subdomains = [top_volume, bottom_volume]

    my_model.temperature = breeder_temperature

    vacuum_TT_rxn = F.SurfaceReactionBC(
        reactant=[H, H],
        gas_pressure=0.0,  # residual pressure from surrounding pipes
        k_r0=membrane_recombo.pre_exp.magnitude / N_A,
        E_kr=membrane_recombo.act_energy.magnitude,
        k_d0=membrane_diss.pre_exp.magnitude / N_A,
        E_kd=membrane_diss.act_energy.magnitude,
        subdomain=bottom_surface,
    )

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=top_surface, value=c_inlet, species=H),
        # vacuum_TT_rxn,
        F.FixedConcentrationBC(subdomain=bottom_surface, value=0.0, species=H),
    ]

    my_model.settings = F.Settings(
        atol=1e-10, rtol=1e-10, transient=False, max_iterations=100
    )

    flux_in = F.SurfaceFlux(surface=top_surface, field=H)
    flux_out = F.SurfaceFlux(surface=bottom_surface, field=H)
    concentration_field_top = F.VTXSpeciesExport(
        filename=f"testing/T_breeder.bp", field=H, subdomain=top_volume
    )
    concentration_field_bottom = F.VTXSpeciesExport(
        filename=f"testing/T_membrane.bp", field=H, subdomain=bottom_volume
    )

    my_model.exports = [
        flux_in,
        flux_out,
        concentration_field_bottom,
        concentration_field_top,
    ]

    return my_model


if __name__ == "__main__":
    c_in = 1e20 / N_A
    my_model = build_festim_model(c_inlet=c_in, penalty_term=1e-7)

    my_model.initialise()
    set_log_level(LogLevel.INFO)
    my_model.run()

    print(f"flux in is {my_model.exports[0].data}")
    print(f"flux out is {my_model.exports[1].data}")

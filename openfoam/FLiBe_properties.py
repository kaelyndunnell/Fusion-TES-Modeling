from fluid_parameters import (
    calculate_reynolds_number,
    plot_reynolds_number_vs_inlet_velocity,
    calculate_initial_k,
    calculate_initial_epsilon,
    calculate_initial_omega,
)
import numpy as np
import festim as F
import h_transport_materials as htm


def calculate_FLiBe_kinematic_viscosity(
    breeder_temperature,
    breeder_density,
    breeder,
    suppress_print=False,
):
    """Calculate the kinematic viscosity of FLiBe at a given temperature and density.
    Used for OpenFOAM simulation.

    Parameters
    ----------
    breeder_temperature : float
        Breeder temperature in K.
    breeder_density : float
        Breeder density in kg/m3.
     breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Kinematic viscosity in m2/s.
    """
    breeder_dynamic_viscosity = 0.000116 * np.exp(
        3755 / breeder_temperature
    )  # Pa.s = kg s / m s2; equation from Williams 2006

    kinematic_viscosity = breeder_dynamic_viscosity / breeder_density  # m2/s

    if not suppress_print:
        print(
            f"Kinematic viscosity of {breeder} at {breeder_temperature}K is {kinematic_viscosity}m2/s."
        )

    return kinematic_viscosity


breeder = "FLiBe"

breeder_temperature = 900  # K from Meschini 2021
FLiBe_density = 2245 - 0.424 * (
    breeder_temperature - 273.15
)  # kg/m3 ; equation from Vidrio 2022

# flow_rate = 815  # kg/s, https://doi.org/10.1016/j.fusengdes.2024.114261
inlet_diameter = 0.13  # m from CAD

k_b = F.k_B  # eV/K, boltzmann constant

flibe_diffusivity = htm.diffusivities.filter(material=htm.FLIBE).mean()
E_D = flibe_diffusivity.act_energy.magnitude  # eV
D_0 = flibe_diffusivity.pre_exp.magnitude  # m2/s

FLiBe_diffusivity = D_0 * np.exp(-E_D / (k_b * breeder_temperature))  # m2/s

inlet_velocity = 3  # m/s

kinematic_viscosity = calculate_FLiBe_kinematic_viscosity(
    breeder_temperature, FLiBe_density, breeder
)

Re = calculate_reynolds_number(
    inlet_velocity, inlet_diameter, kinematic_viscosity, breeder
)

k = calculate_initial_k(inlet_velocity)
epsilon = calculate_initial_epsilon(k, characteristic_length=inlet_diameter)
omega = calculate_initial_omega(k, inlet_diameter)

print(f"Initial turbulence kinetic energy for {breeder}: {k} m2/s2")
print(f"Initial turbulence dissipation rate for {breeder}: {epsilon} m2/s3")
print(f"Initial specific dissipation rate for {breeder}: {omega} 1/s")

# plot_reynolds_number_vs_inlet_velocity(
#     inlet_diameter,
#     kinematic_viscosity,
#     breeder_temperature,
#     breeder,
#     inlet_velocity,
#     flibe=True,
# )

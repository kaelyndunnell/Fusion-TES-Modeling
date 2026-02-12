from fluid_parameters import (
    calculate_inlet_velocity,
    calculate_reynolds_number,
    plot_reynolds_number_vs_inlet_velocity,
    calculate_initial_k,
    calculate_initial_epsilon,
    calculate_initial_omega,
    calculate_schmidt_number,
)
import numpy as np
import festim as F


def calculate_LiPb_kinematic_viscosity(
    breeder_temperature,
    breeder_density,
    breeder,
    suppress_print=False,
):
    """Calculate the kinematic viscosity of LiPb at a given temperature and density.
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
    breeder_dynamic_viscosity = (
        0.0061091
        - 2.2574e-5 * breeder_temperature
        + 3.766e-8 * breeder_temperature**2
        - 2.2887e-11 * breeder_temperature**3
    )  # Pa.s = kg s / m s2; equation from Martelli 2019

    kinematic_viscosity = breeder_dynamic_viscosity / breeder_density  # m2/s

    if not suppress_print:
        print(
            f"Kinematic viscosity of {breeder} at {breeder_temperature}K is {kinematic_viscosity}m2/s."
        )

    return kinematic_viscosity


breeder = "LiPb"

inlet_diameter = 0.13  # m from CAD

breeder_temperature = 603.15  # K from Utili 2023
LiPb_density = (
    10520.35 - 1.19051 * breeder_temperature
)  # kg/m3 ; equation from Martelli 2019

k_b = F.k_B  # eV/K, boltzmann constant
E_D = 19500 * 1.0364e-5  # = 0.202098
LiPb_diffusivity = 4.03e-8 * np.exp(
    -E_D / (k_b * breeder_temperature)
)  # m2/s ; from Utili 2023, 1 J/mol = 1.0364E-5eV

inlet_velocity = 0.5  # m/s, https://doi.org/10.1016/j.jnucmat.2022.153985

kinematic_viscosity = calculate_LiPb_kinematic_viscosity(
    breeder_temperature, LiPb_density, breeder, suppress_print=True
)

Re = calculate_reynolds_number(
    inlet_velocity, inlet_diameter, kinematic_viscosity, breeder, suppress_print=True
)

print(f"Reynolds number for {breeder} at {inlet_velocity}m/s is {Re}.")

Sc = calculate_schmidt_number(
    kinematic_viscosity=kinematic_viscosity,
    diffusivity=LiPb_diffusivity,
    breeder=breeder,
)

k = calculate_initial_k(inlet_velocity)
epsilon = calculate_initial_epsilon(k, characteristic_length=inlet_diameter)
omega = calculate_initial_omega(k, inlet_diameter)

print(f"Density of {breeder} at {breeder_temperature}K is {LiPb_density}kg/m3.")
print(f"Initial turbulence kinetic energy for {breeder}: {k} m2/s2")
print(f"Initial turbulence dissipation rate for {breeder}: {epsilon} m2/s3")
print(f"Initial specific dissipation rate for {breeder}: {omega} 1/s")

plot_reynolds_number_vs_inlet_velocity(
    inlet_diameter, kinematic_viscosity, breeder_temperature, breeder, inlet_velocity
)

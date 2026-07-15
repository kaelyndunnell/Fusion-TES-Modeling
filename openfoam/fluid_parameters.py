import numpy as np
import matplotlib.pyplot as plt


def calculate_inlet_velocity(
    flow_rate,
    inlet_diameter,
    breeder_density,
    breeder,
    suppress_print=False,
):
    """Calculate the inlet velocity of fluid breeder at a given flow rate, inlet diameter, breeder density, and temperature.
    Used for OpenFOAM simulation.

    Parameters
    ----------
    flow_rate : float
        Flow rate in kg/s.
    inlet_diameter : float
        Inlet diameter in m.
    breder_density : float
        Breeder fluid density in kg/m3.
    breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Inlet velocity in m/s.
    """

    inlet_area = np.pi * (inlet_diameter / 2) ** 2  # m^2

    inlet_velocity = flow_rate * breeder_density ** (-1) * inlet_area ** (-1)  # m/s
    if not suppress_print:
        print(
            f"Inlet velocity for {breeder} flow rate of {flow_rate}kg/s is {inlet_velocity}m/s."
        )

    return inlet_velocity


def calculate_reynolds_number(
    inlet_velocity,
    characteristic_length,
    kinematic_viscosity,
    breeder,
    suppress_print=False,
):
    """Calculate the reynolds number of a fluid breeder given its inlet velocity, characteristic length, and kinematic viscosity.
    Used to determine if turbulence is needed for OpenFOAM simulation.

    Parameters
    ----------
    inlet_velocity : float
        Inlet velocity in m/s.
    characteristic_length : float
        Characteristic length in m.
    kinematic_viscosity : float
        Kinematic viscosity in m2/s.
     breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Reynolds number (dimensionless).
    """

    reynolds_number = (inlet_velocity * characteristic_length) / kinematic_viscosity

    if not suppress_print:
        print(f"Reynolds number for {breeder} is {reynolds_number}.")

        if reynolds_number > 3500:
            print(f"Flow is turbulent.")
        else:
            print(f"Flow is laminar.")

    return reynolds_number


def calculate_schmidt_number(
    kinematic_viscosity,
    diffusivity,
    breeder,
    suppress_print=False,
):
    """Calculate the schmidt number of a fluid breeder given its kinematic viscosity and mass diffusivity.
    Used for turbulent diffusion term in FESTIM.

    Parameters
    ----------
    kinematic_viscosity : float
        Kinematic viscosity in m2/s.
    diffusivity : float
        Fluid diffusivity in m2/s.
    breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Schmidt number (dimensionless).
    """

    schmidt_number = kinematic_viscosity / diffusivity

    if not suppress_print:
        print(f"Schmidt number for {breeder} is {schmidt_number}.")

    return schmidt_number


def plot_reynolds_number_vs_inlet_velocity(
    characteristic_length,
    kinematic_viscosity,
    breeder_temperature,
    breeder,
    model_velocity,
    flibe=False,
):
    """Plot the reynolds number of a fluid breeder with varying inlet velocities.
    Assumes constant characteristic length and kinematic viscosity.

    Parameters
    ----------
    characteristic_length : float
        Characteristic length in m.
    kinematic_viscosity : float
        Kinematic viscosity in m2/s.
    breeder_temperature : float
        Breeder temperature in K.
     breeder : str
        Breeder fluid name.
    """
    if flibe:
        inlet_velocities = np.linspace(0, 0.5, 100000)  # m/s

    else:
        inlet_velocities = np.linspace(0, 1e-2, 100000)  # m/s
    Re_numbers = []

    for inlet_velocity in inlet_velocities:  # m/s
        Re_numbers.append(
            calculate_reynolds_number(
                inlet_velocity,
                characteristic_length,
                kinematic_viscosity,
                breeder,
                suppress_print=True,
            )
        )

    plt.plot(inlet_velocities, Re_numbers, "b-")
    plt.axhline(y=3500, color="r", linestyle="--", label="Turbulence Threshold")
    plt.axvline(x=model_velocity, color="g", label="Model Inlet Velocity")
    plt.xlabel("Inlet Velocity (m/s)")
    plt.ylabel("Reynolds Number")
    plt.title(
        f"Reynolds Number vs Inlet Velocity for {breeder} at {breeder_temperature}K."
    )
    plt.legend()
    plt.show()


# reference: https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras-k-epsilon.html


def calculate_initial_k(inlet_velocity):
    """Calculate initial turbulence kinetic energy.

    Parameters
    ----------
    inlet_velocity : float
        Inlet velocity in m/s.

    Returns
    -------
    float
        Initial turbulence kinetic energy in m2/s2.
    """

    # assume initial turbulence is isotropic so that U'2_x = U'2_y = U'2_z
    return (3 / 2) * (0.05 * inlet_velocity) ** 2  # 5% of inlet velocity


def calculate_initial_epsilon(k, characteristic_length):
    """Calculate initial turbulence dissipation rate.

    Parameters
    ----------
    k : float
        Initial turbulence kinetic energy in m2/s2.
    characteristic_length : float
        Characteristic length in m.

    Returns
    -------
    float
        Initial turbulence dissipation rate in m2/s3.
    """
    return 0.09 ** (3 / 4) * k ** (3 / 2) / characteristic_length


def calculate_initial_omega(k, characteristic_length):
    """Calculate initial specific dissipation rate.

    Parameters
    ----------
    k : float
        Initial turbulence kinetic energy in m2/s2.
    characteristic_length : float
        Characteristic length in m.

    Returns
    -------
    float
        Initial specific dissipation rate in 1/s.
    """
    return np.sqrt(k) / (0.09 ** (1 / 4) * characteristic_length)

from fluid_parameters import (
    calculate_reynolds_number,
    calculate_initial_k,
    calculate_initial_omega,
)
from LiPb_properties import calculate_LiPb_kinematic_viscosity
from FLiBe_properties import calculate_FLiBe_kinematic_viscosity
import numpy as np
import festim as F
import os
import shutil
import h_transport_materials as htm


def change_variable_in_openfoam_file(filename, old_value, new_value):
    with open(filename, "r", encoding="utf-8") as f:
        content = f.read()

    content = content.replace(str(old_value), str(new_value))

    with open(filename, "w") as f:
        f.write(content)


# CONSTANT PARAMETERS
inlet_diameter = 9e-3  # m from CAD
k_b = F.k_B  # eV/K, boltzmann constant


### LIPB ###

breeder = "LiPb"
breeder_temperature = 603.15  # K from Utili 2023
LiPb_density = (
    10520.35 - 1.19051 * breeder_temperature
)  # kg/m3 ; equation from Martelli 2019
E_D = 19500 * 1.0364e-5  # = 0.202098
LiPb_diffusivity = 4.03e-8 * np.exp(
    -E_D / (k_b * breeder_temperature)
)  # m2/s ; from Utili 2023, 1 J/mol = 1.0364E-5eV
kinematic_viscosity = calculate_LiPb_kinematic_viscosity(
    breeder_temperature, LiPb_density, breeder, suppress_print=True
)

# velocity parametrization
velocities = np.arange(0.1, 1.1, 0.1)  # m/s

for vel in velocities:
    Re = calculate_reynolds_number(
        vel,
        inlet_diameter,
        kinematic_viscosity,
        breeder,
        suppress_print=True,
    )
    k = calculate_initial_k(vel)
    omega = calculate_initial_omega(k, inlet_diameter)

    openfoam_folder = f"openfoam/velocity_parametrization/lipb/vel_{vel:.1f}"
    os.makedirs(openfoam_folder, exist_ok=True)

    shutil.copytree(
        "openfoam/lipb_simple/0/",
        openfoam_folder + "/0",
    )  # p, nut files are the same as benchmark case
    shutil.copytree(
        "openfoam/lipb_simple/system/",
        openfoam_folder + "/system",
    )
    shutil.copytree(
        "openfoam/lipb_simple/constant/",
        openfoam_folder + "/constant",
    )

    variables_dict = {
        "U": [0.5, vel],
        "k": [0.000937, k],
        "omega": [0.43, omega],
    }

    for name, values in variables_dict.items():
        change_variable_in_openfoam_file(
            filename=openfoam_folder + "/0/" + name,
            old_value=values[0],
            new_value=values[1],
        )

    shutil.copy("meshing/tes_openfoam.msh", openfoam_folder)

### FLIBE ###

breeder = "FLiBe"
breeder_temperature = 900  # K from Meschini 2021
FLiBe_density = 2245 - 0.424 * (
    breeder_temperature - 273.15
)  # kg/m3 ; equation from Vidrio 2022
flibe_diffusivity = htm.diffusivities.filter(material=htm.FLIBE).mean()
E_D = flibe_diffusivity.act_energy.magnitude  # eV
D_0 = flibe_diffusivity.pre_exp.magnitude  # m2/s

FLiBe_diffusivity = D_0 * np.exp(-E_D / (k_b * breeder_temperature))  # m2/s

kinematic_viscosity = calculate_FLiBe_kinematic_viscosity(
    breeder_temperature, FLiBe_density, breeder
)

# velocity parametrization
velocities = np.arange(1, 3.2, 0.2)  # m/s

for vel in velocities:
    Re = calculate_reynolds_number(
        vel,
        inlet_diameter,
        kinematic_viscosity,
        breeder,
        suppress_print=True,
    )
    k = calculate_initial_k(vel)
    omega = calculate_initial_omega(k, inlet_diameter)

    openfoam_folder = f"openfoam/velocity_parametrization/flibe/vel_{vel:.1f}"
    os.makedirs(openfoam_folder, exist_ok=True)

    shutil.copytree(
        "openfoam/flibe_simple/0/",
        openfoam_folder + "/0",
    )  # p, nut files are the same as flibe simple case
    shutil.copytree(
        "openfoam/flibe_simple/system/",
        openfoam_folder + "/system",
    )
    shutil.copytree(
        "openfoam/flibe_simple/constant/",
        openfoam_folder + "/constant",
    )

    variables_dict = {
        "U": [2.7, vel],
        "k": [0.027, k],
        "omega": [2.58, omega],
    }

    for name, values in variables_dict.items():
        change_variable_in_openfoam_file(
            filename=openfoam_folder + "/0/" + name,
            old_value=values[0],
            new_value=values[1],
        )

    shutil.copy("meshing/tes_openfoam.msh", openfoam_folder)

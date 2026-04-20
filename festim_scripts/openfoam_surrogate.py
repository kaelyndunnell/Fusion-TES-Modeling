from autoemulate.simulations.base import Simulator
import torch
from openfoam.parametrize_velocity import change_variable_in_openfoam_file
from openfoam.FLiBe_properties import (
    calculate_FLiBe_kinematic_viscosity,
    calculate_initial_k,
    calculate_initial_omega,
    calculate_reynolds_number,
)
import festim as F
import h_transport_materials as htm
import numpy as np
import os
import shutil
import subprocess


# CONSTANT PARAMETERS
inlet_diameter = 9e-3  # m from CAD
k_b = F.k_B  # eV/K, boltzmann constant

# breeder parameters
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


class OpenFOAMProblem(Simulator):
    def _forward(self, x: torch.Tensor) -> torch.Tensor:
        v_in = x[:, 0]

        # convert to floats
        v_in = v_in.item()

        # CREATE OPENFOAM MODEL
        Re = calculate_reynolds_number(
            v_in, inlet_diameter, kinematic_viscosity, breeder, suppress_print=True
        )
        k = calculate_initial_k(v_in)
        omega = calculate_initial_omega(k, inlet_diameter)

        openfoam_folder = f"openfoam/velocity_parametrization/flibe/vel_{v_in:.1f}"
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
            "U": [2.7, v_in],
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

        # RUN OPENFOAM MODEL
        subprocess.run(
            [
                "openfoam/velocity_parametrization/run_flibe_case.sh",
                "-p " + openfoam_folder,  # case pathway
                "-b " + breeder,  # breeder
                "-v " + str(v_in),  # velocity
            ]
        )

        # solutions file --> TODO: maybe instead of pulling specific timestep for final to pass to festim we can write a script that pulls the max time folder cause this will be the one that is the converged solution
        openfoam_output = openfoam_folder + "/tes.foam"

        y = torch.tensor([openfoam_output]).T

        return y


simulator = OpenFOAMProblem(
    parameters_range={
        "v_in": (1, 3),
    },  # ranges for each variable
    output_names=["foam_file"],
)

n_samples = 2

X = simulator.sample_inputs(n_samples)

Y, _ = simulator.forward_batch(X, allow_failures=False)

print(Y)

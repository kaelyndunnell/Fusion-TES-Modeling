from autoemulate.simulations.base import Simulator
import torch
import festim as F
from main import build_festim_model
from autoemulate import AutoEmulate
import h_transport_materials as htm
import numpy as np
import os
import sys
import shutil
import subprocess

# add openfoam/ to path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, ".."))
sys.path.append(parent_dir)

from openfoam.FLiBe_properties import (
    calculate_FLiBe_kinematic_viscosity,
    calculate_initial_k,
    calculate_initial_omega,
    calculate_reynolds_number,
)
from openfoam.change_variable_openfoam import change_variable_in_openfoam_file


class FestimProblem(Simulator):
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

        openfoam_folder = f"openfoam/velocity_parametrization/{breeder}/vel_{v_in:.1f}"
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
                "openfoam/velocity_parametrization/run_single_case.sh",
                "-p " + openfoam_folder,  # case pathway
                "-b " + breeder,  # breeder
                "-v " + str(v_in),  # velocity
            ]
        )

        # solutions file
        openfoam_output = openfoam_folder

        c_in = x[:, 0]
        residual_pressure = x[:, 1]

        # convert to float
        c_in = c_in.item()
        residual_pressure = residual_pressure.item()

        # FESTIM MODEL
        model = build_festim_model(
            c_inlet=c_in,
            residual_pressure=residual_pressure,
            breeder="flibe",
            openfoam_data_folder=openfoam_output,
            festim_mesh_file="meshing/test_festim.msh",
            results_folder=f"flibe_festim_results/inlet_{c_in}_pressure_{residual_pressure}",
        )

        # solve the model
        model.initialise()
        model.run()

        # extract c at outlet and permeation flux
        c_out = model.exports[0].data
        permeation_flux = model.exports[-1].data

        y = torch.tensor([c_out, permeation_flux]).T

        # ensure the output is a 2D tensor
        if y.ndim == 1:
            y = y.unsqueeze(1)

        return y


# BREEDER/OPENFOAM PARAMETERS
inlet_diameter = 9e-3  # m from CAD
k_b = F.k_B  # eV/K, boltzmann constant

# breeder parameters
breeder = "flibe"
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

# set up model
simulator = FestimProblem(
    parameters_range={
        "v_in": (0.5, 3),
        "c_in": (1e15, 1e25),
        "residual_pressure": (0.0, 0.0),
    },  # ranges for each variable
    output_names=["c_out", "permeation_flux"],
)

# training data
n_samples = 20

X = simulator.sample_inputs(n_samples)
Y, _ = simulator.forward_batch(X, allow_failures=False)

# train the model
ae = AutoEmulate(X, Y, log_level="info")
emulator = [r for r in ae.results if r.model_name == "GaussianProcessRBF"][0]

print(f"Selected model: {emulator.model_name} with id: {emulator.id}")

# plotting
# ae.plot_preds(emulator, output_names=simulator.output_names) # TODO: add plotting visualization that works

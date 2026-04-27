from autoemulate.simulations.base import Simulator
import torch
import festim as F
from main import build_festim_model, findDir
from autoemulate import AutoEmulate
import os
import sys
import h_transport_materials as htm
import numpy as np
import shutil
import subprocess
import matplotlib.pyplot as plt

# add openfoam/ to path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, ".."))
sys.path.append(parent_dir)

from openfoam.LiPb_properties import (
    calculate_LiPb_kinematic_viscosity,
    calculate_initial_k,
    calculate_initial_omega,
    calculate_reynolds_number,
)
from openfoam.change_variable_openfoam import change_variable_in_openfoam_file


def create_folders(
    new_folder_path, bench_folder_path, openfoam_mesh, turbulent_variables_dict
):
    os.makedirs(new_folder_path, exist_ok=True)

    shutil.copytree(
        bench_folder_path + "/0/",
        new_folder_path + "/0",
    )  # p, nut files are the same as lipb simple case
    shutil.copytree(
        bench_folder_path + "/system/",
        new_folder_path + "/system",
    )
    shutil.copytree(
        bench_folder_path + "/constant/",
        new_folder_path + "/constant",
    )

    for name, values in turbulent_variables_dict.items():
        change_variable_in_openfoam_file(
            filename=new_folder_path + "/0/" + name,
            old_value=values[0],
            new_value=values[1],
        )

    shutil.copy(openfoam_mesh, new_folder_path)


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

        variables_dict = {
            "U": [0.5, v_in],
            "k": [0.000937, k],
            "omega": [0.43, omega],
        }

        # TANK PATHS
        TANK_NEW_FOLDER = (
            f"openfoam/velocity_parametrization/lipb/inlet_tank_{v_in:.1f}m_s"
        )
        TANK_BENCH_FOLDER = "openfoam/inlet_tank"
        TANK_MESH = "meshing/inlet_tank.msh"

        # PIPE_PATHS
        PIPE_NEW_FOLDER = f"openfoam/velocity_parametrization/lipb/pipe_{v_in:.1f}m_s"
        PIPE_BENCH_FOLDER = "openfoam/lipb_simple"
        PIPE_MESH = "meshing/tes_openfoam.msh"

        # make tank folders
        create_folders(
            new_folder_path=TANK_NEW_FOLDER,
            bench_folder_path=TANK_BENCH_FOLDER,
            openfoam_mesh=TANK_MESH,
            turbulent_variables_dict=variables_dict,
        )

        # make pipe folders
        create_folders(
            new_folder_path=PIPE_NEW_FOLDER,
            bench_folder_path=PIPE_BENCH_FOLDER,
            openfoam_mesh=PIPE_MESH,
            turbulent_variables_dict=variables_dict,
        )

        # RUN OPENFOAM TANK MODEL
        subprocess.run(
            [
                "openfoam/velocity_parametrization/run_tank_case.sh",
                "-p " + TANK_NEW_FOLDER,  # case pathway
                "-b " + breeder,  # breeder
                "-v " + str(v_in),  # velocity
            ]
        )

        # RUN OPENFOAM MODEL
        subprocess.run(
            [
                "openfoam/velocity_parametrization/run_pipe_case.sh",
                "-p " + PIPE_NEW_FOLDER,  # pipe case pathway
                "-b " + breeder,  # breeder
                "-v " + str(v_in),  # velocity
            ]
        )

        # solutions file
        openfoam_output = PIPE_NEW_FOLDER

        c_in = x[:, 0]
        residual_pressure = x[:, 1]

        # convert to float
        c_in = c_in.item()
        residual_pressure = residual_pressure.item()

        # FESTIM MODEL
        model = build_festim_model(
            c_inlet=c_in,
            residual_pressure=residual_pressure,
            breeder="lipb",
            openfoam_data_folder=openfoam_output,
            festim_mesh_file="meshing/test_festim.msh",
            results_folder=f"lipb_festim_results/inlet_{c_in}_pressure_{residual_pressure}",
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
breeder = "lipb"
breeder_temperature = 603.15  # K from Utili 2023
LiPb_density = (
    10520.35 - 1.19051 * breeder_temperature
)  # kg/m3 ; equation from Martelli 2019

lipb_diffusivity = (
    htm.diffusivities.filter(material=htm.LIPB)
    .filter(exclude=True, isotope="H")
    .filter(exclude=True, isotope="D")
    .mean()
)
E_D = lipb_diffusivity.act_energy.magnitude  # eV
D_0 = lipb_diffusivity.pre_exp.magnitude  # m2/s
LiPb_diffusivity = D_0 * np.exp(-E_D / (k_b * breeder_temperature))  # m2/s

kinematic_viscosity = calculate_LiPb_kinematic_viscosity(
    breeder_temperature, LiPb_density, breeder, suppress_print=True
)

# set up model
simulator = FestimProblem(
    parameters_range={
        "v_in": (0.1, 2),
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
fig = ae.plot(emulator)
fig.savefig("lipb_plot.png")
plt.close(fig)

from autoemulate.simulations.base import Simulator
import torch
import festim as F
from main import build_festim_model
from autoemulate import AutoEmulate
import os
import sys
import h_transport_materials as htm
import numpy as np
import shutil
import subprocess
import matplotlib.pyplot as plt
from autoemulate.core.plotting import create_and_plot_slice

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
        try:
            v_in = x[:, 0]

            # convert to floats
            v_in = v_in.item()
            v_in = round(v_in, 2)

            # # CREATE OPENFOAM MODEL
            # tank_inlet_diameter = 0.0275 * 2

            # Re = calculate_reynolds_number(
            #     v_in, tank_inlet_diameter, kinematic_viscosity, breeder, suppress_print=True
            # )
            # k = calculate_initial_k(v_in)
            # omega = calculate_initial_omega(k, tank_inlet_diameter)

            # tank_variables_dict = {
            #     "U": [0.5, v_in],
            #     "k": [0.000937, k],
            #     "omega": [0.43, omega],
            # }

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

            # # TANK PATHS
            # TANK_NEW_FOLDER = f"{parent_dir}/openfoam/velocity_parametrization/lipb/inlet_tank_{v_in:.2f}m_s"
            # TANK_BENCH_FOLDER = f"{parent_dir}/openfoam/inlet_tank"
            # TANK_MESH = f"{parent_dir}/meshing/inlet_tank.msh"

            # PIPE_PATHS
            PIPE_NEW_FOLDER = (
                f"{parent_dir}/openfoam/velocity_parametrization/lipb/pipe_{v_in:.2f}m_s"
            )
            PIPE_BENCH_FOLDER = f"{parent_dir}/openfoam/lipb_simple"
            PIPE_MESH = f"{parent_dir}/meshing/tes_openfoam.msh"

            if os.path.isdir(PIPE_NEW_FOLDER):
                print("Skipping to FESTIM simulation.")
            else:
                # # make tank folders
                # create_folders(
                #     new_folder_path=TANK_NEW_FOLDER,
                #     bench_folder_path=TANK_BENCH_FOLDER,
                #     openfoam_mesh=TANK_MESH,
                #     turbulent_variables_dict=tank_variables_dict,
                # )

                # make pipe folders
                create_folders(
                    new_folder_path=PIPE_NEW_FOLDER,
                    bench_folder_path=PIPE_BENCH_FOLDER,
                    openfoam_mesh=PIPE_MESH,
                    turbulent_variables_dict=variables_dict,
                )

                # # RUN OPENFOAM TANK MODEL
                # subprocess.run(
                #     [
                #         f"{parent_dir}/openfoam/velocity_parametrization/run_tank_case.sh",
                #         "-p " + TANK_NEW_FOLDER,  # case pathway
                #         "-b " + breeder,  # breeder
                #         "-v " + f"{v_in:.2f}",  # velocity
                #     ]
                # )

                # RUN OPENFOAM MODEL
                subprocess.run(
                    [
                        f"{parent_dir}/openfoam/velocity_parametrization/run_pipe_case.sh",
                        "-p " + PIPE_NEW_FOLDER,  # pipe case pathway
                        "-b " + breeder,  # breeder
                        "-v " + f"{v_in:.2f}",  # velocity
                    ]
                )

            # solutions file
            openfoam_output = PIPE_NEW_FOLDER

            c_in = x[:, 1]
            # residual_pressure = x[:, 2]

            # convert to float
            c_in = 10**(c_in.item())
            # c_in = 10**c_in
            print(c_in)

            # FESTIM MODEL
            model = build_festim_model(
                c_inlet=c_in,
                residual_pressure=0.0,
                breeder="lipb",
                openfoam_data_folder=openfoam_output,
                festim_mesh_file=f"{parent_dir}/meshing/tes_festim.msh",
                results_folder=f"lipb_festim_results/emulator_2/inlet_{c_in}_velocity_{v_in}",
            )

            # solve the model
            model.initialise()
            model.run()

            # extract c at outlet and permeation flux
            c_out = model.exports[0].data
            permeation_flux = model.exports[1].data

            if c_out is None or np.isnan(c_out): # flag error
                raise ValueError(f"Invalid output: {c_out}")
            
            y = torch.tensor([c_out]).T
            return y
        
        except Exception as e: # raise error
            print(f"Simulation failed for inputs {x}: {e}")
            raise 

if __name__ == "__main__":

    # set up model
    simulator = FestimProblem(
        parameters_range={
            "v_in": (0.01, 2.1),
            "c_in": (np.log10(1e15), np.log10(1e25)),
        },  # ranges for each variable
        # output_names=["c_out", "permeation_flux"],
        output_names=["c_out"],
    )

    # training data
    n_samples = 25

    X = simulator.sample_inputs(n_samples)
    Y, _ = simulator.forward_batch(X, allow_failures=True)
    Y = np.log10(Y)

    np.savetxt("lipb_emulators/emulator_2/inputs.csv", X, delimiter=",")
    np.savetxt("lipb_emulators/emulator_2/outputs.csv", Y, delimiter=",")
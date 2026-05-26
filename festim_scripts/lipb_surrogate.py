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
)
from openfoam.change_variable_openfoam import change_variable_in_openfoam_file

N_A = 6.0221e23 # atms/mol

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
            # CREATE MESH

            if os.path.isdir(f"{parent_dir}/meshing/parametric_meshes"):
                print("Meshes will be dumped into meshing/parametric_meshes.")
            else:
                os.mkdir(f"{parent_dir}/meshing/parametric_meshes")
            
            bend_radius = x[:,0]
            length = x[:,1]

            # convert to floats
            bend_radius = round(bend_radius.item(), 2)
            length = round(length.item(), 2)

            # bend_radius = 0.10
            # length = 1.50

            subprocess.run(
                [
                    f"{parent_dir}/meshing/run_salome.sh",
                    "-m", f"{parent_dir}/meshing/parametric_meshes",  # meshing folder
                    "-r", f"{bend_radius}",  # bend radius (m)
                    "-l", f"{length}",  # pipe length (m)
                ]
            )

            v_in = x[:, 2]

            # convert to floats
            v_in = v_in.item()
            v_in = round(v_in, 2)

            # v_in = 1.40

            # BREEDER/OPENFOAM PARAMETERS
            inlet_diameter = 9e-3  # m from CAD
            k = calculate_initial_k(v_in)
            omega = calculate_initial_omega(k, inlet_diameter)

            variables_dict = {
                "U": [0.5, v_in],
                "k": [0.000937, k],
                "omega": [0.43, omega],
            }

            # PIPE_PATHS
            PIPE_NEW_FOLDER = (
                f"{parent_dir}/openfoam/velocity_parametrization/lipb/pipe_{v_in:.2f}m_s_r{bend_radius:.2f}_l{length:.2f}"
            )
            PIPE_BENCH_FOLDER = f"{parent_dir}/openfoam/lipb_simple"
            PIPE_MESH = f"{parent_dir}/meshing/parametric_meshes/one_vol_0.0046_{bend_radius:.2f}_{length:.2f}.unv"

            
            if os.path.isdir(PIPE_NEW_FOLDER):
                print("Skipping to FESTIM simulation.")
            else:
            # make pipe folders
                create_folders(
                    new_folder_path=PIPE_NEW_FOLDER,
                    bench_folder_path=PIPE_BENCH_FOLDER,
                    openfoam_mesh=PIPE_MESH,
                    turbulent_variables_dict=variables_dict,
                )

                # RUN OPENFOAM MODEL
                subprocess.run(
                    [
                        f"{parent_dir}/openfoam/velocity_parametrization/run_pipe_case.sh",
                        "-p " + PIPE_NEW_FOLDER,  # pipe case pathway
                        "-b " + f"lipb",  # breeder
                        "-v " + f"{v_in:.2f}",  # velocity
                        "-m" + PIPE_MESH, # mesh file
                        "-r" + f"{bend_radius:.2f}", # bend radius
                        "-l" + f"{length:.2f}", # length
                    ]
                )

            # solutions file
            openfoam_output = PIPE_NEW_FOLDER

            c_in = x[:, 3]

            # convert to float
            c_in = 10**(c_in.item())
            print(f"Proceeding for inlet concentration of {c_in} T/m3.")

            # FESTIM MODEL
            model = build_festim_model(
                c_inlet=c_in,
                residual_pressure=0.0,
                breeder="lipb",
                openfoam_data_folder=openfoam_output,
                festim_mesh_file=f"{parent_dir}/meshing/parametric_meshes/two_vol_0.0046_{bend_radius:.2f}_{length:.2f}.med", 
                results_folder=f"lipb_festim_results/with_meshes/in_{c_in:.2e}_vel_{v_in:.2e}_bend_{bend_radius:.2f}_len_{length:.2f}",
            )

            # solve the model
            model.initialise()
            model.run()

            # extract c at outlet and permeation flux
            c_out = model.exports[0].data
            permeation_flux = model.exports[1].data

            if c_out is None or np.isnan(c_out): # flag error
                raise ValueError(f"Invalid output: {c_out}")
            
            y = torch.tensor([c_out, permeation_flux]).T
            y = np.log10(y)

            # ensure the output is a 2D tensor
            if y.ndim == 1:
                y = y.unsqueeze(1)

            return y
    
        except Exception as e: # raise error
            print(f"Simulation failed for inputs {x}: {e}")
            raise 

if __name__ == "__main__":

    # set up model
    simulator = FestimProblem(
        parameters_range={
            "bend_radius": (0.02,2),
            "length": (1,8),
            "v_in": (0.01, 3.5),
            "c_in": (np.log10(1e24/N_A), np.log10(1e25/N_A)),
        },  # ranges for each variable
        output_names=["c_out", "permeation_flux"],
    )

    # training data
    n_samples = 1

    X = simulator.sample_inputs(n_samples)
    Y, _ = simulator.forward_batch(X, allow_failures=True)

    results_folder = "lipb_emulators/parametric_mesh"

    if os.path.isdir(results_folder):
        print(f"Results will be dumped into {results_folder}")
    else:
        os.mkdir(results_folder)

    np.savetxt(f"{results_folder}/inputs.csv", X, delimiter=",")
    np.savetxt(f"{results_folder}/outputs.csv", Y, delimiter=",")

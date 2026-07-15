import matplotlib.pyplot as plt
import numpy as np
import os
from autoemulate import AutoEmulate
from autoemulate.core.plotting import create_and_plot_slice
from lipb_surrogate import FestimProblem, calculate_inlet_velocity
from autoemulate.learners import stream
from autoemulate.emulators import GaussianProcessRBF
import torch
from joblib import dump
import sys

# add openfoam/ to path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, ".."))
sys.path.append(parent_dir)

N_A = 6.0221e23  # atms/mol


def make_gp(x_train, y_train, lr=5e-2):
    return GaussianProcessRBF(
        x_train,
        y_train,
        lr=lr,
        standardize_y=False,
    )


## BUILD AND TRAIN INITIAL EMULATOR ##
inlet_diameter = 9.2e-3  # m
breeder_temperature = 723.15  # K from Utili 2023
LiPb_density = (
    10520.35 - 1.19051 * breeder_temperature
)  # kg/m3 ; equation from Martelli 2019
min_velocity = calculate_inlet_velocity(
    flow_rate=0.2,
    inlet_diameter=inlet_diameter,
    breeder_density=LiPb_density,
    breeder="lipb",
    suppress_print=True,
)

simulator = FestimProblem(
    parameters_range={
        "bend_radius": (0.02, 0.2),
        "length": (0.5, 1.1),
        "v_in": (
            min_velocity,
            0.97,
        ),  # from https://doi.org/10.3390/en16073022, 0.2-4.6 kg/s, min vel about 0.3m/s
        "c_in": (np.log10(1e16 / N_A), np.log10(1e25 / N_A)),
    },  # ranges for each variable
    output_names=["c_out", "permeation_flux"],
)

# load data
weak_folder = "lipb_emulators/parametric_mesh/trained"
X = np.loadtxt(f"{parent_dir}/{weak_folder}/inputs.csv", delimiter=",")
Y = np.loadtxt(f"{parent_dir}/{weak_folder}/outputs.csv", delimiter=",")

# train initial emulator
print("training initial emulator...")
emulator = make_gp(X, Y, 0.1)
emulator.fit(X, Y)

strong_folder = weak_folder  # +"/trained" # results folder

x_train = torch.from_numpy(X)
y_train = torch.from_numpy(Y)

## ACTIVE LEARNING ##
print("building active learning...")

# build active learner
learner = stream.Random(
    simulator=simulator,
    emulator=emulator,
    x_train=x_train,
    y_train=y_train,
    p_query=0.3,
    show_progress=True,
)

# in case of failure, patch forward method
original_forward = simulator.forward


def safe_forward(x):
    try:
        return original_forward(x)
    except Exception as e:
        print(f"Simulation failed for inputs {x}: {e}")
        return None


simulator.forward = safe_forward

# patch fit so learner keeps going
original_fit = learner.fit


def safe_fit(*args):
    try:
        return original_fit(*args)
    except Exception as e:
        print(f"Skipping failed simulation: {e}")


learner.fit = safe_fit

# stream samples
print("active learning...")
X_stream = simulator.sample_inputs(500)  # need sufficient amount to see metrics on plot
learner.fit_samples(X_stream)

# save emulator
path = strong_folder
if not os.path.exists(path):
    os.makedirs(path)

# save trained data
x_np = (
    learner.x_train.numpy()
    if hasattr(learner.x_train, "numpy")
    else np.array(learner.x_train)
)
y_np = (
    learner.y_train.numpy()
    if hasattr(learner.y_train, "numpy")
    else np.array(learner.y_train)
)

np.savetxt(os.path.join(path, "trained_inputs.csv"), x_np, delimiter=",")
np.savetxt(os.path.join(path, "trained_outputs.csv"), y_np, delimiter=",")

emulator_filepath = os.path.join(path, "emulator.joblib")
dump(learner.emulator, emulator_filepath)
print("active learning complete!")

## PLOTTING ##
# LEARNER METRICS PLOT
print("plotting learner metrics...")

fig, axs = plt.subplots(
    nrows=len(learner.metrics), ncols=1, sharex=True, figsize=(8, 15)
)
if len(learner.metrics) == 1:
    axs = [axs]

for i, (k, v) in enumerate(learner.metrics.items()):
    axs[i].plot(v, c="k", alpha=0.8)
    axs[i].set_ylabel(k)
axs[-1].set_xlabel("Iterations")

axs[1].set_ylim(0, 1)
plt.savefig(f"{strong_folder}/learner_metrics.png")
print("plotting completed! script finished.")

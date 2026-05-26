import matplotlib.pyplot as plt
import numpy as np
import os
from autoemulate import AutoEmulate
from autoemulate.core.plotting import create_and_plot_slice
from lipb_surrogate import FestimProblem
from autoemulate.learners import stream
from autoemulate.emulators import GaussianProcessRBF
import torch
from joblib import dump


def make_gp(x_train, y_train, lr=5e-2):
    return GaussianProcessRBF(
        x_train, 
        y_train, 
        lr=lr, 
        standardize_y=False,
    )

## BUILD AND TRAIN INITIAL EMULATOR ##
simulator = FestimProblem(
    parameters_range={
        "v_in": (0.01, 2.1),
        "c_in": (np.log10(1e15), np.log10(1e25)),
    },  # ranges for each variable
    output_names=["c_out"],
)

# load data
X = np.loadtxt("/home/kaelyn/Fusion-TES-Modeling/lipb_emulators/emulator_1/inputs_updated.csv", delimiter=",") 
Y = np.loadtxt("/home/kaelyn/Fusion-TES-Modeling/lipb_emulators/emulator_1/outputs_updated.csv", delimiter=",") 

# train initial emulator 
print("training initial emulator...")
emulator = make_gp(X,Y, 0.1)
emulator.fit(X,Y)

# PLOT INITIAL EMULATOR
fig_mean, axs = create_and_plot_slice(
    emulator,
    output_idx=0,
    parameters_range=simulator.parameters_range,
    quantile=0.5,
    param_pair=(0, 1),
)
plt.scatter(X[:, 0], X[:, 1])
plt.suptitle(f"{simulator.output_names[0]}")
plt.savefig("lipb_emulators/trained_emulator_1/weak_emulator.png")

print('initial emulator trained and plotted!')

x_train = torch.from_numpy(X)
y_train = torch.from_numpy(Y).squeeze().reshape(-1,1)

## ACTIVE LEARNING ##
print('running active learning...')

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
X_stream = simulator.sample_inputs(50) # need sufficient amount to see metrics on plot 
learner.fit_samples(X_stream)

# save emulator
path = "lipb_emulators/trained_emulator_1"
if not os.path.exists(path):
    os.makedirs(path)

# save trained data
x_np = learner.x_train.numpy() if hasattr(learner.x_train, 'numpy') else np.array(learner.x_train)
y_np = learner.y_train.numpy() if hasattr(learner.y_train, 'numpy') else np.array(learner.y_train)

np.savetxt(os.path.join(path, "trained_inputs.csv"), x_np, delimiter=",")
np.savetxt(os.path.join(path, "trained_outputs.csv"), y_np, delimiter=",")

emulator_filepath = os.path.join(path, "emulator.joblib")
dump(learner.emulator, emulator_filepath)

## PLOTTING ##
# LEARNER METRICS PLOT
fig, axs = plt.subplots(
    nrows=len(learner.metrics), ncols=1, sharex=True, figsize=(8, 15)
)
for i, (k, v) in enumerate(learner.metrics.items()):
    axs[i].plot(v, c="k", alpha=0.8)
    axs[i].set_ylabel(k)
axs[-1].set_xlabel("Iterations")

axs[1].set_ylim(0, 1)
plt.savefig("lipb_emulators/trained_emulator_1/learner_metrics.png")

# MEAN AND VARIANCE PLOT OF TRAINED EMULATOR
fig_mean, axs = create_and_plot_slice(
    learner.emulator,
    output_idx=0,
    parameters_range=simulator.parameters_range,
    quantile=0.5,
    param_pair=(0, 1),
)
plt.scatter(learner.x_train[:, 0], learner.x_train[:, 1])
plt.suptitle(f"{simulator.output_names[0]}")
plt.savefig("lipb_emulators/trained_emulator_1/strong_emulator.png")

print('trained model trained and plotted!')
print('finishing plotting...')

# COMPARE INITIAL AND TRAINED EMULATORS PLOT

# original simulator
simulator_original = FestimProblem(
    parameters_range={
        "v_in": (0.01, 2.1),
        "c_in": (np.log10(1e15), np.log10(1e25)),
    },  # ranges for each variable
    output_names=["c_out"],
)

# test weak emulator
X_test = simulator_original.sample_inputs(5)
Y_mean_weak, var_weak= emulator.predict_mean_and_variance(X_test)
Y_true, _ = simulator_original.forward_batch(X_test, allow_failures=True)

# test strong emulator
Y_mean_strong, var_strong = learner.emulator.predict_mean_and_variance(X_test)
Y_std_strong = var_strong.sqrt()

# sort based on x index for first column for plotting 
idx0 = np.argsort(X_test[:, 0])
x_plot0 = X_test[idx0, 0]
y_true0 = Y_true.flatten()[idx0]

y_mean_weak0 = Y_mean_weak.flatten()[idx0]
y_mean_strong0 = Y_mean_strong.flatten()[idx0]
y_std_strong0 = Y_std_strong.flatten()[idx0]

# sort based on x index for second column 
idx1 = np.argsort(X_test[:, 1])
x_plot1 = X_test[idx1, 1]
y_true1 = Y_true.flatten()[idx1]

y_mean_weak1 = Y_mean_weak.flatten()[idx1]
y_mean_strong1 = Y_mean_strong.flatten()[idx1]
y_std_strong1 = Y_std_strong.flatten()[idx1]

# plotting
fig, axs = plt.subplots(1,2, figsize=(18,7), sharey=True)

axs[0].plot(x_plot0, y_true0, label='Simulator', alpha=0.5, c='k')
axs[0].plot(x_plot0, y_mean_weak0, label='Initial Emulator')
axs[0].plot(x_plot0, y_mean_strong0, label='Trained Emulator')
axs[0].fill_between(x_plot0, y_mean_strong0 - y_std_strong0, y_mean_strong0 + y_std_strong0, alpha=0.2, label='Confidence')
axs[0].set_xlabel("Inlet velocity (m s$^{-1}$)", fontsize=20)
axs[0].spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("Outlet concentration (T m$^{-3}$)", fontsize=20)
axs[0].tick_params(axis='x', labelsize=20)
axs[0].tick_params(axis='y', labelsize=20)

axs[1].plot(x_plot1, y_true1, label='Simulator', alpha=0.5, c='k')
axs[1].plot(x_plot1, y_mean_weak1, label='Initial Emulator') 
axs[1].plot(x_plot1, y_mean_strong1, label='Trained Emulator') 
axs[1].fill_between(x_plot1, y_mean_strong1 - y_std_strong1, y_mean_strong1 + y_std_strong1, alpha=0.2, label='Confidence')
axs[1].set_xlabel("Inlet concentration (T m$^{-3}$)", fontsize=20)
axs[1].spines[["top", "right"]].set_visible(False)
axs[1].tick_params(axis='x', labelsize=20)
axs[1].tick_params(axis='y', labelsize=20)

plt.legend()
# plt.suptitle(f"Emulator Comparison")
plt.tight_layout()
plt.savefig("lipb_emulators/trained_emulator_1/compare_emulators.png")

print('plotting completed!')
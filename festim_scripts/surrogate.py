from autoemulate.simulations.base import Simulator
import torch
import festim as F
from main import build_festim_model


class FestimProblem(Simulator):
    def _forward(self, x: torch.Tensor) -> torch.Tensor:
        c_in = x[:, 0]
        residual_pressure = x[:, 1]

        # convert to float
        c_in = c_in.item()
        residual_pressure = residual_pressure.item()
        model = build_festim_model(
            c_inlet=c_in,
            residual_pressure=residual_pressure,
            breeder="flibe",
            openfoam_data_file="openfoam/flibe_simple/tes.foam",
            openfoam_final_time=1430,
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


simulator = FestimProblem(
    parameters_range={
        "c_in": (1e15, 1e25),
        "residual_pressure": (0.0, 0.0),
    },  # ranges for each variable
    output_names=["c_out", "permeation_flux"],
)

n_samples = 2

X = simulator.sample_inputs(n_samples)

Y, _ = simulator.forward_batch(X, allow_failures=False)

print(Y)

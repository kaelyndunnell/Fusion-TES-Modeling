# pipe thickness vs. temperature
import numpy as np
from main import compute_efficiency
import matplotlib.pyplot as plt

length = 0.3  # m
height_fluid = 3e-2
temperature = 700  # K
pipe_thickness = 1e-3  # m

concentration_values = np.logspace(16, 22, 10)

efficiency_values = []
for concentration in concentration_values:
    efficiency = compute_efficiency(
        temperature=temperature,
        length=length,
        height_fluid=height_fluid,
        pipe_thickness=pipe_thickness,
        c_inlet=concentration,
    )
    print(efficiency)
    efficiency_values.append(efficiency)

print(efficiency_values)

plt.loglog(concentration_values, efficiency_values)

plt.title("Contour Plot")
plt.ylabel("Temperature (K)")
plt.xlabel("Pipe Thickness (m)")
plt.show()

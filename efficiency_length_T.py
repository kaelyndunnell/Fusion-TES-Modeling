import numpy as np
from main import compute_efficiency
import matplotlib.pyplot as plt

# length vs. temperature
length_values = np.linspace(0.3, 0.5, 10)  # m
temperature_values = np.linspace(500, 800, 10)  # K
pipe_thickness = 4e-3
height_fluid = 3e-2

efficiency_values = []
for length in length_values:
    efficiency_col = []
    for T in temperature_values:
        print(f"length = {length}, T = {T}")
        efficiency = compute_efficiency(
            temperature=T,
            length=length,
            height_fluid=height_fluid,
            pipe_thickness=pipe_thickness,
        )
        efficiency_col.append(efficiency)
    efficiency_values.append(efficiency_col)

# contour plots
XX, YY = np.meshgrid(length_values, temperature_values)

efficiency_values = np.array(efficiency_values).T

CF = plt.contourf(XX, YY, efficiency_values, levels=100)
plt.scatter(XX, YY)

plt.colorbar(CF, label="Efficiency")

plt.title("Contour Plot")
plt.ylabel("Temperature (K)")
plt.xlabel("Length (m)")
plt.savefig("efficiency_length_T.png")
plt.show()

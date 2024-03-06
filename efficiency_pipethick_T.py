# pipe thickness vs. temperature
import numpy as np
from main import compute_efficiency
import matplotlib.pyplot as plt

length = 0.3  # m
height_fluid = 3e-2
temperature_values = np.linspace(500, 800, 6)  # K
pipe_thickness_values = np.linspace(1e-4, 1e-2, 3)  # m

efficiency_values = []
for pipe_thickness in pipe_thickness_values:
    efficiency_col = []
    for T in temperature_values:
        print(f"pipe thickness = {pipe_thickness}, T = {T}")
        efficiency = compute_efficiency(
            temperature=T,
            length=length,
            height_fluid=height_fluid,
            pipe_thickness=pipe_thickness,
        )
        efficiency_col.append(efficiency)
        print(efficiency)
    efficiency_values.append(efficiency_col)

# contour plots
XX, YY = np.meshgrid(pipe_thickness_values, temperature_values)

efficiency_values = np.array(efficiency_values).T

print(efficiency_values)

CF = plt.contourf(XX, YY, efficiency_values, levels=100)
plt.scatter(XX, YY)

plt.colorbar(CF, label="Efficiency")

plt.title("Contour Plot")
plt.ylabel("Temperature (K)")
plt.xlabel("Pipe Thickness (m)")
plt.savefig("efficiency_thickness_T.png")
plt.show()
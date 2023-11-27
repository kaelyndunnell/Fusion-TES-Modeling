# height_fluid vs. temperature parametric study
import numpy as np
from main import compute_efficiency
import matplotlib.pyplot as plt

pipe_thickness = 4e-3 # m
height_fluid_values = np.linspace(pipe_thickness, 1e-2, 5)  # m
temperature_values = np.linspace(500, 800, 5)  # K
length = 0.3  # m

efficiency_values = []
for height_fluid in height_fluid_values:
    efficiency_col = []
    for T in temperature_values:
        print(f"fluid height = {height_fluid}, T = {T}")
        efficiency = compute_efficiency(
            temperature=T,
            length=length,
            height_fluid=height_fluid,
            pipe_thickness=pipe_thickness,
        )
        efficiency_col.append(efficiency)
    efficiency_values.append(efficiency_col)

# contour plots
XX, YY = np.meshgrid(height_fluid_values, temperature_values)

print(efficiency_values)
efficiency_values = np.array(efficiency_values).T

CF = plt.contourf(XX, YY, efficiency_values, levels=100)
plt.scatter(XX, YY)

plt.colorbar(CF, label="Efficiency")

plt.title("Contour Plot")
plt.ylabel("Temperature (K)")
plt.xlabel("Fluid Height (m)")
plt.savefig("efficiency_height_T.png")
plt.show()
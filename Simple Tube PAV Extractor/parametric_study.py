import numpy as np
from main import compute_efficienty
import matplotlib.pyplot as plt

length_values = np.linspace(0.3, 0.5, 10)  # m
temperature_values = np.linspace(500, 800, 10)  # K

efficiency_values = []
for length in length_values:
    for T in temperature_values:
        print(f"length = {length}, T = {T}")
        efficiency = compute_efficienty(temperature=T, length=length)
        efficiency_values.append(efficiency)

# contour plots
XX, YY = np.meshgrid(length_values, temperature_values) 

efficiency_values = np.array(efficiency_values)
efficiency_values = efficiency_values.reshape(len(temperature_values), len(length_values))

CF = plt.contourf(XX, YY, efficiency_values, levels=100)
plt.scatter(XX, YY)

plt.colorbar(CF, label="Efficiency")

plt.title('Contour Plot') 
plt.ylabel('Temperature (K)') 
plt.xlabel('Length (m)')
plt.savefig("efficiency_contour_plot.png")
plt.show()

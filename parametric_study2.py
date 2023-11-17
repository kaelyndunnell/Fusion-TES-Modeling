import numpy as np
from main import compute_efficienty
import matplotlib.pyplot as plt

# length vs. temperature 
length_values = np.linspace(0.3, 0.5, 10)  # m
temperature_values = np.linspace(500, 800, 10)  # K
pipe_thickness = 4e-3
height_fluid=3e-2

efficiency_values = []
for length in length_values:
    for T in temperature_values:
        print(f"length = {length}, T = {T}")
        efficiency = compute_efficienty(temperature=T, length=length, height_fluid=height_fluid, pipe_thickness=pipe_thickness)
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

# height_fluid vs. temperature 
length = 0.3  # m
pipe_thickness = 4e-3
temperature_values = np.linspace(500, 800, 10)  # K
height_fluid_values = np.linspace(1e-4, 1e-2, 10) # m

efficiency_values = []
for height_fluid in height_fluid_values:
    for T in temperature_values:
        print(f"fluid height = {height_fluid}, T = {T}")
        efficiency = compute_efficienty(temperature=T, length=length, height_fluid=height_fluid, pipe_thickness=pipe_thickness)
        efficiency_values.append(efficiency)

# contour plots
XX, YY = np.meshgrid(height_fluid_values, temperature_values) 

efficiency_values = np.array(efficiency_values)
efficiency_values = efficiency_values.reshape(len(temperature_values), len(height_fluid_values))

CF = plt.contourf(XX, YY, efficiency_values, levels=100)
plt.scatter(XX, YY)

plt.colorbar(CF, label="Efficiency")

plt.title('Contour Plot') 
plt.ylabel('Temperature (K)') 
plt.xlabel('Fluid Height (m)')
plt.savefig("efficiency_contour_plot.png")
plt.show()

# pipe thickness vs. temperature 
length = 0.3  # m
height_fluid=3e-2
temperature_values = np.linspace(500, 800, 10)  # K
pipe_thickness_values = np.linspace(1e-4, 1e-2, 10) # m

efficiency_values = []
for pipe_thickness in pipe_thickness_values:
    for T in temperature_values:
        print(f"pipe thickness = {pipe_thickness}, T = {T}")
        efficiency = compute_efficienty(temperature=T, length=length, height_fluid=height_fluid, pipe_thickness=pipe_thickness)
        efficiency_values.append(efficiency)

# contour plots
XX, YY = np.meshgrid(pipe_thickness_values, temperature_values) 

efficiency_values = np.array(efficiency_values)
efficiency_values = efficiency_values.reshape(len(temperature_values), len(pipe_thickness_values))

CF = plt.contourf(XX, YY, efficiency_values, levels=100)
plt.scatter(XX, YY)

plt.colorbar(CF, label="Efficiency")

plt.title('Contour Plot') 
plt.ylabel('Temperature (K)') 
plt.xlabel('Pipe Thickness (m)')
plt.savefig("efficiency_contour_plot.png")
plt.show()
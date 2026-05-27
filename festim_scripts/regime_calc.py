import numpy as np

# calc mass trans coeff

# example case is 1.66 mol/m3
N_A = 6.022e23
c_in = 1.66 * N_A

D_breeder = 4.6e-8
diameter = 9.2e-3  # m
viscosity = 1.1956969781435462e-07  # kinematic
v_in = 2.9  # m/s


Re = v_in * diameter / viscosity
Sc = viscosity / D_breeder  # 0.7
# Sh = 0.3 + (
#     0.62 * Re ** (1 / 2) * Sc ** (1 / 3) * (1 + (0.4 / Sc) ** (2 / 3)) ** (-1 / 4)
# ) * (1 + (Re / 28200) ** (5 / 8))
Sh = 0.023 * Re ** (0.8) * Sc ** (1 / 3)

mass_trans_coeff = Sh * D_breeder / diameter
print(f"mass transfer coefficient is {mass_trans_coeff}")

membrane_thickness = 0.4e-3  # m
D_membrane = 5e-8
K_breeder = 1.6e21
K_membrane = 7.59e22
# mass_trans_coeff = (
#     5.1e-3  # https://www.sciencedirect.com/science/article/abs/pii/092037969190064W
# ) # better agreement with this value but still good agreement with the value calculated above


K_r = 2.7e-25  # recombo coeff

partition_parameter = (
    D_membrane / mass_trans_coeff * K_membrane / K_breeder / membrane_thickness
)
print(f"partition parameter is {partition_parameter}")

######### UROGORRI 2023 MODEL #########

W_0 = K_r * c_in * membrane_thickness / (D_membrane * K_breeder / K_membrane)

print(f"mixed permeation number is {W_0}")

calculated_flux = (  # assuming LDR from Alberghi
    mass_trans_coeff * c_in * partition_parameter / (partition_parameter + 1)
)

# calculated_flux = mass_trans_coeff * c_in  # assuming LLR from Alberghi

# festim_flux = 0.0025916 * N_A  # for v_in = 0.05
festim_flux = 0.0028 * N_A  # for v_in = 2.9
print(f"calculated permeation flux = {calculated_flux}")
print(f"actual flux (from FESTIM) = {festim_flux}")

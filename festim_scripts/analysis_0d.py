import sympy as sp

c_s1 = sp.Symbol("c_s1")
c_s2 = sp.Symbol("c_s2")
# ro cannot be zero
ro = sp.Symbol("ro", nonzero=True)
ri = sp.Symbol("ri", nonzero=True)
D_s = sp.Symbol("D_s", nonzero=True)
K_L = sp.Symbol("K_L", nonzero=True)
K_S = sp.Symbol("K_S", nonzero=True)
K_r = sp.Symbol("K_r", nonzero=True)
K_d = sp.Symbol("K_d", nonzero=True)
c_l1 = sp.Symbol("c_l1", nonzero=True)
c_l2 = sp.Symbol("c_l2", nonzero=True)
P_V = sp.Symbol("P_V")

K_T = sp.Symbol("K_T", nonzero=True)

J = sp.Symbol("J")

flux_through_liquid = sp.Eq(J, K_T * (c_l1 - c_l2))

flux_through_tube = sp.Eq(J, D_s * (c_s1 - c_s2) / (ri * sp.log(ro / ri)))

flux_at_vacuum = sp.Eq(J, K_r * c_s2**2 - K_d * P_V)

conservation_of_chemical_potential = sp.Eq(c_l2 / K_L, c_s1 / K_S)

# solve system of equations
solution = sp.solve(
    [
        flux_through_liquid,
        flux_through_tube,
        flux_at_vacuum,
        conservation_of_chemical_potential,
    ],
    (J, c_s1, c_s2, c_l2),
    dict=True,
)


def replace_values(solution, new_vals):
    new_solution = []
    for sol in solution:
        new_sol = {}
        for symbol, value in sol.items():
            new_val = value.subs(new_vals)
            new_sol[symbol] = new_val
        new_solution.append(new_sol)
    return new_solution


def flux(new_vals):
    new_solution = replace_values(solution, new_vals)
    return new_solution[0][J]


# plot solution with matplotlib
import matplotlib.pyplot as plt

# replace all symbols with values
c_l1_val = 1  # molT/m3
ro_val = 5.00e-3  # m
ri_val = 4.75e-3  # m

new_vals = {
    ro: ro_val,
    ri: ri_val,
    D_s: 1e-9,
    K_L: 1e-3,
    K_S: 2e-3,
    K_r: 1e-3,
    K_d: 1e-3,
    c_l1: c_l1_val,
    P_V: 0,
    K_T: 1e-9,
}
new_values = replace_values(solution, new_vals)
fig, ax = plt.subplots()
ax.plot(
    [0, ri_val, ri_val, ro_val],
    [c_l1_val, new_values[0][c_l2], new_values[0][c_s1], new_values[0][c_s2]],
)
plt.xlabel("r (m)")
plt.ylabel("Concentration (molT/m3)")
plt.show()


def arrhenius_law(pre_exp, activation_energy, temperature, mod=sp):
    R = 8.314  # J/(mol K)
    return pre_exp * mod.exp(-activation_energy / (R * temperature))


def linton_sherwood_correlation(Re, Sc):
    return 0.023 * Re**0.83 * Sc ** (1 / 3)


def mass_transfer_coeff_from_schmidt(Sh, D, ri):
    return Sh * D / (2 * ri)


# data from https://doi.org/10.1080/15361055.2023.2196237
T = 673  # K
Sh = linton_sherwood_correlation(Re=1e5, Sc=100)
D_L = arrhenius_law(8.3e-9, 7.37e3, T)
Q = 1.26e-4  # m3/s
x0 = 0  # m
x1 = 5  # m
ri_val = 4.75e-3  # m
ro_val = 5.00e-3  # m

# from https://doi.org/10.1016/j.fusengdes.2018.11.028
density_of_pbli = 10520.35 - 1.19051 * T  # kg/m3
molar_mass_pbli = 214.1410  # g/mol
# in mol/kg
molar_mass_pbli = 1 / molar_mass_pbli * 1e3

new_vals = {
    ro: ro_val,
    ri: ri_val,
    D_s: arrhenius_law(2.9e-8, 4.2e3, T),
    K_L: arrhenius_law(4.29e-6, 12.8e3, T) * density_of_pbli * molar_mass_pbli,
    # K_L: 1.08e-6 * density_of_pbli * molar_mass_pbli,
    # K_L: arrhenius_law(2.32e-8, 1.35e3, T) * density_of_pbli * molar_mass_pbli,
    K_S: arrhenius_law(0.138, -29.0e3, T),  # this is the correct base value
    # K_S: arrhenius_law(0.126, -35.3e3, T),  # upper bound
    K_r: arrhenius_law(1.3, 127e3, T),
    K_d: 1e-12,  # osef
    P_V: 0,
    K_T: mass_transfer_coeff_from_schmidt(Sh, D_L, ri_val),
}

flux_given_c = lambda c: flux({c_l1: c, **new_vals})
import numpy as np

# solve the ODE with scipy
# dc/dx = - 2 pi R / Q * J(c)
# mol/m = m * m-3 * s * mol/m2/s

from scipy.integrate import solve_ivp

# lambdify the flux function
flux_given_c = sp.lambdify([c_l1], flux({**new_vals}))


def dc_dx(x, c):
    return -2 * np.pi * ri_val / Q * flux_given_c(c)


efficiencies = []
initial_concentrations = np.logspace(-6, 0, num=50)  # molT/m3

for c0 in initial_concentrations:
    sol = solve_ivp(dc_dx, (x0, x1), [c0], t_eval=np.linspace(x0, x1, 100))
    # add colour to curve based on c0
    plt.plot(sol.t, sol.y[0], color=plt.cm.viridis(c0 / initial_concentrations[-1]))
    efficiency = 1 - sol.y[0][-1] / c0
    efficiencies.append(efficiency)
plt.ylim(bottom=0)
plt.xlabel("x (m)")
plt.ylabel("Concentration (molT/m3)")
# colorbar
sm = plt.cm.ScalarMappable(
    cmap=plt.cm.viridis,
    norm=plt.Normalize(vmin=initial_concentrations[0], vmax=initial_concentrations[-1]),
)
sm.set_array([])
plt.colorbar(sm, label="Initial concentration (molT/m3)", ax=plt.gca())
plt.show()

plt.plot(initial_concentrations, flux_given_c(initial_concentrations))
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Concentration (molT/m3)")
plt.ylabel("Flux (mol/m2/s)")
plt.show()

plt.loglog(initial_concentrations, efficiencies)
plt.tight_layout()
plt.xlabel("Initial concentration (molT/m3)")
plt.ylabel("Efficiency")
plt.show()

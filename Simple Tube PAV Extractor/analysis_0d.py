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
c_l1_val = 1
ro_val = 0.15
ri_val = 0.1

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
plt.show()


new_vals = {
    ro: ro_val,
    ri: ri_val,
    D_s: 1e-2,
    K_L: 1e-3,
    K_S: 2e-3,
    K_r: 1e-3,
    K_d: 1e-3,
    P_V: 0,
    K_T: 1e-2,
}

flux_given_c = lambda c: flux({c_l1: c, **new_vals})
import numpy as np

concentration_values = np.logspace(1, 4)
flux_values = [flux_given_c(c) for c in concentration_values]

plt.plot(concentration_values, flux_values)
plt.show()

# solve the ODE with scipy
# dc/dx = - 2 pi R / Q * J(c)

from scipy.integrate import solve_ivp

R = ri_val
Q = 1e-3
c0 = 1
x0 = 0
x1 = 1


# lambdify the flux function
flux_given_c = sp.lambdify([c_l1], flux({**new_vals}))


def dc_dx(x, c):
    return -2 * np.pi * R / Q * flux_given_c(c)


sol = solve_ivp(dc_dx, (x0, x1), [c0], t_eval=np.linspace(x0, x1, 100))
plt.plot(sol.t, sol.y[0])
plt.ylim(bottom=0)
plt.show()

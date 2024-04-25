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

# replace all symbols with values
c_l1_val = 1
ro_val = 0.15
ri_val = 0.1
new_values = {}
for symbol in [c_l2, c_s1, c_s2]:
    new_val = solution[0][symbol].subs(
        {
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
    )
    new_values[symbol] = new_val

# plot solution with matplotlib
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(
    [0, ri_val, ri_val, ro_val],
    [c_l1_val, new_values[c_l2], new_values[c_s1], new_values[c_s2]],
)
plt.show()

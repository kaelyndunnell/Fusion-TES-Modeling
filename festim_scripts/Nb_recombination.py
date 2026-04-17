import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("festim_scripts/Nb_recombo_data.csv")

temp = 10e3 / df.iloc[:, 0]
recombo_coeff = df.iloc[:, 1]  # cm^4/s

plt.figure()
plt.plot(temp, recombo_coeff, "o", label="Data points")

# best fit line
p = np.polyfit(temp, np.log(recombo_coeff), 1)

a = np.exp(p[1])
b = p[0]

x_line = np.linspace(min(temp), max(temp), 30)
y_line = a * np.exp(b * x_line)


plt.plot(x_line, y_line, linestyle="--", linewidth=2)

print(f"Best Fit Line (y={a:.2e} * exp({b:.2e}x)")

plt.yscale("log")

plt.ylabel("Recombination Coefficient (cm4/s)")
plt.xlabel("Temperature (K)")
plt.title("K_r0 vs. Temp from Hayakawa (2003)")
# plt.show()

# calc recombo coeff for T = 900K for our case
festim_temp = 900  # K
recombo_coeff = a * np.exp(b * festim_temp)

print(
    f"Recombination Coefficient for H in Nb at {festim_temp}K is {recombo_coeff:.2e} cm^4/s."
)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_breeder = pd.read_csv("bend_T_breeder.csv")
df_membrane = pd.read_csv("bend_T_membrane.csv")

df_breeder = df_breeder.dropna()  # get rid of nan values
df_membrane = df_membrane.dropna()  # get rid of nan values

# df = pd.read_csv("bend_vel.csv")
# df = df.dropna()

# # U:0 is vel in x direction, arc_length is distance along line which is along pipe diameter in y direction in this case
# velocity = np.sqrt(df["U:0"] ** 2 + df["U:1"] ** 2 + df["U:2"] ** 2)
# distance = df["arc_length"]

distance_b = df_breeder["arc_length"]
distance_m = df_membrane["arc_length"]

concentration_b = df_breeder["T_6"]
concentration_m = df_membrane["T_7"]


fig, axs = plt.subplots(figsize=(6.4, 3))
# axs.plot(distance, velocity, color="steelblue", linewidth=1)
axs.plot(distance_b, concentration_b, color="purple", linewidth=1)
axs.scatter(distance_m, concentration_m, color="purple", s=0.5)
axs.set_xlabel("Arc Length (m)")
# axs.set_ylabel(r"Velocity Magnitude (m s$^{-1}$)")
axs.set_ylabel(r"c$_{\mathrm{T}}$ (mol m$^{-3}$)")
axs.spines["top"].set_visible(False)
axs.spines["right"].set_visible(False)
axs.set_yscale("log")
# axs.set_xticks([])  # remove ticks
# axs.set_yticks([])
axs.set_ylim(1e-3, 1e7)
plt.tight_layout()
plt.savefig("bend_T.pdf")
plt.show()

import h_transport_materials as htm
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("festim_scripts/Nb_recombo_data.csv")

df["1/T"] = df.iloc[:, 0]
df["T"] = 1 / df["1/T"]

df["K_r"] = df.iloc[:, 1]  # cm^4/s


nb_recomb = htm.RecombinationCoeff(
    data_T=list(df["T"].values),
    data_y=list(df["K_r"].values) * htm.ureg["cm^4/particle/s"],
    source="10.1238/Physica.Topical.103a00113",
    material=htm.NIOBIUM,
)


if __name__ == "__main__":
    print(df["T"].values)
    plt.figure()
    htm.plotting.plot(nb_recomb)
    plt.plot(
        df["1/T"],
        (df["K_r"].to_numpy() * htm.ureg["cm^4/particle/s"]).to("m^4/particle/s"),
        "o",
        label="Data points",
    )
    plt.xlabel("1/T (K^-1)")
    plt.ylabel("Recombination Coefficient (m^4/s)")
    plt.yscale("log")
    print(nb_recomb)
    print(nb_recomb.range)
    plt.show()

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("data/two_moment_dt_sweep.csv")

df = df.sort_values("dt_s")

plt.plot(
    df["dt_s"],
    df["SSmax_percent"],
    marker="o",
)

plt.xlabel("Time step dt [s]")
plt.ylabel("Maximum supersaturation [%]")
plt.title("Time-step convergence of SSmax")
plt.tight_layout()
plt.savefig("data/two_moment_dt_convergence_SSmax.png", dpi=200)
plt.close()

plt.plot(
    df["dt_s"],
    df["qc_final"],
    marker="o",
)

plt.xlabel("Time step dt [s]")
plt.ylabel("Final cloud water mixing ratio [kg/kg]")
plt.title("Time-step convergence of final qc")
plt.tight_layout()
plt.savefig("data/two_moment_dt_convergence_qc.png", dpi=200)
plt.close()

print("Saved:")
print("data/two_moment_dt_convergence_SSmax.png")
print("data/two_moment_dt_convergence_qc.png")

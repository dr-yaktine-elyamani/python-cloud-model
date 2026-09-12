import pandas as pd
import matplotlib.pyplot as plt


filename = "data/two_moment_w_N_sweep_dt025.csv"

df = pd.read_csv(filename)


# ------------------------------------------------------------
# Plot 1: SSmax versus aerosol number
# ------------------------------------------------------------

for w in sorted(df["w_m_s"].unique()):

    sub = df[
        df["w_m_s"] == w
    ].sort_values(
        "aerosol_N_m3"
    )

    plt.plot(
        sub["aerosol_N_m3"],
        sub["SSmax_percent"],
        marker="o",
        label=f"w={w} m/s",
    )

plt.xlabel(
    "Aerosol number concentration [m$^{-3}$]"
)

plt.ylabel(
    "Maximum supersaturation [%]"
)

plt.title(
    "SSmax sensitivity to aerosol number (dt=0.25 s)"
)

plt.xscale("log")
plt.legend()
plt.tight_layout()

plt.savefig(
    "data/two_moment_SSmax_vs_N_dt025.png",
    dpi=200,
)

plt.close()


# ------------------------------------------------------------
# Plot 2: activated fraction versus aerosol number
# ------------------------------------------------------------

for w in sorted(df["w_m_s"].unique()):

    sub = df[
        df["w_m_s"] == w
    ].sort_values(
        "aerosol_N_m3"
    )

    plt.plot(
        sub["aerosol_N_m3"],
        sub["activated_fraction"],
        marker="o",
        label=f"w={w} m/s",
    )

plt.xlabel(
    "Aerosol number concentration [m$^{-3}$]"
)

plt.ylabel(
    "Activated fraction"
)

plt.title(
    "Activated fraction sensitivity to aerosol number (dt=0.25 s)"
)

plt.xscale("log")
plt.ylim(0.0, 1.05)
plt.legend()
plt.tight_layout()

plt.savefig(
    "data/two_moment_activated_fraction_vs_N_dt025.png",
    dpi=200,
)

plt.close()


# ------------------------------------------------------------
# Plot 3: activation time versus updraft
# ------------------------------------------------------------

for N in sorted(
    df["aerosol_N_m3"].unique()
):

    sub = df[
        df["aerosol_N_m3"] == N
    ].sort_values(
        "w_m_s"
    )

    plt.plot(
        sub["w_m_s"],
        sub["activation_time_s"],
        marker="o",
        label=f"N={N:.0e} m$^{{-3}}$",
    )

plt.xlabel(
    "Updraft velocity [m s$^{-1}$]"
)

plt.ylabel(
    "Activation time [s]"
)

plt.title(
    "Activation time sensitivity to updraft (dt=0.25 s)"
)

plt.legend()
plt.tight_layout()

plt.savefig(
    "data/two_moment_activation_time_vs_w_dt025.png",
    dpi=200,
)

plt.close()


# ------------------------------------------------------------
# Plot 4: mean radius versus aerosol number
# ------------------------------------------------------------

for w in sorted(df["w_m_s"].unique()):

    sub = df[
        df["w_m_s"] == w
    ].sort_values(
        "aerosol_N_m3"
    )

    plt.plot(
        sub["aerosol_N_m3"],
        sub["mean_radius_um"],
        marker="o",
        label=f"w={w} m/s",
    )

plt.xlabel(
    "Aerosol number concentration [m$^{-3}$]"
)

plt.ylabel(
    "Final mean droplet radius [$\\mu$m]"
)

plt.title(
    "Droplet radius sensitivity to aerosol number (dt=0.25 s)"
)

plt.xscale("log")
plt.legend()
plt.tight_layout()

plt.savefig(
    "data/two_moment_mean_radius_vs_N_dt025.png",
    dpi=200,
)

plt.close()


print("Saved:")
print(
    "data/two_moment_SSmax_vs_N_dt025.png"
)
print(
    "data/two_moment_activated_fraction_vs_N_dt025.png"
)
print(
    "data/two_moment_activation_time_vs_w_dt025.png"
)
print(
    "data/two_moment_mean_radius_vs_N_dt025.png"
)

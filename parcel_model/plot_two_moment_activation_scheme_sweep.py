import os

import pandas as pd
import matplotlib.pyplot as plt


INPUT = (
    "data/"
    "two_moment_activation_scheme_sweep_dt025.csv"
)

OUTDIR = "data"

os.makedirs(
    OUTDIR,
    exist_ok=True,
)

df = pd.read_csv(
    INPUT
)

schemes = [
    "simple_kappa",
    "ARG1998",
]

w_values = sorted(
    df["w_m_s"].unique()
)


# ---------------------------------------------------------
# Helper
# ---------------------------------------------------------

def subset(scheme, w):
    return (
        df[
            (df["activation_scheme"] == scheme)
            & (df["w_m_s"] == w)
        ]
        .sort_values("aerosol_N_m3")
    )


# ---------------------------------------------------------
# 1. SSmax versus aerosol concentration
# ---------------------------------------------------------

plt.figure()

for scheme in schemes:

    for w in w_values:

        d = subset(
            scheme,
            w,
        )

        plt.plot(
            d["aerosol_N_m3"] / 1.0e6,
            d["SSmax_percent"],
            marker="o",
            label=f"{scheme}, w={w:g}",
        )

plt.xlabel(
    "Aerosol number concentration (cm$^{-3}$)"
)

plt.ylabel(
    "Maximum supersaturation (%)"
)

plt.title(
    "Parcel SSmax: activation-scheme comparison"
)

plt.grid(
    True,
    alpha=0.3,
)

plt.legend(
    fontsize=7,
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTDIR,
        "two_moment_scheme_SSmax_vs_N_dt025.png",
    ),
    dpi=200,
)

plt.close()


# ---------------------------------------------------------
# 2. Final activated fraction versus aerosol concentration
# ---------------------------------------------------------

plt.figure()

for scheme in schemes:

    for w in w_values:

        d = subset(
            scheme,
            w,
        )

        plt.plot(
            d["aerosol_N_m3"] / 1.0e6,
            d["activated_fraction"],
            marker="o",
            label=f"{scheme}, w={w:g}",
        )

plt.xlabel(
    "Aerosol number concentration (cm$^{-3}$)"
)

plt.ylabel(
    "Final activated fraction"
)

plt.title(
    "Activated fraction: activation-scheme comparison"
)

plt.ylim(
    0.0,
    1.05,
)

plt.grid(
    True,
    alpha=0.3,
)

plt.legend(
    fontsize=7,
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTDIR,
        "two_moment_scheme_fraction_vs_N_dt025.png",
    ),
    dpi=200,
)

plt.close()


# ---------------------------------------------------------
# 3. 50% activation delay
# ---------------------------------------------------------

plt.figure()

for scheme in schemes:

    for w in w_values:

        d = subset(
            scheme,
            w,
        )

        plt.plot(
            d["aerosol_N_m3"] / 1.0e6,
            d["activation_50pct_delay_s"],
            marker="o",
            label=f"{scheme}, w={w:g}",
        )

plt.xlabel(
    "Aerosol number concentration (cm$^{-3}$)"
)

plt.ylabel(
    "Time from saturation to 50% activation (s)"
)

plt.title(
    "50% activation timing"
)

plt.grid(
    True,
    alpha=0.3,
)

plt.legend(
    fontsize=7,
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTDIR,
        "two_moment_scheme_50pct_delay_vs_N_dt025.png",
    ),
    dpi=200,
)

plt.close()


# ---------------------------------------------------------
# 4. ARG parcel SSmax versus independent ARG prediction
# ---------------------------------------------------------

arg = (
    df[
        df["activation_scheme"] == "ARG1998"
    ]
    .copy()
)


plt.figure()

for w in w_values:

    d = (
        arg[
            arg["w_m_s"] == w
        ]
        .sort_values(
            "aerosol_N_m3"
        )
    )

    plt.plot(
        d["aerosol_N_m3"] / 1.0e6,
        d["SSmax_percent"],
        marker="o",
        label=f"Parcel, w={w:g}",
    )

    plt.plot(
        d["aerosol_N_m3"] / 1.0e6,
        d["ARG_Smax_prediction_percent"],
        marker="x",
        linestyle="--",
        label=f"ARG benchmark, w={w:g}",
    )

plt.xlabel(
    "Aerosol number concentration (cm$^{-3}$)"
)

plt.ylabel(
    "Maximum supersaturation (%)"
)

plt.title(
    "ARG1998 parcel result vs independent ARG benchmark"
)

plt.grid(
    True,
    alpha=0.3,
)

plt.legend(
    fontsize=7,
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTDIR,
        "two_moment_ARG_parcel_vs_benchmark_SSmax_dt025.png",
    ),
    dpi=200,
)

plt.close()


# ---------------------------------------------------------
# 5. ARG parcel fraction versus independent ARG prediction
# ---------------------------------------------------------

plt.figure()

for w in w_values:

    d = (
        arg[
            arg["w_m_s"] == w
        ]
        .sort_values(
            "aerosol_N_m3"
        )
    )

    plt.plot(
        d["aerosol_N_m3"] / 1.0e6,
        d["activated_fraction"],
        marker="o",
        label=f"Parcel, w={w:g}",
    )

    plt.plot(
        d["aerosol_N_m3"] / 1.0e6,
        d["ARG_fraction_prediction"],
        marker="x",
        linestyle="--",
        label=f"ARG benchmark, w={w:g}",
    )

plt.xlabel(
    "Aerosol number concentration (cm$^{-3}$)"
)

plt.ylabel(
    "Activated fraction"
)

plt.title(
    "ARG1998 activated fraction vs independent ARG benchmark"
)

plt.ylim(
    0.0,
    1.05,
)

plt.grid(
    True,
    alpha=0.3,
)

plt.legend(
    fontsize=7,
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTDIR,
        "two_moment_ARG_parcel_vs_benchmark_fraction_dt025.png",
    ),
    dpi=200,
)

plt.close()


print(
    "Saved 5 figures:"
)

print(
    "data/two_moment_scheme_SSmax_vs_N_dt025.png"
)

print(
    "data/two_moment_scheme_fraction_vs_N_dt025.png"
)

print(
    "data/two_moment_scheme_50pct_delay_vs_N_dt025.png"
)

print(
    "data/two_moment_ARG_parcel_vs_benchmark_SSmax_dt025.png"
)

print(
    "data/two_moment_ARG_parcel_vs_benchmark_fraction_dt025.png"
)

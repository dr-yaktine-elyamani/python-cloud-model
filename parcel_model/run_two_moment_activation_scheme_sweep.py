import csv
import os

from parcel_model.run_two_moment_warm import run


# ---------------------------------------------------------
# Settings
# ---------------------------------------------------------

w_values = [
    0.2,
    0.5,
    1.0,
    2.0,
    3.0,
]

N_values = [
    1.0e7,
    5.0e7,
    1.0e8,
    2.0e8,
]

schemes = [
    "simple_kappa",
    "ARG1998",
]

dt = 0.25

aerosol_radius = 0.05e-6
aerosol_kappa = 0.3
aerosol_sigma = 1.4


# ---------------------------------------------------------
# Output
# ---------------------------------------------------------

os.makedirs(
    "data",
    exist_ok=True,
)

outfile = (
    "data/"
    "two_moment_activation_scheme_sweep_dt025.csv"
)


# ---------------------------------------------------------
# Run sweep
# ---------------------------------------------------------

rows = []

for scheme in schemes:

    for w in w_values:

        for aerosol_N in N_values:

            print()
            print(
                f"Running "
                f"scheme={scheme}, "
                f"w={w}, "
                f"N={aerosol_N:.3e}"
            )

            result = run(
                w=w,
                aerosol_N=aerosol_N,
                aerosol_radius=aerosol_radius,
                aerosol_kappa=aerosol_kappa,
                aerosol_sigma=aerosol_sigma,
                activation_scheme=scheme,
                dt=dt,
                verbose=False,
            )

            row = {
                "activation_scheme":
                    scheme,

                "w_m_s":
                    w,

                "aerosol_N_m3":
                    aerosol_N,

                "dt_s":
                    dt,

                "saturation_time_s":
                    result[
                        "saturation_time_s"
                    ],

                "activation_time_s":
                    result[
                        "activation_time_s"
                    ],

                "activation_1pct_time_s":
                    result[
                        "activation_1pct_time_s"
                    ],

                "activation_1pct_delay_s":
                    result[
                        "activation_1pct_delay_s"
                    ],

                "activation_50pct_time_s":
                    result[
                        "activation_50pct_time_s"
                    ],

                "activation_50pct_delay_s":
                    result[
                        "activation_50pct_delay_s"
                    ],

                "activation_90pct_time_s":
                    result[
                        "activation_90pct_time_s"
                    ],

                "activation_90pct_delay_s":
                    result[
                        "activation_90pct_delay_s"
                    ],

                "S_at_1pct_activation":
                    result[
                        "S_at_1pct_activation"
                    ],

                "S_at_50pct_activation":
                    result[
                        "S_at_50pct_activation"
                    ],

                "S_at_90pct_activation":
                    result[
                        "S_at_90pct_activation"
                    ],

                "SSmax":
                    result["SSmax"],

                "SSmax_percent":
                    result[
                        "SSmax_percent"
                    ],

                "activated_fraction":
                    result[
                        "activated_fraction"
                    ],

                "fraction_at_SSmax":
                    result[
                        "fraction_at_SSmax"
                    ],

                "Nc_final":
                    result["Nc_final"],

                "qc_final":
                    result["qc_final"],

                "mean_radius_um":
                    result[
                        "mean_radius_um"
                    ],

                "ARG_Sm_percent":
                    result[
                        "ARG_Sm_at_saturation_percent"
                    ],

                "ARG_Smax_prediction_percent":
                    result[
                        "ARG_Smax_prediction_percent"
                    ],

                "ARG_fraction_prediction":
                    result[
                        "ARG_fraction_prediction"
                    ],

                "water_error":
                    result["water_error"],
            }

            rows.append(
                row
            )

            print(
                f"  saturation="
                f"{row['saturation_time_s']} s"
            )

            print(
                f"  1% delay="
                f"{row['activation_1pct_delay_s']} s"
            )

            print(
                f"  50% delay="
                f"{row['activation_50pct_delay_s']} s"
            )

            print(
                f"  90% delay="
                f"{row['activation_90pct_delay_s']} s"
            )

            print(
                f"  SSmax="
                f"{row['SSmax_percent']:.4f} %"
            )

            print(
                f"  final fraction="
                f"{row['activated_fraction']:.6f}"
            )

            if scheme == "ARG1998":

                print(
                    f"  ARG predicted Smax="
                    f"{row['ARG_Smax_prediction_percent']:.4f} %"
                )

                print(
                    f"  ARG predicted fraction="
                    f"{row['ARG_fraction_prediction']:.6f}"
                )


# ---------------------------------------------------------
# Save CSV
# ---------------------------------------------------------

fieldnames = list(
    rows[0].keys()
)

with open(
    outfile,
    "w",
    newline="",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    writer.writerows(
        rows
    )


print()
print(
    "Saved:"
)
print(
    outfile
)

print()
print(
    f"Total simulations = {len(rows)}"
)

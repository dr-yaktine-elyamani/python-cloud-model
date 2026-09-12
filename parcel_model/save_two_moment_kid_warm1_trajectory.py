"""Save the experimental two-moment KiD warm1 parcel trajectory.

The purpose of this script is to save parcel-model fields at the same
30-second output interval used by the KiD warm1 NetCDF output.

This enables a diagnostic KiD-versus-parcel trajectory comparison.

This is not a validation experiment.
"""

import csv
import math
import os

from parcel_model.two_moment_warm import (
    TwoMomentConfig,
    initialize_two_moment_state,
    apply_activation_step,
    apply_condensation_step,
    liquid_supersaturation,
)

from parcel_model.aerosol import AerosolPopulation
from parcel_model.activation import check_activation

from parcel_model.abdul_razzak_1998 import (
    ARGAerosolMode,
    ammonium_sulfate_B,
    curvature_A,
    mode_critical_supersaturation,
    activated_fraction_ARG,
)

from parcel_model.run_two_moment_warm import (
    maxwell_G_liquid,
)


# ---------------------------------------------------------
# KiD warm1-matched settings
# ---------------------------------------------------------

T0 = 297.6651916503906
P0 = 99724.30419921875
QV0 = 0.014959459193050861

Z0 = 25.0

AEROSOL_N = 50e6
AEROSOL_RADIUS = 0.05e-6
AEROSOL_KAPPA = 0.3
AEROSOL_SIGMA = 1.4

DT = 0.25
T_END = 3600.0

OUTPUT_INTERVAL = 30.0

OUTPUT_FILE = (
    "data/python_two_moment_kid_warm1_trajectory.csv"
)


def kid_warm1_updraft(t):
    """KiD warm1 sinusoidal vertical velocity."""

    if 0.0 <= t < 600.0:

        return (
            2.0
            * math.sin(
                math.pi
                * t
                / 600.0
            )
        )

    return 0.0


def main():

    os.makedirs(
        "data",
        exist_ok=True,
    )

    config = TwoMomentConfig(
        dt=DT,
        t_end=T_END,
        T0=T0,
        p0=P0,
        qv0=QV0,
        aerosol_N=AEROSOL_N,
        aerosol_radius=AEROSOL_RADIUS,
        aerosol_kappa=AEROSOL_KAPPA,
        aerosol_sigma=AEROSOL_SIGMA,
    )

    state = initialize_two_moment_state(
        config
    )

    aerosol = AerosolPopulation(
        name="kid_warm1_aerosol",
        N=config.aerosol_N,
        radius=config.aerosol_radius,
        kappa=config.aerosol_kappa,
        rho_p=config.aerosol_density,
    )

    arg_B = ammonium_sulfate_B(
        soluble_mass_fraction=1.0
    )

    arg_mode = ARGAerosolMode(
        N=config.aerosol_N,
        am=config.aerosol_radius,
        sigma=config.aerosol_sigma,
        B=arg_B,
    )

    g = 9.81
    cp = 1004.0
    Rd = 287.05
    Lv = 2.5e6

    t = 0.0
    z = Z0

    next_output = OUTPUT_INTERVAL

    rows = []

    # -------------------------------------------------
    # Main integration
    # -------------------------------------------------

    while t < config.t_end:

        step_dt = min(
            DT,
            config.t_end - t,
        )

        config.dt = step_dt

        forcing_time = (
            t
            + 0.5 * step_dt
        )

        w_current = (
            kid_warm1_updraft(
                forcing_time
            )
        )

        rho_air = (
            state.p
            / (Rd * state.T)
        )

        config.rho_air = rho_air

        # ---------------------------------------------
        # Dynamics
        # ---------------------------------------------

        dz = (
            w_current
            * step_dt
        )

        z += dz

        state.p -= (
            rho_air
            * g
            * dz
        )

        state.T -= (
            g / cp
        ) * dz

        # ---------------------------------------------
        # Supersaturation before microphysics
        # ---------------------------------------------

        state.S = (
            liquid_supersaturation(
                state.qv,
                state.T,
                state.p,
            )
        )

        S_before_activation = (
            state.S
        )

        # ---------------------------------------------
        # Lognormal Köhler-threshold activation
        # ---------------------------------------------

        A_arg = curvature_A(
            state.T
        )

        Sc = (
            mode_critical_supersaturation(
                arg_mode,
                A_arg,
            )
        )

        if S_before_activation <= 0.0:

            target_fraction = 0.0

        else:

            target_fraction = (
                activated_fraction_ARG(
                    mode=arg_mode,
                    Smax=S_before_activation,
                    Sm=Sc,
                )
            )

        qc_before_microphysics = (
            state.qc
        )

        apply_activation_step(
            state,
            config,
            activated_fraction=(
                target_fraction
            ),
        )

        # ---------------------------------------------
        # Condensational growth
        # ---------------------------------------------

        config.G_liquid = (
            maxwell_G_liquid(
                state.T,
                state.p,
            )
        )

        apply_condensation_step(
            state,
            config,
        )

        # ---------------------------------------------
        # Latent heating
        # ---------------------------------------------

        delta_qc_total = (
            state.qc
            - qc_before_microphysics
        )

        if delta_qc_total != 0.0:

            state.T += (
                Lv / cp
            ) * delta_qc_total

            state.S = (
                liquid_supersaturation(
                    state.qv,
                    state.T,
                    state.p,
                )
            )

        t += step_dt

        # ---------------------------------------------
        # Save every 30 seconds
        # ---------------------------------------------

        if (
            t + 1.0e-9
            >= next_output
        ):

            if config.aerosol_N > 0.0:

                activated_fraction = (
                    state.Nc
                    / config.aerosol_N
                )

            else:

                activated_fraction = 0.0

            rows.append(
                {
                    "time_s":
                        t,

                    "parcel_height_m":
                        z,

                    "w_m_s":
                        kid_warm1_updraft(t),

                    "temperature_K":
                        state.T,

                    "pressure_Pa":
                        state.p,

                    "qv_kg_kg":
                        state.qv,

                    "S_fraction":
                        state.S,

                    "S_percent":
                        100.0 * state.S,

                    "qc_kg_kg":
                        state.qc,

                    "Nc_m3":
                        state.Nc,

                    "activated_fraction":
                        activated_fraction,

                    "mean_radius_um":
                        state.mean_radius
                        * 1.0e6,

                    "G_liquid_m2_s":
                        config.G_liquid,
                }
            )

            next_output += (
                OUTPUT_INTERVAL
            )

    # -------------------------------------------------
    # Write CSV
    # -------------------------------------------------

    fieldnames = list(
        rows[0].keys()
    )

    with open(
        OUTPUT_FILE,
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
        "Python warm1 trajectory saved"
    )

    print(
        "-----------------------------"
    )

    print(
        f"Saved: {OUTPUT_FILE}"
    )

    print(
        f"Rows: {len(rows)}"
    )

    print()

    print(
        "Selected times"
    )

    print(
        "--------------"
    )

    for target in (
        300.0,
        390.0,
        420.0,
        480.0,
        570.0,
        600.0,
    ):

        row = min(
            rows,
            key=lambda r: abs(
                r["time_s"]
                - target
            ),
        )

        print()

        print(
            f"time = "
            f"{row['time_s']:.1f} s"
        )

        print(
            f"z = "
            f"{row['parcel_height_m']:.2f} m"
        )

        print(
            f"w = "
            f"{row['w_m_s']:.4f} m/s"
        )

        print(
            f"T = "
            f"{row['temperature_K']:.3f} K"
        )

        print(
            f"qv = "
            f"{row['qv_kg_kg']:.6e}"
        )

        print(
            f"S = "
            f"{row['S_percent']:.4f} %"
        )

        print(
            f"qc = "
            f"{row['qc_kg_kg']:.6e}"
        )

        print(
            f"Nc = "
            f"{row['Nc_m3']:.6e}"
        )

        print(
            f"activated fraction = "
            f"{row['activated_fraction']:.6f}"
        )


if __name__ == "__main__":

    main()

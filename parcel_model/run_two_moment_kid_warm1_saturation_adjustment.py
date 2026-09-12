"""
KiD warm1 parcel experiment using saturation adjustment.

This is a diagnostic experiment.

It keeps the same:
- KiD warm1 initial thermodynamic state
- sinusoidal vertical-velocity forcing
- lognormal Kohler activation

but replaces explicit Maxwell condensational growth with a
saturation-adjustment calculation inspired by the condensation
treatment in Thompson09.

This is NOT a full reproduction of Thompson09 and is NOT a
validation of the parcel model.
"""

import math

from parcel_model.two_moment_warm import (
    TwoMomentConfig,
    initialize_two_moment_state,
    apply_activation_step,
    liquid_supersaturation,
)

from parcel_model.abdul_razzak_1998 import (
    ARGAerosolMode,
    ammonium_sulfate_B,
    curvature_A,
    mode_critical_supersaturation,
    activated_fraction_ARG,
)


# ---------------------------------------------------------
# Constants
# ---------------------------------------------------------

RD = 287.05
CP = 1004.0
G = 9.81
LV = 2.5e6

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


# ---------------------------------------------------------
# KiD warm1 forcing
# ---------------------------------------------------------

def kid_warm1_w(t):
    """
    KiD warm1 vertical velocity.

    w = 2 sin(pi t / 600), t < 600 s
    w = 0 afterwards.
    """

    if t < 600.0:
        return 2.0 * math.sin(
            math.pi * t / 600.0
        )

    return 0.0


# ---------------------------------------------------------
# Saturation adjustment
# ---------------------------------------------------------

def saturation_adjustment(
    state,
    config,
):
    """
    Adjust vapour/cloud water toward liquid saturation.

    Diagnostic saturation-adjustment experiment.

    The adjustment includes latent heating and solves for
    the condensed-water increment using Newton iterations.

    This is inspired by the logic used in Thompson09,
    but is not intended to be a line-by-line reproduction.
    """

    if state.S <= 0.0:
        return state

    qv_initial = state.qv
    qc_initial = state.qc
    T_initial = state.T

    # Initial estimate of condensed mass.
    dq = 0.0

    for _ in range(10):

        # Trial thermodynamic state after condensation.
        qv_trial = max(
            qv_initial - dq,
            1.0e-12,
        )

        T_trial = (
            T_initial
            + LV * dq / CP
        )

        S_trial = liquid_supersaturation(
            qv_trial,
            T_trial,
            state.p,
        )

        # Close enough to saturation.
        if abs(S_trial) < 1.0e-10:
            break

        # Numerical derivative dS/ddq.
        eps = 1.0e-9

        dq_plus = dq + eps

        qv_plus = max(
            qv_initial - dq_plus,
            1.0e-12,
        )

        T_plus = (
            T_initial
            + LV * dq_plus / CP
        )

        S_plus = liquid_supersaturation(
            qv_plus,
            T_plus,
            state.p,
        )

        derivative = (
            S_plus - S_trial
        ) / eps

        if abs(derivative) < 1.0e-14:
            break

        dq_new = (
            dq
            - S_trial / derivative
        )

        # Condensation cannot exceed available vapour.
        dq = min(
            max(dq_new, 0.0),
            qv_initial,
        )

    state.qv = (
        qv_initial - dq
    )

    state.qc = (
        qc_initial + dq
    )

    state.T = (
        T_initial
        + LV * dq / CP
    )

    state.S = liquid_supersaturation(
        state.qv,
        state.T,
        state.p,
    )

    return state


# ---------------------------------------------------------
# Main experiment
# ---------------------------------------------------------

def run():

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

    arg_B = ammonium_sulfate_B(
        soluble_mass_fraction=1.0
    )

    arg_mode = ARGAerosolMode(
        N=config.aerosol_N,
        am=config.aerosol_radius,
        sigma=config.aerosol_sigma,
        B=arg_B,
    )

    t = 0.0
    z = Z0

    saturation_time = None
    activation_time = None

    max_S = state.S
    time_of_max_S = 0.0

    max_qc = state.qc
    time_of_max_qc = 0.0

    initial_total_water = (
        state.qv + state.qc
    )

    selected_times = {
        300.0,
        390.0,
        420.0,
        480.0,
        570.0,
        600.0,
    }

    print()
    print(
        "Python warm1 saturation-adjustment experiment"
    )
    print(
        "---------------------------------------------"
    )

    while t < T_END - 0.5 * DT:

        # Midpoint forcing.
        t_mid = t + 0.5 * DT

        w = kid_warm1_w(
            t_mid
        )

        rho_air = (
            state.p
            / (RD * state.T)
        )

        config.rho_air = rho_air

        # Parcel displacement.
        dz = w * DT

        z += dz

        # Dry adiabatic pressure change.
        dp = (
            -rho_air
            * G
            * dz
        )

        state.p += dp

        # Dry adiabatic cooling.
        state.T += (
            -G / CP
            * dz
        )

        # Supersaturation before microphysics.
        state.S = liquid_supersaturation(
            state.qv,
            state.T,
            state.p,
        )

        S_before_activation = state.S

        if (
            saturation_time is None
            and S_before_activation >= 0.0
        ):
            saturation_time = (
                t + DT
            )

        # -------------------------------------------------
        # Lognormal Kohler threshold activation
        # -------------------------------------------------

        A_arg = curvature_A(
            state.T
        )

        Sc = mode_critical_supersaturation(
            arg_mode,
            A_arg,
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

        Nc_before = state.Nc

        apply_activation_step(
            state,
            config,
            target_fraction,
        )

        if (
            activation_time is None
            and state.Nc > Nc_before
        ):
            activation_time = (
                t + DT
            )

        # -------------------------------------------------
        # Saturation adjustment instead of Maxwell growth
        # -------------------------------------------------

        saturation_adjustment(
            state,
            config,
        )

        # Update diagnostics.
        if state.S > max_S:

            max_S = state.S
            time_of_max_S = (
                t + DT
            )

        if state.qc > max_qc:

            max_qc = state.qc
            time_of_max_qc = (
                t + DT
            )

        t += DT

        # Selected diagnostic times.
        if any(
            abs(t - target) < 0.5 * DT
            for target in selected_times
        ):

            print()
            print(
                f"time = {t:.1f} s"
            )

            print(
                f"z = {z:.2f} m"
            )

            print(
                f"w = {kid_warm1_w(t):.4f} m/s"
            )

            print(
                f"T = {state.T:.3f} K"
            )

            print(
                f"qv = {state.qv:.6e}"
            )

            print(
                f"S = {100.0 * state.S:.6f} %"
            )

            print(
                f"qc = {state.qc:.6e}"
            )

            print(
                f"Nc = {state.Nc:.6e}"
            )

            if config.aerosol_N > 0.0:

                fraction = (
                    state.Nc
                    / config.aerosol_N
                )

            else:

                fraction = 0.0

            print(
                "activated fraction = "
                f"{fraction:.6f}"
            )

    final_total_water = (
        state.qv
        + state.qc
    )

    water_error = (
        final_total_water
        - initial_total_water
    )

    print()
    print(
        "Final diagnostics"
    )
    print(
        "-----------------"
    )

    print(
        "comparison_status = "
        "diagnostic experiment; not validation"
    )

    print(
        "condensation_scheme = "
        "saturation_adjustment"
    )

    print(
        "activation_scheme = "
        "lognormal_kohler"
    )

    print(
        "saturation_time_s =",
        saturation_time,
    )

    print(
        "activation_time_s =",
        activation_time,
    )

    print(
        "SSmax_percent =",
        100.0 * max_S,
    )

    print(
        "time_of_SSmax_s =",
        time_of_max_S,
    )

    print(
        "max_qc =",
        max_qc,
    )

    print(
        "time_of_max_qc_s =",
        time_of_max_qc,
    )

    print(
        "Nc_final =",
        state.Nc,
    )

    print(
        "qc_final =",
        state.qc,
    )

    print(
        "qv_final =",
        state.qv,
    )

    print(
        "T_final_K =",
        state.T,
    )

    print(
        "final_height_m =",
        z,
    )

    print(
        "water_error =",
        water_error,
    )


if __name__ == "__main__":
    run()

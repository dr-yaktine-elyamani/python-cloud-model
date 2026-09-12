"""
Diagnostic KiD warm1 fixed-environment experiment.

Purpose
-------
This experiment isolates the effect of the thermodynamic framework used
for the KiD warm1 comparison.

It is NOT a replacement for the Lagrangian parcel model and it is NOT
a validation experiment.

KiD warm1 uses:
    L_FIX_THETA = True
    L_PUPDATE   = False

Therefore this diagnostic follows the warm1 vertical trajectory while
diagnosing temperature and pressure from a fixed RICO-like environmental
theta/Exner profile rather than carrying an interactive parcel
temperature/pressure history.

Microphysics is intentionally kept simple:
    - saturation adjustment for positive supersaturation
    - no rain/autoconversion/sedimentation
    - no accumulated latent-heating feedback on environmental theta

The purpose is to determine how much of the KiD/Python difference comes
from the thermodynamic framework before changing microphysical
parameters.
"""

import math

from parcel_model.two_moment_warm import liquid_supersaturation


# ---------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------

DT = 0.25
T_END = 600.0

Z0 = 25.0

P0_REF = 100000.0
P_SURF = 100000.0

RD = 287.05
CP = 1004.0
G = 9.81

R_ON_CP = RD / CP

# Approximate first KiD-output water-vapour value used in the previous
# matched warm1 experiments.
QV0 = 0.014959459193050861


# ---------------------------------------------------------------------
# Warm1 forcing
# ---------------------------------------------------------------------

def warm1_w(t):
    """
    KiD warm1 vertical-velocity forcing.

    w = 2 sin(pi t / 600) for t < 600 s
    w = 0 afterwards
    """
    if t < 600.0:
        return 2.0 * math.sin(math.pi * t / 600.0)
    return 0.0


# ---------------------------------------------------------------------
# Fixed RICO theta profile
# ---------------------------------------------------------------------

def rico_theta(z):
    """
    RICO potential-temperature profile used by KiD when l_fix_theta=True.

    Control points:
        z =    0 m : theta = 297.90 K
        z =  740 m : theta = 297.90 K
        z = 3260 m : theta = 312.66 K
    """

    if z <= 740.0:
        return 297.90

    if z >= 3260.0:
        return 312.66

    fraction = (z - 740.0) / (3260.0 - 740.0)

    return 297.90 + fraction * (312.66 - 297.90)


# ---------------------------------------------------------------------
# Fixed hydrostatic Exner profile
# ---------------------------------------------------------------------

def build_environment(
    z_max=3000.0,
    dz=1.0,
):
    """
    Construct a simple fixed hydrostatic Exner profile from the fixed
    RICO theta profile.

    This follows the KiD hydrostatic relationship structurally:

        dPi/dz = -g / (cp * theta)

    The profile is diagnostic and is not claimed to reproduce every
    detail of the KiD grid construction.
    """

    z_values = [0.0]
    theta_values = [rico_theta(0.0)]

    exner_surface = (P_SURF / P0_REF) ** R_ON_CP
    exner_values = [exner_surface]

    z = 0.0

    while z < z_max:

        z_new = min(z + dz, z_max)

        theta_old = rico_theta(z)
        theta_new = rico_theta(z_new)
        theta_mean = 0.5 * (theta_old + theta_new)

        dpi = (
            G
            * (z_new - z)
            / (CP * theta_mean)
        )

        exner_new = exner_values[-1] - dpi

        z_values.append(z_new)
        theta_values.append(theta_new)
        exner_values.append(exner_new)

        z = z_new

    return z_values, theta_values, exner_values


Z_ENV, THETA_ENV, EXNER_ENV = build_environment()


def interpolate(x, xs, ys):
    """Simple linear interpolation."""

    if x <= xs[0]:
        return ys[0]

    if x >= xs[-1]:
        return ys[-1]

    i = int(x)

    if i >= len(xs) - 1:
        i = len(xs) - 2

    x0 = xs[i]
    x1 = xs[i + 1]

    y0 = ys[i]
    y1 = ys[i + 1]

    if x1 == x0:
        return y0

    f = (x - x0) / (x1 - x0)

    return y0 + f * (y1 - y0)


def environmental_state(z):
    """
    Return fixed-environment theta, Exner, T and pressure at height z.
    """

    theta = rico_theta(z)

    exner = interpolate(
        z,
        Z_ENV,
        EXNER_ENV,
    )

    T = theta * exner

    p = P0_REF * exner ** (1.0 / R_ON_CP)

    return theta, exner, T, p


# ---------------------------------------------------------------------
# Saturation adjustment with fixed environmental temperature
# ---------------------------------------------------------------------

def saturation_adjust_fixed_T(
    qv,
    qc,
    T,
    p,
):
    """
    Condense excess vapour to saturation at fixed environmental T.

    No latent-heating feedback is retained because this experiment is
    designed specifically to isolate the fixed-theta thermodynamic
    framework.

    Total qv + qc is conserved.
    """

    S = liquid_supersaturation(
        qv,
        T,
        p,
    )

    if S <= 0.0:
        return qv, qc, S

    qv_initial = qv

    # Find saturated qv by bisection.
    low = 0.0
    high = qv_initial

    for _ in range(60):

        mid = 0.5 * (low + high)

        S_mid = liquid_supersaturation(
            mid,
            T,
            p,
        )

        if S_mid > 0.0:
            high = mid
        else:
            low = mid

    qv_sat = 0.5 * (low + high)

    dq = max(
        0.0,
        qv_initial - qv_sat,
    )

    qv = qv_initial - dq
    qc = qc + dq

    S = liquid_supersaturation(
        qv,
        T,
        p,
    )

    return qv, qc, S


# ---------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------

def run():

    t = 0.0
    z = Z0

    qv = QV0
    qc = 0.0

    initial_total_water = qv + qc

    saturation_time = None
    saturation_height = None

    max_qc = 0.0
    time_of_max_qc = 0.0

    selected_times = {
        300.0,
        390.0,
        420.0,
        480.0,
        570.0,
        600.0,
    }

    print(
        "Python warm1 fixed-environment diagnostic"
    )
    print(
        "-----------------------------------------"
    )

    while t < T_END - 0.5 * DT:

        t_mid = t + 0.5 * DT
        w = warm1_w(t_mid)

        z += w * DT
        t += DT

        theta, exner, T, p = environmental_state(z)

        S_before = liquid_supersaturation(
            qv,
            T,
            p,
        )

        if (
            saturation_time is None
            and S_before >= 0.0
        ):
            saturation_time = t
            saturation_height = z

        qv, qc, S = saturation_adjust_fixed_T(
            qv,
            qc,
            T,
            p,
        )

        if qc > max_qc:
            max_qc = qc
            time_of_max_qc = t

        for target in selected_times:
            if abs(t - target) < 0.5 * DT:

                print()
                print(
                    f"time = {t:.1f} s"
                )
                print(
                    f"z = {z:.2f} m"
                )
                print(
                    f"w = {warm1_w(t):.4f} m/s"
                )
                print(
                    f"theta = {theta:.3f} K"
                )
                print(
                    f"exner = {exner:.8f}"
                )
                print(
                    f"T = {T:.3f} K"
                )
                print(
                    f"p = {p:.2f} Pa"
                )
                print(
                    f"qv = {qv:.6e}"
                )
                print(
                    f"S = {100.0*S:.6f} %"
                )
                print(
                    f"qc = {qc:.6e}"
                )

                break

    theta, exner, T, p = environmental_state(z)

    final_total_water = qv + qc

    water_error = (
        final_total_water
        - initial_total_water
    )

    print()
    print("Final diagnostics")
    print("-----------------")

    print(
        "comparison_status = "
        "fixed-environment diagnostic; not validation"
    )
    print(
        "thermodynamic_framework = "
        "fixed RICO theta + fixed hydrostatic Exner"
    )
    print(
        "condensation_scheme = "
        "fixed-temperature saturation adjustment"
    )

    print(
        "saturation_time_s =",
        saturation_time,
    )
    print(
        "saturation_height_m =",
        saturation_height,
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
        "qc_final =",
        qc,
    )
    print(
        "qv_final =",
        qv,
    )
    print(
        "T_final_K =",
        T,
    )
    print(
        "p_final_Pa =",
        p,
    )
    print(
        "theta_final_K =",
        theta,
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

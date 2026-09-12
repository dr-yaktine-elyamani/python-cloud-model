"""KiD warm1-matched experimental parcel run.

This runner is intended for a diagnostic comparison with the KiD-A
GCSS warm1 case (icase=101).

Matched / approximately matched settings
----------------------------------------
KiD forcing:
    w(t) = 2 sin(pi t / 600), for t < 600 s
    w(t) = 0,                for t >= 600 s

Aerosol:
    N     = 50e6 m-3
    rd    = 0.05e-6 m
    sigma = 1.4

Initial thermodynamic state:
    Taken from the first KiD output level at z=25 m, t=30 s:
        T  = 297.6651916503906 K
        p  = 997.2430419921875 mb
        qv = 0.014959459193050861 kg/kg

The thermodynamic state is therefore an approximate matched
initialisation, not an exact reconstruction of KiD at t=0.

Important:
KiD is a 1-D Eulerian cloud model with Thompson09 microphysics,
advection, divergence, collision-coalescence and sedimentation.
This parcel model is much simpler. Results are diagnostic
comparisons, not model validation.
"""

import math

from parcel_model.two_moment_warm import (
    TwoMomentConfig,
    initialize_two_moment_state,
    apply_activation_step,
    apply_condensation_step,
    liquid_supersaturation,
    total_water,
)

from parcel_model.aerosol import AerosolPopulation
from parcel_model.activation import check_activation

from parcel_model.abdul_razzak_1998 import (
    ARGAerosolMode,
    ammonium_sulfate_B,
    curvature_A,
    mode_critical_supersaturation,
    activated_fraction_ARG,
    smax_ARG_single,
)

from parcel_model.run_two_moment_warm import (
    maxwell_G_liquid,
    arg_growth_coefficient,
)


# ---------------------------------------------------------
# KiD warm1 reference settings
# ---------------------------------------------------------

KID_Z0 = 25.0

KID_T0 = 297.6651916503906

KID_P0 = (
    997.2430419921875
    * 100.0
)

KID_QV0 = 0.014959459193050861

KID_W_PEAK = 2.0

KID_FORCING_DURATION = 600.0


def kid_warm1_updraft(t):
    """KiD warm1 prescribed vertical velocity [m/s]."""

    if 0.0 <= t < KID_FORCING_DURATION:

        return (
            KID_W_PEAK
            * math.sin(
                math.pi
                * t
                / KID_FORCING_DURATION
            )
        )

    return 0.0


def run(
    aerosol_N=50e6,
    aerosol_radius=0.05e-6,
    aerosol_kappa=0.3,
    aerosol_sigma=1.4,
    activation_scheme="lognormal_kohler",
    arg_soluble_mass_fraction=1.0,
    dt=0.25,
    t_end=3600.0,
    verbose=True,
):

    # -------------------------------------------------
    # Input validation
    # -------------------------------------------------

    if aerosol_N < 0.0:
        raise ValueError(
            "aerosol_N must be nonnegative"
        )

    if aerosol_radius <= 0.0:
        raise ValueError(
            "aerosol_radius must be positive"
        )

    if aerosol_kappa < 0.0:
        raise ValueError(
            "aerosol_kappa must be nonnegative"
        )

    if aerosol_sigma <= 0.0:
        raise ValueError(
            "aerosol_sigma must be positive"
        )

    if not (
        0.0
        <= arg_soluble_mass_fraction
        <= 1.0
    ):
        raise ValueError(
            "arg_soluble_mass_fraction must "
            "be between 0 and 1"
        )

    if dt <= 0.0:
        raise ValueError(
            "dt must be positive"
        )

    if t_end <= 0.0:
        raise ValueError(
            "t_end must be positive"
        )

    if activation_scheme not in (
        "simple_kappa",
        "lognormal_kohler",
        "ARG1998",
    ):
        raise ValueError(
            "activation_scheme must be "
            "'simple_kappa', "
            "'lognormal_kohler', "
            "or legacy alias 'ARG1998'"
        )

    # -------------------------------------------------
    # Backward-compatible activation naming
    # -------------------------------------------------

    activation_scheme_input = (
        activation_scheme
    )

    if activation_scheme == "ARG1998":

        parcel_activation_scheme = (
            "lognormal_kohler"
        )

    else:

        parcel_activation_scheme = (
            activation_scheme
        )

    # -------------------------------------------------
    # Matched KiD warm1 configuration
    # -------------------------------------------------

    config = TwoMomentConfig(
        dt=dt,
        t_end=t_end,

        T0=KID_T0,
        p0=KID_P0,
        qv0=KID_QV0,

        aerosol_N=aerosol_N,
        aerosol_radius=aerosol_radius,
        aerosol_kappa=aerosol_kappa,
        aerosol_sigma=aerosol_sigma,
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

    # -------------------------------------------------
    # ARG aerosol mode
    # -------------------------------------------------

    arg_B = ammonium_sulfate_B(
        soluble_mass_fraction=(
            arg_soluble_mass_fraction
        )
    )

    arg_mode = ARGAerosolMode(
        N=config.aerosol_N,
        am=config.aerosol_radius,
        sigma=config.aerosol_sigma,
        B=arg_B,
    )

    # -------------------------------------------------
    # Constants
    # -------------------------------------------------

    g = 9.81
    cp = 1004.0
    Rd = 287.05
    Lv = 2.5e6

    # -------------------------------------------------
    # Initial parcel position/time
    # -------------------------------------------------

    t = 0.0
    z = KID_Z0

    initial_water = total_water(
        state
    )

    initial_T = state.T
    initial_p = state.p
    initial_qv = state.qv

    # -------------------------------------------------
    # Saturation / activation diagnostics
    # -------------------------------------------------

    saturation_time = None
    saturation_height = None
    w_at_saturation = None

    activation_time = None
    activation_height = None

    S_at_activation = None
    Sc_at_activation = None

    # -------------------------------------------------
    # Activation thresholds
    # -------------------------------------------------

    activation_1pct_time = None
    activation_50pct_time = None
    activation_90pct_time = None

    activation_1pct_height = None
    activation_50pct_height = None
    activation_90pct_height = None

    S_at_1pct_activation = None
    S_at_50pct_activation = None
    S_at_90pct_activation = None

    activation_1pct_delay = None
    activation_50pct_delay = None
    activation_90pct_delay = None

    # -------------------------------------------------
    # Supersaturation diagnostics
    # -------------------------------------------------

    SSmax = state.S

    time_of_SSmax = 0.0
    height_of_SSmax = z
    w_at_SSmax = kid_warm1_updraft(0.0)

    fraction_at_SSmax = (
        state.Nc / config.aerosol_N
        if config.aerosol_N > 0.0
        else 0.0
    )

    # -------------------------------------------------
    # Maximum cloud water
    # -------------------------------------------------

    max_qc = state.qc
    time_of_max_qc = 0.0
    height_of_max_qc = z

    # -------------------------------------------------
    # ARG1998 standalone analytical benchmark
    # -------------------------------------------------

    ARG_Smax_prediction = None
    ARG_Sm_at_saturation = None
    ARG_fraction_prediction = None
    ARG_G_at_saturation = None

    ARG_G_representative_radius = (
        config.initial_droplet_radius
    )

    next_print = 0.0

    # -------------------------------------------------
    # Main time loop
    # -------------------------------------------------

    while t < config.t_end:

        step_dt = min(
            dt,
            config.t_end - t,
        )

        if step_dt <= 0.0:
            break

        config.dt = step_dt

        # Use velocity at the midpoint of the step.
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
        # Parcel displacement
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
        # Pre-microphysics supersaturation
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

        diagnostic_time = (
            t + step_dt
        )

        # ---------------------------------------------
        # First saturation crossing
        # ---------------------------------------------

        if (
            saturation_time is None
            and S_before_activation >= 0.0
        ):

            saturation_time = (
                diagnostic_time
            )

            saturation_height = z

            w_at_saturation = (
                w_current
            )

            A_arg = curvature_A(
                state.T
            )

            ARG_Sm_at_saturation = (
                mode_critical_supersaturation(
                    arg_mode,
                    A_arg,
                )
            )

            ARG_G_at_saturation = (
                arg_growth_coefficient(
                    state.T,
                    state.p,
                    ARG_G_representative_radius,
                )
            )

            arg_smax_result = (
                smax_ARG_single(
                    mode=arg_mode,
                    T=state.T,
                    p=state.p,
                    ps=(
                        __import__(
                            "parcel_model.thermodynamics",
                            fromlist=["esat_water"],
                        ).esat_water(
                            state.T
                        )
                    ),
                    V=w_current,
                    G=ARG_G_at_saturation,
                )
            )

            ARG_Smax_prediction = (
                arg_smax_result["Smax"]
            )

            ARG_fraction_prediction = (
                activated_fraction_ARG(
                    mode=arg_mode,
                    Smax=ARG_Smax_prediction,
                    Sm=ARG_Sm_at_saturation,
                )
            )

        # ---------------------------------------------
        # Pre-microphysics SSmax
        # ---------------------------------------------

        if S_before_activation > SSmax:

            SSmax = (
                S_before_activation
            )

            time_of_SSmax = (
                diagnostic_time
            )

            height_of_SSmax = z

            w_at_SSmax = (
                w_current
            )

            fraction_at_SSmax = (
                state.Nc
                / config.aerosol_N
                if config.aerosol_N > 0.0
                else 0.0
            )

        # ---------------------------------------------
        # Activation scheme
        # ---------------------------------------------

        if (
            parcel_activation_scheme
            == "simple_kappa"
        ):

            activated, Sc = (
                check_activation(
                    S=S_before_activation,
                    aerosol=aerosol,
                    T=state.T,
                )
            )

            target_fraction = (
                aerosol.activated_fraction
            )

        else:

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

        # ---------------------------------------------
        # Activation source
        # ---------------------------------------------

        Nc_before = state.Nc

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
        # Current activated fraction
        # ---------------------------------------------

        if config.aerosol_N > 0.0:

            current_fraction = (
                state.Nc
                / config.aerosol_N
            )

        else:

            current_fraction = 0.0

        # ---------------------------------------------
        # First non-zero activation
        # ---------------------------------------------

        if (
            activation_time is None
            and Nc_before <= 0.0
            and state.Nc > 0.0
        ):

            activation_time = (
                diagnostic_time
            )

            activation_height = z

            S_at_activation = (
                S_before_activation
            )

            Sc_at_activation = Sc

        # ---------------------------------------------
        # 1% activation
        # ---------------------------------------------

        if (
            activation_1pct_time is None
            and current_fraction >= 0.01
        ):

            activation_1pct_time = (
                diagnostic_time
            )

            activation_1pct_height = z

            S_at_1pct_activation = (
                S_before_activation
            )

            if saturation_time is not None:

                activation_1pct_delay = (
                    activation_1pct_time
                    - saturation_time
                )

        # ---------------------------------------------
        # 50% activation
        # ---------------------------------------------

        if (
            activation_50pct_time is None
            and current_fraction >= 0.50
        ):

            activation_50pct_time = (
                diagnostic_time
            )

            activation_50pct_height = z

            S_at_50pct_activation = (
                S_before_activation
            )

            if saturation_time is not None:

                activation_50pct_delay = (
                    activation_50pct_time
                    - saturation_time
                )

        # ---------------------------------------------
        # 90% activation
        # ---------------------------------------------

        if (
            activation_90pct_time is None
            and current_fraction >= 0.90
        ):

            activation_90pct_time = (
                diagnostic_time
            )

            activation_90pct_height = z

            S_at_90pct_activation = (
                S_before_activation
            )

            if saturation_time is not None:

                activation_90pct_delay = (
                    activation_90pct_time
                    - saturation_time
                )

        # ---------------------------------------------
        # Dynamic Maxwell growth coefficient
        # ---------------------------------------------

        config.G_liquid = (
            maxwell_G_liquid(
                state.T,
                state.p,
            )
        )

        # ---------------------------------------------
        # Condensation / evaporation
        # ---------------------------------------------

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

        # ---------------------------------------------
        # Post-microphysics SSmax
        # ---------------------------------------------

        if state.S > SSmax:

            SSmax = state.S

            time_of_SSmax = (
                diagnostic_time
            )

            height_of_SSmax = z

            w_at_SSmax = (
                w_current
            )

            fraction_at_SSmax = (
                state.Nc
                / config.aerosol_N
                if config.aerosol_N > 0.0
                else 0.0
            )

        # ---------------------------------------------
        # Maximum parcel cloud-water mixing ratio
        # ---------------------------------------------

        if state.qc > max_qc:

            max_qc = state.qc

            time_of_max_qc = (
                diagnostic_time
            )

            height_of_max_qc = z

        t += step_dt

        # ---------------------------------------------
        # Screen output
        # ---------------------------------------------

        if (
            verbose
            and t >= next_print
        ):

            print(
                f"t={t:7.1f} s  "
                f"z={z:8.2f} m  "
                f"w={w_current:6.3f} m/s  "
                f"T={state.T:7.2f} K  "
                f"S={100.0 * state.S:8.4f} %  "
                f"Nc={state.Nc:12.4e}  "
                f"qc={state.qc:10.4e}"
            )

            next_print += 100.0

    # -------------------------------------------------
    # Final diagnostics
    # -------------------------------------------------

    if config.aerosol_N > 0.0:

        activated_fraction = (
            state.Nc
            / config.aerosol_N
        )

    else:

        activated_fraction = 0.0

    if (
        saturation_time is not None
        and activation_time is not None
    ):

        time_from_saturation_to_activation = (
            activation_time
            - saturation_time
        )

    else:

        time_from_saturation_to_activation = (
            None
        )

    final_water = total_water(
        state
    )

    water_error = (
        final_water
        - initial_water
    )

    # -------------------------------------------------
    # Results
    # -------------------------------------------------

    result = {

        "case":
            "KiD warm1 matched parcel experiment",

        "comparison_status":
            "diagnostic comparison; not validation",

        "activation_scheme":
            activation_scheme_input,

        "parcel_activation_scheme":
            parcel_activation_scheme,

        "ARG1998_role":
            "standalone analytical benchmark",

        "forcing":
            "w=2*sin(pi*t/600) for t<600 s; w=0 afterwards",

        "w_peak_m_s":
            KID_W_PEAK,

        "forcing_duration_s":
            KID_FORCING_DURATION,

        "initial_height_m":
            KID_Z0,

        "initial_T_K":
            initial_T,

        "initial_p_Pa":
            initial_p,

        "initial_qv":
            initial_qv,

        "aerosol_N_m3":
            aerosol_N,

        "aerosol_radius_m":
            aerosol_radius,

        "aerosol_kappa":
            aerosol_kappa,

        "aerosol_sigma":
            aerosol_sigma,

        "dt_s":
            dt,

        "requested_t_end_s":
            t_end,

        # Saturation
        "saturation_time_s":
            saturation_time,

        "saturation_height_m":
            saturation_height,

        "w_at_saturation_m_s":
            w_at_saturation,

        # Activation
        "activation_time_s":
            activation_time,

        "activation_height_m":
            activation_height,

        "time_to_activation_s":
            time_from_saturation_to_activation,

        "S_at_activation":
            S_at_activation,

        "Sc_at_activation":
            Sc_at_activation,

        # 1%
        "activation_1pct_time_s":
            activation_1pct_time,

        "activation_1pct_height_m":
            activation_1pct_height,

        "activation_1pct_delay_s":
            activation_1pct_delay,

        "S_at_1pct_activation":
            S_at_1pct_activation,

        # 50%
        "activation_50pct_time_s":
            activation_50pct_time,

        "activation_50pct_height_m":
            activation_50pct_height,

        "activation_50pct_delay_s":
            activation_50pct_delay,

        "S_at_50pct_activation":
            S_at_50pct_activation,

        # 90%
        "activation_90pct_time_s":
            activation_90pct_time,

        "activation_90pct_height_m":
            activation_90pct_height,

        "activation_90pct_delay_s":
            activation_90pct_delay,

        "S_at_90pct_activation":
            S_at_90pct_activation,

        # Supersaturation
        "SSmax":
            SSmax,

        "SSmax_percent":
            100.0 * SSmax,

        "time_of_SSmax_s":
            time_of_SSmax,

        "height_of_SSmax_m":
            height_of_SSmax,

        "w_at_SSmax_m_s":
            w_at_SSmax,

        "fraction_at_SSmax":
            fraction_at_SSmax,

        # ARG benchmark
        "ARG_Sm_at_saturation":
            ARG_Sm_at_saturation,

        "ARG_Sm_at_saturation_percent":
            (
                100.0
                * ARG_Sm_at_saturation
                if ARG_Sm_at_saturation is not None
                else None
            ),

        "ARG_Smax_prediction":
            ARG_Smax_prediction,

        "ARG_Smax_prediction_percent":
            (
                100.0
                * ARG_Smax_prediction
                if ARG_Smax_prediction is not None
                else None
            ),

        "ARG_fraction_prediction":
            ARG_fraction_prediction,

        "ARG_G_at_saturation":
            ARG_G_at_saturation,

        "ARG_G_representative_radius_m":
            ARG_G_representative_radius,

        "ARG_G_radius_is_provisional":
            True,

        # Cloud water
        "max_qc":
            max_qc,

        "time_of_max_qc_s":
            time_of_max_qc,

        "height_of_max_qc_m":
            height_of_max_qc,

        # Final state
        "activated_fraction":
            activated_fraction,

        "Nc_final":
            state.Nc,

        "qc_final":
            state.qc,

        "qv_final":
            state.qv,

        "mean_radius_um":
            state.mean_radius
            * 1.0e6,

        "T_final_K":
            state.T,

        "p_final_Pa":
            state.p,

        "final_height_m":
            z,

        "height_gain_m":
            z - KID_Z0,

        "simulation_time_s":
            t,

        "water_error":
            water_error,
    }

    # -------------------------------------------------
    # Print final diagnostics
    # -------------------------------------------------

    if verbose:

        print()
        print(
            "KiD warm1 matched parcel diagnostics"
        )
        print(
            "------------------------------------"
        )

        for key, value in result.items():

            print(
                f"{key} = {value}"
            )

    return result


if __name__ == "__main__":

    run()

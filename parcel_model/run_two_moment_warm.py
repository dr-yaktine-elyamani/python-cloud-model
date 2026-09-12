"""Experimental warm-cloud two-moment parcel runner.

This is a provisional warm-cloud model.

Parcel activation options
-------------------------

1. simple_kappa
   Simplified monodisperse kappa-Kohler activation.

2. lognormal_kohler
   Lognormal critical-supersaturation threshold activation.
   The aerosol critical-supersaturation distribution is based
   on the single-mode formulation used by Abdul-Razzak,
   Ghan & Rivera-Carpio (1998).

   This is NOT the complete analytical ARG1998 activation
   parameterization. It uses the instantaneous parcel
   supersaturation as a cumulative activation threshold.

For backward compatibility, the old name "ARG1998" is still
accepted as an alias for "lognormal_kohler".

Separately, the runner computes a standalone ARG1998
analytical benchmark at first saturation:
    - ARG_Smax_prediction
    - ARG_fraction_prediction

The model is not yet a validated aerosol-resolved
microphysics scheme.
"""

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
from parcel_model.thermodynamics import esat_water

from parcel_model.abdul_razzak_1998 import (
    ARGAerosolMode,
    ammonium_sulfate_B,
    curvature_A,
    mode_critical_supersaturation,
    activated_fraction_ARG,
    growth_coefficient_ARG1998,
    smax_ARG_single,
)


Rv = 461.5
rho_l = 1000.0
Lv = 2.5e6

P_REF = 101325.0


def diffusivity_water_vapour(T, P):
    """Water-vapour diffusivity in air [m2/s]."""

    D0 = 2.11e-5

    return (
        D0
        * (T / 273.15) ** 1.94
        * (P_REF / P)
    )


def thermal_conductivity_air(T):
    """Approximate thermal conductivity of air [W/m/K]."""

    k0 = 0.024

    return (
        k0
        * (T / 273.15) ** 0.9
    )


def maxwell_G_liquid(T, P):
    """Simplified Maxwell liquid growth coefficient [m2/s]."""

    esw = esat_water(T)

    D = diffusivity_water_vapour(T, P)
    k = thermal_conductivity_air(T)

    A = (
        rho_l
        * Rv
        * T
        / (D * esw)
    )

    B = (
        rho_l
        * Lv**2
        / (k * Rv * T**2)
    )

    return 1.0 / (A + B)


def arg_growth_coefficient(
    T,
    P,
    representative_radius,
):
    """
    ARG-style condensational growth coefficient.

    The kinetic corrections in the current implementation
    require a representative wet-droplet radius.

    The choice of representative radius is currently
    provisional and is retained only for the standalone
    ARG analytical benchmark.
    """

    ps = esat_water(T)

    Dv = diffusivity_water_vapour(
        T,
        P,
    )

    Ka = thermal_conductivity_air(T)

    return growth_coefficient_ARG1998(
        r=representative_radius,
        T=T,
        ps=ps,
        Dv=Dv,
        Ka=Ka,
    )


def run(
    w=1.0,
    aerosol_N=100e6,
    aerosol_radius=0.05e-6,
    aerosol_kappa=0.3,
    aerosol_sigma=1.4,
    activation_scheme="simple_kappa",
    arg_soluble_mass_fraction=1.0,
    dt=0.5,
    t_end=3600.0,
    verbose=True,
):

    # -------------------------------------------------
    # Input validation
    # -------------------------------------------------

    if w < 0.0:
        raise ValueError(
            "w must be nonnegative"
        )

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
    # Backward-compatible naming
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
    # Model configuration
    # -------------------------------------------------

    config = TwoMomentConfig(
        dt=dt,
        t_end=t_end,
        aerosol_N=aerosol_N,
        aerosol_radius=aerosol_radius,
        aerosol_kappa=aerosol_kappa,
        aerosol_sigma=aerosol_sigma,
    )

    state = initialize_two_moment_state(
        config
    )

    aerosol = AerosolPopulation(
        name="warm_cloud_aerosol",
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

    T_min = 273.15

    # -------------------------------------------------
    # Initial state
    # -------------------------------------------------

    t = 0.0
    z = 0.0

    initial_water = total_water(
        state
    )

    saturation_time = None

    activation_time = None
    S_at_activation = None
    Sc_at_activation = None

    # -------------------------------------------------
    # Activation-threshold diagnostics
    # -------------------------------------------------

    activation_1pct_time = None
    activation_50pct_time = None
    activation_90pct_time = None

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

    fraction_at_SSmax = (
        state.Nc / config.aerosol_N
        if config.aerosol_N > 0.0
        else 0.0
    )

    # -------------------------------------------------
    # Standalone ARG1998 analytical diagnostics
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
    # Time loop
    # -------------------------------------------------

    while (
        t < config.t_end
        and state.T > T_min
    ):

        step_dt = min(
            dt,
            config.t_end - t,
        )

        if step_dt <= 0.0:
            break

        config.dt = step_dt

        rho_air = (
            state.p
            / (Rd * state.T)
        )

        config.rho_air = rho_air

        # ---------------------------------------------
        # Parcel ascent
        # ---------------------------------------------

        dz = w * step_dt

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

        state.S = liquid_supersaturation(
            state.qv,
            state.T,
            state.p,
        )

        S_before_activation = (
            state.S
        )

        # ---------------------------------------------
        # First saturation crossing
        # ---------------------------------------------

        if (
            saturation_time is None
            and S_before_activation >= 0.0
        ):

            saturation_time = (
                t + step_dt
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

            # -----------------------------------------
            # Provisional representative radius used
            # only for the standalone ARG benchmark.
            # -----------------------------------------

            ARG_G_at_saturation = (
                arg_growth_coefficient(
                    state.T,
                    state.p,
                    ARG_G_representative_radius,
                )
            )

            ps = esat_water(
                state.T
            )

            arg_smax_result = (
                smax_ARG_single(
                    mode=arg_mode,
                    T=state.T,
                    p=state.p,
                    ps=ps,
                    V=w,
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
                t + step_dt
            )

            fraction_at_SSmax = (
                state.Nc / config.aerosol_N
                if config.aerosol_N > 0.0
                else 0.0
            )

        # ---------------------------------------------
        # Parcel activation scheme
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

            # -----------------------------------------
            # Lognormal Kohler-threshold activation.
            #
            # This uses the ARG single-mode critical
            # supersaturation distribution, but it is
            # NOT the complete analytical ARG1998
            # parameterization.
            #
            # The current parcel supersaturation is
            # treated as the maximum threshold reached
            # so far. apply_activation_step is
            # irreversible, so activation is cumulative.
            # -----------------------------------------

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
        # Activated fraction after activation
        # ---------------------------------------------

        if config.aerosol_N > 0.0:

            current_fraction = (
                state.Nc
                / config.aerosol_N
            )

        else:

            current_fraction = 0.0

        # ---------------------------------------------
        # First nonzero activation
        # ---------------------------------------------

        if activation_time is None:

            if (
                Nc_before <= 0.0
                and state.Nc > 0.0
            ):

                activation_time = (
                    t + step_dt
                )

                S_at_activation = (
                    S_before_activation
                )

                Sc_at_activation = Sc

        # ---------------------------------------------
        # 1% activation diagnostic
        # ---------------------------------------------

        if (
            activation_1pct_time is None
            and current_fraction >= 0.01
        ):

            activation_1pct_time = (
                t + step_dt
            )

            S_at_1pct_activation = (
                S_before_activation
            )

            if saturation_time is not None:

                activation_1pct_delay = (
                    activation_1pct_time
                    - saturation_time
                )

        # ---------------------------------------------
        # 50% activation diagnostic
        # ---------------------------------------------

        if (
            activation_50pct_time is None
            and current_fraction >= 0.50
        ):

            activation_50pct_time = (
                t + step_dt
            )

            S_at_50pct_activation = (
                S_before_activation
            )

            if saturation_time is not None:

                activation_50pct_delay = (
                    activation_50pct_time
                    - saturation_time
                )

        # ---------------------------------------------
        # 90% activation diagnostic
        # ---------------------------------------------

        if (
            activation_90pct_time is None
            and current_fraction >= 0.90
        ):

            activation_90pct_time = (
                t + step_dt
            )

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
                t + step_dt
            )

            fraction_at_SSmax = (
                state.Nc / config.aerosol_N
                if config.aerosol_N > 0.0
                else 0.0
            )

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
                f"T={state.T:7.2f} K  "
                f"S={100.0 * state.S:8.4f} %  "
                f"Nc={state.Nc:12.4e}  "
                f"qc={state.qc:10.4e}  "
                f"G={config.G_liquid:10.4e}"
            )

            next_print += 100.0

    # -------------------------------------------------
    # Final activated fraction
    # -------------------------------------------------

    if config.aerosol_N > 0.0:

        activated_fraction = (
            state.Nc
            / config.aerosol_N
        )

    else:

        activated_fraction = 0.0

    # -------------------------------------------------
    # Time from saturation to first activation
    # -------------------------------------------------

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

    # -------------------------------------------------
    # Water conservation
    # -------------------------------------------------

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

        # ---------------------------------------------
        # Scheme names
        # ---------------------------------------------

        # Preserve the original input so old tests and
        # sweep files remain backward compatible.
        "activation_scheme":
            activation_scheme_input,

        # Scientifically clearer parcel-scheme name.
        "parcel_activation_scheme":
            parcel_activation_scheme,

        "ARG1998_role":
            "standalone analytical benchmark",

        "w_m_s":
            w,

        "aerosol_N_m3":
            aerosol_N,

        "aerosol_radius_m":
            aerosol_radius,

        "aerosol_kappa":
            aerosol_kappa,

        "aerosol_sigma":
            aerosol_sigma,

        "ARG_B":
            arg_B,

        "ARG_soluble_mass_fraction":
            arg_soluble_mass_fraction,

        "dt_s":
            dt,

        "requested_t_end_s":
            t_end,

        # ---------------------------------------------
        # Saturation / first activation
        # ---------------------------------------------

        "saturation_time_s":
            saturation_time,

        "activation_time_s":
            activation_time,

        "time_to_activation_s":
            time_from_saturation_to_activation,

        "time_from_saturation_to_activation_s":
            time_from_saturation_to_activation,

        "S_at_activation":
            S_at_activation,

        "Sc_at_activation":
            Sc_at_activation,

        # ---------------------------------------------
        # Threshold diagnostics
        # ---------------------------------------------

        "activation_1pct_time_s":
            activation_1pct_time,

        "activation_1pct_delay_s":
            activation_1pct_delay,

        "S_at_1pct_activation":
            S_at_1pct_activation,

        "activation_50pct_time_s":
            activation_50pct_time,

        "activation_50pct_delay_s":
            activation_50pct_delay,

        "S_at_50pct_activation":
            S_at_50pct_activation,

        "activation_90pct_time_s":
            activation_90pct_time,

        "activation_90pct_delay_s":
            activation_90pct_delay,

        "S_at_90pct_activation":
            S_at_90pct_activation,

        # ---------------------------------------------
        # Main parcel diagnostics
        # ---------------------------------------------

        "activated_fraction":
            activated_fraction,

        "SSmax":
            SSmax,

        "SSmax_percent":
            100.0 * SSmax,

        "time_of_SSmax_s":
            time_of_SSmax,

        "fraction_at_SSmax":
            fraction_at_SSmax,

        # ---------------------------------------------
        # Standalone ARG1998 analytical benchmark
        # ---------------------------------------------

        "ARG_Sm_at_saturation":
            ARG_Sm_at_saturation,

        "ARG_Sm_at_saturation_percent":
            (
                100.0
                * ARG_Sm_at_saturation
                if ARG_Sm_at_saturation
                is not None
                else None
            ),

        "ARG_Smax_prediction":
            ARG_Smax_prediction,

        "ARG_Smax_prediction_percent":
            (
                100.0
                * ARG_Smax_prediction
                if ARG_Smax_prediction
                is not None
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

        # ---------------------------------------------
        # Final state
        # ---------------------------------------------

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

        "height_gain_m":
            z,

        "simulation_time_s":
            t,

        "water_error":
            water_error,

        "stop_reason":
            (
                "warm-cloud temperature limit reached"
                if state.T <= T_min
                else
                "requested simulation time reached"
            ),
    }

    # -------------------------------------------------
    # Print diagnostics
    # -------------------------------------------------

    if verbose:

        print()

        print(
            "Two-moment warm-cloud diagnostics"
        )

        print(
            "---------------------------------"
        )

        for key, value in result.items():

            print(
                f"{key} = {value}"
            )

    return result


if __name__ == "__main__":
    run()

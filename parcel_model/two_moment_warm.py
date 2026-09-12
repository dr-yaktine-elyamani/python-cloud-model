"""Experimental warm-cloud two-moment model.

This module is under development. It is not yet a validated
two-moment microphysics scheme.
"""

import math
from dataclasses import dataclass


@dataclass
class TwoMomentConfig:
    # Numerical settings
    dt: float = 0.5
    t_end: float = 3600.0

    # Initial thermodynamic conditions
    T0: float = 283.15
    p0: float = 90000.0
    qv0: float = 0.00817

    # Aerosol population
    aerosol_N: float = 100e6
    aerosol_radius: float = 0.05e-6
    aerosol_kappa: float = 0.3
    aerosol_sigma: float = 1.4
    aerosol_density: float = 1770.0

    # Physical parameters
    rho_water: float = 1000.0
    rho_air: float = 1.0
    G_liquid: float = 8.0e-12

    # Initial cloud droplet size
    initial_droplet_radius: float = 1.0e-6


@dataclass
class TwoMomentState:
    # Thermodynamic state
    T: float
    p: float
    qv: float

    # Prognostic cloud moments
    Nc: float = 0.0
    qc: float = 0.0

    # Diagnostic quantities
    S: float = 0.0
    Sc: float = 0.0
    mean_radius: float = 0.0

    # Activation diagnostics
    SSmax: float = 0.0
    activation_time: float = -1.0
    activated_fraction: float = 0.0


def mean_radius_from_moments(
    Nc: float,
    qc: float,
    rho_air: float,
    rho_water: float,
) -> float:
    """
    Diagnose the equivalent-volume droplet radius.

    Nc: droplet number concentration [m^-3]
    qc: cloud-water mixing ratio [kg/kg dry air]
    rho_air: dry-air density [kg/m^3]
    rho_water: liquid-water density [kg/m^3]

    Assumes a monodisperse population.
    """
    if Nc <= 0.0 or qc <= 0.0:
        return 0.0

    return (
        3.0 * rho_air * qc
        / (4.0 * math.pi * rho_water * Nc)
    ) ** (1.0 / 3.0)
def vapour_pressure_from_mixing_ratio(qv: float, p: float) -> float:
    """Convert water-vapour mixing ratio [kg/kg] to pressure [Pa]."""
    epsilon = 0.622
    return qv * p / (epsilon + qv)


def liquid_supersaturation(qv: float, T: float, p: float) -> float:
    """Calculate supersaturation over liquid water."""
    from parcel_model.thermodynamics import Sw

    e = vapour_pressure_from_mixing_ratio(qv, p)
    return float(Sw(e, T))
def activation_number_tendency(
    Nc: float,
    aerosol_N: float,
    activated_fraction: float,
    dt: float,
) -> float:
    """
    Calculate the activation source for cloud droplet number.

    Nc and aerosol_N are number concentrations [m^-3].
    Returns dNc/dt [m^-3 s^-1].

    This is a provisional activation closure, not a complete
    prognostic aerosol or two-moment microphysics scheme.
    """
    if dt <= 0.0:
        raise ValueError("dt must be positive")

    target_Nc = aerosol_N * min(max(activated_fraction, 0.0), 1.0)

    # Only newly activated droplets contribute to the source.
    new_activation = max(target_Nc - Nc, 0.0)

    return new_activation / dt
def activation_mass_tendency(
    activation_number_rate: float,
    initial_radius: float,
    rho_water: float,
    rho_air: float,
) -> float:
    """
    Cloud-water mass source associated with newly activated droplets.

    activation_number_rate: [m^-3 s^-1]
    initial_radius: [m]
    rho_water: [kg/m^3]
    rho_air: [kg/m^3]

    Returns cloud-water mixing-ratio tendency [kg/kg/s].

    Assumes newly activated droplets enter the cloud category
    with the prescribed initial radius.
    """
    if activation_number_rate <= 0.0:
        return 0.0

    if initial_radius <= 0.0 or rho_water <= 0.0 or rho_air <= 0.0:
        raise ValueError("Radii and densities must be positive")

    mass_per_droplet = (
        (4.0 / 3.0)
        * math.pi
        * rho_water
        * initial_radius**3
    )

    return activation_number_rate * mass_per_droplet / rho_air
def condensation_tendency(
    Nc: float,
    qc: float,
    S: float,
    G_liquid: float,
    rho_air: float,
    rho_water: float,
) -> float:
    """
    Maxwell-type bulk condensation/evaporation tendency.

    Returns dqc/dt [kg/kg/s].
    Positive: condensation.
    Negative: evaporation.

    Monodisperse closure: radius is diagnosed from Nc and qc.
    """
    if Nc <= 0.0 or qc <= 0.0:
        return 0.0

    r = mean_radius_from_moments(
        Nc, qc, rho_air, rho_water
    )

    if r <= 0.0:
        return 0.0

    dr_dt = G_liquid * S / r

    return (
        Nc
        * 4.0
        * math.pi
        * rho_water
        * r**2
        * dr_dt
        / rho_air
    )
def apply_condensation_step(
    state: TwoMomentState,
    config: TwoMomentConfig,
) -> TwoMomentState:
    """
    Apply one conservative condensation/evaporation step.

    Transfers water between vapour and cloud liquid.
    This is a provisional fixed-temperature, fixed-pressure step.
    """
    if config.dt <= 0.0:
        raise ValueError("dt must be positive")

    rate = condensation_tendency(
        Nc=state.Nc,
        qc=state.qc,
        S=state.S,
        G_liquid=config.G_liquid,
        rho_air=config.rho_air,
        rho_water=config.rho_water,
    )

    # Requested liquid-water change over the timestep.
    delta_qc = rate * config.dt

    # Condensation cannot consume more vapour than available.
    # Evaporation cannot remove more cloud water than available.
    delta_qc = min(delta_qc, state.qv)
    delta_qc = max(delta_qc, -state.qc)

    state.qv -= delta_qc
    state.qc += delta_qc

    # Recalculate diagnostics from the updated prognostic variables.
    state.S = liquid_supersaturation(state.qv, state.T, state.p)
    state.mean_radius = mean_radius_from_moments(
        state.Nc,
        state.qc,
        config.rho_air,
        config.rho_water,
    )

        # Update activation diagnostics.
    state.SSmax = max(state.SSmax, state.S)

    if config.aerosol_N > 0.0:
        state.activated_fraction = min(
            max(state.Nc / config.aerosol_N, 0.0),
            1.0,
        )
    else:
        state.activated_fraction = 0.0
    
        # Update activation diagnostics.
    state.SSmax = max(state.SSmax, state.S)

    if config.aerosol_N > 0.0:
        state.activated_fraction = min(
            max(state.Nc / config.aerosol_N, 0.0),
            1.0,
        )
    else:
        state.activated_fraction = 0.0

    return state

    # Update activation diagnostics.
    state.SSmax = max(state.SSmax, state.S)

    if config.aerosol_N > 0.0:
        state.activated_fraction = min(
            max(state.Nc / config.aerosol_N, 0.0),
            1.0,
        )
    else:
        state.activated_fraction = 0.0
def initialize_two_moment_state(
    config: TwoMomentConfig,
) -> TwoMomentState:
    """Initialize the warm-cloud two-moment state."""

    state = TwoMomentState(
        T=config.T0,
        p=config.p0,
        qv=config.qv0,
        Nc=0.0,
        qc=0.0,
    )

    state.S = liquid_supersaturation(
        state.qv, state.T, state.p
    )

    state.mean_radius = mean_radius_from_moments(
        state.Nc,
        state.qc,
        config.rho_air,
        config.rho_water,
    )

    return state
def total_water(state: TwoMomentState) -> float:
    """Total vapour plus cloud liquid mixing ratio [kg/kg]."""
    return state.qv + state.qc


def validate_state(
    state: TwoMomentState,
    config: TwoMomentConfig,
) -> None:
    """Check basic physical and numerical invariants."""
    values = (
        state.T, state.p, state.qv, state.Nc, state.qc,
        state.S, state.mean_radius,
    )

    if not all(math.isfinite(value) for value in values):
        raise ValueError("Non-finite value in two-moment state")

    if state.T <= 0.0 or state.p <= 0.0:
        raise ValueError("Temperature and pressure must be positive")

    if state.qv < 0.0 or state.qc < 0.0 or state.Nc < 0.0:
        raise ValueError("Negative water or droplet number")

    if state.Nc > config.aerosol_N * (1.0 + 1.0e-12):
        raise ValueError("Cloud droplet number exceeds aerosol population")

    if state.Nc == 0.0 and state.qc > 0.0:
        raise ValueError("Cloud water exists without cloud droplets")
def apply_activation_step(
    state: TwoMomentState,
    config: TwoMomentConfig,
    activated_fraction: float,
) -> TwoMomentState:
    """
    Apply a provisional activation source conservatively.

    Newly activated droplets enter with a prescribed wet radius.
    The associated liquid mass is transferred from water vapour.
    This is not yet a complete aerosol-resolved activation scheme.
    """

    if config.dt <= 0.0:
        raise ValueError("dt must be positive")

    # Calculate number activation tendency.
    number_rate = activation_number_tendency(
        Nc=state.Nc,
        aerosol_N=config.aerosol_N,
        activated_fraction=activated_fraction,
        dt=config.dt,
    )

    # Calculate the liquid-water mass associated
    # with the newly activated droplets.
    mass_rate = activation_mass_tendency(
        activation_number_rate=number_rate,
        initial_radius=config.initial_droplet_radius,
        rho_water=config.rho_water,
        rho_air=config.rho_air,
    )

    # Convert tendencies to changes over the timestep.
    delta_Nc = number_rate * config.dt
    delta_qc = mass_rate * config.dt

    # Activation cannot consume more vapour than is available.
    if delta_qc > state.qv and delta_qc > 0.0:
        factor = state.qv / delta_qc
        delta_Nc *= factor
        delta_qc = state.qv

    # Update prognostic variables.
    state.Nc += delta_Nc
    state.qc += delta_qc
    state.qv -= delta_qc

    # Recalculate supersaturation.
    state.S = liquid_supersaturation(
        state.qv,
        state.T,
        state.p,
    )

    # Recalculate mean droplet radius.
    state.mean_radius = mean_radius_from_moments(
        state.Nc,
        state.qc,
        config.rho_air,
        config.rho_water,
    )

    # Maximum supersaturation reached so far.
    state.SSmax = max(
        state.SSmax,
        state.S,
    )

    # Fraction of the aerosol population activated.
    if config.aerosol_N > 0.0:
        state.activated_fraction = min(
            max(
                state.Nc / config.aerosol_N,
                0.0,
            ),
            1.0,
        )
    else:
        state.activated_fraction = 0.0

    return state


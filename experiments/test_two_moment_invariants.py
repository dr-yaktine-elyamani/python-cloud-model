def apply_microphysics_step(
    state: TwoMomentState,
    config: TwoMomentConfig,
    activated_fraction: float,
) -> TwoMomentState:
    """
    Apply one provisional warm-cloud microphysics step.

    1. Activate new droplets and transfer their initial mass.
    2. Apply condensation or evaporation.
    3. Check physical invariants and water conservation.

    This is a development prototype, not a validated full scheme.
    """
    water_before = total_water(state)

    apply_activation_step(state, config, activated_fraction)
    apply_condensation_step(state, config)

    validate_state(state, config)

    water_after = total_water(state)
    error = water_after - water_before

    if abs(error) > 1.0e-12:
        raise ValueError(
            f"Water conservation failed: error={error}"
        )

    return state

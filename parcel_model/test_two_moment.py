import unittest

from parcel_model.two_moment_warm import (
    TwoMomentConfig,
    initialize_two_moment_state,
    apply_activation_step,
    apply_condensation_step,
    total_water,
)


class TestTwoMomentActivation(unittest.TestCase):

    def test_activation_increases_droplet_number(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        Nc_before = state.Nc

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        self.assertGreater(state.Nc, Nc_before)

    def test_activation_conserves_water(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        water_before = total_water(state)

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        water_after = total_water(state)

        self.assertAlmostEqual(
            water_before,
            water_after,
            places=14,
        )

    def test_droplet_number_does_not_exceed_aerosol(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        apply_activation_step(
            state,
            config,
            activated_fraction=1.0,
        )

        self.assertLessEqual(
            state.Nc,
            config.aerosol_N,
        )

    def test_condensation_conserves_water(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        water_before = total_water(state)

        apply_condensation_step(
            state,
            config,
        )

        water_after = total_water(state)

        self.assertAlmostEqual(
            water_before,
            water_after,
            places=14,
        )

    def test_activated_fraction_diagnostic(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        self.assertAlmostEqual(
            state.activated_fraction,
            0.5,
            places=12,
        )

    def test_ssmax_diagnostic(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        state.S = 0.01
        state.SSmax = 0.0

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        self.assertGreaterEqual(
            state.SSmax,
            state.S,
        )


if __name__ == "__main__":
    unittest.main()

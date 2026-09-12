import unittest

from parcel_model.two_moment_warm import (
    TwoMomentConfig,
    initialize_two_moment_state,
    apply_activation_step,
    apply_condensation_step,
    total_water,
)


class TestTwoMomentDiagnostics(unittest.TestCase):

    def test_initial_state_has_no_cloud_water(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        self.assertEqual(state.Nc, 0.0)
        self.assertEqual(state.qc, 0.0)

    def test_activation_fraction_is_consistent(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        expected = state.Nc / config.aerosol_N

        self.assertAlmostEqual(
            state.activated_fraction,
            expected,
            places=12,
        )

    def test_activation_conserves_total_water(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        water_before = total_water(state)

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        self.assertAlmostEqual(
            total_water(state),
            water_before,
            places=14,
        )

    def test_condensation_conserves_total_water(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        water_before = total_water(state)

        apply_condensation_step(state, config)

        self.assertAlmostEqual(
            total_water(state),
            water_before,
            places=14,
        )

    def test_droplet_number_is_bounded(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        apply_activation_step(
            state,
            config,
            activated_fraction=1.0,
        )

        self.assertGreaterEqual(state.Nc, 0.0)
        self.assertLessEqual(state.Nc, config.aerosol_N)

    def test_ssmax_retains_previous_peak(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        state.SSmax = 0.02

        apply_activation_step(
            state,
            config,
            activated_fraction=0.5,
        )

        self.assertGreaterEqual(state.SSmax, 0.02)

    def test_no_activation_when_target_is_zero(self):
        config = TwoMomentConfig()
        state = initialize_two_moment_state(config)

        apply_activation_step(
            state,
            config,
            activated_fraction=0.0,
        )

        self.assertEqual(state.Nc, 0.0)
        self.assertEqual(state.qc, 0.0)

    def test_activation_does_not_create_water(self):
        config = TwoMomentConfig(qv0=0.0)
        state = initialize_two_moment_state(config)

        apply_activation_step(
            state,
            config,
            activated_fraction=1.0,
        )

        self.assertEqual(state.qv, 0.0)
        self.assertEqual(state.qc, 0.0)
        self.assertEqual(state.Nc, 0.0)


if __name__ == "__main__":
    unittest.main()

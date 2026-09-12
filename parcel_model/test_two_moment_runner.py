import unittest

from parcel_model.run_two_moment_warm import (
    run,
    maxwell_G_liquid,
)


class TestTwoMomentWarmRunner(unittest.TestCase):

    def test_runner_returns_dict(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertIsInstance(
            result,
            dict,
        )


    def test_activation_occurs(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertIsNotNone(
            result["activation_time_s"]
        )

        self.assertGreater(
            result["Nc_final"],
            0.0,
        )


    def test_activation_threshold_consistency(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertIsNotNone(
            result["S_at_activation"]
        )

        self.assertIsNotNone(
            result["Sc_at_activation"]
        )

        self.assertGreaterEqual(
            result["S_at_activation"],
            result["Sc_at_activation"],
        )


    def test_SSmax_not_below_activation_S(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertGreaterEqual(
            result["SSmax"],
            result["S_at_activation"],
        )


    def test_SSmax_time_within_simulation(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertGreaterEqual(
            result["time_of_SSmax_s"],
            0.0,
        )

        self.assertLessEqual(
            result["time_of_SSmax_s"],
            result["simulation_time_s"],
        )


    def test_activated_fraction_bounded(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertGreaterEqual(
            result["activated_fraction"],
            0.0,
        )

        self.assertLessEqual(
            result["activated_fraction"],
            1.0,
        )


    def test_fraction_at_SSmax_bounded(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertGreaterEqual(
            result["fraction_at_SSmax"],
            0.0,
        )

        self.assertLessEqual(
            result["fraction_at_SSmax"],
            1.0,
        )


    def test_water_species_nonnegative(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertGreaterEqual(
            result["qv_final"],
            0.0,
        )

        self.assertGreaterEqual(
            result["qc_final"],
            0.0,
        )


    def test_final_temperature_is_near_warm_limit(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertGreater(
            result["T_final_K"],
            272.0,
        )

        self.assertLessEqual(
            result["T_final_K"],
            273.15,
        )


    def test_maxwell_G_positive(self):

        G = maxwell_G_liquid(
            T=283.15,
            P=90000.0,
        )

        self.assertGreater(
            G,
            0.0,
        )


    def test_maxwell_G_changes_with_temperature(self):

        G1 = maxwell_G_liquid(
            T=283.15,
            P=90000.0,
        )

        G2 = maxwell_G_liquid(
            T=273.15,
            P=90000.0,
        )

        self.assertNotEqual(
            G1,
            G2,
        )


    def test_saturation_occurs_before_activation(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertIsNotNone(
            result["saturation_time_s"]
        )

        self.assertIsNotNone(
            result["activation_time_s"]
        )

        self.assertLessEqual(
            result["saturation_time_s"],
            result["activation_time_s"],
        )


    def test_time_to_activation_is_nonnegative(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertIsNotNone(
            result["time_to_activation_s"]
        )

        self.assertGreaterEqual(
            result["time_to_activation_s"],
            0.0,
        )


    def test_time_to_activation_definition(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        expected = (
            result["activation_time_s"]
            - result["saturation_time_s"]
        )

        self.assertAlmostEqual(
            result["time_to_activation_s"],
            expected,
            places=12,
        )


    def test_time_to_activation_alias_consistency(self):

        result = run(
            w=1.0,
            dt=0.5,
            verbose=False,
        )

        self.assertAlmostEqual(
            result["time_to_activation_s"],
            result[
                "time_from_saturation_to_activation_s"
            ],
            places=12,
        )


if __name__ == "__main__":

    unittest.main()

import unittest

from parcel_model.run_two_moment_warm import run


class TestActivationThresholdDiagnostics(
    unittest.TestCase
):

    def test_arg_threshold_times_are_ordered(self):

        result = run(
            w=1.0,
            aerosol_N=100e6,
            aerosol_radius=0.05e-6,
            aerosol_sigma=1.4,
            activation_scheme="ARG1998",
            dt=0.25,
            verbose=False,
        )

        t1 = result["activation_1pct_time_s"]
        t50 = result["activation_50pct_time_s"]
        t90 = result["activation_90pct_time_s"]

        self.assertIsNotNone(t1)
        self.assertIsNotNone(t50)
        self.assertIsNotNone(t90)

        self.assertLessEqual(t1, t50)
        self.assertLessEqual(t50, t90)

    def test_simple_threshold_times_are_ordered(self):

        result = run(
            w=1.0,
            aerosol_N=100e6,
            aerosol_radius=0.05e-6,
            activation_scheme="simple_kappa",
            dt=0.25,
            verbose=False,
        )

        t1 = result["activation_1pct_time_s"]
        t50 = result["activation_50pct_time_s"]
        t90 = result["activation_90pct_time_s"]

        self.assertIsNotNone(t1)
        self.assertIsNotNone(t50)
        self.assertIsNotNone(t90)

        self.assertLessEqual(t1, t50)
        self.assertLessEqual(t50, t90)

    def test_arg_threshold_delays_are_nonnegative(self):

        result = run(
            w=1.0,
            aerosol_N=100e6,
            aerosol_radius=0.05e-6,
            aerosol_sigma=1.4,
            activation_scheme="ARG1998",
            dt=0.25,
            verbose=False,
        )

        for key in [
            "activation_1pct_delay_s",
            "activation_50pct_delay_s",
            "activation_90pct_delay_s",
        ]:
            self.assertIsNotNone(result[key])
            self.assertGreaterEqual(
                result[key],
                0.0,
            )

    def test_unreached_threshold_is_none(self):

        result = run(
            w=0.2,
            aerosol_N=200e6,
            aerosol_radius=0.05e-6,
            activation_scheme="simple_kappa",
            dt=0.25,
            verbose=False,
        )

        self.assertLess(
            result["activated_fraction"],
            0.50,
        )

        self.assertIsNone(
            result["activation_50pct_time_s"]
        )

        self.assertIsNone(
            result["activation_90pct_time_s"]
        )


if __name__ == "__main__":
    unittest.main()

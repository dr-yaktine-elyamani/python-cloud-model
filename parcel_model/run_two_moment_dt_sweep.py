"""Time-step convergence study for the experimental warm-cloud model."""

import csv
import os

from parcel_model.run_two_moment_warm import run


def main():
    os.makedirs("data", exist_ok=True)

    dt_values = [2.0, 1.0, 0.5, 0.25, 0.1]

    results = []

    for dt in dt_values:
        print(f"Running dt = {dt} s")

        result = run(
            w=1.0,
            aerosol_N=100e6,
            dt=dt,
            verbose=False,
        )

        results.append(result)

        print(
            f"  activation={result['activation_time_s']} s, "
            f"SSmax={result['SSmax_percent']:.6f} %, "
            f"fraction={result['activated_fraction']:.6f}, "
            f"qc={result['qc_final']:.6e}"
        )

    outfile = "data/two_moment_dt_sweep.csv"

    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(results[0].keys()),
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved {len(results)} simulations to {outfile}")


if __name__ == "__main__":
    main()

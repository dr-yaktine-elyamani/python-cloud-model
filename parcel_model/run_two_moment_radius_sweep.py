"""Dry-aerosol-radius sensitivity sweep.

Experimental warm-cloud two-moment parcel model.
"""

import csv
import os

from parcel_model.run_two_moment_warm import run


def main():

    os.makedirs("data", exist_ok=True)

    # Dry aerosol radii [m]
    radius_values = [
        0.02e-6,
        0.03e-6,
        0.05e-6,
        0.07e-6,
        0.10e-6,
    ]

    # Keep these fixed for this experiment.
    w = 1.0
    aerosol_N = 100e6
    aerosol_kappa = 0.3
    dt = 0.25

    results = []

    for radius in radius_values:

        print(
            f"Running radius="
            f"{radius * 1.0e6:.3f} um"
        )

        result = run(
            w=w,
            aerosol_N=aerosol_N,
            aerosol_radius=radius,
            aerosol_kappa=aerosol_kappa,
            dt=dt,
            verbose=False,
        )

        results.append(result)

        Sc_percent = (
            100.0 * result["Sc_at_activation"]
            if result["Sc_at_activation"] is not None
            else None
        )

        print(
            f"  saturation="
            f"{result['saturation_time_s']} s, "
            f"activation="
            f"{result['activation_time_s']} s, "
            f"time_to_activation="
            f"{result['time_to_activation_s']} s, "
            f"Sc="
            f"{Sc_percent:.4f} %, "
            f"SSmax="
            f"{result['SSmax_percent']:.4f} %, "
            f"fraction="
            f"{result['activated_fraction']:.4f}"
        )

    outfile = (
        "data/two_moment_radius_sweep_dt025.csv"
    )

    with open(
        outfile,
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(
                results[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(results)

    print()

    print(
        f"Saved {len(results)} simulations "
        f"to {outfile}"
    )


if __name__ == "__main__":
    main()

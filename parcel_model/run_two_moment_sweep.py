"""Experimental updraft and aerosol-number sensitivity sweep."""

import csv
import os

from parcel_model.run_two_moment_warm import run


def main():

    os.makedirs("data", exist_ok=True)

    w_values = [
        0.2,
        0.5,
        1.0,
        2.0,
        3.0,
    ]

    aerosol_N_values = [
        10e6,
        50e6,
        100e6,
        200e6,
    ]

    dt = 0.25

    results = []

    for w in w_values:

        for aerosol_N in aerosol_N_values:

            print(
                f"Running w={w} m/s, "
                f"N={aerosol_N:.2e} m^-3, "
                f"dt={dt} s"
            )

            result = run(
                w=w,
                aerosol_N=aerosol_N,
                dt=dt,
                verbose=False,
            )

            results.append(result)

            print(
                f"  activation="
                f"{result['activation_time_s']} s, "
                f"SSmax="
                f"{result['SSmax_percent']:.4f} %, "
                f"fraction="
                f"{result['activated_fraction']:.4f}, "
                f"r="
                f"{result['mean_radius_um']:.2f} um"
            )

    outfile = (
        "data/two_moment_w_N_sweep_dt025.csv"
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

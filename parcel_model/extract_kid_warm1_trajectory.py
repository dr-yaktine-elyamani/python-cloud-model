"""Extract KiD warm1 fields along the approximate parcel trajectory.

The parcel trajectory follows the KiD warm1 forcing:

    w(t) = 2 sin(pi t / 600),  t < 600 s
    w(t) = 0,                  t >= 600 s

Starting height:
    z0 = 25 m

For each KiD output time, this script:
1. computes the approximate parcel height,
2. finds the nearest KiD vertical grid level,
3. extracts thermodynamic and cloud variables,
4. writes a CSV for comparison with the Python parcel model.

This is a diagnostic comparison, not validation.
"""

import csv
import math

from netCDF4 import Dataset


KID_FILE = (
    "/home/dryaktine/KiD-A/"
    "output/warm1_output.nc"
)

OUTPUT_CSV = (
    "data/kid_warm1_parcel_trajectory.csv"
)

Z0 = 25.0
WMAX = 2.0
TFORCING = 600.0


def parcel_height(t):
    """Analytical height from KiD warm1 sinusoidal forcing."""

    if t <= TFORCING:

        return (
            Z0
            + (
                WMAX
                * TFORCING
                / math.pi
            )
            * (
                1.0
                - math.cos(
                    math.pi
                    * t
                    / TFORCING
                )
            )
        )

    return (
        Z0
        + (
            2.0
            * WMAX
            * TFORCING
            / math.pi
        )
    )


def nearest_index(values, target):
    """Return index of value nearest to target."""

    return min(
        range(len(values)),
        key=lambda i: abs(
            float(values[i])
            - target
        ),
    )


with Dataset(KID_FILE) as ds:

    time = ds.variables["time"][:]
    zgrid = ds.variables["z"][:]

    pressure = ds.variables["pressure"]
    temperature = ds.variables["temperature"]
    vapour = ds.variables["vapour"]
    RH = ds.variables["RH"]
    w = ds.variables["w"]
    cloud_mass = ds.variables["cloud_mass"]
    rain_mass = ds.variables["rain_mass"]

    rows = []

    for it in range(len(time)):

        t = float(
            time[it]
        )

        z_parcel = (
            parcel_height(t)
        )

        iz = nearest_index(
            zgrid,
            z_parcel,
        )

        z_kid = float(
            zgrid[iz]
        )

        row = {
            "time_s":
                t,

            "parcel_height_m":
                z_parcel,

            "kid_nearest_height_m":
                z_kid,

            "height_difference_m":
                z_kid
                - z_parcel,

            "w_m_s":
                float(
                    w[iz, it]
                ),

            "temperature_K":
                float(
                    temperature[iz, it]
                ),

            "pressure_mb":
                float(
                    pressure[iz, it]
                ),

            "qv_kg_kg":
                float(
                    vapour[iz, it]
                ),

            "RH_percent":
                float(
                    RH[iz, it]
                ),

            "cloud_mass_kg_kg":
                float(
                    cloud_mass[iz, it]
                ),

            "rain_mass_kg_kg":
                float(
                    rain_mass[iz, it]
                ),
        }

        rows.append(
            row
        )


fieldnames = list(
    rows[0].keys()
)

with open(
    OUTPUT_CSV,
    "w",
    newline="",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    writer.writerows(
        rows
    )


print()
print(
    "KiD warm1 trajectory extraction complete"
)
print(
    "-----------------------------------------"
)
print(
    f"Saved: {OUTPUT_CSV}"
)
print(
    f"Rows: {len(rows)}"
)

print()
print(
    "Selected times"
)
print(
    "--------------"
)

for target_time in (
    300.0,
    390.0,
    420.0,
    480.0,
    570.0,
    600.0,
):

    nearest_row = min(
        rows,
        key=lambda r: abs(
            r["time_s"]
            - target_time
        ),
    )

    print()
    print(
        f"time = "
        f"{nearest_row['time_s']:.1f} s"
    )

    print(
        f"parcel z = "
        f"{nearest_row['parcel_height_m']:.2f} m"
    )

    print(
        f"KiD z = "
        f"{nearest_row['kid_nearest_height_m']:.2f} m"
    )

    print(
        f"w = "
        f"{nearest_row['w_m_s']:.4f} m/s"
    )

    print(
        f"T = "
        f"{nearest_row['temperature_K']:.3f} K"
    )

    print(
        f"qv = "
        f"{nearest_row['qv_kg_kg']:.6e}"
    )

    print(
        f"RH = "
        f"{nearest_row['RH_percent']:.4f} %"
    )

    print(
        f"cloud = "
        f"{nearest_row['cloud_mass_kg_kg']:.6e}"
    )

    print(
        f"rain = "
        f"{nearest_row['rain_mass_kg_kg']:.6e}"
    )

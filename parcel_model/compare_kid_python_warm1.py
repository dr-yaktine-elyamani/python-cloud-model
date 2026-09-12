"""Compare KiD warm1 with the experimental Python parcel trajectory."""

import csv
import os

import matplotlib.pyplot as plt


KID_FILE = "data/kid_warm1_parcel_trajectory.csv"

PYTHON_FILE = (
    "data/python_two_moment_kid_warm1_trajectory.csv"
)

OUTPUT_CSV = (
    "data/kid_python_warm1_comparison.csv"
)


def read_csv(filename):
    with open(filename, newline="") as f:
        return list(csv.DictReader(f))


kid = read_csv(KID_FILE)
python = read_csv(PYTHON_FILE)


kid_by_time = {
    round(float(row["time_s"]), 6): row
    for row in kid
}

python_by_time = {
    round(float(row["time_s"]), 6): row
    for row in python
}


common_times = sorted(
    set(kid_by_time)
    & set(python_by_time)
)


rows = []


for t in common_times:

    k = kid_by_time[t]
    p = python_by_time[t]

    kid_T = float(
        k["temperature_K"]
    )

    py_T = float(
        p["temperature_K"]
    )

    kid_qv = float(
        k["qv_kg_kg"]
    )

    py_qv = float(
        p["qv_kg_kg"]
    )

    kid_cloud = float(
        k["cloud_mass_kg_kg"]
    )

    kid_rain = float(
        k["rain_mass_kg_kg"]
    )

    py_qc = float(
        p["qc_kg_kg"]
    )

    kid_total_condensed = (
        kid_cloud
        + kid_rain
    )

    kid_total_water = (
        kid_qv
        + kid_cloud
        + kid_rain
    )

    py_total_water = (
        py_qv
        + py_qc
    )

    if kid_cloud > 0.0:

        qc_ratio = (
            py_qc
            / kid_cloud
        )

    else:

        qc_ratio = None

    rows.append(
        {
            "time_s":
                t,

            "kid_height_m":
                float(
                    k["kid_nearest_height_m"]
                ),

            "python_height_m":
                float(
                    p["parcel_height_m"]
                ),

            "kid_T_K":
                kid_T,

            "python_T_K":
                py_T,

            "delta_T_python_minus_kid_K":
                py_T - kid_T,

            "kid_qv":
                kid_qv,

            "python_qv":
                py_qv,

            "delta_qv_python_minus_kid":
                py_qv - kid_qv,

            "kid_cloud":
                kid_cloud,

            "kid_rain":
                kid_rain,

            "kid_total_condensed":
                kid_total_condensed,

            "python_qc":
                py_qc,

            "python_qc_over_kid_cloud":
                qc_ratio,

            "kid_total_water":
                kid_total_water,

            "python_total_water":
                py_total_water,

            "total_water_difference":
                py_total_water
                - kid_total_water,

            "python_S_percent":
                float(
                    p["S_percent"]
                ),

            "kid_RH_percent":
                float(
                    k["RH_percent"]
                ),
        }
    )


# ---------------------------------------------------------
# Save combined CSV
# ---------------------------------------------------------

os.makedirs(
    "data",
    exist_ok=True,
)

with open(
    OUTPUT_CSV,
    "w",
    newline="",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys(),
    )

    writer.writeheader()

    writer.writerows(
        rows
    )


# ---------------------------------------------------------
# Extract plotting arrays
# ---------------------------------------------------------

time = [
    row["time_s"]
    for row in rows
]

kid_T = [
    row["kid_T_K"]
    for row in rows
]

py_T = [
    row["python_T_K"]
    for row in rows
]

kid_qv = [
    row["kid_qv"]
    for row in rows
]

py_qv = [
    row["python_qv"]
    for row in rows
]

kid_cloud = [
    row["kid_cloud"]
    for row in rows
]

kid_condensed = [
    row["kid_total_condensed"]
    for row in rows
]

py_qc = [
    row["python_qc"]
    for row in rows
]


# ---------------------------------------------------------
# Temperature plot
# ---------------------------------------------------------

plt.figure()

plt.plot(
    time,
    kid_T,
    label="KiD trajectory",
)

plt.plot(
    time,
    py_T,
    label="Python parcel",
)

plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    "Temperature (K)"
)

plt.title(
    "KiD warm1: temperature comparison"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    "data/kid_python_warm1_temperature.png",
    dpi=200,
)

plt.close()


# ---------------------------------------------------------
# Vapour plot
# ---------------------------------------------------------

plt.figure()

plt.plot(
    time,
    kid_qv,
    label="KiD trajectory",
)

plt.plot(
    time,
    py_qv,
    label="Python parcel",
)

plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    "Water-vapour mixing ratio (kg/kg)"
)

plt.title(
    "KiD warm1: water-vapour comparison"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    "data/kid_python_warm1_qv.png",
    dpi=200,
)

plt.close()


# ---------------------------------------------------------
# Condensed-water plot
# ---------------------------------------------------------

plt.figure()

plt.plot(
    time,
    kid_cloud,
    label="KiD cloud mass",
)

plt.plot(
    time,
    kid_condensed,
    label="KiD cloud + rain",
)

plt.plot(
    time,
    py_qc,
    label="Python cloud water",
)

plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    "Mixing ratio (kg/kg)"
)

plt.title(
    "KiD warm1: condensed-water comparison"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    "data/kid_python_warm1_condensed_water.png",
    dpi=200,
)

plt.close()


# ---------------------------------------------------------
# Print selected comparison points
# ---------------------------------------------------------

print()

print(
    "KiD versus Python warm1 comparison"
)

print(
    "----------------------------------"
)

print(
    f"Saved: {OUTPUT_CSV}"
)

print()

for target in (
    300.0,
    390.0,
    420.0,
    480.0,
    570.0,
    600.0,
):

    row = min(
        rows,
        key=lambda r: abs(
            r["time_s"]
            - target
        ),
    )

    print(
        f"time = {row['time_s']:.0f} s"
    )

    print(
        "  delta T = "
        f"{row['delta_T_python_minus_kid_K']:.6f} K"
    )

    print(
        "  delta qv = "
        f"{row['delta_qv_python_minus_kid']:.6e}"
    )

    print(
        "  KiD cloud = "
        f"{row['kid_cloud']:.6e}"
    )

    print(
        "  KiD rain = "
        f"{row['kid_rain']:.6e}"
    )

    print(
        "  Python qc = "
        f"{row['python_qc']:.6e}"
    )

    if (
        row["python_qc_over_kid_cloud"]
        is not None
    ):

        print(
            "  Python qc / KiD cloud = "
            f"{100.0 * row['python_qc_over_kid_cloud']:.2f} %"
        )

    print(
        "  KiD total water = "
        f"{row['kid_total_water']:.6e}"
    )

    print(
        "  Python total water = "
        f"{row['python_total_water']:.6e}"
    )

    print()


print(
    "Plots saved:"
)

print(
    "data/kid_python_warm1_temperature.png"
)

print(
    "data/kid_python_warm1_qv.png"
)

print(
    "data/kid_python_warm1_condensed_water.png"
)

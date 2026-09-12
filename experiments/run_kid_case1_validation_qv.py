import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from experiments.run_kid_case1 import (
    KiDCase1Config,
    initialize_state,
    kid_case1_updraft,
    update_activation,
    update_liquid_growth,
    update_warm_rain,
)

OUTDIR = "validation_results/qv_sensitivity"
os.makedirs(OUTDIR, exist_ok=True)

qv_tests = {
    "RH100": 0.00860,
    "RH101": 0.00868,
    "Original_RH116": 0.01000,
}

def run_case(label, qv_value):
    config = KiDCase1Config()
    config.initial_qv = qv_value

    state = initialize_state(config)

    history = {
        "time": [], "cloud_mass": [], "rain_mass": [],
        "surface_rain_rate": [], "S": [], "Sc": [],
        "droplet_radius": [], "dr_dt": []
    }

    times = np.arange(0, config.t_end + config.dt, config.dt)

    for t in times:
        w = kid_case1_updraft(t, config.w1, config.t1, config.t2)

        state = update_activation(state, config)
        state = update_liquid_growth(state, config, w)

        if config.enable_rain:
            state = update_warm_rain(state, config)

        history["time"].append(t)
        history["cloud_mass"].append(state.get("cloud_mass", 0.0))
        history["rain_mass"].append(state.get("rain_mass", 0.0))
        history["surface_rain_rate"].append(state.get("surface_rain_rate", 0.0))
        history["S"].append(state.get("S", 0.0))
        history["Sc"].append(state.get("Sc", 0.0))
        history["droplet_radius"].append(state.get("droplet_radius", 0.0))
        history["dr_dt"].append(state.get("dr_dt", 0.0))

    df = pd.DataFrame(history)
    csv_path = f"{OUTDIR}/kid_case1_{label}.csv"
    df.to_csv(csv_path, index=False)

    return df

results = {}

for label, qv in qv_tests.items():
    print(f"Running {label}: initial_qv = {qv}")
    results[label] = run_case(label, qv)

plt.figure(figsize=(8, 6))

for label, df in results.items():
    plt.plot(df["time"], df["cloud_mass"], label=label)

plt.xlabel("Time (s)")
plt.ylabel("Cloud mass (kg/kg)")
plt.title("Initial Moisture Sensitivity Validation")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

fig_path = f"{OUTDIR}/initial_moisture_sensitivity_validation.png"
plt.savefig(fig_path, dpi=300)

print("Saved figure:", fig_path)
print("Saved CSV files in:", OUTDIR)

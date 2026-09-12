"""Plot aerosol-radius and kappa sensitivity sweeps."""

import os

import pandas as pd
import matplotlib.pyplot as plt


def main():

    os.makedirs("data", exist_ok=True)

    # --------------------------------------------------------
    # Read data
    # --------------------------------------------------------

    radius_df = pd.read_csv(
        "data/two_moment_radius_sweep_dt025.csv"
    )

    kappa_df = pd.read_csv(
        "data/two_moment_kappa_sweep_dt025.csv"
    )

    # Convert radius from m to micrometres
    radius_df["radius_um"] = (
        radius_df["aerosol_radius_m"]
        * 1.0e6
    )

    # Critical supersaturation in percent
    radius_df["Sc_percent"] = (
        100.0
        * radius_df["Sc_at_activation"]
    )

    kappa_df["Sc_percent"] = (
        100.0
        * kappa_df["Sc_at_activation"]
    )

    # --------------------------------------------------------
    # 1. Radius vs critical supersaturation
    # --------------------------------------------------------

    plt.figure()

    plt.plot(
        radius_df["radius_um"],
        radius_df["Sc_percent"],
        marker="o",
    )

    plt.xlabel(
        "Dry aerosol radius (µm)"
    )

    plt.ylabel(
        "Critical supersaturation Sc (%)"
    )

    plt.title(
        "Critical supersaturation vs aerosol radius\n"
        "(w = 1 m/s, N = 1e8 m⁻³, κ = 0.3, dt = 0.25 s)"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        "data/two_moment_Sc_vs_radius_dt025.png",
        dpi=200,
    )

    plt.close()

    # --------------------------------------------------------
    # 2. Radius vs time to activation
    # --------------------------------------------------------

    plt.figure()

    plt.plot(
        radius_df["radius_um"],
        radius_df["time_to_activation_s"],
        marker="o",
    )

    plt.xlabel(
        "Dry aerosol radius (µm)"
    )

    plt.ylabel(
        "Time from saturation to activation (s)"
    )

    plt.title(
        "Activation delay vs aerosol radius\n"
        "(w = 1 m/s, N = 1e8 m⁻³, κ = 0.3, dt = 0.25 s)"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        "data/two_moment_activation_delay_vs_radius_dt025.png",
        dpi=200,
    )

    plt.close()

    # --------------------------------------------------------
    # 3. Radius vs SSmax
    # --------------------------------------------------------

    plt.figure()

    plt.plot(
        radius_df["radius_um"],
        radius_df["SSmax_percent"],
        marker="o",
    )

    plt.xlabel(
        "Dry aerosol radius (µm)"
    )

    plt.ylabel(
        "Maximum supersaturation (%)"
    )

    plt.title(
        "SSmax vs aerosol radius\n"
        "(w = 1 m/s, N = 1e8 m⁻³, κ = 0.3, dt = 0.25 s)"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        "data/two_moment_SSmax_vs_radius_dt025.png",
        dpi=200,
    )

    plt.close()

    # --------------------------------------------------------
    # 4. Radius vs activated fraction
    # --------------------------------------------------------

    plt.figure()

    plt.plot(
        radius_df["radius_um"],
        radius_df["activated_fraction"],
        marker="o",
    )

    plt.xlabel(
        "Dry aerosol radius (µm)"
    )

    plt.ylabel(
        "Activated fraction"
    )

    plt.title(
        "Activated fraction vs aerosol radius\n"
        "(w = 1 m/s, N = 1e8 m⁻³, κ = 0.3, dt = 0.25 s)"
    )

    plt.grid(True)

    plt.ylim(
        0.0,
        1.05,
    )

    plt.tight_layout()

    plt.savefig(
        "data/two_moment_activated_fraction_vs_radius_dt025.png",
        dpi=200,
    )

    plt.close()

    # --------------------------------------------------------
    # 5. Kappa vs critical supersaturation
    # --------------------------------------------------------

    plt.figure()

    plt.plot(
        kappa_df["aerosol_kappa"],
        kappa_df["Sc_percent"],
        marker="o",
    )

    plt.xlabel(
        "Aerosol hygroscopicity κ"
    )

    plt.ylabel(
        "Critical supersaturation Sc (%)"
    )

    plt.title(
        "Critical supersaturation vs κ\n"
        "(w = 1 m/s, N = 1e8 m⁻³, radius = 0.05 µm, dt = 0.25 s)"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        "data/two_moment_Sc_vs_kappa_dt025.png",
        dpi=200,
    )

    plt.close()

    # --------------------------------------------------------
    # 6. Kappa vs time to activation
    # --------------------------------------------------------

    plt.figure()

    plt.plot(
        kappa_df["aerosol_kappa"],
        kappa_df["time_to_activation_s"],
        marker="o",
    )

    plt.xlabel(
        "Aerosol hygroscopicity κ"
    )

    plt.ylabel(
        "Time from saturation to activation (s)"
    )

    plt.title(
        "Activation delay vs κ\n"
        "(w = 1 m/s, N = 1e8 m⁻³, radius = 0.05 µm, dt = 0.25 s)"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        "data/two_moment_activation_delay_vs_kappa_dt025.png",
        dpi=200,
    )

    plt.close()

    # --------------------------------------------------------
    # 7. Kappa vs SSmax
    # --------------------------------------------------------

    plt.figure()

    plt.plot(
        kappa_df["aerosol_kappa"],
        kappa_df["SSmax_percent"],
        marker="o",
    )

    plt.xlabel(
        "Aerosol hygroscopicity κ"
    )

    plt.ylabel(
        "Maximum supersaturation (%)"
    )

    plt.title(
        "SSmax vs κ\n"
        "(w = 1 m/s, N = 1e8 m⁻³, radius = 0.05 µm, dt = 0.25 s)"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        "data/two_moment_SSmax_vs_kappa_dt025.png",
        dpi=200,
    )

    plt.close()

    # --------------------------------------------------------
    # 8. Kappa vs activated fraction
    # --------------------------------------------------------

    plt.figure()

    plt.plot(
        kappa_df["aerosol_kappa"],
        kappa_df["activated_fraction"],
        marker="o",
    )

    plt.xlabel(
        "Aerosol hygroscopicity κ"
    )

    plt.ylabel(
        "Activated fraction"
    )

    plt.title(
        "Activated fraction vs κ\n"
        "(w = 1 m/s, N = 1e8 m⁻³, radius = 0.05 µm, dt = 0.25 s)"
    )

    plt.grid(True)

    plt.ylim(
        0.0,
        1.05,
    )

    plt.tight_layout()

    plt.savefig(
        "data/two_moment_activated_fraction_vs_kappa_dt025.png",
        dpi=200,
    )

    plt.close()

    print(
        "Saved aerosol sensitivity plots:"
    )

    print(
        "  data/two_moment_Sc_vs_radius_dt025.png"
    )

    print(
        "  data/two_moment_activation_delay_vs_radius_dt025.png"
    )

    print(
        "  data/two_moment_SSmax_vs_radius_dt025.png"
    )

    print(
        "  data/two_moment_activated_fraction_vs_radius_dt025.png"
    )

    print(
        "  data/two_moment_Sc_vs_kappa_dt025.png"
    )

    print(
        "  data/two_moment_activation_delay_vs_kappa_dt025.png"
    )

    print(
        "  data/two_moment_SSmax_vs_kappa_dt025.png"
    )

    print(
        "  data/two_moment_activated_fraction_vs_kappa_dt025.png"
    )


if __name__ == "__main__":
    main()

# Python Cloud Model – Validation and KiD Comparison Summary

## 1. Objective

This validation work investigates the behaviour of the Python cloud microphysics model and compares its activation, warm-cloud, and mixed-phase behaviour with KiD-A.

The main diagnostics considered are:

- saturation time
- activation time
- time from saturation to activation
- activated fraction
- maximum supersaturation (SSmax)
- cloud droplet number concentration
- cloud water
- ice number concentration
- frozen-water mass
- depositional growth
- sedimentation effects

The aim is not only to identify differences between the Python model and KiD, but to determine which physical and numerical mechanisms produce those differences.

## 2. Warm-Cloud Activation Diagnostics

The Python model now diagnoses activation using the aerosol activation calculation rather than relying only on qualitative cloud-water behaviour.

Diagnostics include:

- first activation time
- saturation-to-activation delay
- fraction activated at first activation
- times to 1%, 10%, 50%, and 90% activation
- SSmax
- maximum cloud droplet number
- cloud-water onset
- maximum cloud water

## 3. Timestep Sensitivity

The clean fixed-cooling experiment was tested with:

- dt = 2.0 s
- dt = 1.0 s
- dt = 0.5 s
- dt = 0.25 s
- dt = 0.1 s

The main activation diagnostics converge closely as the timestep decreases.

For example:

- saturation time approaches approximately 743 s
- first activation occurs at approximately 777 s
- saturation-to-activation delay is approximately 34 s
- SSmax approaches approximately 1.363%
- maximum cloud water remains approximately 9.07e-4 kg/kg

The fraction present at the first discrete activation timestep is timestep-sensitive, so threshold-based diagnostics (1%, 10%, 50%, 90%) provide a more robust description of activation.

## 4. Updraft Sensitivity

Two forcing approaches were investigated.

### Fixed-cooling experiment

Changing updraft speed while retaining an imposed fixed cooling rate changed parcel height but did not materially change:

- saturation time
- activation time
- SSmax
- maximum cloud water

This demonstrated that updraft was not physically coupled to supersaturation production in that configuration.

### Adiabatic-proxy experiment

When cooling was coupled to vertical motion, a physically meaningful updraft sensitivity emerged.

Increasing updraft speed produced:

- earlier saturation
- earlier activation
- shorter saturation-to-activation delay
- larger SSmax
- larger cloud-water production

This identifies the thermodynamic coupling between vertical velocity and cooling as essential for interpreting updraft sensitivity.

## 5. Mixed-Phase KiD Comparison

A mixed-phase Python experiment was compared with the KiD mixed1 case using Thompson09 microphysics.

A same-layer comparison was made near:

- z = 400 m

The biological-IN concentration was varied to investigate ice-number sensitivity.

For bio_N = 2500 m^-3:

Python maximum ice number:
- approximately 2465 m^-3

KiD maximum ice number near 400 m:
- approximately 2450 m^-3

Therefore, ice-number concentrations can be closely matched while substantial differences remain in frozen-water mass.

## 6. Frozen-Phase Partitioning

KiD includes frozen condensate categories for:

- cloud ice
- snow
- graupel

In the mixed1 comparison at z ≈ 400 m over 0–3600 s, graupel mass remained zero, so the relevant total frozen condensate is cloud ice + snow.

At z ≈ 400 m, snow represents an important fraction of the frozen condensate.

The maximum total frozen mixing ratio in KiD is therefore larger than the cloud-ice mixing ratio alone.

This demonstrated that comparing Python qice only with KiD cloud ice alone is not a complete frozen-water comparison.

## 7. Sedimentation Sensitivity

KiD mixed1 was run with sedimentation enabled and disabled.

At t = 3600 s:

Sedimentation ON:
- frozen mixing ratio at 400 m = 1.4432e-7 kg/kg
- frozen column path = 2.0718e-4

Sedimentation OFF:
- frozen mixing ratio at 400 m = 6.2966e-7 kg/kg
- frozen column path = 9.9382e-4

OFF / ON ratios:

- frozen water at 400 m ≈ 4.36
- frozen column path ≈ 4.80

This demonstrates that sedimentation is a major sink of frozen condensate in KiD and explains an important part of the difference from the non-sedimenting single-parcel Python model.

## 8. Depositional-Growth Comparison

KiD Thompson09 process diagnostics were examined using:

- pri_inu
- pri_ide
- prs_ide
- prs_sde

Python deposition rates were converted from kg m^-3 s^-1 to kg kg^-1 s^-1 before comparison.

Over 0–3600 s:

Python integrated net deposition:
- 1.5246454e-5 kg/kg

KiD integrated frozen deposition at z ≈ 400 m:
- 2.906967e-6 kg/kg

Python / KiD integrated deposition ratio:
- approximately 5.24

The Python model therefore produces substantially larger cumulative depositional growth over the comparison period.

The instantaneous difference is strongly time-dependent. During some periods Python deposition is only about twice the KiD value, while during low-deposition phases in KiD the difference becomes much larger.

## 9. Current Interpretation

The mixed-phase discrepancy cannot be attributed to a single mechanism.

The evidence currently indicates contributions from:

1. stronger and more persistent cumulative depositional growth in the Python model;
2. partitioning of frozen condensate between ice and snow in Thompson09;
3. sedimentation and removal of frozen condensate in KiD;
4. structural differences between a simplified single-parcel Maxwell-growth model and the bulk Thompson09 microphysics scheme.

Therefore, agreement in ice number concentration does not imply agreement in frozen-water mass.

## 10. Remaining Validation Work

The principal validation experiments, publication figures, and quantitative summary tables are now complete and have undergone a final review of units, experiment provenance, and terminology.

Remaining work is focused on manuscript preparation and interpretation:

- integrate the validated quantitative results into the Methods and Results sections;
- develop the Discussion around the identified physical and numerical sources of model disagreement;
- prepare final figure captions and table descriptions that clearly distinguish baseline and sensitivity experiments;
- perform a final manuscript-level consistency check before submission.

## 11. Working Research Direction

A possible research direction emerging from the validation is:

> Which physical and numerical processes control differences in cloud activation and mixed-phase condensate evolution between a simplified parcel microphysics model and KiD/Thompson09?

The current experiments suggest that process-level diagnostics are essential for interpreting model disagreement rather than relying only on final cloud-water or ice-water amounts.

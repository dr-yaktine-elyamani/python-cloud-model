# Manuscript Draft

## Proposed Title

To be determined after the Results and Discussion are structured.

## Research Question

Which physical and numerical processes control differences in cloud activation and mixed-phase condensate evolution between a simplified parcel microphysics model and KiD/Thompson09?

## 1. Introduction

## 2. Methods

### 2.1 Simplified Python parcel microphysics model
A simplified parcel microphysics model was used to investigate warm-cloud activation and condensational growth under controlled thermodynamic forcing. Water-vapour supersaturation was calculated from the vapour pressure and saturation vapour pressure over liquid water, with the latter evaluated using a Bolton-style formulation. Aerosol activation was represented using a κ-Köhler critical supersaturation criterion. For an aerosol particle of dry diameter Dp and hygroscopicity parameter κ, the critical supersaturation was calculated from the Kelvin curvature term. Below the critical supersaturation no particles were activated; between Sc and 2Sc, the activated fraction increased linearly from zero to unity; and at S ≥ 2Sc, complete activation was assumed.

Condensational droplet growth was represented using a simplified Maxwell-type growth law, dr/dt = (G/r)S, where r is the droplet radius, S is the water supersaturation, and G is the liquid-water growth coefficient. The growth coefficient was set to 8.0 × 10^-12 m^2 s^-1. Radius growth was converted to a bulk cloud-water mixing-ratio tendency using the diagnosed droplet number concentration, a water density of 1000 kg m^-3, and an air density of 1.0 kg m^-3. Condensation was permitted only when S > 0 and activated droplets were present. For the controlled warm-cloud benchmark, negative radius growth was suppressed, so droplet evaporation was not represented by this growth routine.
### 2.2 Warm-cloud activation diagnostics
Warm-cloud activation was quantified using a set of time-resolved diagnostics designed to distinguish the onset of saturation from the progressive activation of the aerosol population. Saturation time was defined as the first model time at which water supersaturation became non-negative (S >= 0), while first activation time was defined as the first time at which the diagnosed droplet number concentration became greater than zero. The delay from saturation to first activation was calculated as the difference between these two times. Because the first non-zero activated fraction can be sensitive to the discrete model timestep, additional threshold-based diagnostics were defined as the first times at which 1%, 10%, 50%, and 90% of the aerosol population had activated (t1, t10, t50, and t90, respectively).

Additional bulk diagnostics were used to characterize the warm-cloud response. Maximum supersaturation (SSmax) was defined as the maximum water supersaturation attained during the simulation and was reported as a percentage. The maximum droplet number concentration (Nc,max) and maximum cloud-water mixing ratio (qc,max) were also recorded. Cloud-water onset was defined as the first time at which the cloud-water mixing ratio reached 1.0 × 10^-7 kg kg^-1. Together, these diagnostics allowed activation timing, activation completeness, supersaturation evolution, and condensate development to be compared consistently across the sensitivity experiments.
### 2.3 Abdul-Razzak et al. (1998) activation benchmark
The aerosol-activation formulation was benchmarked against the single-aerosol parameterization of Abdul-Razzak et al. (1998) using the conditions associated with their Figure 5 test case. The benchmark used a temperature of 283.15 K, pressure of 800 hPa, updraft velocity of 5 m s^-1, aerosol number concentration of 200 cm^-3, median dry radius of 0.01 µm, and geometric standard deviation of 2.5. Ammonium sulfate was used as the aerosol composition. The calculation evaluated the critical supersaturation, the dimensionless eta and zeta parameters, the predicted maximum supersaturation, and the resulting activated fraction.

For the benchmark calculation, water surface tension was evaluated using the IAPWS correlation, while water-vapour diffusivity followed the temperature- and pressure-dependent scaling used elsewhere in the parcel-model framework. The model-calculated critical supersaturation, eta, and zeta were compared with the corresponding quantities reported for the Abdul-Razzak et al. (1998) benchmark. In addition, the published benchmark values of critical supersaturation, eta, and zeta were inserted directly into the Abdul-Razzak algebra to provide comparison values for maximum supersaturation and activated fraction. This separated evaluation of the implemented thermophysical closure from verification of the parameterization algebra, without fitting the model to the published benchmark values.
### 2.4 Aerosol and timestep sensitivity experiments
Aerosol sensitivities were evaluated using a controlled warm-cloud baseline with an initial temperature of 283.15 K, pressure of 99398.334 Pa, water-vapour mixing ratio of 0.0073923 kg kg^-1, timestep of 0.5 s, total simulation time of 3600 s, and a prescribed cooling rate of 0.001 K s^-1. Vertical velocity was set to zero so that the activation response could be isolated from dynamical forcing. The baseline aerosol population used a number concentration of 50 cm^-3, dry radius of 0.05 µm, and hygroscopicity parameter κ = 0.3.

Three one-at-a-time aerosol sensitivity experiments were performed. Aerosol number concentration was varied between 25, 50, 100, and 200 cm^-3; dry particle radius was varied between 0.03, 0.05, 0.07, and 0.10 µm; and κ was varied between 0.1, 0.3, 0.6, and 1.0. In each experiment, only the parameter under investigation was varied while the remaining aerosol and thermodynamic parameters were held at their baseline values. The activation diagnostics defined above, including threshold activation times, maximum supersaturation, maximum droplet number concentration, and maximum cloud-water mixing ratio, were evaluated for every simulation.

Numerical sensitivity was assessed using a timestep-convergence experiment based on the same controlled warm-cloud configuration. Timesteps of 2.0, 1.0, 0.5, 0.25, and 0.1 s were tested while all physical parameters and forcing conditions were held fixed. Convergence was evaluated using saturation and activation timing, the threshold-based activation diagnostics, maximum supersaturation, and maximum cloud-water mixing ratio. This timestep analysis was used to assess numerical sensitivity of the current explicit sequential integration structure; no alternative time-integration solver was implemented for solver-to-solver comparison.
### 2.5 KiD/Thompson09 mixed-phase comparisonNumerical sensitivity was assessed using a timestep-convergence experiment based on the same controlled warm-cloud configuration. Timesteps of 2.0, 1.0, 0.5, 0.25, and 0.1 s were tested while all physical parameters and forcing conditions were held fixed. Convergence was evaluated using saturation and activation timing, the threshold-based activation diagnostics, maximum supersaturation, and maximum cloud-water mixing ratio. This timestep analysis was used to assess numerical sensitivity of the current explicit sequential integration structure; no alternative time-integration solver was implemented for solver-to-solver comparison.
The mixed-phase comparison used output from the simplified Python parcel model and the KiD testbed configured with the Thompson09 microphysics scheme. The comparison was restricted to the first 3600 s and to the KiD model level at 400 m. To isolate differences in condensate evolution from differences in ice number concentration, a matched-number Python sensitivity experiment with a biological ice-nucleating particle concentration of 2500 m^-3 was used. This configuration produced a maximum ice number concentration close to that simulated by KiD/Thompson09 over the comparison period.

Python ice mass was diagnosed internally as mass per unit volume (kg m^-3) and was converted to an ice mixing ratio (kg kg^-1) using the ideal-gas estimate of air density, rho_air = P/(Rd T), with P = 101325 Pa and Rd = 287.05 J kg^-1 K^-1. The converted Python ice mixing ratio was compared with KiD cloud-ice mixing ratio and with the sum of KiD cloud ice and snow. Ice number concentration was also compared directly between the two models. This allowed ice-number agreement and frozen-condensate disagreement to be evaluated separately.
### 2.6 Mixed-phase process diagnostics
Mixed-phase process attribution was examined using integrated depositional growth and sedimentation-sensitivity diagnostics. Positive depositional growth was integrated over the 0–3600 s comparison period for the matched-number Python experiment and for KiD/Thompson09 at 400 m. This provided a direct measure of the cumulative vapour-to-ice depositional tendency represented by each model over the common analysis interval.

The role of sedimentation was assessed separately within KiD/Thompson09 by comparing simulations with sedimentation enabled and disabled. Frozen condensate at 400 m was evaluated for both configurations, allowing the effect of sedimentation on the retained frozen-water mass to be isolated from the depositional-growth comparison. This diagnostic was used to determine how much of the Python–KiD condensate difference could be attributed to sedimentation rather than to deposition alone.
## 3. Results

### 3.1 Activation formulation and numerical validation
The Abdul-Razzak et al. (1998) benchmark showed close agreement between the implemented activation parameterization and the published benchmark quantities. The calculated critical supersaturation was 0.01731 compared with 0.01762 for the benchmark, corresponding to a relative difference of -1.76%. The calculated eta and zeta parameters differed from the published values by +1.68% and -0.12%, respectively. Using the complete implemented parameterization, the predicted maximum supersaturation was 1.767%, compared with 1.764% obtained from the published-parameter algebra, a difference of 0.16%. The corresponding activated fractions were 0.506 and 0.500, respectively, representing a relative difference of 1.12%.
The timestep-convergence experiment showed that the principal warm-cloud diagnostics were only weakly sensitive to timestep over the tested range. Reducing the timestep from 2.0 to 0.1 s changed SSmax from 1.3658% to 1.3628%, a relative difference of 0.224%. The 50% activation time changed from 794.0 to 794.3 s, corresponding to a relative difference of -0.038%, while the maximum cloud-water mixing ratio changed by only 0.00083%. These results indicate convergence of the principal bulk and threshold-based activation diagnostics over the tested timestep range.
### 3.2 Aerosol controls on warm-cloud activation
Aerosol number concentration exerted its strongest influence on the maximum supersaturation rather than on the onset of activation. Increasing aerosol number concentration from 25 to 200 cm^-3 reduced SSmax from 1.931% to 0.680%, while first activation remained at 777 s across all four simulations. The 50% activation time changed only slightly, from 794.0 to 794.5 s, whereas the maximum droplet number concentration scaled directly with the prescribed aerosol number concentration. Cloud-water onset occurred progressively earlier as aerosol number increased, from 801.0 s at 25 cm^-3 to 787.0 s at 200 cm^-3, while qc,max changed only modestly from 9.14 × 10^-4 to 9.01 × 10^-4 kg kg^-1.

Dry particle radius produced a substantially stronger response in activation timing. Increasing the dry radius from 0.03 to 0.10 µm advanced first activation from 816.0 to 755.0 s and reduced t50 from 854.0 to 761.0 s. The delay between saturation and first activation therefore decreased from 73.5 to 12.5 s. Over the same range, SSmax decreased from 1.400% to 1.357%. In contrast, qc,max remained nearly unchanged at approximately 9.07 × 10^-4 kg kg^-1.
### 3.3 Mixed-phase comparison with KiD/Thompson09

### 3.4 Deposition, frozen-phase partitioning, and sedimentation

## 4. Discussion

### 4.1 Physical controls on activation differences

### 4.2 Numerical sensitivity and model structure

### 4.3 Why similar ice number does not imply similar frozen-water mass

### 4.4 Limitations of the simplified parcel framework

## 5. Conclusions

## References

ALGORITHM
- [x] Use DataFrames.jl to perform data retrieval, selection and plot
- [ ] Time efficient code: write on file every x times, with x reasonably large
- [x] Implement automatic safecopy of data
- [x] Employ symmetry: many calculations are just pure repetitions
- [x] Use normalized K
- [x] Re-use previous value in step for fast convergence
- [x] Restrict to MBZ calculation AND CHECK MULTIPLIERS!
- [x] Insert graphs in report
- [x] Investigate superconducting phase
- [ ] Finite size scaling on the boundary
- [x] Extract and plot Fermi surface
- [x] Fake rigid shift of hopping
- [x] Implement bands renormalization in superconducting phase
- [x] Move RenormalizeBands => RenormalizeBands (d-wave part of g)
- [ ] Implement triplet simulations
- [ ] Correct free energy and compute AF free energy
- [ ] New structure for AF simulations
- [ ] Discard record-g mode and implement single mode
- [ ] Test new general plot function

HM: AF
- [ ] heatmaps
- [x] scan
- [ ] interesting U scan: vary temperature, vary doping
- [ ] write free energy module

HM: SU-SINGLET
- [ ] heatmaps s-wave (negative U too)
- [ ] scan s-wave (negative U too)
- [ ] write free energy module

EHM: AF
- [ ] half-filling in depth analysis
- [ ] heatmaps: keep doping null
- [ ] scan: keep doping null
- [ ] write free energy module
- [ ] plot the metallic bands

EHM: FAKEAF
- [ ] do it later

EHM: SU-SINGLET
- [ ] sanity check: compare optimized weights run vs dummy run
- [ ] filtered run: re-compute the NaN points => write the filtered run module
- [ ] heatmaps+RMPs s+s*-wave (negative U too)
- [ ] heatmaps+RMPs d-wave
- [ ] heatmaps+RMPs s+s*+d-wave

EHM: FAKESU-SINGLET
- [ ] do it later

EHM: SU-TRIPLET

EHM: FAKESU-TRIPLET
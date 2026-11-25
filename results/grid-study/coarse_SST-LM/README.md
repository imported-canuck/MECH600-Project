# Coarse grid SST-LM run

NACA 0012 at Reynolds number 130,000 and 10 degrees incidence. Wall-resolved SST-LM with 0.16% nominal inlet turbulence, a 0.20 m chord and the 25-chord outer-extent setting. The OpenFOAM mesh matches the archived coarse SST mesh byte for byte. The solver ran a 1,500-iteration pilot and resumed to iteration 7,000.

- [Convergence history](report-7000/convergence.png)
- [Surface diagnostics](report-7000/surface-diagnostics.png)
- [Velocity field](report-7000/flow-field.png)
- [Reattachment history](report-7000/reattachment-history.csv)
- [Final summary](report-7000/summary.csv)
- [Pilot summary](report-1500/summary.csv)
- [OpenFOAM case](cases/low_lm_016/baseline.foam)

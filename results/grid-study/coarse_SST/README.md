# Coarse grid SST run

NACA 0012 at Reynolds number 130,000 and 10 degrees incidence. Wall-resolved SST, nominal inlet turbulence of 0.16%, 0.20 m chord, and farfield radius of 25 chords. OpenFOAM v2312 used 4 MPI ranks and the same startup and final schemes as the medium run. The solver ran a 1,500-iteration pilot and resumed to iteration 7,000.

| Quantity | Result |
| --- | --- |
| Cells | 93312 |
| Upper-surface separation x/c | 0.0151224 |
| Upper-surface reattachment x/c | 0.0718209 |
| Reattachment span across last three checkpoints | 0.000169654c |
| Maximum upper-surface y+ | 0.154945 |
| Mean CL / CD over last 400 iterations | 0.907892 / 0.0315306 |
| Peak-to-peak CL / CD over last 400 iterations | 0.00164985 / 0.000194422 |
| Force plateau criterion met | False |
| Largest final initial residual | 0.00333477 |
| Solver wall-clock seconds | 637.149 |

- [Convergence history](report-7000/convergence.png)
- [Surface diagnostics](report-7000/surface-diagnostics.png)
- [Velocity field](report-7000/flow-field.png)
- [ParaView midspan image](report-7000/paraview-velocity.png)
- [Pilot report](report-1500/summary.csv)
- [Reconstructed OpenFOAM case](cases/low_sst_016/baseline.foam)

The report folders also contain CSV histories. Open the case marker in ParaView and select the saved iteration. To load the ParaView state, select the local case marker when prompted for the input file.

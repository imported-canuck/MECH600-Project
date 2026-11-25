# Medium kkLOmega on the blockMesh C-grid

Wall-resolved NACA 0012 at Reynolds number 130,000, 10 degrees incidence, and 0.16% nominal inlet turbulence. This run reused the verified 209,952-cell medium mesh and includes the 1,500-iteration pilot and the 7,000-iteration endpoint.

| Metric at 7,000 iterations | Value |
| --- | ---: |
| First separation x/c | 0.01300 |
| First reattachment x/c | 0.05037 |
| Reattachment span over the last three checkpoints | 0.030 mm |
| Mean lift coefficient over the final force window | 0.83522 |
| Mean drag coefficient over the final force window | 0.10104 |
| Peak upper wall y+ | 0.091 |
| Peak lower wall y+ | 0.090 |
| Force plateau criterion | Missed |
| Solver wall time in seconds | 1787 |

Both kkLOmega runs use upwind convection and 0.15 relaxation for `kt`, `kl`, and `omega`. Those settings were chosen for stability and also differ from the previous SST and SST-LM turbulence numerics. The plateau criterion and checkpoint span describe convergence separately from reaching the iteration cap.

- [Convergence history](report-7000/convergence.png)
- [Surface diagnostics](report-7000/surface-diagnostics.png)
- [Velocity field](report-7000/flow-field.png)
- [Reattachment history](report-7000/reattachment-history.csv)
- [Final summary](report-7000/summary.csv)
- [Pilot summary](report-1500/summary.csv)
- [OpenFOAM case](cases/low_kkl_016/baseline.foam)
- [ParaView velocity](report-7000/paraview-velocity.png)
- [ParaView leading edge](report-7000/paraview-leading-edge.png)

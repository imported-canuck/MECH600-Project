# MECH 600 Step 4 baseline

Fresh OpenFOAM v2312 NACA 0012 simulation on four ranks using the medium, 209,952-cell grid. The wall-resolved k–ω SST case used a 0.20 m chord, 10 m/s freestream, 10° incidence, chord Reynolds number 130,000, 25-chord farfield radius and 0.16% nominal inlet turbulence. The run started from uniform fields, with 500 startup iterations followed by the final schemes.

The cluster regenerated the grid using source commit 80740db. Node coordinates rounded to 0.1 nm, physical tags and element connectivity matched the Windows grid. The Gmsh files had different byte hashes. The cluster scratch quota interrupted the iteration-1,500 field write. The job restarted from the intact iteration-1,250 checkpoint in cluster home storage and then continued from the recovered iteration-1,500 fields. Slurm and OpenFOAM logs record these stages.

## Iteration 7,000

| Saved iteration | Upper-surface separation x/c | Upper-surface reattachment x/c |
| --- | --- | --- |
| 6500 | 0.0149609 | 0.0690205 |
| 6750 | 0.0148226 | 0.0691907 |
| 7000 | 0.014704 | 0.0693744 |

The last three checkpoints span 0.000354c (0.071 mm). The iteration target for reattachment was a span of at most 0.002c (0.4 mm). The positions come from linearly interpolated zero crossings of signed wall skin friction. The tunnel comparison point is approximately x/c = 0.10, as reported alongside measured mean velocity vectors. [Hata et al. (2002), Figure 4 and §3.2](https://www.jstage.jst.go.jp/article/jjsass/50/582/50_582_293/_pdf)

- Pressure initial residual went from 0.00556247 at iteration 501.0 in the final scheme to 0.0032716 at iteration 7000.0. The largest final initial residual was 0.0032716. Normalized flux imbalance was 3.25665e-08.
- Mean C_L was 0.965789 and mean C_D was 0.0300168 over the final 400 iterations. Their peak-to-peak ranges were 0.00444828 and 0.000258146, respectively. The run exceeded the project's force tolerances of 0.0001 and 0.00001.
- Maximum y+ was 0.104852 on the upper surface and 0.105059 on the lower surface. The modeled upstream turbulence intensity at the near-airfoil probe was 0.150526%.
- The accumulated solver stages used 1632.45 wall-clock seconds and 1.81384 core-hours.

## Mesh audit

OpenFOAM passed the basic mesh checks. The extended check flagged 30,144 cells with determinant below 0.001. Maximum non-orthogonality was 44.863°, maximum skewness was 1.228 and maximum aspect ratio was 448.77. See the OpenFOAM mesh logs and [mesh-audit.png](mesh-audit.png).

## Files

- [Convergence](report-7000/convergence.png), [surface diagnostics](report-7000/surface-diagnostics.png), and [velocity overview](report-7000/flow-field.png) show the final solver data. The corresponding CSVs record the force, residual, surface and checkpoint histories.
- [ParaView velocity](report-7000/paraview-velocity.png) and [leading-edge detail](report-7000/paraview-leading-edge.png) show the midspan slice of the reconstructed velocity field. Open `cases/low_sst_016/baseline.foam` in ParaView and select iteration 7,000. When loading the saved ParaView state, point the reader to that local file. Iteration numbers represent steady SIMPLE solver steps.
- The OpenFOAM case includes the mesh, reconstructed fields and dictionaries. The earlier report records the iteration-1,500 pilot. The `logs` folder records the cluster runs and post-processing.

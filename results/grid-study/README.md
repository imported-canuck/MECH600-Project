# NACA 0012 grid study

The three archived SST runs used wall-resolved SST at Reynolds number 130,000, 10 degrees incidence, and 0.16% nominal inlet turbulence. Each case used a 1,500-iteration pilot and a 7,000-iteration endpoint.

| Grid | Cells | Separation x/c | Reattachment x/c | Late reattachment span (mm) | Main figures |
| --- | ---: | ---: | ---: | ---: | --- |
| [coarse](coarse_SST/README.md) | 93,312 | 0.01512 | 0.07182 | 0.034 | [Convergence](coarse_SST/report-7000/convergence.png), [velocity](coarse_SST/report-7000/flow-field.png), [surface](coarse_SST/report-7000/surface-diagnostics.png) |
| [medium](medium_SST/README.md) | 209,952 | 0.01470 | 0.06937 | 0.071 | [Convergence](medium_SST/report-7000/convergence.png), [velocity](medium_SST/report-7000/flow-field.png), [surface](medium_SST/report-7000/surface-diagnostics.png) |
| [fine](fine_SST/README.md) | 472,392 | 0.02013 | 0.06279 | 1.214 | [Convergence](fine_SST/report-7000/convergence.png), [velocity](fine_SST/report-7000/flow-field.png), [surface](fine_SST/report-7000/surface-diagnostics.png) |

The late span covers iterations 6,500, 6,750 and 7,000. The iteration target for reattachment is 0.4 mm.

**Reattachment GCI status**  A late reattachment span exceeds 0.002 chord

- [Three-grid comparison figure](comparison/grid-comparison.png)
- [Comparison table and GCI diagnostics](comparison/grid-comparison.csv)
- [Hata et al. reference, Figure 4 and Section 3.2](https://www.jstage.jst.go.jp/article/jjsass/50/582/50_582_293/_pdf#page=3)

The individual folders contain both iteration reports, reconstructed OpenFOAM fields, mesh audits, logs, and ParaView images.

## SST and SST-LM model comparison

The wall-resolved SST-LM cases use the same meshes as the corresponding SST cases. The comparison uses the archived 7,000-iteration reports.

- [Coarse SST-LM](coarse_SST-LM/README.md)
- [Medium SST-LM](medium_SST-LM/README.md)
- [Fine SST-LM](fine_SST-LM/README.md)
- [Matched-mesh comparison figure](comparison/model-comparison.png)
- [Matched-mesh comparison table](comparison/model-comparison.csv)

At 7,000 iterations, SST-LM gives reattachment closer to Hata's approximately 0.10c on all three matched meshes. Its coarse and medium cases reattach around 0.116c and its fine case around 0.0714c. The coarse, medium and fine reattachment spans across iterations 6,500, 6,750 and 7,000 are 0.152 mm, 0.414 mm and 3.91 mm. The medium and fine cases exceed the 0.4 mm iteration target, so these endpoints provide a provisional model comparison.

## Medium k-kl-omega

The 0.16% nominal inlet turbulence case reached 7,000 iterations on the original medium grid. First reattachment was at 0.04430c and its span over the last three checkpoints was 1.520 mm. The run missed the force plateau criterion.

Both k-kl-omega cases used upwind convection and 0.15 relaxation for the turbulence fields. The medium comparison records these settings alongside the SST and SST-LM results.

- [Original-grid case and reports](medium_KKL/README.md)
- [Comparison with the blockMesh case](../../results2/comparison/medium-mesh-comparison.png)

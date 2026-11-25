# Medium NACA 0012 blockMesh reruns

The SST and SST-LM cases share a 12-block C-grid based on the high-re-v2 meshing approach. The mesh retains the original medium case's 209,952 cells, 1.317 µm target full first-cell height, 0.20 m chord, 25-chord domain setting and single spanwise cell. The blocks concentrate cells near the airfoil and expand the inner wake strip toward the outlet. The outer boundary shape and cell distribution differ from the original Gmsh-format grid.

Both cases use the original medium flow settings of Reynolds number 130,000, 10 degrees incidence and 0.16% nominal inlet turbulence. Both jobs completed the 1,500-iteration pilot and the 7,000-iteration final run.

| Case | Slurm job | Destination for results |
| --- | ---: | --- |
| SST | 1257 | [Case and reports](medium_SST/README.md) |
| SST-LM | 1258 | [Case and reports](medium_SST-LM/README.md) |

## Mesh audit

OpenFOAM v2312 produced 209,952 cells and passed standard mesh checks. Maximum non-orthogonality was 57.02 degrees, maximum skewness was 1.321, and maximum aspect ratio was 480.6. The only failed extended check flagged 33,154 cells below the determinant threshold of 0.001. The first-cell height is a design input. Wall y+ is evaluated from the solved flow.

- [Mesh generator](../mesh/generate_blockmesh.py)
- [Generated blockMeshDict](mesh/system/blockMeshDict)
- [Standard checkMesh log](mesh/log.checkMesh.basic)
- [Extended checkMesh log](mesh/log.checkMesh)
- [Shared OpenFOAM mesh hashes](mesh/polyMesh.sha256)
- [McGill meshing tutorial](https://mcgillteaching.github.io/MECH501_tutorial/openfoam/meshing.html#which-one-for-your-project-case)

The scripts folder contains the cluster build and run scripts. On the cluster, the generator and scripts were placed together in the run directory. The build script calls blockMesh and checkMesh. The run script takes the shared run directory, the archived medium SST baseline directory and the model name.

## Comparison with the original medium grid

The two reruns use the same audited OpenFOAM mesh. The original and blockMesh comparisons use 209,952 cells and the iteration-7,000 reports. Hata et al. measured reattachment at approximately 0.10c.

| Model | Original reattachment x/c | blockMesh reattachment x/c | blockMesh displacement from Hata | blockMesh late span |
| --- | ---: | ---: | ---: | ---: |
| SST | 0.06937 | 0.07312 | 5.38 mm upstream | 0.056 mm |
| SST-LM | 0.11664 | 0.14518 | 9.04 mm downstream | 0.494 mm |

The late span covers iterations 6,500, 6,750, and 7,000. Original SST and SST-LM spans were 0.071 mm and 0.414 mm. The blockMesh SST span meets the 0.4 mm target and SST-LM is slightly above it. The blockMesh runs had peak wall y+ of 0.104 and 0.100, respectively. Both first separation points fall inside Hata's measured range of 0.011c to 0.025c. Both force histories missed the plateau criterion.

SST-LM gave the closer reattachment location on the original medium grid. SST gives the closer location on the blockMesh grid. The difference reflects a change in the mesh layout, including cell placement and outer-domain shape.

- [Medium-mesh comparison plot](comparison/medium-mesh-comparison.png)
- [Medium-mesh comparison table](comparison/medium-mesh-comparison.csv)
- [Comparison script](scripts/compare-medium.py)
- [Hata et al. Figure 4 and Section 3.2](https://www.jstage.jst.go.jp/article/jjsass/50/582/50_582_293/_pdf#page=3)

## Medium k-kl-omega

Both meshes reached 7,000 iterations at 0.16% nominal inlet turbulence. The table records first reattachment because both cases predicted multiple separated regions.

| Mesh | First reattachment x/c | Distance upstream of Hata | Late reattachment span |
| --- | ---: | ---: | ---: |
| Original C-grid | 0.04430 | 11.14 mm | 1.520 mm |
| blockMesh C-grid | 0.05037 | 9.93 mm | 0.030 mm |

The span uses iterations 6,500, 6,750 and 7,000. Both force histories missed the plateau criterion. The k-kl-omega cases use upwind convection and 0.15 relaxation for the turbulence equations, so the comparison includes changes in the model and in its numerical settings.

- [blockMesh case and reports](medium_KKL/README.md)
- [Original-grid case and reports](../results/grid-study/medium_KKL/README.md)
- [Comparison across all three models](comparison/medium-mesh-comparison.png)
- [Comparison data](comparison/medium-mesh-comparison.csv)

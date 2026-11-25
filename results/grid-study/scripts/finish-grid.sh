#!/bin/bash -l
set -euo pipefail
run=$1
case_path="$run/cases/low_sst_016"
iteration=$(python3 -c 'import json,sys; s=json.load(open(sys.argv[1])); assert s["status"]=="finished" and s["exit_code"]==0; print(s["latest_iteration"])' "$case_path/run-status.json")
[[ "$iteration" = 7000 ]]
for field in U p k omega nut phi wallShearStress yPlus; do
  test -s "$case_path/$iteration/$field"
done
if [[ -s "$case_path/constant/polyMesh/sets/underdeterminedCells" ]]; then
  ( module load openfoam/v2312 gmsh/4.13.1 && python3 "$run/plot-mesh-audit.py" "$run" ) > "$run/logs/mesh-visual.log" 2>&1
fi
for checkpoint in 1500 "$iteration"; do
  report="$run/report-$checkpoint"
  test -s "$report/summary.csv"
  ( module load paraview/6.0.0 && OMP_NUM_THREADS=1 mpirun -np 1 --host localhost pvbatch --force-offscreen-rendering "$run/render-paraview.py" "$case_path" "$checkpoint" "$report" ) > "$run/logs/paraview-$checkpoint.log" 2>&1
  tail -1 "$run/logs/paraview-$checkpoint.log"
done

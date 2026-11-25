#!/bin/bash -l
set -euo pipefail
run=$1
step=$2
case_path="$run/cases/low_sst_016"
if [[ ! -s "$case_path/$step/p" ]]; then
  grep -aq "^Time = $((step + 5))$" "$case_path/log.solve"
  for rank in "$case_path"/processor[0-9]*; do
    for field in U p k omega nut wallShearStress yPlus; do test -s "$rank/$step/$field"; done
  done
  module load openfoam/v2312
  reconstructPar -case "$case_path" -time "$step" > "$run/logs/reconstruct-$step.log" 2>&1
fi
python3 "$run/read-checkpoint.py" "$run" "$step"

#!/bin/bash -l
set -euo pipefail
umask 077
run=$1
bash -l "$run/run-baseline.sbatch" "$run" --mesh-backend gmsh --initialization uniform --warmup 500 --write-interval 250 --iterations 1500
if ! bash -l "$run/collect-baseline.sh" "$run" > "$run/logs/pilot-report.log" 2>&1; then
  echo 'Pilot report failed; continuing the solver from the saved checkpoint.'
fi
bash -l "$run/run-baseline.sbatch" "$run" --resume --iterations 7000
bash -l "$run/collect-baseline.sh" "$run" > "$run/logs/final-report.log" 2>&1

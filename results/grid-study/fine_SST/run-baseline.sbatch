#!/bin/bash -l
set -eo pipefail
run_root=$1
shift
module load openfoam/v2312 gmsh/4.13.1
command -v gmshToFoam checkMesh simpleFoam
cd "$run_root/code"
exec python3 study.py run low_sst_016 --work-dir "$run_root/cases" --flavor openfoam-com --ranks "$SLURM_NTASKS" "$@"

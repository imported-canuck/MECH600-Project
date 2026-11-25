#!/usr/bin/env bash
#SBATCH --job-name=naca0012
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

project=${NACA_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}
work_dir=${NACA_WORK_DIR:-${SCRATCH:-$project}/naca0012-runs}
python=${NACA_PYTHON:-python3}
launcher=${NACA_MPI_LAUNCHER:-mpirun}
target=${1:?Usage: sbatch slurm/run.sbatch CASE_OR_ARRAY_GROUP [study.py run options]}
shift
cd "$project"

# An optional site-specific file can load modules or source an OpenFOAM bashrc.
if [[ -n ${NACA_FOAM_ENV:-} ]]; then
    set +u
    source "$NACA_FOAM_ENV"
    set -u
fi

if [[ -n ${SLURM_ARRAY_TASK_ID:-} ]]; then
    listing=$("$python" study.py list --group "$target" --names)
    mapfile -t cases <<< "$listing"
    if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#cases[@]} )); then
        printf 'Array index %s is out of range for group %s (%s cases).\n' \
            "$SLURM_ARRAY_TASK_ID" "$target" "${#cases[@]}" >&2
        exit 2
    fi
    target=${cases[SLURM_ARRAY_TASK_ID]}
fi

exec "$python" study.py run "$target" --work-dir "$work_dir" \
    --ranks "${SLURM_NTASKS:-1}" --launcher "$launcher" "$@"

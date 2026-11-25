#!/bin/bash -l
# Run and archive one model using the common medium blockMesh grid.
set -eo pipefail
umask 077

if [[ $# -ne 3 || ( "$3" != "SST" && "$3" != "SST-LM" ) ]]; then
    echo "Usage: sbatch run-medium.sbatch RUN_ROOT BASELINE_ROOT SST|SST-LM" >&2
    exit 2
fi
run=$(realpath "$1")
baseline=$(realpath "$2")
model=$3
name="medium_$model"
output="$run/$name"
test -s "$run/mesh/mesh-status.json"
test -d "$output"
if [[ -e "$output/log.warmup" || -e "$output/log.solve" ]]; then
    echo "Solver history already exists at $output" >&2
    exit 2
fi

module load openfoam/v2312
export FOAM_FILEHANDLER=uncollated
export OMP_NUM_THREADS=1
work=$(mktemp -d "${SLURM_TMPDIR:-/tmp}/mech600-results2-${SLURM_JOB_ID}.XXXXXX")
printf '%s\n' "$work" > "$output/node-workdir"
hostname -s > "$output/node"
mkdir -p "$work/logs"
cp -R --no-preserve=all -- "$baseline/code" "$work/code"
cp -- "$run/generate_blockmesh.py" "$work/code/mesh/"
cp -- "$run/solve-on-mesh.py" "$work/"
case_id=low_sst_016
if [[ "$model" == SST-LM ]]; then case_id=low_lm_016; fi
for helper in collect-baseline.sh plot-step4.py render-paraview.py; do
    sed -e "s/low_sst_016/$case_id/g" \
        -e "s/wall-resolved SST/wall-resolved $model/g" \
        -e "s/ | SST | / | $model | /g" \
        -e 's| gmsh/4\.13\.1||g' \
        "$baseline/$helper" > "$work/$helper"
done

for step in 1500 7000; do
    python3 "$work/solve-on-mesh.py" "$baseline" "$run/mesh" "$work" "$output" "$model" "$step"
    bash "$work/collect-baseline.sh" "$work" > "$output/report-$step.log" 2>&1
    cp -- "$work/report-$step/summary.csv" "$output/summary-$step.csv"
done
case_dir="$work/cases/$case_id"
for field in U p k omega wallShearStress yPlus; do test -s "$case_dir/7000/$field"; done
if [[ "$model" == SST-LM ]]; then
    for field in gammaInt ReThetat; do test -s "$case_dir/7000/$field"; done
fi
touch "$case_dir/baseline.foam"
if ! (module load paraview/6.0.0 && mpirun -np 1 --host localhost \
    pvbatch --force-offscreen-rendering "$work/render-paraview.py" "$case_dir" 7000 "$work/report-7000"
) > "$output/paraview.log" 2>&1; then
    echo "ParaView rendering failed. The diagnostic plots and fields are available."
fi
cp -- "$work/report-7000/reattachment-history.csv" "$output/reattachment-history.csv"

# Copy only the final full field while retaining both visual and CSV reports.
archive="$output/$name.tar.gz"
tar --dereference -czf "$archive.part" -C "$work" \
    code logs report-1500 report-7000 collect-baseline.sh plot-step4.py render-paraview.py solve-on-mesh.py \
    "cases/$case_id/"{0,constant,system,7000,postProcessing,case.json,mesh.json,run-status.json,baseline.foam,log.blockMesh,log.warmup,log.solve,log.checkMesh,log.checkMesh.basic}
tar -tzf "$archive.part" > "$output/manifest.txt"
mv -- "$archive.part" "$archive"
echo "Verified final archive: $archive"

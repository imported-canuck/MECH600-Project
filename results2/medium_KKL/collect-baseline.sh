#!/bin/bash -l
set -eo pipefail
run_root=$1
case_path="$run_root/cases/low_kkl_016"
module load openfoam/v2312 >/dev/null 2>&1
iteration=$(python3 -c 'import json,sys; s=json.load(open(sys.argv[1])); assert s["status"]=="finished", s; print(s["latest_iteration"])' "$case_path/run-status.json")
report="$run_root/report-$iteration"
mkdir -p "$report"
reconstructPar -case "$case_path" -newTimes > "$run_root/logs/reconstruct-$iteration.log" 2>&1
python3 "$run_root/code/analysis/analyze.py" "$case_path" --out "$report"
python3 "$run_root/plot-step4.py" "$run_root" "$report"
python3 - "$report" <<'PY'
import csv,sys
from pathlib import Path
report=Path(sys.argv[1])
with (report/'summary.csv').open() as stream:
    summary=next(csv.DictReader(stream))
for key in ('final_iteration','cl','cd','cl_range','cd_range','force_plateau','last_max_initial_residual','normalized_flux_imbalance','x_sep_over_c','x_reatt_over_c','y_plus_max_upper','y_plus_max_lower','probe_tu_pct','wall_seconds','core_hours'):
    print(f'{key}={summary.get(key)}')
print((report/'reattachment-history.csv').read_text().strip())
PY

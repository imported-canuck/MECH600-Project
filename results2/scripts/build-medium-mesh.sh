#!/bin/bash -l
# Build and check the one shared mesh used by both medium-grid simulations.
set -eo pipefail
umask 077

if [[ $# -ne 2 ]]; then
    echo "Usage: bash build-medium-mesh.sh RUN_ROOT BASELINE_ROOT" >&2
    exit 2
fi
run=$(realpath "$1")
baseline=$(realpath "$2")
if [[ -e "$run/mesh" ]]; then
    echo "A checked mesh already exists at $run/mesh" >&2
    exit 2
fi

module load openfoam/v2312
export FOAM_FILEHANDLER=uncollated

# Keep a failed mesh and its logs available without publishing it for solving.
mesh=$(mktemp -d "$run/mesh-build.XXXXXX")
printf 'Mesh work directory: %s\n' "$mesh"
mkdir -p "$mesh/system" "$mesh/constant"
cp -- "$baseline/cases/low_sst_016/system/"{controlDict,fvSchemes,fvSolution} "$mesh/system/"
python3 "$run/generate_blockmesh.py" --output "$mesh/system/blockMeshDict" \
    --metadata "$mesh/mesh.json"
blockMesh -case "$mesh" > "$mesh/log.blockMesh" 2>&1

python3 - "$baseline/code" "$mesh" <<'PY'
import json
import os
from pathlib import Path
import re
import sys

sys.path.insert(0, sys.argv[1])
import study

mesh = Path(sys.argv[2])
audit = study.check_imported_mesh(mesh, os.environ.copy())
report = (mesh / "log.checkMesh").read_text()
measured_cells = int(re.search(r"(?m)^\s*cells:\s*(\d+)", report)[1])
spec = json.loads((mesh / "mesh.json").read_text())
if measured_cells != spec["cells"] or measured_cells != 209_952:
    raise RuntimeError(f"Unexpected medium grid cell count: {measured_cells}")

# blockMesh already assigns the patch types. Inspect them without editing them.
boundary = (mesh / "constant/polyMesh/boundary").read_text()
expected = {"airfoil": "wall", "farfield": "patch", "frontAndBack": "empty"}
for name, kind in expected.items():
    block = re.search(rf"\b{name}\s*\{{([^}}]+)\}}", boundary)
    if not block or not re.search(rf"\btype\s+{kind}\s*;", block[1]):
        raise RuntimeError(f"Expected {kind} patch {name}")
study.save_json(mesh / "mesh-status.json", {
    "mesh_checked": True,
    "mesh_geometry_audit": audit,
})
print(f"Shared medium mesh: {measured_cells:,} cells. Geometry audit: {audit}.")
PY

mv -T -- "$mesh" "$run/mesh"
echo "Checked medium mesh available at $run/mesh"

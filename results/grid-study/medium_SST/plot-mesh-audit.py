import runpy
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(sys.argv[1])
sys.path.insert(0, str(root / "code" / "analysis"))
from surface import foam_list

case = root / "cases" / "low_sst_016"
mesh = case / "constant" / "polyMesh"

def read_array(path, dtype):
    count, text = foam_list(path)
    return count, np.fromstring(text.replace("(", " ").replace(")", " "), sep=" ", dtype=dtype)

point_count, points = read_array(mesh / "points", float)
points = points.reshape(point_count, 3)
face_count, faces = read_array(mesh / "faces", np.int64)
faces = faces.reshape(face_count, 5)[:, 1:]
_, owners = read_array(mesh / "owner", np.int64)
flag_count, flags = read_array(mesh / "sets" / "underdeterminedCells", np.int64)
assert len(flags) == flag_count
front = np.all(np.abs(points[faces, 2]) < 1e-10, axis=1)
front_centers = points[faces[front], :2].mean(axis=1) / 0.2
locations = front_centers[np.isin(owners[front], flags)]
assert len(locations) == flag_count
shape = runpy.run_path(str(root / "code" / "mesh" / "generate.py"))["thickness"]
x = np.linspace(0, 1, 500)
y = np.array([shape(v) for v in x])

fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout="constrained")
fig.suptitle(f"Medium C-grid | OpenFOAM found {flag_count:,} cells with determinant < 0.001", fontweight="bold")
for ax, title, xlim, ylim in ((axes[0], "Full domain", (-26, 32), (-27, 33)),
                               (axes[1], "Leading-edge region", (-0.02, 0.22), (-0.085, 0.085))):
    ax.scatter(*locations.T, s=3, linewidths=0, alpha=0.28, color="#ba4d15")
    ax.fill_between(x, -y, y, color="#172d3d")
    ax.set(xlabel="$x/c$", ylabel="$y/c$", xlim=xlim, ylim=ylim, title=title, aspect="equal")
    ax.grid(alpha=0.15)
fig.savefig(root / "mesh-audit.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print(f"Plotted {len(locations):,} checkMesh cell locations")

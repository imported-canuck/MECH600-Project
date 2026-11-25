import csv
import json
import re
import sys
from pathlib import Path
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np

root = Path(sys.argv[1])
case = root / "cases" / "low_kkl_016"
report = Path(sys.argv[2])
sys.path.insert(0, str(root / "code" / "analysis"))
import surface


def rows(name):
    with (report / name).open(newline="") as source:
        return list(csv.DictReader(source))


def numbers(records, column):
    return np.array([float(row[column]) if row.get(column) else np.nan for row in records])


config = json.loads((case / "case.json").read_text())
status = json.loads((case / "run-status.json").read_text())
checkpoint = int(status["latest_iteration"])
summary = rows("summary.csv")[0]
surface_rows = rows("low_kkl_016_surface.csv")
force = rows("low_kkl_016_force_history.csv")
residual = rows("low_kkl_016_residuals.csv")
label = f"NACA 0012 | wall-resolved k-kl-omega | iteration {checkpoint:,}"

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.18, "savefig.dpi": 170})

# The bundled surface analyzer normally reads only the latest iteration.
# Reuse it for each reconstructed checkpoint to track the same wall-shear crossing.
checkpoints = []
for folder in sorted((p for p in case.iterdir() if p.is_dir() and p.name.isdigit()),
                     key=lambda p: int(p.name)):
    if int(folder.name) > checkpoint:
        continue
    if not (folder / "p").is_file() or surface.find_field(folder, "wallShearStress") is None:
        continue
    with patch.object(surface, "latest_time", return_value=folder):
        data, _ = surface.surface_rows(case, config)
    metrics = surface.surface_metrics(data, config["wall_treatment"])
    checkpoints.append({"iteration": int(folder.name),
                        "x_sep_over_c": metrics.get("x_sep_over_c"),
                        "x_reatt_over_c": metrics.get("x_reatt_over_c")})
with (report / "reattachment-history.csv").open("w", newline="") as target:
    writer = csv.DictWriter(target, fieldnames=["iteration", "x_sep_over_c", "x_reatt_over_c"])
    writer.writeheader()
    writer.writerows(checkpoints)

figure, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
figure.suptitle(label + "\nConvergence history", fontweight="bold")
ax = axes[0, 0]
for field in ("p", "Ux", "Uy", "kt", "kl", "omega"):
    values = numbers(residual, field)
    ok = np.isfinite(values) & (values > 0)
    ax.semilogy(numbers(residual, "time")[ok], values[ok], label=field, lw=0.9)
ax.legend(ncol=3, fontsize=8, frameon=False)
ax.set(xlabel="SIMPLE iteration", ylabel="Initial residual", title="Equation residuals")
for ax, field, title, color in ((axes[0, 1], "cl", "Lift coefficient", "#076595"),
                                (axes[1, 0], "cd", "Drag coefficient", "#c35d17")):
    final_forces = [row for row in force if float(row["time"]) > config["warmup"]]
    ax.plot(numbers(final_forces, "time"), numbers(final_forces, field), lw=1.2, color=color)
    ax.axvspan(max(config["warmup"], checkpoint - 400), checkpoint,
               color="#94a3b8", alpha=0.22, label="Final 400 iterations")
    ax.set(xlabel="SIMPLE iteration", ylabel="$C_L$" if field == "cl" else "$C_D$", title=title)
    ax.legend(frameon=False, fontsize=8)
ax = axes[1, 1]
valid = [row for row in checkpoints if row["x_reatt_over_c"] is not None]
if valid:
    ax.plot(numbers(valid, "iteration"), numbers(valid, "x_reatt_over_c"), "o-", color="#167461", lw=1.2)
else:
    ax.text(0.5, 0.5, "No closed upper-surface bubble found", transform=ax.transAxes, ha="center")
ax.set(xlabel="Saved iteration", ylabel="$x_r/c$", title="Upper-surface reattachment")
figure.savefig(report / "convergence.png", bbox_inches="tight")
plt.close(figure)

figure, axes = plt.subplots(3, 1, figsize=(9, 8), layout="constrained")
figure.suptitle(label + "\nSurface diagnostics", fontweight="bold")
for side, color in (("upper", "#076595"), ("lower", "#c35d17")):
    points = [row for row in surface_rows if row["side"] == side]
    for ax, field, ylabel in zip(axes, ("cp", "cf", "y_plus"), ("$C_p$", "$C_f$", "$y^+$")):
        ax.plot(numbers(points, "x_over_c"), numbers(points, field), color=color,
                label=side.capitalize(), lw=1.2)
        ax.set(xlabel="$x/c$", ylabel=ylabel)
axes[0].set(xlim=(0, 1), title="Surface pressure")
axes[0].invert_yaxis()
axes[0].legend(frameon=False)
axes[1].set(xlim=(0, 0.18), title="Skin friction near the leading edge")
axes[1].axhline(0, color="0.4", lw=0.8)
for key, text, color in (("x_sep_over_c", "Separation", "#93408d"),
                         ("x_reatt_over_c", "Reattachment", "#167461")):
    if summary.get(key):
        axes[1].axvline(float(summary[key]), color=color, ls="--", lw=1, label=text)
axes[1].legend(frameon=False, ncol=4, fontsize=8)
axes[2].set(xlim=(0, 1), title="Wall distance in viscous units")
figure.savefig(report / "surface-diagnostics.png", bbox_inches="tight")
plt.close(figure)

# Each hexahedron owns exactly one front-plane quad. Color that quad with
# the solved cell velocity so the plot preserves the airfoil and the mesh.
def foam_array(path, dtype):
    count, content = surface.foam_list(path)
    cleaned = content.replace("(", " ").replace(")", " ")
    return count, np.fromstring(cleaned, sep=" ", dtype=dtype)

mesh = case / "constant" / "polyMesh"
n, coordinates = foam_array(mesh / "points", float)
coordinates = coordinates.reshape(n, 3)
n, raw_faces = foam_array(mesh / "faces", np.int64)
raw_faces = raw_faces.reshape(n, 5)
assert np.all(raw_faces[:, 0] == 4), "Expected quadrilateral faces"
faces = raw_faces[:, 1:]
_, owners = foam_array(mesh / "owner", np.int64)
front_faces = np.all(np.abs(coordinates[faces, 2]) < 1e-10, axis=1)
front_owners = owners[front_faces]
quads = coordinates[faces[front_faces], :2] / config["chord"]

velocity_text = (case / str(checkpoint) / "U").read_text()
match = re.search(r"\binternalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;",
                  velocity_text, re.S)
assert match is not None, "Expected a nonuniform internal velocity field"
count, vectors = match.groups()
velocity = np.fromstring(vectors.replace("(", " ").replace(")", " "), sep=" ").reshape(int(count), 3)
assert len(front_owners) == len(velocity) == len(np.unique(front_owners))
front_u = velocity[front_owners] / config["speed"]

figure, axes = plt.subplots(2, 1, figsize=(10, 7), layout="constrained")
figure.suptitle(label + "\nSolved velocity field", fontweight="bold")
views = (
    (axes[0], np.linalg.norm(front_u, axis=1), (-0.18, 1.65), (-0.20, 0.32),
     "viridis", Normalize(0, 1.8), "Speed around the airfoil and wake", "$|U|/U_\u221e$"),
    (axes[1], front_u[:, 0], (-0.008, 0.165), (0, 0.105),
     "RdBu", TwoSlopeNorm(vmin=-0.15, vcenter=0, vmax=1.6),
     "Chordwise velocity near the leading edge", "$U_x/U_\u221e$"),
)
low = quads.min(axis=1)
high = quads.max(axis=1)
for ax, values, xlim, ylim, palette, scale, title, bar_label in views:
    visible = ((high[:, 0] >= xlim[0]) & (low[:, 0] <= xlim[1]) &
               (high[:, 1] >= ylim[0]) & (low[:, 1] <= ylim[1]))
    colors = PolyCollection(quads[visible], array=values[visible], cmap=palette, norm=scale,
                            edgecolors="face", linewidths=0.18)
    ax.add_collection(colors)
    ax.set(xlim=xlim, ylim=ylim, xlabel="$x/c$", ylabel="$y/c$", title=title, aspect="equal")
    ax.grid(False)
    figure.colorbar(colors, ax=ax, label=bar_label, shrink=0.85, extend="both")
for key, color in (("x_sep_over_c", "#93408d"), ("x_reatt_over_c", "#167461")):
    if summary.get(key):
        axes[1].axvline(float(summary[key]), color=color, ls="--", lw=1.2)
figure.savefig(report / "flow-field.png", bbox_inches="tight")
plt.close(figure)
print(f"Plotted {checkpoint} with {len(checkpoints)} reconstructed checkpoints.")

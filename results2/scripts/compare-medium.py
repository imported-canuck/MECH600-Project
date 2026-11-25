"""Compare the completed medium SST, SST-LM and k-kl-omega cases at iteration 7,000."""

import argparse
import csv
from dataclasses import dataclass, field
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODELS = {
    "SST": ("kOmegaSST", "#ac492a", "o", -0.17),
    "SST-LM": ("kOmegaSSTLM", "#176578", "s", 0.0),
    "KKL": ("kkLOmega", "#5b4a8e", "^", 0.17),
}
LABELS = {"SST": "SST", "SST-LM": "SST-LM", "KKL": "KKL (upwind)"}
CHECKPOINTS = (6500, 6750, 7000)
ENDPOINT = CHECKPOINTS[-1]
HATA_SEPARATION = (0.011, 0.025)
HATA_REATTACHMENT = 0.10
CHORD_MM = 200.0
TARGET_SPAN_MM = 0.4


def read_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    value = row.get(key)
    if value is None or not value.strip():
        return None
    result = float(value)
    return result if math.isfinite(result) else None


@dataclass
class Report:
    summary: dict
    separation: float | None
    reattachment: float | None
    late_range: tuple | None
    profile: list = field(default_factory=list)

    def row(self, layout, model):
        lo, hi = self.late_range or (None, None)
        span = (hi - lo) * CHORD_MM if lo is not None else None
        error = None if self.reattachment is None else (
            self.reattachment - HATA_REATTACHMENT
        ) * CHORD_MM
        return {
            "mesh_layout": layout,
            "model": model,
            "iterations": ENDPOINT,
            "turbulence_convection": "upwind" if model == "KKL" else "limitedLinear",
            "turbulence_relaxation": 0.15 if model == "KKL" else 0.35,
            "cells": number(self.summary, "cells"),
            "separation_x_over_c": self.separation,
            "reattachment_x_over_c": self.reattachment,
            "hata_separation_min_over_c": HATA_SEPARATION[0],
            "hata_separation_max_over_c": HATA_SEPARATION[1],
            "hata_reattachment_approx_over_c": HATA_REATTACHMENT,
            "reattachment_minus_hata_mm": error,
            "absolute_reattachment_error_mm": abs(error) if error is not None else None,
            "reattachment_late_min_over_c": lo,
            "reattachment_late_max_over_c": hi,
            "reattachment_late_span_mm": span,
            "within_0p4mm_span": span <= TARGET_SPAN_MM if span is not None else None,
            "force_plateau": self.summary.get("force_plateau"),
            **{key: number(self.summary, key) for key in (
                "cl", "cd", "cl_range", "cd_range", "y_plus_max_upper", "y_plus_max_lower",
                "nprocs", "wall_seconds", "core_hours", "last_max_initial_residual",
            )},
        }


def load_report(root, model):
    folder = root / f"medium_{model}" / f"report-{ENDPOINT}"
    path = folder / "summary.csv"
    if not path.is_file():
        return None
    rows = read_rows(path)
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one case in {path}")
    summary = rows[0]
    if summary.get("model") != MODELS[model][0] or summary.get("level") != "medium":
        raise ValueError(f"Unexpected model or grid in {path}")
    if (summary.get("run_status") != "finished" or number(summary, "exit_code") != 0
            or number(summary, "final_iteration") != ENDPOINT
            or number(summary, "surface_time") not in (None, ENDPOINT)):
        raise ValueError(f"Expected a completed {ENDPOINT}-iteration case in {path}")

    xr = number(summary, "x_reatt_over_c")
    late_range = None
    history = folder / "reattachment-history.csv"
    if history.is_file() and xr is not None:
        positions = {number(row, "iteration"): number(row, "x_reatt_over_c")
                     for row in read_rows(history)}
        late = [positions.get(step) for step in CHECKPOINTS]
        if late[-1] is not None and not math.isclose(late[-1], xr, abs_tol=1e-8):
            raise ValueError(f"The last checkpoint disagrees with {path}")
        if all(position is not None for position in late):
            late_range = (min(late), max(late))

    profile = []
    surface = folder / f"{summary['id']}_surface.csv"
    if surface.is_file():
        for row in read_rows(surface):
            if row.get("side") == "upper" and number(row, "time") == ENDPOINT:
                x, cf = number(row, "x_over_c"), number(row, "cf")
                if x is not None and cf is not None:
                    profile.append((x, cf))
    return Report(summary, number(summary, "x_sep_over_c"), xr, late_range, sorted(profile))


def save_figure(data, path):
    figure, (original, block, endpoints) = plt.subplots(
        1, 3, figsize=(13, 4.3), layout="constrained", width_ratios=(1, 1, 1.05)
    )
    figure.suptitle("NACA 0012 medium mesh comparison  |  Iteration 7,000")
    original.sharey(block)
    for axis, (layout, reports) in zip((original, block), data.items()):
        axis.axhline(0, color="#686868", linewidth=0.85)
        axis.axvspan(*HATA_SEPARATION, color="#737373", alpha=0.15)
        axis.axvline(HATA_REATTACHMENT, color="#535353", linestyle=":")
        for model, report in reports.items():
            if report is not None and report.profile:
                x, cf = zip(*report.profile)
                axis.plot(x, cf, color=MODELS[model][1], lw=1.6, label=LABELS[model])
        axis.set(xlim=(0, 0.18), xlabel="Chordwise position, $x/c$", title=layout)
        if axis.lines[2:]:
            axis.legend(loc="upper right", frameon=False, fontsize=9)
    original.set_ylabel("Upper-surface skin friction, $C_f$")
    block.tick_params(labelleft=False)

    endpoints.axhline(HATA_REATTACHMENT, color="#535353", linestyle=":", label="Hata ≈ 0.10")
    for model, (_, color, marker, offset) in MODELS.items():
        points = []
        for index, reports in enumerate(data.values()):
            report = reports[model]
            if report is None or report.reattachment is None:
                continue
            x = index + offset
            points.append((x, report.reattachment))
            if report.late_range is not None:
                lo, hi = report.late_range
                endpoints.vlines(x, lo, hi, color=color, lw=1.5)
                endpoints.hlines([lo, hi], x - 0.025, x + 0.025, color=color, lw=1.0)
        if points:
            endpoints.plot(*zip(*points), marker=marker, linestyle="none", color=color, label=LABELS[model])
    endpoints.set(xlim=(-0.35, 1.35), xticks=range(2), xticklabels=list(data),
                  ylabel="Reattachment position, $x_r/c$", title="First reattachment and late ranges")
    endpoints.legend(frameon=False, fontsize=8)
    for axis in (original, block, endpoints):
        axis.set_axisbelow(True)
        axis.grid(alpha=0.22)
        axis.spines[["top", "right"]].set_visible(False)
    figure.savefig(path, dpi=190)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path, help="Original results/grid-study folder")
    parser.add_argument("rerun", type=Path, help="New results2 folder")
    parser.add_argument("output", type=Path, help="Comparison output folder")
    args = parser.parse_args()

    data = {
        layout: {model: load_report(root, model) for model in MODELS}
        for layout, root in (("Original C-grid", args.original), ("12-block C-grid", args.rerun))
    }
    rows = []
    reference = None
    for layout, reports in data.items():
        for model, report in reports.items():
            if report is None:
                continue
            conditions = tuple(number(report.summary, key) for key in (
                "cells", "reynolds", "alpha", "radius", "tu", "first_cell_height_over_c",
            )) + (report.summary.get("wall_treatment"),)
            if reference is None:
                reference = conditions
            if conditions != reference:
                raise ValueError(f"Conditions differ from the original case in {layout} {model}")
            rows.append(report.row(layout, model))
    if not rows:
        parser.error("No completed medium reports were found")

    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "medium-mesh-comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    save_figure(data, args.output / "medium-mesh-comparison.png")


if __name__ == "__main__":
    main()

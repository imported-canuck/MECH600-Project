"""Compare matching 7,000-iteration SST and SST-LM airfoil reports."""

import argparse
import csv
from dataclasses import dataclass, field
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


LEVELS = ("coarse", "medium", "fine")
MODELS = ("SST", "SST-LM")
OPENFOAM_MODELS = {"SST": "kOmegaSST", "SST-LM": "kOmegaSSTLM"}
COLORS = {"SST": "#ac492a", "SST-LM": "#176578"}
MARKERS = {"SST": "o", "SST-LM": "s"}
CHECKPOINTS = (6500, 6750, 7000)
ENDPOINT = CHECKPOINTS[-1]
TARGET_SPAN_OVER_C = 0.002
CHORD_MM = 200
HATA_SEPARATION = (0.011, 0.025)
HATA_REATTACHMENT = 0.10


def read_rows(path):
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    value = row.get(key)
    if value is None or not value.strip():
        return None
    value = float(value)
    return value if math.isfinite(value) else None


@dataclass
class Report:
    summary: dict = field(default_factory=dict)
    separation: float | None = None
    reattachment: float | None = None
    late_range: tuple | None = None
    profile: list = field(default_factory=list)

    @property
    def status(self):
        if not self.summary:
            return "Report unavailable"
        if self.reattachment is None:
            return "Reattachment unavailable"
        if self.late_range is None:
            return "Insufficient late checkpoints"
        lo, hi = self.late_range
        return "Within late-span target" if hi - lo <= TARGET_SPAN_OVER_C else "Above late-span target"


def load_report(root, level, model):
    folder = root / f"{level}_{model}" / f"report-{ENDPOINT}"
    rows = read_rows(folder / "summary.csv")
    if not rows:
        return Report()
    if len(rows) != 1:
        raise ValueError(f"Expected one case in {folder / 'summary.csv'}")
    summary = rows[0]
    if summary.get("model") != OPENFOAM_MODELS[model] or summary.get("level") != level:
        raise ValueError(f"Model or mesh level does not match {folder}")
    if (summary.get("run_status") != "finished" or number(summary, "exit_code") != 0
            or number(summary, "final_iteration") != ENDPOINT
            or number(summary, "surface_time") not in (None, ENDPOINT)):
        raise ValueError(f"Expected a completed {ENDPOINT}-iteration report in {folder}")

    peaks = {number(row, "iteration"): number(row, "x_reatt_over_c")
             for row in read_rows(folder / "reattachment-history.csv")}
    late = [peaks.get(step) for step in CHECKPOINTS]
    reattachment = number(summary, "x_reatt_over_c")
    late_range = None
    if reattachment is not None and all(position is not None for position in late):
        if not math.isclose(late[-1], reattachment, abs_tol=1e-8):
            raise ValueError(f"Final reattachment differs between the summary and history in {folder}")
        late_range = (min(late), max(late))

    profile = []
    if level == "medium":
        for row in read_rows(folder / f"{summary['id']}_surface.csv"):
            if row.get("side") == "upper" and number(row, "time") == ENDPOINT:
                position, cf = number(row, "x_over_c"), number(row, "cf")
                if position is not None and cf is not None:
                    profile.append((position, cf))
        profile.sort()
    return Report(summary, number(summary, "x_sep_over_c"), reattachment, late_range, profile)


def difference(a, b):
    return a - b if a is not None and b is not None else None


def comparison_row(level, models):
    sst, lm = (models[model] for model in MODELS)
    if sst.summary and lm.summary:
        conditions = ("cells", "reynolds", "alpha", "radius", "tu", "first_cell_height_over_c")
        for key in conditions:
            if number(sst.summary, key) != number(lm.summary, key):
                raise ValueError(f"{key} differs in the matched {level} reports")
        if sst.summary.get("wall_treatment") != lm.summary.get("wall_treatment"):
            raise ValueError(f"Wall treatment differs in the matched {level} reports")

    cells = number(sst.summary or lm.summary, "cells")
    row = {
        "grid": level,
        "cells": int(cells) if cells is not None else None,
        "iteration": ENDPOINT,
        "hata_separation_min_over_c": HATA_SEPARATION[0],
        "hata_separation_max_over_c": HATA_SEPARATION[1],
        "hata_approx_reattachment_over_c": HATA_REATTACHMENT,
    }
    for model, prefix in (("SST", "sst"), ("SST-LM", "sst_lm")):
        report = models[model]
        lo, hi = report.late_range or (None, None)
        span = difference(hi, lo)
        row.update({
            f"{prefix}_separation_over_c": report.separation,
            f"{prefix}_reattachment_over_c": report.reattachment,
            f"{prefix}_reattachment_minus_hata_over_c": difference(report.reattachment, HATA_REATTACHMENT),
            f"{prefix}_late_min_over_c": lo,
            f"{prefix}_late_max_over_c": hi,
            f"{prefix}_late_span_mm": span * CHORD_MM if span is not None else None,
            f"{prefix}_reattachment_status": report.status,
        })
        for key in ("y_plus_max_upper", "y_plus_max_lower", "cl", "cd", "core_hours"):
            row[f"{prefix}_{key}"] = number(report.summary, key)
    row["sst_lm_minus_sst_separation_over_c"] = difference(lm.separation, sst.separation)
    row["sst_lm_minus_sst_reattachment_over_c"] = difference(lm.reattachment, sst.reattachment)
    return row


def save_plot(data, path):
    fig, (wall, convergence) = plt.subplots(1, 2, figsize=(10, 4.3), layout="constrained")
    fig.suptitle("NACA 0012  |  SST and SST-LM at iteration 7,000")
    wall.axhline(0, color="#666666", linewidth=0.85)
    wall.axvspan(*HATA_SEPARATION, color="#737373", alpha=0.13, label="Hata separation interval")
    wall.axvline(HATA_REATTACHMENT, color="#535353", linestyle=":", label="Hata reattachment")
    for model, report in data["medium"].items():
        if report.profile:
            x, cf = zip(*report.profile)
            wall.plot(x, cf, color=COLORS[model], linewidth=1.7, label=model)
    wall.set(xlim=(0, 0.15), xlabel="Chordwise location, $x/c$", ylabel="Skin-friction coefficient, $C_f$",
             title="Upper surface on the medium grid")
    wall.legend(loc="upper right", framealpha=0.95, edgecolor="none", fontsize=8)

    convergence.axhline(HATA_REATTACHMENT, color="#535353", linestyle=":", label="Hata $x_r/c$ ≈ 0.10")
    for model, offset in (("SST", -0.1), ("SST-LM", 0.1)):
        color = COLORS[model]
        plotted = []
        for index, level in enumerate(LEVELS):
            report = data[level][model]
            if report.reattachment is None:
                continue
            x = index + offset
            plotted.append((x, report.reattachment))
            if report.late_range is not None:
                lo, hi = report.late_range
                convergence.vlines(x, lo, hi, color=color, linewidth=1.4)
                convergence.hlines([lo, hi], x - 0.035, x + 0.035, color=color, linewidth=1.1)
        if plotted:
            x, xr = zip(*plotted)
            convergence.plot(x, xr, color=color, marker=MARKERS[model], linestyle="none", label=model)
    convergence.set(xlim=(-0.4, 2.4), xticks=range(3), xticklabels=[level.title() for level in LEVELS],
                    ylabel="Reattachment, $x_r/c$", title="Reattachment and late checkpoint range")
    convergence.legend(frameon=False, fontsize=8)
    for axis in (wall, convergence):
        axis.grid(alpha=0.25)
        axis.set_axisbelow(True)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study_root", type=Path, help="Folder containing the named grid and model results")
    parser.add_argument("output_dir", type=Path, help="Destination for the comparison CSV and PNG")
    args = parser.parse_args()

    data = {level: {model: load_report(args.study_root, level, model) for model in MODELS} for level in LEVELS}
    if not any(report.summary for models in data.values() for report in models.values()):
        parser.error("No completed model reports were found")
    rows = [comparison_row(level, data[level]) for level in LEVELS]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "model-comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    save_plot(data, args.output_dir / "model-comparison.png")


if __name__ == "__main__":
    main()

"""Compare the upper-surface separation bubble on three MECH 600 grids."""

import argparse
import csv
import math
from pathlib import Path
from statistics import geometric_mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


LEVELS = ("coarse", "medium", "fine")
COLORS = dict(zip(LEVELS, ("#b45625", "#246c96", "#25866c")))
CONDITIONS = ("model", "wall_treatment", "reynolds", "alpha", "radius", "tu")
HATA_SEPARATION = (0.011, 0.025)
HATA_REATTACHMENT = 0.10
ITERATION_TOLERANCE = 0.002  # Fraction of chord, or 0.4 mm for this case.
DRIFT_DIAGNOSTIC = 0.10  # Conservative comparison to a mesh-to-mesh shift.


def read_csv(path):
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def number(row, column):
    value = row.get(column)
    if value is None or not value.strip():
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def read_report(root, level, iteration):
    run = root / f"{level}_SST"
    if not run.is_dir():
        run = root / level
    report = run / f"report-{iteration}"
    summaries = read_csv(report / "summary.csv")
    if len(summaries) > 1:
        raise ValueError(f"Expected one case in {report / 'summary.csv'}")
    summary = summaries[0] if summaries else {}
    case_id = summary.get("id", "low_sst_016")
    time = number(summary, "surface_time")
    surface = []
    for row in read_csv(report / f"{case_id}_surface.csv"):
        if row["side"] == "upper" and number(row, "time") == time:
            x = number(row, "x_over_c")
            cf = number(row, "cf")
            if x is not None and cf is not None:
                surface.append((x, cf))
    surface.sort()

    history = sorted(
        (int(float(row["iteration"])), number(row, "x_reatt_over_c"))
        for row in read_csv(report / "reattachment-history.csv")
    )
    late = [(step, x) for step, x in history if step <= iteration][-3:]
    late_positions = []
    if len(late) == 3 and late[-1][0] == iteration and all(x is not None for _, x in late):
        late_positions = [x for _, x in late]

    xs = number(summary, "x_sep_over_c")
    xr = number(summary, "x_reatt_over_c")
    if xs is not None and xr is not None and not 0 <= xs < xr <= 1:
        raise ValueError(f"Reattachment must be downstream of separation in {report}")
    if xs is not None and xr is not None:
        bubble_status = "closed bubble"
    elif xs is not None:
        bubble_status = "separation without detected reattachment"
    elif not summary:
        bubble_status = "report unavailable"
    elif not surface:
        bubble_status = "wall-shear data unavailable"
    elif min(cf for _, cf in surface) >= 0:
        bubble_status = "forward shear across sampled upper surface"
    else:
        bubble_status = "no qualifying leading-edge separation bubble"

    cells = number(summary, "cells")
    area = number(summary, "fluid_area")
    h = math.sqrt(area / cells) if area and area > 0 and cells and cells > 0 else None
    y_plus = [number(summary, f"y_plus_max_{side}") for side in ("upper", "lower")]
    residuals = read_csv(report / f"{case_id}_residuals.csv")
    residual = next((row for row in reversed(residuals) if number(row, "time") == iteration), {})
    result = {
        "grid": level,
        "cells": int(cells) if cells is not None else None,
        "representative_h_m": h,
        "fluid_area_m2": area,
        "final_iteration": int(float(summary["final_iteration"])) if summary.get("final_iteration") else None,
        "x_sep_over_c": xs,
        "x_reatt_over_c": xr,
        "bubble_length_over_c": xr - xs if xs is not None and xr is not None else None,
        "bubble_status": bubble_status,
        "separation_intervals": summary.get("upper_separation_intervals"),
        "late_checkpoints": ", ".join(str(step) for step, _ in late),
        "late_reatt_min_over_c": min(late_positions, default=None),
        "late_reatt_max_over_c": max(late_positions, default=None),
        "late_reatt_span_over_c": max(late_positions) - min(late_positions) if late_positions else None,
        "max_y_plus": max((y for y in y_plus if y is not None), default=None),
        "cl_mean": number(summary, "cl"),
        "cl_range": number(summary, "cl_range"),
        "cd_mean": number(summary, "cd"),
        "cd_range": number(summary, "cd_range"),
        "force_plateau": summary.get("force_plateau"),
        "pressure_initial_residual": number(residual, "p"),
        "normalized_flux_imbalance": number(summary, "normalized_flux_imbalance"),
        "core_hours": number(summary, "core_hours"),
        "mesh_geometry_audit": summary.get("mesh_geometry_audit"),
        "run_status": summary.get("run_status"),
    }
    return result, summary, surface


def finer_ratio(finer, coarser):
    hf = finer["representative_h_m"]
    hc = coarser["representative_h_m"]
    return hc / hf if hf is not None and hc is not None else None


def reattachment_gci(records, summaries, iteration):
    coarse, medium, fine = records
    r21 = finer_ratio(fine, medium)
    r32 = finer_ratio(medium, coarse)
    r = geometric_mean((r21, r32)) if r21 and r32 else None
    result = dict.fromkeys(("gci_status", "common_r", "observed_order_p", "richardson_reatt_over_c",
                            "fine_gci_abs_over_c", "fine_gci_percent", "late_span_fraction_of_grid_shift"))
    x3, x2, x1 = [row["x_reatt_over_c"] for row in records]
    span3, span2, span1 = [row["late_reatt_span_over_c"] for row in records]

    if not all(summaries):
        status = "Need all three grid reports"
    elif any(s.get("run_status") != "finished" or number(s, "exit_code") != 0 for s in summaries):
        status = "Need three successful solves"
    elif any(row["final_iteration"] != iteration for row in records):
        status = "A report does not match the requested iteration"
    elif any(s.get("mesh_ok") != "True" for s in summaries):
        status = "A grid failed the basic mesh check"
    elif any(len({s.get(key) for s in summaries}) != 1 or not summaries[0].get(key) for key in CONDITIONS):
        status = "Grid cases have incomplete or different physical settings"
    elif any(s.get("level") != level for s, level in zip(summaries, LEVELS)):
        status = "Grid levels do not match the requested family"
    elif r is None or min(r21, r32) < 1.3 or abs(r21 / r32 - 1) > 0.01:
        status = "Need nearly equal refinement ratios of at least 1.3"
    elif any(row["bubble_status"] != "closed bubble" for row in records):
        status = "Need a closed separation bubble on each grid"
    elif any(span is None for span in (span3, span2, span1)):
        status = "Need three late reattachment samples on each grid"
    elif max(span3, span2, span1) > ITERATION_TOLERANCE:
        status = "A late reattachment span exceeds 0.002 chord"
    else:
        e21 = x2 - x1
        e32 = x3 - x2
        if e21 * e32 <= 0 or abs(e32) <= abs(e21):
            status = "Grid differences do not shrink monotonically"
        else:
            error_ratio = abs(e32 / e21)
            gci = 1.25 * abs(e21) / (error_ratio - 1)
            span_fraction = max(max(span2, span1) / abs(e21), max(span3, span2) / abs(e32))
            result.update(common_r=r, observed_order_p=math.log(error_ratio) / math.log(r),
                          richardson_reatt_over_c=x1 - e21 / (error_ratio - 1),
                          fine_gci_abs_over_c=gci, fine_gci_percent=100 * gci / abs(x1),
                          late_span_fraction_of_grid_shift=span_fraction)
            status = "Computed for reattachment"
            if span_fraction > DRIFT_DIAGNOSTIC:
                status = "Provisional: late drift exceeds the conservative 10% diagnostic"
    result["gci_status"] = status
    return result


def comparison_plot(output, records, surfaces, iteration):
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), layout="constrained")
    fig.suptitle(f"NACA 0012 · wall-resolved SST · iteration {iteration:,}", fontweight="bold")
    cf_axes, reattach_axes = axes
    cf_axes.axhline(0, color="0.4", linewidth=0.85)
    cf_axes.axvspan(*HATA_SEPARATION, color="0.4", alpha=0.12, label="Hata separation interval")
    cf_axes.axvline(HATA_REATTACHMENT, color="0.4", linestyle="--", linewidth=1,
                   label="Hata reattachment (approx.)")
    reattach_axes.axhline(HATA_REATTACHMENT, color="0.4", linestyle="--", linewidth=1,
                         label="Hata reattachment (approx.)")

    # Include the spacing of a completed grid even when its solution has no reattachment point.
    sizes_mm = [1000 * row["representative_h_m"] for row in records if row["representative_h_m"] is not None]
    if sizes_mm:
        left, right = min(sizes_mm), max(sizes_mm)
        margin = 0.16 * max(right - left, 0.15 * right)
        reattach_axes.set_xlim(left - margin, right + margin)

    for index, (row, surface) in enumerate(zip(records, surfaces)):
        level = row["grid"]
        color = COLORS[level]
        in_view = [(x, cf) for x, cf in surface if 0 <= x <= 0.15]
        if in_view:
            x, cf = zip(*in_view)
            cf_axes.plot(x, cf, color=color, linewidth=1.6, label=level.capitalize())
        xr = row["x_reatt_over_c"]
        h = row["representative_h_m"]
        if xr is None:
            if row["run_status"] == "finished" and row["final_iteration"] == iteration:
                message = f"{level.capitalize()}\nReattachment unavailable"
                if h is None:
                    reattach_axes.text(0.22 + 0.28 * index, 0.04, message, transform=reattach_axes.transAxes,
                                       ha="center", va="bottom", color=color, fontsize=8)
                else:
                    reattach_axes.text(1000 * h, 0.04, message, transform=reattach_axes.get_xaxis_transform(),
                                       ha="center", va="bottom", color=color, fontsize=8)
            continue
        cf_axes.plot(xr, 0, "o", color=color, markersize=4)
        if h is None:
            continue
        h_mm = 1000 * h
        reattach_axes.plot(h_mm, xr, "o", color=color, label=level.capitalize())
        low, high = row["late_reatt_min_over_c"], row["late_reatt_max_over_c"]
        if low is not None:
            reattach_axes.vlines(h_mm, low, high, color=color, linewidth=2)

    cf_axes.set(xlim=(0, 0.15), xlabel=r"Distance from leading edge, $x/c$",
                ylabel=r"Skin friction, $C_f$", title="Upper-surface skin friction near the leading edge")
    reattach_axes.set(xlabel=r"Representative cell width $\sqrt{A/N}$ (mm)",
                      ylabel=r"Reattachment location, $x_r/c$",
                      title="Reattachment by mesh (whiskers span the last three checkpoints)")
    for axis in axes:
        axis.grid(alpha=0.22)
        axis.legend(frameon=False, fontsize=9, ncol=2)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study_root", type=Path, help="Folder containing coarse, medium and fine result folders")
    parser.add_argument("output_dir", type=Path, help="Destination for the comparison CSV and PNG")
    parser.add_argument("--iteration", type=int, default=7000, help="Iteration shared by the grid reports")
    args = parser.parse_args()
    if not args.study_root.is_dir() or args.iteration < 1:
        parser.error("Provide an existing study root and a positive iteration")

    records, summaries, surfaces = zip(*(read_report(args.study_root, level, args.iteration) for level in LEVELS))
    if not any(summaries):
        parser.error(f"No report-{args.iteration} summaries were found under {args.study_root}")
    gci = reattachment_gci(records, summaries, args.iteration)
    for index, row in enumerate(records):
        ratio = finer_ratio(row, records[index - 1]) if index else None
        row["ratio_to_next_coarser_grid"] = ratio
        row.update(gci if index == 2 else dict.fromkeys(gci))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_csv = args.output_dir / "grid-comparison.csv"
    with output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    comparison_plot(args.output_dir / "grid-comparison.png", records, surfaces, args.iteration)
    print(f"Reattachment GCI: {gci['gci_status']}")
    print(f"Comparison files written to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()

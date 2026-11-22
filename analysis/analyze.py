#!/usr/bin/env python3
"""Summarize the NACA 0012 case matrix with Python's standard library."""

import argparse
from bisect import bisect_left
import csv
import json
import math
from pathlib import Path
import re
from statistics import fmean
import sys

from surface import NUMBER, VECTOR, surface_metrics, surface_rows


STUDIES = Path(__file__).resolve().parents[1] / "studies"


def json_file(path):
    return json.loads(path.read_text(encoding="utf-8-sig")) if path.is_file() else {}


def study_analysis():
    plans = {}
    for path in sorted(STUDIES.glob("*/study.json")):
        study = json_file(path)
        group = study["defaults"]["group"]
        if group in plans:
            raise ValueError(f"{path}: duplicate analysis group {group!r}.")
        plans[group] = study["analysis"]
    return plans


def csv_file(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        if fields:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def segment_files(case, object_name, names):
    paths = []
    parent = case / "postProcessing" / object_name
    if parent.is_dir():
        for start in parent.iterdir():
            if start.is_dir() and re.fullmatch(NUMBER, start.name):
                path = next((start / name for name in names if (start / name).is_file()), None)
                if path:
                    paths.append((float(start.name), path))
    return sorted(paths)


def merge_segments(segments, read_segment):
    history = {}
    for start, path in segments:
        rows = read_segment(path)
        if rows:
            # A restart supersedes the older tail, even if the rerun ends earlier.
            history = {time: row for time, row in history.items() if time < start}
            history.update((row["time"], row) for row in rows)
    return [history[time] for time in sorted(history)]


def coeff_segment(path):
    rows = []
    keys = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.lstrip().startswith("#"):
            entries = line.lstrip()[1:].strip().lower().split()
            if entries and entries[0] == "time" and "cl" in entries and "cd" in entries:
                keys = [entries.index(key) for key in ("time", "cl", "cd")]
            continue
        if not keys:
            continue
        parts = line.split()
        try:
            values = [float(parts[i]) for i in keys]
        except (ValueError, IndexError):
            continue
        if all(map(math.isfinite, values)):
            rows.append(dict(zip(("time", "cl", "cd"), values)))
    return rows


def force_segment(path):
    data = path.read_text(encoding="utf-8", errors="replace").splitlines()
    offset = None
    columns = ()
    component_names = [f"{part}_{axis}" for part in ("pressure", "viscous") for axis in "xyz"]
    for line in data:
        if line.lstrip().startswith("#"):
            headings = line.lstrip()[1:].strip().lower().split()
            if re.search(r"\bforces\s*\(pressure\s+viscous", line, re.I):
                offset = int(bool(re.search(r"\bCofR\b", line, re.I)))
            elif all(label in headings for label in component_names):
                columns = tuple(headings.index(label) for label in component_names)
            elif re.search(r"\btotal\s+pressure\s+viscous\b", line, re.I):
                offset = 1
        elif columns:
            parts = line.replace("(", "").replace(")", "").split()
            try:
                time = float(parts[0])
                components = tuple(float(parts[i]) for i in columns)
            except (ValueError, IndexError):
                continue
            yield {"time": time, "pressure": components[:3], "viscous": components[3:]}
        elif offset is not None:
            time = re.match(rf"\s*({NUMBER})", line)
            groups = VECTOR.findall(line)
            if time and len(groups) >= offset + 2:
                yield {"time": float(time[1]),
                       "pressure": tuple(map(float, groups[offset])),
                       "viscous": tuple(map(float, groups[offset + 1]))}


def probe_segment(path, vector):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        time = re.match(rf"\s*({NUMBER})\s+", line)
        if not time:
            continue
        values = VECTOR.findall(line[time.end():]) if vector else line[time.end():].split()
        if values:
            try:
                value = math.hypot(*map(float, values[0])) if vector else float(values[0])
            except ValueError:
                continue
            if math.isfinite(value) and abs(value) < 1e20:
                rows.append({"time": float(time[1]), "value": value})
    return rows


def upstream_metrics(case):
    fields = {}
    for field in ("U", "k", "omega"):
        segments = segment_files(case, "upstream", [field])
        history = merge_segments(segments, lambda path: probe_segment(path, field == "U"))
        fields[field] = {row["time"]: row["value"] for row in history}
    times = fields["U"].keys() & fields["k"].keys() & fields["omega"].keys()
    if not times:
        return {}
    latest = max(times)
    u, k, omega = (fields[field][latest] for field in ("U", "k", "omega"))
    if u <= 0 or k < 0:
        return {}
    return {"probe_time": latest, "probe_speed": u, "probe_k": k, "probe_omega": omega,
            "probe_tu_pct": 100 * math.sqrt(2 * k / 3) / u}


def flux_metrics(case):
    values = {}
    for name in ("netFlux", "absoluteFlux"):
        segments = segment_files(case, name, ["surfaceFieldValue.dat", f"{name}.dat"])
        history = merge_segments(segments, lambda path: probe_segment(path, vector=False))
        values[name] = {row["time"]: row["value"] for row in history}
    times = values["netFlux"].keys() & values["absoluteFlux"].keys()
    if not times:
        return {}
    time = max(times)
    net, total = (values[key][time] for key in ("netFlux", "absoluteFlux"))
    if total <= 0:
        return {}
    return {"flux_time": time, "net_boundary_flux": net,
            "abs_boundary_flux": total, "normalized_flux_imbalance": abs(net) / total}


def mesh_metrics(case):
    path = case / "log.checkMesh"
    if not path.is_file():
        return {}
    log = path.read_text(encoding="utf-8", errors="replace")
    fields = {
        "cells": r"(?m)^\s*cells:\s*(\d+)",
        "max_aspect": rf"Max aspect ratio\s*=\s*({NUMBER})",
        "max_nonorthogonality": rf"non-orthogonality\s+Max:\s*({NUMBER})",
        "max_skewness": rf"Max skewness\s*=\s*({NUMBER})",
    }
    basic_path = case / "log.checkMesh.basic"
    basic = basic_path.read_text(encoding="utf-8", errors="replace") if basic_path.is_file() else log
    measured = {"mesh_ok": "Mesh OK." in basic and "Failed " not in basic,
                "mesh_all_geometry_ok": "Mesh OK." in log and "Failed " not in log}
    for name, pattern in fields.items():
        match = re.search(pattern, log, re.I)
        if match:
            measured[name] = int(match[1]) if name == "cells" else float(match[1])
    return measured


def solver_log(path):
    if not path.is_file():
        return [], {}
    time = None
    residuals = {}
    clock = None
    ranks = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(rf"\s*Time\s*=\s*({NUMBER})\s*$", line)
        if match:
            new_time = float(match[1])
            if time is not None and new_time <= time:
                residuals = {t: data for t, data in residuals.items() if t < new_time}
            time = new_time
            residuals[time] = {}
        if time is not None:
            solved = re.search(rf"Solving for\s+([\w.]+),\s*Initial residual\s*=\s*({NUMBER})", line)
            if solved:
                row = residuals[time]
                row.setdefault(solved[1], float(solved[2]))
            continuity = re.search(rf"time step continuity errors.*?global\s*=\s*({NUMBER})", line)
            if continuity:
                residuals[time]["continuity_global"] = float(continuity[1])
        timing = re.search(rf"ClockTime\s*=\s*({NUMBER})", line)
        if timing:
            clock = float(timing[1])
        proc = re.match(r"\s*nProcs\s*:\s*(\d+)", line)
        if proc:
            ranks = int(proc[1])
    return [{"time": time, **residuals[time]} for time in sorted(residuals)], {
        "clock_seconds": clock, "ranks": ranks,
    }


def report_runs(case, meta, window, cl_tol, cd_tol, output):
    name = meta["id"]
    status = json_file(case / "run-status.json")
    row = {key: meta.get(key) for key in
           ("id", "group", "model", "level", "wall_treatment", "reynolds", "alpha", "radius", "tu")}
    row.update(mesh_metrics(case))
    generated = json_file(case / "mesh.json")
    row["generated_cells"] = generated.get("cells")
    for key in ("fluid_area", "first_cell_height_over_c", "nose_streamwise_spacing_over_c",
                "trailing_edge_streamwise_spacing_over_c"):
        row[key] = generated.get(key)
    row.update(upstream_metrics(case))
    row.update(flux_metrics(case))
    notes = []

    history = merge_segments(
        segment_files(case, "forceCoeffs", ["coefficient.dat", "forceCoeffs.dat"]), coeff_segment)
    residuals, log_run = solver_log(case / "log.solve")
    warm, log_warm = solver_log(case / "log.warmup")
    stage_start = meta.get("warmup", 0) or 0
    if warm:
        stage_start = max(stage_start, warm[-1]["time"])
    nprocs = status.get("nprocs") or log_run.get("ranks") or meta.get("nprocs")
    wall = status.get("wall_seconds")
    if wall is None and log_run.get("clock_seconds") is not None:
        wall = log_run["clock_seconds"] + (log_warm.get("clock_seconds") or 0)
    row.update({
        "run_status": status.get("status"),
        "mesh_geometry_audit": status.get("mesh_geometry_audit"),
        "init_wall_seconds": status.get("init_wall_seconds"),
        "exit_code": status.get("exit_code"),
        "nprocs": nprocs,
        "wall_seconds": wall,
        "core_hours": wall * nprocs / 3600 if wall is not None and nprocs else None,
        "logged_solve_steps": len(residuals),
        "logged_warmup_steps": len(warm),
    })
    if residuals:
        csv_file(output / f"{name}_residuals.csv", residuals)
        row["final_iteration"] = residuals[-1]["time"]
        row["last_log_continuity_global"] = residuals[-1].get("continuity_global")
        last_residuals = [value for key, value in residuals[-1].items()
                          if key not in ("time", "continuity_global")]
        row["last_max_initial_residual"] = max(last_residuals, default=None)

    actual_history = [entry for entry in history if entry["time"] > stage_start]
    if actual_history:
        csv_file(output / f"{name}_force_history.csv", history)
        final_time = actual_history[-1]["time"]
        threshold = final_time - window
        entries = [item for item in actual_history if item["time"] >= threshold]
        sizes = len(entries)
        coverage = final_time - entries[0]["time"]
        row.update({"window_samples": sizes, "window_span": coverage, "last_force_time": final_time})
        for key in ("cl", "cd"):
            data = [entry[key] for entry in entries]
            row[key] = fmean(data)
            row[f"{key}_range"] = max(data) - min(data)
        row["force_plateau"] = (
            actual_history[0]["time"] <= threshold
            and sizes >= 3 and row["cl_range"] <= cl_tol and row["cd_range"] <= cd_tol
        )
        if not row["force_plateau"]:
            notes.append("Force plateau criterion was not met.")
        forces = merge_segments(
            segment_files(case, "forces", ["force.dat", "forces.dat"]),
            lambda path: list(force_segment(path)),
        )
        force_window = [item for item in forces if item["time"] > stage_start and item["time"] >= threshold]
        if force_window:
            angle = math.radians(meta["alpha"])
            direction = (math.cos(angle), math.sin(angle))
            norm = 0.5 * meta.get("rho", 1.0) * meta["speed"]**2 * meta["chord"] * meta["span"]
            for component in ("pressure", "viscous"):
                row[f"cd_{component}"] = fmean(
                    sum(a * b for a, b in zip(direction, item[component])) / norm
                    for item in force_window
                )
    else:
        notes.append("No final-scheme lift and drag history found.")
    try:
        surface, messages = surface_rows(case, meta)
        if surface:
            csv_file(output / f"{name}_surface.csv", surface)
            row.update(surface_metrics(surface, meta["wall_treatment"]))
        notes.extend(messages)
    except (OSError, ValueError, TypeError, IndexError) as error:
        surface = []
        notes.append(f"Surface data skipped ({error}).")
    row["notes"] = " ".join(notes)
    return row, history, surface


def reference_rows(reference):
    alpha = reference["alpha"]
    result = []
    for grit, (a, b) in reference["points_by_grit"].items():
        fraction = (alpha - a[0]) / (b[0] - a[0])
        result.append({"grit": int(grit), "alpha": alpha,
                       "cl": a[1] + fraction * (b[1] - a[1]),
                       "cd": a[2] + fraction * (b[2] - a[2]),
                       "source": reference["source"]})
    return result


def write_references(output, plans):
    force_targets = {}
    for group, plan in plans.items():
        if force := plan.get("force_reference"):
            name = force["name"]
            rows = reference_rows(force)
            csv_file(output / f"{name}_interpolated.csv", rows)
            target = next(row for row in rows if row["grit"] == force["primary_grit"])
            force_targets[group] = (name, target)
        if bubble := plan.get("bubble_reference"):
            csv_file(output / f"{bubble['name']}_bubble_reference.csv", [bubble["measurements"]])
    return force_targets


def gci_table(by_name, case_ids):
    records = [by_name.get(name, {}) for name in case_ids]
    fine, medium, coarse = records
    cells = [item.get("cells") for item in records]
    areas = [item.get("fluid_area") for item in records]
    r21 = r32 = None
    if all(cells):
        areas = areas if all(areas) else [1, 1, 1]
        h1, h2, h3 = [math.sqrt(area / n) for area, n in zip(areas, cells)]
        r21, r32 = h2 / h1, h3 / h2
    result = []
    for key in ("cl", "cd"):
        values = [item.get(key) for item in records]
        row = {"quantity": key, "fine": values[0], "medium": values[1], "coarse": values[2],
               "n_fine": cells[0], "n_medium": cells[1], "n_coarse": cells[2],
               "r21": r21, "r32": r32, "status": "missing grids or coefficients"}
        result.append(row)
        if not all(cells) or any(x is None for x in values):
            continue
        e21, e32 = values[1] - values[0], values[2] - values[1]
        row["grid_value_span"] = max(values) - min(values)
        if any(record.get("exit_code") not in (None, 0) for record in records):
            row["status"] = "a baseline solve failed"
        elif any(record.get("run_status") not in (None, "finished") for record in records):
            row["status"] = "baseline solve still in progress or incomplete"
        elif any(record.get("mesh_ok") is False for record in records):
            row["status"] = "a baseline mesh did not pass checkMesh"
        elif len({item.get("radius") for item in records}) != 1:
            row["status"] = "grid family must share the same outer extent"
        elif not all(record.get("force_plateau") for record in records):
            row["status"] = "force plateau not established on all grids"
        elif min(r21, r32) < 1.3 or abs(r21 / r32 - 1) > 0.01:
            row["status"] = "refinement ratios differ or are below 1.3"
        elif e21 * e32 <= 0 or abs(e32) <= abs(e21):
            row["status"] = "nonmonotone or nonshrinking coefficients"
        elif all(record.get(f"{key}_range") is not None for record in records) and (
                max(fine[f"{key}_range"], medium[f"{key}_range"]) > 0.1 * abs(e21)
                or max(medium[f"{key}_range"], coarse[f"{key}_range"]) > 0.1 * abs(e32)):
            row["status"] = "iterative range exceeds 10% of a grid difference"
        elif abs(values[0]) < 1e-12:
            row["status"] = "fine coefficient too close to zero for a relative GCI"
        else:
            r = math.sqrt(r21 * r32)
            ratio = abs(e32 / e21)
            order = math.log(ratio) / math.log(r)
            absolute = 1.25 * abs(e21) / (ratio - 1)
            row.update({"p": order, "extrapolated": (ratio * values[0] - values[1]) / (ratio - 1),
                        "gci_abs": absolute, "gci_pct": 100 * absolute / abs(values[0]),
                        "status": "computed for this grid family"})
    return result


def contrasts(by_name, pairs):
    result = []
    for name, (a, b) in pairs.items():
        left, right = by_name.get(a, {}), by_name.get(b, {})
        row = {"comparison": name, "case_a": a, "case_b": b}
        for field in ("cl", "cd", "x_sep_over_c", "x_reatt_over_c", "probe_tu_pct", "probe_k", "probe_omega"):
            if left.get(field) is not None and right.get(field) is not None:
                row[f"delta_{field}"] = left[field] - right[field]
        if row.get("delta_cd") is not None:
            row["delta_cd_counts"] = row["delta_cd"] / 1e-4
        row["both_force_plateaus"] = bool(left.get("force_plateau") and right.get("force_plateau"))
        result.append(row)
    return result


def cp_errors(surface, reference, group):
    errors = []
    for side in ("upper", "lower"):
        points = [item for item in surface if item["side"] == side]
        coords = [item["x_over_c"] for item in points]
        for item in reference:
            if item["side"] != side or item["group"] != group:
                continue
            x = float(item["x_over_c"])
            i = bisect_left(coords, x)
            if i < len(coords) and coords[i] == x:
                value = points[i]["cp"]
            elif 0 < i < len(coords):
                fraction = (x - coords[i - 1]) / (coords[i] - coords[i - 1])
                value = points[i - 1]["cp"] + fraction * (points[i]["cp"] - points[i - 1]["cp"])
            else:
                continue
            reference_value = float(item["cp"])
            errors.append({"group": group, "side": side, "x_over_c": x, "cp_reference": reference_value,
                           "cp_cfd": value, "cp_error": value - reference_value,
                           "source": item.get("source", ""), "source_url": item.get("source_url", "")})
    return errors


def make_plots(output, histories, surfaces, by_name, reference, plans):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is unavailable. CSV files were written.", file=sys.stderr)
        return
    for group, plan in plans.items():
        plot_study(plt, output, histories, surfaces, by_name, reference, group, plan)


def plot_study(plt, output, histories, surfaces, by_name, reference, group, plan):
    if plan.get("plot_force_history"):
        figure, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
        for name, history in histories.items():
            if not history or by_name[name]["group"] != group:
                continue
            for axis, key in zip(axes, ("cl", "cd")):
                axis.plot([entry["time"] for entry in history], [entry[key] for entry in history], label=name)
                axis.set_ylabel(r"$C_L$" if key == "cl" else r"$C_D$")
        axes[-1].set_xlabel("Steady solver iteration")
        if axes[0].lines:
            axes[0].legend(fontsize=7, ncols=2)
            figure.savefig(output / f"{group}_force_history.png", dpi=180)
        plt.close(figure)

    case_ids = plan.get("gci_cases", [])
    grids = [by_name.get(name, {}) for name in case_ids]
    if grids and all(row.get("cells") and row.get("cl") is not None for row in grids):
        h_ratio = [math.sqrt(grids[0]["cells"] / row["cells"]) for row in grids]
        figure, axes = plt.subplots(1, 2, figsize=(8, 3.5), constrained_layout=True)
        for axis, result, label in zip(axes, gci_table(by_name, case_ids), (r"$C_L$", r"$C_D$")):
            key = result["quantity"]
            centers = [grid[key] for grid in grids]
            axis.plot(h_ratio, centers, "o-", label="Computed")
            if result.get("extrapolated") is not None:
                axis.scatter(0, result["extrapolated"], marker="x", label="Richardson")
                axis.errorbar(1, centers[0], yerr=result["gci_abs"], fmt="none", capsize=4,
                              label="Fine-grid GCI")
            axis.set(xlabel=r"$h/h_\mathrm{fine}$", ylabel=label)
            axis.ticklabel_format(axis="y", useOffset=False)
            axis.legend(fontsize=7)
        figure.savefig(output / f"{group}_grid_convergence.png", dpi=180)
        plt.close(figure)

    if wall_cases := plan.get("wall_resolution_cases"):
        figure, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
        for name in wall_cases:
            rows = surfaces.get(name, [])
            for side, axis in zip(("upper", "lower"), axes):
                points = [row for row in rows if row["side"] == side and (row["y_plus"] or 0) > 0]
                if points:
                    axis.semilogy([row["x_over_c"] for row in points],
                                  [row["y_plus"] for row in points], label=name)
        if axes[0].lines:
            has_wall_functions = any(surfaces.get(name) and by_name[name]["wall_treatment"] == "function"
                                     for name in wall_cases)
            for side, axis in zip(("Upper", "Lower"), axes):
                axis.axhline(1, color="grey", linewidth=0.7, linestyle="--")
                if has_wall_functions:
                    axis.axhspan(30, 100, color="grey", alpha=0.15)
                axis.set(ylabel=side + r" $y^+$", xlim=(0, 1))
            axes[0].legend(fontsize=7, ncols=2)
            axes[-1].set_xlabel(r"$x/c$")
            figure.savefig(output / f"{group}_y_plus.png", dpi=180)
        plt.close(figure)

    surface_plot = plan["surface_plot"]
    figure, (cp_ax, cf_ax) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
    for name, rows in surfaces.items():
        if not rows or name not in surface_plot["cases"]:
            continue
        for side in surface_plot["sides"]:
            linestyle = "-" if side == "upper" else "--"
            points = [item for item in rows if item["side"] == side]
            curve, = cp_ax.plot([item["x_over_c"] for item in points], [item["cp"] for item in points],
                                linestyle=linestyle, label=name if side == "upper" else "_nolegend_")
            if any(item["cf"] is not None for item in points):
                cf_ax.plot([item["x_over_c"] for item in points], [item["cf"] for item in points],
                           color=curve.get_color(), linestyle=linestyle)
    if cp_ax.lines:
        points = [item for item in reference if item["group"] == group]
        sources = list(dict.fromkeys(item.get("source") or "Measured" for item in points))
        for i, source in enumerate(sources):
            selected = [item for item in points if (item.get("source") or "Measured") == source]
            symbols = ("o", "s", "x", "^")
            cp_ax.scatter([float(item["x_over_c"]) for item in selected],
                          [float(item["cp"]) for item in selected], color="black", s=10,
                          marker=symbols[i % len(symbols)], label=source)
        cp_ax.set(ylabel=r"$C_p$", xlim=surface_plot["x_limits"])
        cp_ax.invert_yaxis()
        cp_ax.legend(fontsize=7)
        cf_ax.axhline(0, color="black", linewidth=0.6)
        if bubble := plan.get("bubble_reference"):
            measurements = bubble["measurements"]
            cf_ax.axvspan(measurements["sep_min_x_over_c"], measurements["sep_max_x_over_c"],
                          color="grey", alpha=0.2, label=f"{bubble['label']} separation")
            cf_ax.axvline(measurements["reatt_approx_x_over_c"], color="grey",
                          linestyle=":", label=f"{bubble['label']} reattachment")
            cf_ax.legend(fontsize=7)
        cf_ax.set(xlabel=r"$x/c$", ylabel=r"$C_f$")
        figure.savefig(output / f"{group}_surface.png", dpi=180)
    plt.close(figure)


def main():
    plans = study_analysis()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", type=Path, help="Directory containing generated cases, or one case.")
    parser.add_argument("--out", type=Path, default=Path("analysis-results"))
    parser.add_argument("--group", choices=plans, help="analyze only cases in this study group")
    parser.add_argument("--window", type=float, default=400, help="Final iteration span for force statistics.")
    parser.add_argument("--cl-tol", type=float, default=1e-4, help="Maximum peak-to-peak lift drift.")
    parser.add_argument("--cd-tol", type=float, default=1e-5, help="Maximum peak-to-peak drag drift.")
    parser.add_argument("--cp", type=Path, action="append", default=[],
                        help="Optional reference CSV with group,side,x_over_c,cp,source columns.")
    parser.add_argument("--plots", action="store_true", help="Generate figures if matplotlib is installed.")
    args = parser.parse_args()
    if min(args.window, args.cl_tol, args.cd_tol) <= 0:
        parser.error("the window and both tolerances must be positive")
    case_dirs = [args.cases] if (args.cases / "case.json").is_file() else [
        path.parent for path in sorted(args.cases.glob("*/case.json"))
    ]
    cases = []
    for case in case_dirs:
        meta = json_file(case / "case.json")
        if args.group is None or meta["group"] == args.group:
            cases.append((case, meta))
    if not cases:
        group = f" for group {args.group!r}" if args.group else ""
        parser.error(f"no cases with case.json{group} found in {args.cases}")
    groups = {meta["group"] for case, meta in cases}
    plans = {group: plan for group, plan in plans.items() if group in groups}
    args.out.mkdir(parents=True, exist_ok=True)

    reference, rows, histories, surfaces = [], [], {}, {}
    cp_summary = []
    for path in args.cp:
        with path.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            if not {"group", "side", "x_over_c", "cp"} <= set(reader.fieldnames or []):
                parser.error(f"{path} requires group, side, x_over_c and cp columns")
            reference.extend(reader)
    force_targets = write_references(args.out, plans)
    for case, meta in cases:
        row, history, surface = report_runs(
            case, meta, args.window, args.cl_tol, args.cd_tol, args.out)
        name = row["id"]
        histories[name], surfaces[name] = history, surface
        if row["group"] in force_targets and row.get("cl") is not None:
            reference_name, target = force_targets[row["group"]]
            for key in ("cl", "cd"):
                row[f"{key}_error_vs_{reference_name}_pct"] = 100 * (row[key] - target[key]) / target[key]
            row[f"cd_error_vs_{reference_name}_counts"] = (row["cd"] - target["cd"]) / 1e-4
        errors = cp_errors(surface, reference, row["group"])
        if errors:
            csv_file(args.out / f"{name}_cp_comparison.csv", errors)
            row["cp_reference_points"] = len(errors)
            sources = list(dict.fromkeys(item["source"] for item in errors))
            for source in sources:
                selected = [item for item in errors if item["source"] == source]
                rmse = math.sqrt(fmean(item["cp_error"] ** 2 for item in selected))
                cp_summary.append({"id": name, "source": source, "matched_points": len(selected), "cp_rmse": rmse})
                if len(sources) == 1:
                    row["cp_rmse"] = rmse
        rows.append(row)
    by_name = {row["id"]: row for row in rows}
    gci_rows, contrast_rows = [], []
    for plan in plans.values():
        if case_ids := plan.get("gci_cases"):
            gci_rows.extend(gci_table(by_name, case_ids))
        contrast_rows.extend(contrasts(by_name, plan["contrasts"]))
    csv_file(args.out / "summary.csv", rows)
    csv_file(args.out / "gci.csv", gci_rows)
    csv_file(args.out / "contrasts.csv", contrast_rows)
    if cp_summary:
        csv_file(args.out / "cp_errors_by_source.csv", cp_summary)
    if args.plots:
        make_plots(args.out, histories, surfaces, by_name, reference, plans)
    print(f"Analyzed {len(rows)} cases. Reports: {args.out.resolve()}")


if __name__ == "__main__":
    main()

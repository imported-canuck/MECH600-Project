#!/usr/bin/env python3
"""Write a medium, wall-resolved NACA 0012 blockMeshDict for the low-Re study."""

# The 12-block layout is adapted from highre/mesh/create_mesh.py on the
# teammate's high-re-v2 branch, commit 1389627. OpenFOAM subdivides the blocks.

import argparse
import json
import math
from pathlib import Path


CHORD = 0.20  # m
SPAN = 0.002  # m
RADIUS = 25.0  # chord lengths
SPLIT = 0.30  # divide forward and aft blocks here
INNER_BAND = 0.10  # chord lengths
OUTLET_BAND = 0.75  # widen the inner wake band toward the outlet
FARFIELD_ANGLE = math.radians(25)
FIRST_HEIGHT = 6.5835497453624065e-6  # full first cell in chord lengths
OUTLET_FIRST_HEIGHT = 0.005  # relax wall-scale spacing at the wake outlet
LEADING_SPACING = 4.0e-4  # first surface cell measured along the airfoil
TRAILING_SPACING = 1.5540899958921273e-4
FRONT, AFT, WAKE = 260, 100, 126
INNER, OUTER = 120, 96

# Edges point in the positive local coordinate directions of OpenFOAM's hex.
I_EDGES = ((0, 1), (3, 2), (7, 6), (4, 5))
J_EDGES = ((0, 3), (1, 2), (5, 6), (4, 7))


def thickness(x):
    """Upper surface of the NASA TMR sharp-trailing-edge profile, in chords."""
    if x in (0, 1):
        return 0.0
    return 0.594689181 * (
        0.298222773 * math.sqrt(x) - 0.127125232 * x - 0.357907906 * x**2
        + 0.291984971 * x**3 - 0.105174606 * x**4
    )


def profile(start, stop, samples=401):
    """Sample the curve with extra spline points near the nose and tail."""
    first, last = math.acos(1 - 2 * start), math.acos(1 - 2 * stop)
    xs = [start]
    for i in range(1, samples - 1):
        theta = first + (last - first) * i / (samples - 1)
        xs.append((1 - math.cos(theta)) / 2)
    xs.append(stop)
    return [(x, thickness(x)) for x in xs]


def inner_front_point(x, y):
    """Offset the inner-band edge along the nose normal."""
    if x == 0:
        return -INNER_BAND, 0.0
    # The next block uses an upward offset. Meet it at the same point.
    if x == SPLIT:
        return x, y + INNER_BAND
    slope = 0.594689181 * (
        0.298222773 / (2 * math.sqrt(x)) - 0.127125232
        - 2 * 0.357907906 * x + 3 * 0.291984971 * x**2 - 4 * 0.105174606 * x**3
    )
    normal_scale = INNER_BAND / math.hypot(slope, 1)
    return x - slope * normal_scale, y + normal_scale


def geometric_growth(first, length, cells):
    """Find the per-cell growth factor whose cell widths add to length."""
    if not 0 < first * cells <= length:
        raise ValueError("Requested first cell is too large for increasing widths")
    lower, upper = 1.0, 2.0
    for _ in range(60):
        factor = (lower + upper) / 2
        total = first * sum(factor**i for i in range(cells))
        if total < length:
            lower = factor
        else:
            upper = factor
    return (lower + upper) / 2


def expansion(first, length, cells):
    """Convert per-cell growth to OpenFOAM's last/first size ratio."""
    return geometric_growth(first, length, cells) ** (cells - 1)


def coordinates(values):
    return "(" + " ".join(f"{value:.12g}" for value in values) + ")"


def write_dictionary(output):
    """Describe block corners, shared edge grading, curves and boundary faces."""
    front_curve = profile(0, SPLIT)
    aft_curve = profile(SPLIT, 1)
    inner_front = [inner_front_point(x, y) for x, y in front_curve]
    inner_aft = [(x, y + INNER_BAND) for x, y in aft_curve]

    far_y = RADIUS * math.cos(FARFIELD_ANGLE)
    # Corners on the nose, tail and wake centerline are common to both halves.
    points = {
        "nose": (0.0, 0.0),
        "tail": (1.0, 0.0),
        "wake": (1 + RADIUS, 0.0),
        "inner_nose": inner_front[0],
        "far_nose": (SPLIT - RADIUS, 0.0),
    }
    upper_points = {
        "split": front_curve[-1],
        "inner_split": inner_front[-1],
        "inner_tail": inner_aft[-1],
        "inner_wake": (1 + RADIUS, OUTLET_BAND),
        "far_split": (SPLIT - RADIUS * math.sin(FARFIELD_ANGLE), far_y),
        "far_tail": (1 + RADIUS * math.sin(FARFIELD_ANGLE), far_y),
        "far_wake": (1 + RADIUS, far_y),
    }
    mirror = {name: name for name in points}
    for name, (x, y) in upper_points.items():
        upper, lower = f"upper_{name}", f"lower_{name}"
        points[upper] = x, y
        points[lower] = x, -y
        mirror[upper] = lower

    # Each upper block lists its four corners counterclockwise, then its counts
    # along the airfoil/wake and away from it.
    upper_blocks = [
        (("nose", "upper_split", "upper_inner_split", "inner_nose"), FRONT, INNER),
        (("upper_split", "tail", "upper_inner_tail", "upper_inner_split"), AFT, INNER),
        (("tail", "wake", "upper_inner_wake", "upper_inner_tail"), WAKE, INNER),
        (("inner_nose", "upper_inner_split", "upper_far_split", "far_nose"), FRONT, OUTER),
        (("upper_inner_split", "upper_inner_tail", "upper_far_tail", "upper_far_split"), AFT, OUTER),
        (("upper_inner_tail", "upper_inner_wake", "upper_far_wake", "upper_far_tail"), WAKE, OUTER),
    ]
    blocks = list(upper_blocks)
    for (a, b, c, d), along, outward in upper_blocks:
        # Reflection flips orientation. Swap the two directions to restore it.
        blocks.append(((mirror[a], mirror[d], mirror[c], mirror[b]), outward, along))

    front = expansion(LEADING_SPACING, sum(math.dist(a, b) for a, b in zip(front_curve, front_curve[1:])), FRONT)
    aft = 1 / expansion(TRAILING_SPACING, sum(math.dist(a, b) for a, b in zip(aft_curve, aft_curve[1:])), AFT)
    wake = expansion(TRAILING_SPACING, RADIUS, WAKE)
    far_wake_first = TRAILING_SPACING * math.dist(points["upper_far_split"], points["upper_far_tail"]) / (1 - SPLIT)
    far_wake = expansion(far_wake_first, math.dist(points["upper_far_tail"], points["upper_far_wake"]), WAKE)

    # Define each edge's distribution once and share it between touching blocks.
    upper_grading = {
        ("nose", "upper_split"): front,
        ("inner_nose", "upper_inner_split"): front,
        ("far_nose", "upper_far_split"): front**0.3,
        ("upper_split", "tail"): aft,
        ("upper_inner_split", "upper_inner_tail"): aft,
        ("upper_far_split", "upper_far_tail"): aft,
        ("tail", "wake"): wake,
        ("upper_inner_tail", "upper_inner_wake"): wake,
        ("upper_far_tail", "upper_far_wake"): far_wake,
    }
    for surface, inner, far, first_height in (
        ("nose", "inner_nose", "far_nose", FIRST_HEIGHT),
        ("upper_split", "upper_inner_split", "upper_far_split", FIRST_HEIGHT),
        ("tail", "upper_inner_tail", "upper_far_tail", FIRST_HEIGHT),
        ("wake", "upper_inner_wake", "upper_far_wake", OUTLET_FIRST_HEIGHT),
    ):
        inner_growth = geometric_growth(first_height, math.dist(points[surface], points[inner]), INNER)
        upper_grading[surface, inner] = inner_growth ** (INNER - 1)
        first_outer = first_height * inner_growth**INNER
        upper_grading[inner, far] = expansion(first_outer, math.dist(points[inner], points[far]), OUTER)

    grading = {}
    for (a, b), ratio in upper_grading.items():
        for edge in ((a, b), (mirror[a], mirror[b])):
            grading[edge] = ratio
            grading[edge[::-1]] = 1 / ratio

    point_ids = {name: i for i, name in enumerate(points)}
    offset_id = len(points)
    text = [
        "FoamFile", "{", "    version 2.0;", "    format ascii;",
        "    class dictionary;", "    object blockMeshDict;", "}", "",
        f"scale {CHORD};", "", "vertices", "(",
    ]
    for z in (0.0, SPAN / CHORD):
        for name, point in points.items():
            text.append(f"    {coordinates((*point, z))} // {name}")
    text.extend((");", "", "blocks", "("))
    edge_counts = {}
    for corners, ni, nj in blocks:
        bottom = [point_ids[name] for name in corners]
        top = [i + offset_id for i in bottom]
        ratios = []
        for edges, count in ((I_EDGES, ni), (J_EDGES, nj)):
            for a, b in edges:
                edge = corners[a % 4], corners[b % 4]
                key = frozenset(edge)
                if edge_counts.setdefault(key, count) != count:
                    raise ValueError(f"Adjacent blocks use different cell counts on {edge}")
                ratios.append(grading[edge])
        ratios.extend([1] * 4)
        text.append(f"    hex {coordinates(bottom + top)} ({ni} {nj} 1) edgeGrading {coordinates(ratios)}")
    text.extend((");", "", "edges", "("))

    # Describe curves on both faces of the single-cell-thick mesh.
    midpoint = (math.pi + math.pi / 2 + FARFIELD_ANGLE) / 2
    arc_point = SPLIT + RADIUS * math.cos(midpoint), RADIUS * math.sin(midpoint)
    for layer, z in enumerate((0.0, SPAN / CHORD)):
        for side in (1, -1):
            for kind, a, b, curve in (
                ("spline", "nose", "upper_split", front_curve[1:-1]),
                ("spline", "upper_split", "tail", aft_curve[1:-1]),
                ("spline", "inner_nose", "upper_inner_split", inner_front[1:-1]),
                ("spline", "upper_inner_split", "upper_inner_tail", inner_aft[1:-1]),
                ("arc", "far_nose", "upper_far_split", [arc_point]),
            ):
                if side < 0:
                    a, b = mirror[a], mirror[b]
                start, end = point_ids[a] + layer * offset_id, point_ids[b] + layer * offset_id
                controls = " ".join(coordinates((x, side * y, z)) for x, y in curve)
                text.append(f"    {kind} {start} {end} " + (f"({controls})" if kind == "spline" else controls))
    text.extend((");", "", "boundary", "("))

    # Boundary loops follow each block's outside orientation. Shared internal
    # edges are omitted so blockMesh joins the upper and lower halves.
    airfoil = ("nose", "upper_split", "tail", "lower_split", "nose")
    farfield = (
        "far_nose", "lower_far_split", "lower_far_tail", "lower_far_wake",
        "lower_inner_wake", "wake", "upper_inner_wake", "upper_far_wake",
        "upper_far_tail", "upper_far_split", "far_nose",
    )
    for name, patch_type, loop in (("airfoil", "wall", airfoil), ("farfield", "patch", farfield)):
        text.extend((f"    {name}", "    {", f"        type {patch_type};", "        faces", "        ("))
        for a, b in zip(loop, loop[1:]):
            i, j = point_ids[a], point_ids[b]
            text.append(f"            {coordinates((i, j, j + offset_id, i + offset_id))}")
        text.extend(("        );", "    }"))
    text.extend((");", "", "defaultPatch", "{", "    name frontAndBack;", "    type empty;", "}", ""))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(text), encoding="utf-8")
    return sum(ni * nj for _, ni, nj in blocks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Path to system/blockMeshDict")
    parser.add_argument("--metadata", type=Path, help="Optional mesh.json report path")
    args = parser.parse_args()
    cells = write_dictionary(args.output)
    if args.metadata:
        data = {
            "backend": "blockMesh", "cells": cells, "chord": CHORD, "span": SPAN,
            "outer_extent_in_chords": RADIUS, "level": "medium", "wall_treatment": "resolved",
            "alpha_deg": 10.0, "transition_mesh": True,
            "first_cell_height_over_c": FIRST_HEIGHT, "first_cell_height_m": FIRST_HEIGHT * CHORD,
            "wall_faces_per_side": FRONT + AFT, "normal_cells": INNER + OUTER, "wake_cells": WAKE,
            "wake_outlet_band_over_c": OUTLET_BAND,
            "wake_outlet_first_height_over_c": OUTLET_FIRST_HEIGHT,
        }
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} for {cells:,} cells and a {FIRST_HEIGHT * CHORD * 1e6:.3f} µm target first cell")


if __name__ == "__main__":
    main()

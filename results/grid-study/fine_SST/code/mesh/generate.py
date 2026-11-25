#!/usr/bin/env python3
"""Generate a one-layer, body-fitted Gmsh 2.2 C-grid for the TMR NACA 0012."""

import argparse
import json
import math
from pathlib import Path


EXTENDED_CHORD = 1.008930411365
LEVELS = {"coarse": 1, "medium": 1.5, "fine": 2.25}
PATCHES = {1: "airfoil", 2: "farfield", 3: "frontAndBack", 4: "fluid"}


def thickness(x):
    if x == 0 or x == 1:
        return 0.0
    s = EXTENDED_CHORD * x
    return (0.6 / EXTENDED_CHORD) * (
        0.2969 * math.sqrt(s) - 0.1260 * s - 0.3516 * s**2
        + 0.2843 * s**3 - 0.1015 * s**4
    )


def wall_normal(x):
    if x == 0:
        return -1.0, 0.0
    s = EXTENDED_CHORD * x
    slope = 0.6 * (
        0.2969 / (2 * math.sqrt(s)) - 0.1260 - 2 * 0.3516 * s
        + 3 * 0.2843 * s**2 - 4 * 0.1015 * s**3
    )
    length = math.hypot(1, slope)
    return -slope / length, 1 / length


def inverse_monotone(function, target, lower, upper):
    for _ in range(45):
        midpoint = (lower + upper) / 2
        if function(midpoint) < target:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2


def body_nodes(count, transition):
    # A smooth extra point density on x/c = 0..0.15 resolves the short bubble.
    def leading_edge_weight(x):
        return x - math.log1p(math.exp(50 * (x - 0.18))) / 50

    origin = leading_edge_weight(0)
    integral = leading_edge_weight(1) - origin

    def distribution(x):
        cosine = 2 * math.asin(math.sqrt(x)) / math.pi
        if transition:
            near_nose = (leading_edge_weight(x) - origin) / integral
            return 0.35 * cosine + 0.65 * near_nose
        # The surface tangent is nearly horizontal at maximum thickness. Reduce
        # the tangential cell width there to tame the first-layer aspect ratio.
        near_maximum_thickness = 0.5 * (1 + math.erf((x - 0.2973) / 0.025))
        return 0.8 * cosine + 0.2 * near_maximum_thickness

    interior = [inverse_monotone(distribution, i / count, 0, 1) for i in range(1, count)]
    return [0.0, *interior, 1.0]


def geometric_distance(first, log_ratio, steps):
    if abs(log_ratio) < 1e-12:
        return first * steps
    return first * math.expm1(log_ratio * steps) / math.expm1(log_ratio)


def graded_edges(marks, end, first, factor, first_count=None):
    """Keep node locations common before each domain-study boundary."""
    edges = [0.0]
    target_ratio = 1.11
    for boundary in [*(mark for mark in marks if mark < end), end]:
        start = edges[-1]
        length = boundary - start
        nominal_steps = math.log1p(length * (target_ratio - 1) / first) / math.log(target_ratio)
        base_count = max(4, 4 * round(nominal_steps / 4))
        if start == 0 and first_count is not None:
            base_count = first_count
        count = round(factor * base_count)
        log_ratio = inverse_monotone(
            lambda ratio: geometric_distance(first, ratio, base_count), length, -6, 600 / base_count
        )
        if abs(geometric_distance(first, log_ratio, base_count) - length) > 1e-9 * length:
            raise ValueError(f"First spacing {first:g} is incompatible with interval {length:g}")
        edges.extend(
            start + geometric_distance(first, log_ratio, i * base_count / count)
            for i in range(1, count)
        )
        edges.append(boundary)
        first *= math.exp(log_ratio * base_count)
    return edges


def replace_wall_layers(edges, height):
    """Change normal layers inside 0.05c and keep the outer grid identical."""
    interface_index = edges.index(0.05)
    outer_spacing = edges[interface_index + 1] - 0.05
    if not 0 < height < 0.05:
        raise ValueError("A replacement first height must be between zero and 0.05 chord")

    nominal_steps = 0.05 / height
    if outer_spacing != height:
        nominal_steps = 0.05 * math.log(outer_spacing / height) / (outer_spacing - height)
    count = max(2, round(nominal_steps))
    ratio = inverse_monotone(lambda q: geometric_distance(height, q, count), 0.05, -6, 600 / count)
    return [0.0, *(geometric_distance(height, ratio, i) for i in range(1, count)), *edges[interface_index:]]


def section_points(x_wall, x_wake):
    points = []
    for x in x_wall:
        if x == 0:
            outer_normal = (-1.0, 0.0)
        else:
            outer_normal = (-math.sqrt(1 - x), math.sqrt(x))
        points.append((x, thickness(x), *wall_normal(x), *outer_normal))
    te_nx, te_ny = wall_normal(1.0)
    for x in x_wake[1:]:
        taper = math.exp(-x / 0.20)
        points.append((1 + x, 0.0, taper * te_nx, 1 + taper * (te_ny - 1), 0.0, 1.0))
    return points


def coordinates(point, distance, sign, alpha):
    x, y, nx, ny, outer_nx, outer_ny = point
    wake_distance = max(0, x - 1)
    if wake_distance:
        # Open the fine wake layers downstream instead of carrying the wall spacing
        # all the way to the exit. Begin with zero growth at the trailing edge.
        opening = wake_distance + 0.20 * math.expm1(-wake_distance / 0.20)
        far_wake = opening * (wake_distance / (wake_distance + 4))**3
        distance = math.sqrt(distance * (distance + 0.10 * opening)) + (
            0.011 * far_wake * math.asinh(math.sqrt(distance / 0.001))
        )
    inner_distance = distance * 0.15 / (0.15 + distance)
    outer_distance = distance - inner_distance
    if wake_distance:
        # Match the sqrt(1 - x/c) radial fan on the airfoil side of the TE.
        outer_nx = math.sqrt(wake_distance / (1 + 25 * wake_distance))
    x += inner_distance * nx + outer_distance * outer_nx
    y = sign * (y + inner_distance * ny + outer_distance * outer_ny)
    downstream = max(0, x - 1)
    # Shear only downstream of the trailing edge to align the wake with the inflow.
    y += math.tan(alpha) * (downstream + 0.4 * math.expm1(-downstream / 0.4))
    return x, y


def section_nodes(points, radial, wall_count, alpha):
    xy = []
    sections = {}
    for side in (1, -1):
        rows = []
        for j, distance in enumerate(radial):
            row = []
            for i, point in enumerate(points):
                if side == -1 and (i == 0 or (j == 0 and i >= wall_count)):
                    row.append(sections[1][j][i])
                else:
                    row.append(len(xy) + 1)
                    xy.append(coordinates(point, distance, side, alpha))
            rows.append(row)
        sections[side] = rows
    return xy, sections


def signed_corner_areas(quad, xy):
    corners = [xy[node - 1] for node in quad]
    for index in range(4):
        a = corners[index - 1]
        b = corners[index]
        c = corners[(index + 1) % 4]
        yield (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])


def quad_area(quad, xy):
    x0, y0 = xy[quad[0] - 1]
    x1, y1 = xy[quad[1] - 1]
    x2, y2 = xy[quad[2] - 1]
    x3, y3 = xy[quad[3] - 1]
    return 0.5 * ((x1 - x3) * (y2 - y0) - (y1 - y3) * (x2 - x0))


def connectivity(xy, sections, wall_count):
    volume = []
    patches = {1: [], 2: [], 3: []}
    node_count = len(xy)
    min_corner = math.inf
    for side, rows in sections.items():
        last_row = len(rows) - 2
        last_col = len(rows[0]) - 2
        for j, (inner, outer) in enumerate(zip(rows, rows[1:])):
            for i in range(len(inner) - 1):
                edges = [(inner[i], inner[i + 1]), (inner[i + 1], outer[i + 1]),
                         (outer[i + 1], outer[i]), (outer[i], inner[i])]
                quad = tuple(edge[0] for edge in edges)
                if side == -1:
                    quad = quad[::-1]
                corner_area = min(signed_corner_areas(quad, xy))
                if corner_area <= 0:
                    raise ValueError(f"Non-convex C-grid: side={side}, row={j}, column={i}")
                min_corner = min(min_corner, corner_area)
                top = tuple(node + node_count for node in quad)
                volume.append(quad + top)
                patches[3].extend((quad[::-1], top))

                def add_boundary(tag, edge):
                    a, b = edge if side == 1 else edge[::-1]
                    patches[tag].append((a, b, b + node_count, a + node_count))

                if j == 0 and i < wall_count:
                    add_boundary(1, edges[0])
                if i == last_col:
                    add_boundary(2, edges[1])
                if j == last_row:
                    add_boundary(2, edges[2])
    return patches, volume, min_corner


def write_ascii(path, coordinates_3d, patches, volume):
    with path.open("w", encoding="ascii", newline="\n") as output:
        output.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n4\n')
        for tag, name in PATCHES.items():
            output.write(f'{3 if tag == 4 else 2} {tag} "{name}"\n')
        output.write(f'$EndPhysicalNames\n$Nodes\n{len(coordinates_3d)}\n')
        for index, (x, y, z) in enumerate(coordinates_3d, 1):
            output.write(f'{index} {x:.16g} {y:.16g} {z:.16g}\n')
        total = len(volume) + sum(len(faces) for faces in patches.values())
        output.write(f'$EndNodes\n$Elements\n{total}\n')
        index = 1
        for tag, cells in (*patches.items(), (4, volume)):
            kind = 5 if tag == 4 else 3
            for cell in cells:
                output.write(f'{index} {kind} 2 {tag} {tag} {" ".join(map(str, cell))}\n')
                index += 1
        output.write('$EndElements\n')


def write_gmsh(gmsh, path, coordinates_3d, patches, volume):
    gmsh.initialize()
    try:
        gmsh.model.add("naca0012")
        for tag, name in PATCHES.items():
            dimension = 3 if tag == 4 else 2
            gmsh.model.addDiscreteEntity(dimension, tag)
            gmsh.model.addPhysicalGroup(dimension, [tag], tag)
            gmsh.model.setPhysicalName(dimension, tag, name)
        gmsh.model.mesh.addNodes(3, 4, range(1, len(coordinates_3d) + 1),
                                 [value for point in coordinates_3d for value in point])
        for tag, cells in (*patches.items(), (4, volume)):
            gmsh.model.mesh.addElementsByType(tag, 5 if tag == 4 else 3, [],
                                              [node for cell in cells for node in cell])
        gmsh.model.mesh.reclassifyNodes()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


def options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--chord", required=True, type=float, help="Chord in metres")
    parser.add_argument("--radius", required=True, type=float, help="Outer extent in chords")
    parser.add_argument("--span", type=float, help="Extrusion depth in metres, default 0.01c")
    parser.add_argument("--alpha", type=float, default=10.0, help="Wake alignment angle in degrees")
    parser.add_argument("--level", required=True, choices=LEVELS)
    parser.add_argument("--wall", required=True, choices=("resolved", "function"))
    parser.add_argument("--transition", action="store_true", help="Increase leading-edge resolution")
    parser.add_argument("--first-height", type=float, help="Full first-cell height divided by chord")
    parser.add_argument("--backend", choices=("auto", "gmsh", "ascii"), default="auto",
                        help="Use Gmsh Python API when available, or write the same ASCII MSH2 directly")
    return parser.parse_args()


def main():
    args = options()
    span = args.span if args.span is not None else 0.01 * args.chord
    if not all(math.isfinite(value) for value in (args.chord, span, args.radius, args.alpha)):
        raise SystemExit("Chord, span, radius, and incidence must be finite")
    if min(args.chord, span) <= 0 or args.radius < 2 or not -45 <= args.alpha <= 45:
        raise SystemExit("Require positive chord/span, radius >= 2c, and incidence in [-45, 45] degrees")
    if args.output.suffix.lower() != ".msh":
        raise SystemExit("Output file extension must be .msh")
    factor = LEVELS[args.level]
    wall_count = round(factor * (240 if args.transition else 128))
    wall = body_nodes(wall_count, args.transition)
    base_wall = body_nodes(240 if args.transition else 128, args.transition)
    first_wake = 1 - base_wall[-2]
    landmarks = [25, 50, 100, 250, 500]
    wake = graded_edges([2, 10, *landmarks], args.radius, first_wake, factor)
    first_resolved = 1e-5 if args.transition else 3e-6
    radial = graded_edges([0.05, *landmarks], args.radius, first_resolved, factor, first_count=80)
    if args.first_height is not None or args.wall == "function":
        radial = replace_wall_layers(radial, args.first_height if args.first_height is not None else 4e-4)

    xy, sections = section_nodes(section_points(wall, wake), radial, wall_count, math.radians(args.alpha))
    patches, volume, min_corner = connectivity(xy, sections, wall_count)
    coords = [(args.chord * x, args.chord * y, z) for z in (0.0, span) for x, y in xy]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    gmsh = None
    if args.backend != "ascii":
        try:
            import gmsh
        except (ImportError, OSError, TypeError) as error:
            if args.backend == "gmsh":
                raise SystemExit(f"Gmsh Python API is unavailable: {error}") from error
    if gmsh is None:
        write_ascii(args.output, coords, patches, volume)
    else:
        write_gmsh(gmsh, args.output, coords, patches, volume)

    areas = [quad_area(cell, xy) * args.chord**2 for cell in volume]
    metadata = {
        "backend": "gmsh" if gmsh else "ascii",
        "cells": len(volume), "nodes": len(coords), "fluid_area": math.fsum(areas),
        "min_cell_volume_m3": min(areas) * span,
        "chord": args.chord, "span": span, "outer_extent_in_chords": args.radius,
        "level": args.level, "wall_treatment": args.wall, "alpha_deg": args.alpha,
        "transition_mesh": args.transition,
        "wall_faces_per_side": wall_count, "normal_cells": len(radial) - 1,
        "wake_cells": len(wake) - 1,
        "first_cell_height_over_c": radial[1], "first_cell_height_m": radial[1] * args.chord,
        "nose_streamwise_spacing_over_c": wall[1],
        "max_streamwise_spacing_before_x_c_0p15": max(
            right - left for left, right in zip(wall, wall[1:]) if right < 0.15
        ),
        "trailing_edge_streamwise_spacing_over_c": 1 - wall[-2],
        "min_corner_cross_product_m2": min_corner * args.chord**2,
        "boundary_faces": {PATCHES[tag]: len(faces) for tag, faces in patches.items()},
    }
    if args.metadata:
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

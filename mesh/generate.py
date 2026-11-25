#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path
import gmsh


# The NASA TMR geometry extends the original profile to a sharp trailing edge
# and rescales it so the chord still runs from x/c = 0 to x/c = 1.
EXTENDED_CHORD = 1.008930411365
REFINEMENT = {"coarse": 1, "medium": 1.5, "fine": 2.25}
DOMAIN_CHECK_RADII = (25, 50, 100, 250, 500)
WALL_BAND = 0.05  # Outer edge of the near-wall band, in chords.

UPPER, LOWER = 1, -1
AIRFOIL, FARFIELD, FRONT_AND_BACK, FLUID = 1, 2, 3, 4
PHYSICAL_NAMES = {
    AIRFOIL: "airfoil",
    FARFIELD: "farfield",
    FRONT_AND_BACK: "frontAndBack",
    FLUID: "fluid",
}


def thickness(x):
    """Return the upper surface height y/c at the chord position x/c."""
    if x == 0 or x == 1:
        return 0.0
    profile_x = EXTENDED_CHORD * x
    return (0.6 / EXTENDED_CHORD) * (
        0.2969 * math.sqrt(profile_x)
        - 0.1260 * profile_x
        - 0.3516 * profile_x**2
        + 0.2843 * profile_x**3
        - 0.1015 * profile_x**4
    )


def wall_normal(x):
    """Return the outward unit normal to the upper surface at x/c."""
    if x == 0:
        return -1.0, 0.0  # The normal points straight upstream at the nose.
    profile_x = EXTENDED_CHORD * x
    # Differentiate thickness to get the slope. Then rotate the tangent
    # (1, slope) by 90 degrees and scale the result to unit length.
    slope = 0.6 * (
        0.2969 / (2 * math.sqrt(profile_x))
        - 0.1260
        - 2 * 0.3516 * profile_x
        + 3 * 0.2843 * profile_x**2
        - 4 * 0.1015 * profile_x**3
    )
    length = math.hypot(1, slope)
    return -slope / length, 1 / length


def invert_increasing(function, target, lower, upper):
    """Find where an increasing function reaches target, using bisection."""
    for _ in range(45):
        midpoint = (lower + upper) / 2
        if function(midpoint) < target:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2


def airfoil_stations(cell_count, transition):
    """Place points along the chord, with extra points where spacing matters."""
    def integrated_nose_density(x):
        # This curve rises quickly near the nose and levels off farther aft.
        return x - math.log1p(math.exp(50 * (x - 0.18))) / 50

    nose_start = integrated_nose_density(0)
    nose_total = integrated_nose_density(1) - nose_start

    def fraction_of_cells(x):
        # Cosine spacing concentrates points at both ends of the airfoil.
        cosine_fraction = 2 * math.asin(math.sqrt(x)) / math.pi
        if transition:
            # Put extra points in the first 15% of the chord for the bubble.
            nose_fraction = (integrated_nose_density(x) - nose_start) / nose_total
            return 0.35 * cosine_fraction + 0.65 * nose_fraction
        # Add points around maximum thickness to reduce cell aspect ratios.
        thickness_fraction = 0.5 * (1 + math.erf((x - 0.2973) / 0.025))
        return 0.8 * cosine_fraction + 0.2 * thickness_fraction

    # Place one point at each fraction of the requested cell count.
    stations = [0.0]
    for index in range(1, cell_count):
        stations.append(invert_increasing(fraction_of_cells, index / cell_count, 0, 1))
    stations.append(1.0)
    return stations


def graded_length(first_spacing, log_growth, cell_count):
    """Sum the widths of cells with a constant growth ratio."""
    if abs(log_growth) < 1e-12:
        return first_spacing * cell_count
    # expm1 keeps the sum accurate when the growth ratio is close to one.
    return first_spacing * math.expm1(log_growth * cell_count) / math.expm1(log_growth)


def graded_edges(landmarks, extent, first_spacing, refinement, first_count=None):
    """Place cell edges along the wake or out toward the farfield."""
    edges = [0.0]
    interval_ends = [landmark for landmark in landmarks if landmark < extent]
    interval_ends.append(extent)
    preferred_growth = 1.11

    # Keeping the same landmarks lets a larger domain extend the existing grid.
    for end in interval_ends:
        start = edges[-1]
        length = end - start
        nominal_count = math.log1p(length * (preferred_growth - 1) / first_spacing)
        nominal_count /= math.log(preferred_growth)
        # Multiples of four give whole cell counts at both refinement factors.
        base_count = max(4, 4 * round(nominal_count / 4))
        if start == 0 and first_count is not None:
            base_count = first_count

        # Fit the growth ratio to hit the end exactly, then sample the same
        # spacing curve more densely for the medium and fine grids.
        log_growth = invert_increasing(
            lambda growth: graded_length(first_spacing, growth, base_count),
            length, -6, 600 / base_count,
        )
        fitted_length = graded_length(first_spacing, log_growth, base_count)
        if abs(fitted_length - length) > 1e-9 * length:
            raise ValueError(f"First spacing {first_spacing:g} is incompatible with interval {length:g}")

        cell_count = round(refinement * base_count)
        for index in range(1, cell_count):
            base_index = index * base_count / cell_count
            edges.append(start + graded_length(first_spacing, log_growth, base_index))
        edges.append(end)
        # Start the next interval where the geometric growth would continue.
        first_spacing *= math.exp(log_growth * base_count)
    return edges


def replace_wall_layers(edges, first_height):
    """Set the first-cell height inside the wall band, keeping the outer grid."""
    if not 0 < first_height < WALL_BAND:
        raise ValueError("A replacement first height must be between zero and 0.05 chord")

    interface = edges.index(WALL_BAND)
    outer_spacing = edges[interface + 1] - WALL_BAND
    # Estimate how many layers connect the new first height to the next
    # outer spacing. Refit the growth ratio after rounding to a whole count.
    nominal_count = WALL_BAND / first_height
    if outer_spacing != first_height:
        nominal_count = WALL_BAND * math.log(outer_spacing / first_height) / (outer_spacing - first_height)
    cell_count = max(2, round(nominal_count))
    log_growth = invert_increasing(
        lambda growth: graded_length(first_height, growth, cell_count),
        WALL_BAND, -6, 600 / cell_count,
    )
    inner_edges = [graded_length(first_height, log_growth, i) for i in range(1, cell_count)]
    return [0.0, *inner_edges, *edges[interface:]]


def c_grid_stations(wall_positions, wake_offsets):
    """Store each starting point and the directions of its inner and outer rows."""
    stations = []
    for x in wall_positions:
        surface_point = (x, thickness(x))
        # Outer rows fan backward at the nose and outward on the sides.
        farfield_direction = (-math.sqrt(1 - x), math.sqrt(x))
        stations.append((surface_point, wall_normal(x), farfield_direction))

    trailing_normal_x, trailing_normal_y = wall_normal(1.0)
    for offset in wake_offsets[1:]:  # The airfoil already includes the trailing edge.
        # Turn the inner rows smoothly from the wall normal toward vertical.
        trailing_edge_weight = math.exp(-offset / 0.20)
        inner_direction = (
            trailing_edge_weight * trailing_normal_x,
            1 + trailing_edge_weight * (trailing_normal_y - 1),
        )
        stations.append(((1 + offset, 0.0), inner_direction, (0.0, 1.0)))
    return stations


def c_grid_point(station, offset, side, wake_slope):
    """Find the grid point at the given station and nominal outward offset."""
    (x, y), (inner_x, inner_y), (outer_x, outer_y) = station
    wake_distance = max(0, x - 1)
    if wake_distance:
        # Let thin wake layers widen downstream. The smooth opening starts
        # at zero at the trailing edge, where the wall and wake rows meet.
        opening = wake_distance + 0.20 * math.expm1(-wake_distance / 0.20)
        farther_opening = opening * (wake_distance / (wake_distance + 4))**3
        offset = math.sqrt(offset * (offset + 0.10 * opening)) + (
            0.011 * farther_opening * math.asinh(math.sqrt(offset / 0.001))
        )
        # Continue the outer fan on the wake side of the trailing edge.
        outer_x = math.sqrt(wake_distance / (1 + 25 * wake_distance))

    # Near points follow the inner direction. The outer direction gradually
    # takes over as the distance from the surface or wake increases.
    inner_offset = offset * 0.15 / (0.15 + offset)
    outer_offset = offset - inner_offset
    x += inner_offset * inner_x + outer_offset * outer_x
    y = side * (y + inner_offset * inner_y + outer_offset * outer_y)

    # Keep the airfoil fixed and turn the downstream wake toward the inflow.
    downstream = max(0, x - 1)
    y += wake_slope * (downstream + 0.4 * math.expm1(-downstream / 0.4))
    return x, y


def plane_grid(stations, row_offsets, wall_cells, alpha):
    """Create 2D points and assign each row its Gmsh node numbers."""
    points = []
    sections = {}
    wake_slope = math.tan(alpha)
    for side in (UPPER, LOWER):
        rows = []
        for row_index, offset in enumerate(row_offsets):
            row = []
            for column, station in enumerate(stations):
                # Both halves share the upstream seam and the wake centreline.
                # Reusing those node numbers joins the two halves of the mesh.
                on_seam = column == 0 or (row_index == 0 and column >= wall_cells)
                if side == LOWER and on_seam:
                    row.append(sections[UPPER][row_index][column])
                else:
                    points.append(c_grid_point(station, offset, side, wake_slope))
                    row.append(len(points))  # Gmsh numbers start at 1.
            rows.append(row)
        sections[side] = rows
    return points, sections


def corner_cross_products(quad, points):
    """Positive turns at all four corners identify a valid convex cell."""
    corners = [points[node - 1] for node in quad]
    for index, (corner_x, corner_y) in enumerate(corners):
        previous_x, previous_y = corners[index - 1]
        next_x, next_y = corners[(index + 1) % 4]
        yield (
            (corner_x - previous_x) * (next_y - corner_y)
            - (corner_y - previous_y) * (next_x - corner_x)
        )


def quad_area(quad, points):
    """Calculate the area from the cross product of the two diagonals."""
    x0, y0 = points[quad[0] - 1]
    x1, y1 = points[quad[1] - 1]
    x2, y2 = points[quad[2] - 1]
    x3, y3 = points[quad[3] - 1]
    return 0.5 * ((x1 - x3) * (y2 - y0) - (y1 - y3) * (x2 - x0))


def extruded_face(start, end, side, plane_node_count):
    """Turn a 2D boundary edge into a face, with an outward orientation."""
    if side == LOWER:
        start, end = end, start
    return start, end, end + plane_node_count, start + plane_node_count


def mesh_elements(points, sections, wall_cells):
    """Extrude each 2D quadrilateral into one hexahedron and label its boundaries."""
    cells = []
    patches = {AIRFOIL: [], FARFIELD: [], FRONT_AND_BACK: []}
    plane_node_count = len(points)
    minimum_corner = math.inf

    for side, rows in sections.items():
        last_cell_row = len(rows) - 2
        last_cell_column = len(rows[0]) - 2
        for row_index, (inner_row, outer_row) in enumerate(zip(rows, rows[1:])):
            for column in range(len(inner_row) - 1):
                inner_left, inner_right = inner_row[column:column + 2]
                outer_left, outer_right = outer_row[column:column + 2]
                front = (inner_left, inner_right, outer_right, outer_left)
                if side == LOWER:
                    front = front[::-1]  # Match the winding on the upper half.
                corner = min(corner_cross_products(front, points))
                if corner <= 0:
                    raise ValueError(f"Non-convex C-grid: side={side}, row={row_index}, column={column}")
                minimum_corner = min(minimum_corner, corner)

                # The back plane is a copy of the front with new node numbers.
                back = tuple(node + plane_node_count for node in front)
                cells.append(front + back)
                patches[FRONT_AND_BACK].extend((front[::-1], back))

                if row_index == 0 and column < wall_cells:
                    face = extruded_face(inner_left, inner_right, side, plane_node_count)
                    patches[AIRFOIL].append(face)
                if column == last_cell_column:  # Downstream boundary.
                    face = extruded_face(inner_right, outer_right, side, plane_node_count)
                    patches[FARFIELD].append(face)
                if row_index == last_cell_row:  # Outermost curved boundary.
                    face = extruded_face(outer_right, outer_left, side, plane_node_count)
                    patches[FARFIELD].append(face)
    return patches, cells, minimum_corner


def write_ascii(path, points, patches, cells):
    """Write a Gmsh 2.2 text file using only Python's standard library."""
    with path.open("w", encoding="ascii", newline="\n") as output:
        output.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n4\n")
        for tag, name in PHYSICAL_NAMES.items():
            dimension = 3 if tag == FLUID else 2
            output.write(f'{dimension} {tag} "{name}"\n')

        output.write(f"$EndPhysicalNames\n$Nodes\n{len(points)}\n")
        for node, (x, y, z) in enumerate(points, 1):
            output.write(f"{node} {x:.16g} {y:.16g} {z:.16g}\n")

        face_count = sum(len(faces) for faces in patches.values())
        output.write(f"$EndNodes\n$Elements\n{len(cells) + face_count}\n")
        element_number = 1
        for tag, elements in (*patches.items(), (FLUID, cells)):
            # Gmsh element types 5 and 3 mean hexahedron and quadrilateral.
            element_type = 5 if tag == FLUID else 3
            for element in elements:
                nodes = " ".join(map(str, element))
                output.write(f"{element_number} {element_type} 2 {tag} {tag} {nodes}\n")
                element_number += 1
        output.write("$EndElements\n")


def write_gmsh(path, points, patches, cells):
    """Export the points and cells through the Gmsh Python API."""
    gmsh.initialize()
    try:
        gmsh.model.add("naca0012")
        for tag, name in PHYSICAL_NAMES.items():
            dimension = 3 if tag == FLUID else 2
            gmsh.model.addDiscreteEntity(dimension, tag)
            gmsh.model.addPhysicalGroup(dimension, [tag], tag)
            gmsh.model.setPhysicalName(dimension, tag, name)

        mesh = gmsh.model.mesh
        coordinates = [coordinate for point in points for coordinate in point]
        mesh.addNodes(3, FLUID, range(1, len(points) + 1), coordinates)
        for tag, elements in (*patches.items(), (FLUID, cells)):
            element_type = 5 if tag == FLUID else 3
            node_numbers = [node for element in elements for node in element]
            mesh.addElementsByType(tag, element_type, [], node_numbers)
        mesh.reclassifyNodes()  # Associate each node with its boundary or volume.
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


def options():
    """Read the settings passed by study.py or entered on the command line."""
    parser = argparse.ArgumentParser(description="Build a one-cell-thick NACA 0012 C-grid.")
    parser.add_argument("--output", required=True, type=Path, help="Output mesh path (.msh)")
    parser.add_argument("--metadata", type=Path, help="Optional JSON report path")
    parser.add_argument("--chord", required=True, type=float, help="Chord in metres")
    parser.add_argument("--radius", required=True, type=float, help="Outer extent in chords")
    parser.add_argument("--span", type=float, help="Extrusion depth in metres, default 0.01c")
    parser.add_argument("--alpha", type=float, default=10.0, help="Wake alignment angle in degrees")
    parser.add_argument("--level", required=True, choices=REFINEMENT)
    parser.add_argument("--wall", required=True, choices=("resolved", "function"))
    parser.add_argument("--transition", action="store_true", help="Increase leading-edge resolution")
    parser.add_argument("--first-height", type=float, help="Full first-cell height divided by chord")
    parser.add_argument(
        "--backend", choices=("auto", "gmsh", "ascii"), default="auto",
        help="Use the Gmsh API (auto or gmsh), or write ASCII MSH2 directly (ascii)",
    )
    return parser.parse_args()


def main():
    """Build the mesh from the requested settings and write its report."""
    args = options()
    span = args.span if args.span is not None else 0.01 * args.chord
    if not all(math.isfinite(value) for value in (args.chord, span, args.radius, args.alpha)):
        raise SystemExit("Chord, span, radius, and incidence must be finite")
    if min(args.chord, span) <= 0 or args.radius < 2 or not -45 <= args.alpha <= 45:
        raise SystemExit("Require positive chord/span, radius >= 2c, and incidence in [-45, 45] degrees")
    if args.output.suffix.lower() != ".msh":
        raise SystemExit("Output file extension must be .msh")

    # Space the stations along the airfoil and wake, and the rows out to the
    # farfield. Use the base airfoil spacing to start the wake at every level.
    refinement = REFINEMENT[args.level]
    base_wall_cells = 240 if args.transition else 128
    wall_cells = round(refinement * base_wall_cells)
    wall_positions = airfoil_stations(wall_cells, args.transition)
    base_positions = airfoil_stations(base_wall_cells, args.transition)
    wake_offsets = graded_edges(
        (2, 10, *DOMAIN_CHECK_RADII), args.radius, 1 - base_positions[-2], refinement,
    )
    first_resolved_spacing = 1e-5 if args.transition else 3e-6
    row_offsets = graded_edges(
        (WALL_BAND, *DOMAIN_CHECK_RADII), args.radius, first_resolved_spacing,
        refinement, first_count=80,
    )
    if args.first_height is not None or args.wall == "function":
        first_height = args.first_height if args.first_height is not None else 4e-4
        row_offsets = replace_wall_layers(row_offsets, first_height)

    stations = c_grid_stations(wall_positions, wake_offsets)
    points, sections = plane_grid(stations, row_offsets, wall_cells, math.radians(args.alpha))
    patches, cells, minimum_corner = mesh_elements(points, sections, wall_cells)
    points_3d = [(args.chord * x, args.chord * y, z) for z in (0.0, span) for x, y in points]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    backend = "gmsh" if args.backend == "auto" else args.backend
    if backend == "ascii":
        write_ascii(args.output, points_3d, patches, cells)
    else:
        write_gmsh(args.output, points_3d, patches, cells)

    # Keep the mesh report used by study.py and the analysis scripts.
    cell_areas = [quad_area(cell, points) * args.chord**2 for cell in cells]
    bubble_region_spacing = max(
        right - left for left, right in zip(wall_positions, wall_positions[1:]) if right < 0.15
    )
    metadata = {
        "backend": backend,
        "cells": len(cells),
        "nodes": len(points_3d),
        "fluid_area": math.fsum(cell_areas),
        "min_cell_volume_m3": min(cell_areas) * span,
        "chord": args.chord,
        "span": span,
        "outer_extent_in_chords": args.radius,
        "level": args.level,
        "wall_treatment": args.wall,
        "alpha_deg": args.alpha,
        "transition_mesh": args.transition,
        "wall_faces_per_side": wall_cells,
        "normal_cells": len(row_offsets) - 1,
        "wake_cells": len(wake_offsets) - 1,
        "first_cell_height_over_c": row_offsets[1],
        "first_cell_height_m": row_offsets[1] * args.chord,
        "nose_streamwise_spacing_over_c": wall_positions[1],
        "max_streamwise_spacing_before_x_c_0p15": bubble_region_spacing,
        "trailing_edge_streamwise_spacing_over_c": 1 - wall_positions[-2],
        "min_corner_cross_product_m2": minimum_corner * args.chord**2,
        "boundary_faces": {PHYSICAL_NAMES[tag]: len(faces) for tag, faces in patches.items()},
    }
    report = json.dumps(metadata, indent=2)
    if args.metadata:
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(report + "\n", encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()

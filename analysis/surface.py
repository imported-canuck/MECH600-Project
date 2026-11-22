"""Extract airfoil patch diagnostics from the generated ASCII OpenFOAM cases."""

from itertools import islice
import math
import re
from statistics import median


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
COMMENTS = re.compile(r"/\*.*?\*/|//[^\n]*", re.S)
VECTOR = re.compile(rf"\(\s*({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*\)")


def stripped(path):
    content = path.read_bytes()
    if re.search(rb"\bformat\s+binary\s*;", content[:1000]):
        raise ValueError(f"{path.name} is binary. Use ASCII/reconstruct before analysis.")
    return COMMENTS.sub("", content.decode("utf-8"))


def block(text, name):
    start = re.search(r'(?m)(?<![\w:])"?' + re.escape(name) + r'"?\s*\{', text)
    if not start:
        raise ValueError(f"Missing {name} block")
    depth = 1
    for token in re.finditer(r"[{}]", text[start.end():]):
        depth += 1 if token[0] == "{" else -1
        if not depth:
            return text[start.end():start.end() + token.start()]
    raise ValueError(f"Unclosed {name} block")


def foam_list(path):
    data = stripped(path)
    header = re.search(r"\bFoamFile\s*\{[^}]*\}", data, re.S)
    if header is None:
        raise ValueError(f"Missing FoamFile header in {path.name}")
    body = re.match(r"\s*(\d+)\s*\(", data[header.end():])
    if body is None:
        raise ValueError(f"Missing counted list in {path.name}")
    return int(body[1]), data[header.end() + body.end():]


def patch_values(path, count, vector=False, owners=None):
    text = stripped(path)
    patch = block(block(text, "boundaryField"), "airfoil")
    uniform = re.search(r"\bvalue\s+uniform\s+([^;]+);", patch)
    if uniform:
        if vector:
            match = VECTOR.fullmatch(uniform[1].strip())
            if match is None:
                raise ValueError(f"Unexpected vector in {path.name}")
            value = tuple(map(float, match.groups()))
        else:
            value = float(uniform[1])
        return [value] * count

    match = re.search(r"\bvalue\s+nonuniform\s+\S+\s+(\d+)\s*\((.*?)\)\s*;", patch, re.S)
    if match is None and owners is not None and re.search(r"\btype\s+zeroGradient\s*;", patch):
        # OpenFOAM omits patch values for zeroGradient pressure. The boundary
        # value is the adjacent cell value indexed by the global face owner list.
        internal = re.search(r"\binternalField\s+uniform\s+([^;]+);", text)
        if internal:
            return [float(internal[1])] * count
        internal = re.search(r"\binternalField\s+nonuniform\s+\S+\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
        if internal:
            values = [float(item) for item in internal[2].split()]
            if len(values) == int(internal[1]) and len(owners) == count and max(owners) < len(values):
                return [values[owner] for owner in owners]
    if match is None or int(match[1]) != count:
        raise ValueError(f"Expected {count} airfoil face values in {path.name}")
    if vector:
        values = [tuple(map(float, item)) for item in VECTOR.findall(match[2])]
    else:
        values = [float(item) for item in match[2].split()]
    if len(values) != count:
        raise ValueError(f"Read {len(values)} face values in {path.name}; expected {count}")
    return values


def airfoil_faces(mesh_dir):
    boundary = block(stripped(mesh_dir / "boundary"), "airfoil")
    start = int(re.search(r"\bstartFace\s+(\d+)\s*;", boundary)[1])
    count = int(re.search(r"\bnFaces\s+(\d+)\s*;", boundary)[1])
    _, vertices = foam_list(mesh_dir / "points")
    points = [tuple(map(float, row)) for row in VECTOR.findall(vertices)]
    _, faces = foam_list(mesh_dir / "faces")
    face_rows = re.finditer(r"\d+\s*\(([\d\s]+)\)", faces)
    coordinates = [
        [points[int(index)] for index in face[1].split()]
        for face in islice(face_rows, start, start + count)
    ]
    if len(coordinates) != count or any(len(face) != 4 for face in coordinates):
        raise ValueError("Expected a complete airfoil patch of extruded quadrilaterals")
    owners_path = mesh_dir / "owner"
    owners = None
    if owners_path.exists():
        _, owner_body = foam_list(owners_path)
        owners = [int(entry[0]) for entry in islice(re.finditer(r"\d+", owner_body), start, start + count)]
    return coordinates, owners


def latest_time(case):
    times = [path for path in case.iterdir()
             if path.is_dir() and re.fullmatch(NUMBER, path.name)
             and float(path.name) > 0 and (path / "p").is_file()]
    return max(times, key=lambda path: float(path.name), default=None)


def find_field(time, field):
    path = time / field
    if path.is_file():
        return path
    matches = sorted(time.glob(f"*:{field}"))
    return matches[0] if matches else None


def surface_rows(case, meta):
    """Return the reconstructed latest-time patch as one row per cell face."""
    time = latest_time(case)
    if time is None:
        return [], ["No reconstructed solution time containing p."]
    faces, owners = airfoil_faces(case / "constant" / "polyMesh")
    count = len(faces)
    p = patch_values(time / "p", count, owners=owners)
    missing = []
    optional = {}
    for name, vector in [("yPlus", False), ("wallShearStress", True)]:
        field = find_field(time, name)
        if not field:
            missing.append(f"No {name} field at time {time.name}.")
            optional[name] = [None] * count
            continue
        try:
            optional[name] = patch_values(field, count, vector=vector)
        except ValueError as error:
            missing.append(str(error))
            optional[name] = [None] * count

    chord = meta["chord"]
    inv_half_q = 2 / meta["speed"] ** 2
    rows = []
    for face, pressure, y_plus, shear in zip(
            faces, p, optional["yPlus"], optional["wallShearStress"]):
        x = sum(point[0] for point in face) / 4
        y = sum(point[1] for point in face) / 4
        a = min(face, key=lambda point: point[0])
        b = max(face, key=lambda point: point[0])
        dx, dy = b[0] - a[0], b[1] - a[1]
        arc = math.hypot(dx, dy)
        if not arc:
            raise ValueError("Degenerate face found on the airfoil")
        projection = (shear[0] * dx + shear[1] * dy) / arc if shear is not None else None
        rows.append({
            "time": float(time.name),
            "side": "upper" if y >= 0 else "lower",
            "x_over_c": x / chord,
            "y_over_c": y / chord,
            "arc_over_c": arc / chord,
            "cp": pressure * inv_half_q,
            "y_plus": y_plus,
            "raw_cf": projection * inv_half_q if projection is not None else None,
        })
    # Infer the traction sign from the attached mid-chord pressure side.
    pressure_side = [row["raw_cf"] for row in rows
                     if row["side"] == "lower" and 0.25 < row["x_over_c"] < 0.65
                     and row["raw_cf"] is not None]
    sign = -1 if pressure_side and median(pressure_side) < 0 else 1
    for row in rows:
        raw = row.pop("raw_cf")
        row["cf"] = sign * raw if raw is not None else None
    return sorted(rows, key=lambda row: (row["side"], row["x_over_c"])), missing


def surface_metrics(rows, wall_treatment):
    if not rows:
        return {}
    values = {"surface_time": rows[0]["time"]}
    wall_function = wall_treatment == "function"
    low, high = (30.0, 100.0) if wall_function else (0.0, 1.0)
    for side in ("upper", "lower"):
        surface = [row for row in rows if row["side"] == side]
        measured = [row for row in surface if row["y_plus"] is not None]
        if measured:
            values[f"y_plus_max_{side}"] = max(row["y_plus"] for row in measured)
            in_range = sum(row["arc_over_c"] for row in measured if low <= row["y_plus"] <= high)
            values[f"y_plus_in_range_{side}_pct"] = 100 * in_range / sum(
                row["arc_over_c"] for row in measured)
        attached = [row for row in measured if row["cf"] is not None and row["cf"] > 0]
        if attached:
            in_range = sum(row["arc_over_c"] for row in attached if low <= row["y_plus"] <= high)
            values[f"attached_y_plus_in_range_{side}_pct"] = 100 * in_range / sum(
                row["arc_over_c"] for row in attached)

    upper_surface = [row for row in rows if row["side"] == "upper"]
    if upper_surface:
        values["cp_min_upper"] = min(row["cp"] for row in upper_surface)
    upper = [row for row in upper_surface if row["cf"] is not None]
    sep = None
    intervals = []
    for left, right in zip(upper, upper[1:]):
        a, b = left["cf"], right["cf"]
        if a >= 0 > b or a < 0 <= b:
            root = left["x_over_c"] - a * (right["x_over_c"] - left["x_over_c"]) / (b - a)
            if a >= 0 > b:
                sep = root
            elif sep is not None:
                intervals.append((sep, root))
                sep = None
    if sep is not None:
        intervals.append((sep, None))
    bubble = next((pair for pair in intervals if 0.001 < pair[0] < 0.15), None)
    if bubble:
        values["x_sep_over_c"], values["x_reatt_over_c"] = bubble
    if intervals:
        values["upper_separation_intervals"] = "|".join(
            f"{a:.6g} to {b:.6g}" if b is not None else f"{a:.6g} to TE" for a, b in intervals
        )
    return values

#!/usr/bin/env python3
"""Prepare and run the NACA 0012 case matrix in an OpenFOAM shell."""

import argparse
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time


PROJECT = Path(__file__).resolve().parent
SOLVERS = {"openfoam-com": "simpleFoam", "foundation": "foamRun"}
SMOKE = dict(group="smoke", chord=1.0, speed=52.08, reynolds=6_000_000,
             alpha=0.0, radius=25.0, level="coarse", model="kOmegaSST",
             wall_treatment="resolved", tu=None, length_scale=None,
             iterations=100, transition_mesh=False)


def load_cases():
    cases = {"smoke_alpha0": SMOKE}
    groups = {"smoke": ("smoke_alpha0",)}
    for path in sorted((PROJECT / "studies").glob("*/study.json")):
        study = json.loads(path.read_text(encoding="utf-8-sig"))
        definitions = study["cases"]
        for name, overrides in definitions.items():
            if name in cases:
                raise ValueError(f"{path}: duplicate case {name!r}.")
            cases[name] = study["defaults"] | overrides
        for name, members in study["groups"].items():
            if name in groups or name == "all":
                raise ValueError(f"{path}: duplicate or reserved group {name!r}.")
            if set(members) - definitions.keys():
                raise ValueError(f"{path}: group {name!r} must use cases defined in that file.")
            groups[name] = tuple(members)
    groups["all"] = tuple(name for name in cases if name != "smoke_alpha0")
    return cases, groups


CASES, GROUPS = load_cases()


def save_json(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def select_flavor(requested, require_environment):
    version = os.environ.get("WM_PROJECT_VERSION", "")
    detected = None
    if re.fullmatch(r"v?\d{4}(?:[.\-].*)?", version):
        detected = "openfoam-com"
    elif re.fullmatch(r"\d{1,2}", version):
        detected = "foundation"
    elif requested == "auto" and require_environment:
        found = [name for name, exe in SOLVERS.items() if shutil.which(exe)]
        if len(found) == 1:
            detected = found[0]
        else:
            raise ValueError("Load OpenFOAM first or select --flavor explicitly.")
    if requested != "auto" and detected and requested != detected:
        raise ValueError(f"Loaded OpenFOAM {version} does not match --flavor {requested}.")
    return detected or ("openfoam-com" if requested == "auto" else requested)


def case_spec(name, args, flavor):
    spec = dict(CASES[name])
    for key in ("radius", "first_height", "tu", "length_scale"):
        value = getattr(args, key)
        if value is not None:
            if key in ("tu", "length_scale") and spec["tu"] is None:
                raise ValueError(f"--{key.replace('_', '-')} requires a case with prescribed inlet turbulence.")
            spec[key] = value
    group = spec["group"]
    iterations = args.iterations or spec["iterations"]
    default_warmup = min(75, iterations * 3 // 4) if group == "smoke" else min(500, iterations // 4)
    warmup = args.warmup if args.warmup is not None else default_warmup
    if not (0 <= warmup < iterations):
        raise ValueError("--warmup must be smaller than --iterations and at least zero.")
    if not 0 < spec["radius"] or (spec["tu"] is not None and not 0 < spec["tu"]):
        raise ValueError("Radius and any specified turbulence intensity must be positive.")
    spec.update(id=name, span=0.01 * spec["chord"], nu=spec["chord"] * spec["speed"] / spec["reynolds"],
                flavor=flavor, nprocs=args.ranks,
                iterations=iterations, warmup=warmup, write_interval=args.write_interval or 250,
                mesh_backend=args.mesh_backend or "auto", initialization=args.initialization or "potential")
    return spec


def prepare(work_dir, spec):
    from foam.dictionaries import write_case

    target = work_dir / spec["id"]
    if target.exists():
        raise FileExistsError(f"{target} already exists. Use a different --work-dir for a fresh case.")
    with tempfile.TemporaryDirectory(prefix=".prepare-", dir=work_dir) as temp_dir:
        case = Path(temp_dir) / spec["id"]
        case.mkdir()
        write_case(case, spec, flavor=spec["flavor"], nprocs=spec["nprocs"],
                   end_time=spec["iterations"], write_interval=spec["write_interval"])
        mesh_command = [sys.executable, str(PROJECT / "mesh" / "generate.py"),
                        "--output", str(case / "airfoil.msh"), "--metadata", str(case / "mesh.json"),
                        "--chord", str(spec["chord"]), "--radius", str(spec["radius"]),
                        "--level", spec["level"], "--wall", spec["wall_treatment"],
                        "--alpha", str(spec["alpha"]), "--backend", spec["mesh_backend"]]
        if spec["transition_mesh"]:
            mesh_command.append("--transition")
        if "first_height" in spec:
            mesh_command.extend(("--first-height", str(spec["first_height"])))
        run_logged(mesh_command, case, "log.mesh", env=os.environ.copy())
        save_json(case / "case.json", spec)
        case.rename(target)
    return target


@contextmanager
def case_lock(work_dir, name):
    """Lock a single case without leaving stale locks after a killed Slurm job."""
    try:
        import fcntl
    except ImportError:
        # Windows can prepare cases. The final rename also protects existing cases.
        yield
        return
    with (work_dir / f".{name}.lock").open("a", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another job is already working on {name}.") from exc
        yield


def run_logged(command, case, name, env):
    log = case / name
    print(f"[{case.name}] {' '.join(command)} -> {name}", flush=True)
    with log.open("a", encoding="utf-8") as output:
        output.write(f"\n# {' '.join(command)}\n")
        output.flush()
        result = subprocess.run(command, cwd=case, env=env, stdout=output, stderr=subprocess.STDOUT)
    if result.returncode:
        lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
        raise RuntimeError(f"{command[0]} exited with {result.returncode}. See {log}.\n" + "\n".join(lines[-16:]))


def set_patches(path):
    contents = path.read_text(encoding="utf-8")
    expected = {"airfoil": "wall", "farfield": "patch", "frontAndBack": "empty"}
    found = set()

    def update(match):
        name, body = match.groups()
        if name not in expected:
            if re.search(r"\bnFaces\s+(?!0\s*;)\d+\s*;", body):
                raise RuntimeError(f"Gmsh imported unexpected nonempty boundary patch {name}.")
            return match[0]
        if not re.search(r"\bnFaces\s+[1-9]\d*\s*;", body):
            raise RuntimeError(f"Gmsh patch {name} has no faces.")
        adjusted, number = re.subn(r"(\btype\s+)\w+(\s*;)", rf"\g<1>{expected[name]}\2", body)
        if number != 1:
            raise RuntimeError(f"Could not set patch type for {name} in {path}.")
        found.add(name)
        return match[0].replace(body, adjusted)

    contents = re.sub(r"^[ \t]*(\w+)\s*\{([^{}]*)\}", update, contents, flags=re.MULTILINE)
    if found != set(expected):
        raise RuntimeError(f"Missing imported patches in {path}: {', '.join(sorted(expected.keys() - found))}")
    path.write_text(contents, encoding="utf-8", newline="\n")


def checkpoint(case, ranks, fields):
    workspaces = [case / f"processor{rank}" for rank in range(ranks)] if ranks > 1 else [case]
    time_sets = []
    for directory in workspaces:
        if not directory.is_dir():
            raise RuntimeError(f"Missing MPI directory {directory}. Use a fresh run directory to rebuild decomposition.")
        time_sets.append({int(path.name) for path in directory.iterdir()
                          if path.is_dir() and path.name.isdecimal()
                          and all((path / field).is_file() for field in fields)})
    common = set.intersection(*time_sets)
    if not common:
        raise RuntimeError("No complete checkpoint was found across the requested ranks.")
    return max(common)


def check_imported_mesh(case, environment):
    basic_log = case / "log.checkMesh.basic"
    basic_log.write_text("", encoding="utf-8")
    run_logged(["checkMesh", "-allTopology"], case, basic_log.name, environment)
    basic = basic_log.read_text(encoding="utf-8", errors="replace")
    if not re.search(r"\bMesh OK\.", basic):
        raise RuntimeError(f"Standard mesh checks did not pass. Inspect {basic_log}.")
    if not re.search(r"Mesh has\s+2\s+solution\b", basic):
        raise RuntimeError(f"The mesh does not have exactly two solution directions. Inspect {basic_log}.")

    audit_log = case / "log.checkMesh"
    audit_log.write_text("", encoding="utf-8")
    run_logged(["checkMesh", "-allTopology", "-allGeometry"], case, audit_log.name, environment)
    audit = audit_log.read_text(encoding="utf-8", errors="replace")
    if re.search(r"\bMesh OK\.", audit):
        return "passed"
    if re.search(r"Failed\s+1\s+mesh checks?\.", audit) and re.search(r"Cells with small determinant", audit):
        print(f"[{case.name}] Extended geometry audit flagged only the cell determinant. "
              f"Check stretched layers in {audit_log}.", flush=True)
        return "cell_determinant_only"
    raise RuntimeError(f"Extended mesh checks did not pass. Inspect {audit_log}.")


def set_stage(case, start, end, interval, stage):
    schemes = case / "system" / "fvSchemes"
    shutil.copyfile(schemes.with_suffix(f".{stage}"), schemes)
    control = case / "system" / "controlDict"
    # Make the end of each stage a regular runTime write. Include start in the
    # divisor so both absolute- and start-relative time origins hit the endpoint.
    values = {"startFrom": "startTime", "startTime": start, "endTime": end,
              "writeInterval": math.gcd(interval, start, end)}
    text = control.read_text(encoding="utf-8")
    for key, value in values.items():
        text, count = re.subn(rf"(?m)^({key}[ \t]+)[^;\n]+;", rf"\g<1>{value};", text)
        if count != 1:
            raise RuntimeError(f"Expected one {key} entry in {control}.")
    control.write_text(text, encoding="utf-8", newline="\n")


def run_case(case, spec, args):
    ranks = spec["nprocs"]
    solver = SOLVERS[spec["flavor"]]
    environment = dict(os.environ, FOAM_FILEHANDLER="uncollated", OMP_NUM_THREADS="1")
    status_file = case / "run-status.json"
    status = json.loads(status_file.read_text(encoding="utf-8")) if status_file.exists() else {}
    wall_seconds = float(status.get("wall_seconds", 0))
    history = bool(status.get("solve_started")) or (case / "log.solve").exists() or (case / "log.warmup").exists()
    if history and not args.resume and not args.mesh_only:
        raise RuntimeError(f"{case} has solver history. Use --resume to extend it or a fresh --work-dir.")
    if args.mesh_only and status.get("mesh_checked"):
        print(f"[{case.name}] Mesh already imported and checked.")
        return

    def record(state, exit_code=None):
        status.update(status=state, exit_code=exit_code, wall_seconds=wall_seconds, nprocs=ranks,
                      openfoam_version=environment.get("WM_PROJECT_VERSION"))
        save_json(status_file, status)

    try:
        record("preparing")
        if not status.get("mesh_checked"):
            run_logged(["gmshToFoam", "airfoil.msh"], case, "log.gmshToFoam", environment)
            set_patches(case / "constant" / "polyMesh" / "boundary")
            status["mesh_geometry_audit"] = check_imported_mesh(case, environment)
            status["mesh_checked"] = True
        if args.mesh_only:
            record("mesh_checked", 0)
            return

        if not history and spec["initialization"] == "potential" and not status.get("potential_initialized"):
            set_stage(case, 0, spec["iterations"], spec["write_interval"], "final")
            potential = ["potentialFoam", "-initialiseUBCs", "-writep"]
            if spec["flavor"] == "openfoam-com":
                potential.extend(("-writephi", "-noFunctionObjects"))
            record("potential_initialization")
            begin = time.monotonic()
            run_logged(potential, case, "log.potentialFoam", environment)
            status["init_wall_seconds"] = time.monotonic() - begin
            status["potential_initialized"] = True
            record("potential_initialization")

        if ranks > 1 and not status.get("decomposed"):
            run_logged(["decomposePar", "-force"], case, "log.decomposePar", environment)
            status["decomposed"] = True
        fields = [p.name for p in (case / "0").iterdir() if p.is_file() and p.name.lower() != "phi"]
        start = checkpoint(case, ranks, fields)
        if start > spec["iterations"]:
            raise ValueError(f"The checkpoint at {start} is beyond --iterations {spec['iterations']}.")
        save_json(case / "case.json", spec)
        launch = [solver]
        if ranks > 1:
            if args.launcher == "srun":
                mpi_flags = [f"--mpi={args.mpi_type}"] if args.mpi_type else []
                launch = ["srun", "--ntasks", str(ranks), *mpi_flags, solver, "-parallel"]
            else:
                launch = ["mpirun", "-np", str(ranks), solver, "-parallel"]

        for stage, end, log_name in (("startup", spec["warmup"], "log.warmup"),
                                     ("final", spec["iterations"], "log.solve")):
            if end <= start:
                continue
            set_stage(case, start, end, spec["write_interval"], stage)
            status["solve_started"] = True
            record(stage)
            begin = time.monotonic()
            try:
                run_logged(launch, case, log_name, environment)
            finally:
                wall_seconds += time.monotonic() - begin
            start = checkpoint(case, ranks, fields)
            if start != end:
                raise RuntimeError(f"Expected checkpoint {end} after {stage}. Latest complete checkpoint is {start}.")
            status["latest_iteration"] = start
            record(stage)

        if ranks > 1:
            run_logged(["reconstructPar", "-latestTime"], case, "log.reconstructPar", environment)
        (case / "airfoil.foam").touch()
        status.pop("error", None)
        record("finished", 0)
    except (OSError, RuntimeError, ValueError) as exc:
        status["error"] = str(exc).splitlines()[0]
        record("failed", 1)
        raise


def positive(value):
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("expected a positive finite value")
    return parsed


def positive_int(value):
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return parsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="list the matrix or selected Slurm array group")
    listing.add_argument("--group", choices=GROUPS, default="all")
    listing.add_argument("--names", action="store_true", help="emit only case IDs")
    for command in ("prepare", "run"):
        sub = commands.add_parser(command, help="generate case inputs" if command == "prepare" else "generate and run an OpenFOAM case")
        sub.add_argument("case", choices=CASES, nargs="+" if command == "prepare" else None)
        sub.add_argument("--work-dir", type=Path, default=Path("runs"))
        sub.add_argument("--flavor", choices=["auto", *SOLVERS], default="auto")
        sub.add_argument("--ranks", type=positive_int, default=1)
        sub.add_argument("--iterations", type=positive_int, help="target total steady iterations, including startup")
        sub.add_argument("--warmup", type=int, help="upwind startup iterations (0 disables it)")
        sub.add_argument("--write-interval", type=positive_int, help="checkpoint interval, default 250")
        sub.add_argument("--mesh-backend", choices=("auto", "gmsh", "ascii"),
                         help="use the Gmsh Python API when installed (auto), require it (gmsh), or write MSH2 directly")
        sub.add_argument("--initialization", choices=("potential", "uniform"),
                         help="initialize velocity and pressure with potentialFoam by default")
        sub.add_argument("--radius", type=positive, help="override farfield radius in chords")
        sub.add_argument("--first-height", type=positive, help="override full first-cell height divided by chord")
        sub.add_argument("--tu", type=positive, help="override prescribed inlet turbulence intensity in percent")
        sub.add_argument("--length-scale", type=positive, help="override prescribed inlet length scale in metres")
        if command == "run":
            sub.add_argument("--resume", action="store_true", help="continue from the latest complete checkpoint")
            sub.add_argument("--mesh-only", action="store_true", help="import and validate the mesh without solving")
            sub.add_argument("--launcher", choices=["mpirun", "srun"], default="mpirun")
            sub.add_argument("--mpi-type", help="Slurm MPI plugin, e.g. pmix, when --launcher srun")
    args = parser.parse_args()
    if args.command == "list":
        for name in GROUPS[args.group]:
            spec = CASES[name]
            print(name if args.names else f"{name:24s} {spec['model']:16s} {spec['level']:6s} "
                                         f"{spec['wall_treatment']:8s} R={spec['radius']:g}c")
        return

    running = args.command == "run"
    flavor = select_flavor(args.flavor, require_environment=running and not getattr(args, "mesh_only", False))
    if running:
        required = ["gmshToFoam", "checkMesh"]
        if not args.mesh_only:
            required.append(SOLVERS[flavor])
            if args.ranks > 1:
                required.extend(("decomposePar", "reconstructPar", args.launcher))
        missing = [name for name in required if not shutil.which(name)]
        if missing:
            raise RuntimeError(f"Load the cluster OpenFOAM environment. Missing commands: {', '.join(missing)}")
    work_dir = args.work_dir.expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    for name in ([args.case] if running else args.case):
        with case_lock(work_dir, name):
            target = work_dir / name
            spec = case_spec(name, args, flavor)
            if running and target.exists():
                previous = json.loads((target / "case.json").read_text(encoding="utf-8"))
                options = {"iterations": args.iterations, "warmup": args.warmup,
                           "radius": args.radius, "first_height": args.first_height,
                           "tu": args.tu, "length_scale": args.length_scale,
                           "mesh_backend": args.mesh_backend, "write_interval": args.write_interval,
                           "initialization": args.initialization}
                for key, value in options.items():
                    if value is None and key in previous:
                        spec[key] = previous[key]
                for key in previous.keys() | spec.keys():
                    if key not in {"iterations", "write_interval"} and previous.get(key) != spec.get(key):
                        raise ValueError(f"Prepared case differs in {key}. Match its settings or use a fresh --work-dir.")
                if spec["iterations"] <= spec["warmup"]:
                    raise ValueError("Target iterations must exceed the prepared warmup length.")
            else:
                target = prepare(work_dir, spec)
            if running:
                if not args.mesh_only and not args.resume and spec["initialization"] == "potential" and not shutil.which("potentialFoam"):
                    raise RuntimeError("potentialFoam is missing. Load the complete OpenFOAM environment. "
                                       "To use a uniform start, prepare in a fresh --work-dir with --initialization uniform.")
                run_case(target, spec, args)
            print(target, flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError) as error:
        sys.exit(f"Error: {error}")

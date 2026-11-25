"""Run one low-Re model on the shared, prechecked medium blockMesh grid."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("model", choices=("SST", "SST-LM"))
    parser.add_argument("iteration", type=int, choices=(1500, 7000))
    args = parser.parse_args()

    sys.path.insert(0, str(args.work / "code"))
    import study
    from foam.dictionaries import write_case

    case_id = "low_sst_016" if args.model == "SST" else "low_lm_016"
    model = "kOmegaSST" if args.model == "SST" else "kOmegaSSTLM"
    case = args.work / "cases" / case_id
    resume = args.iteration == 7000
    if resume:
        spec = json.loads((case / "case.json").read_text())
    else:
        baseline_spec = args.baseline / "cases/low_sst_016/case.json"
        spec = json.loads(baseline_spec.read_text())
        ranks = int(os.environ["SLURM_NTASKS"])
        if ranks != spec["nprocs"]:
            raise RuntimeError(f"Use {spec['nprocs']} MPI ranks, matching the prior medium run")
        spec.update(
            id=case_id, model=model, mesh_backend="blockMesh",
            iterations=args.iteration,
        )
        status = json.loads((args.mesh / "mesh-status.json").read_text())
        if status.get("mesh_checked") is not True:
            raise RuntimeError("The common mesh was not checked")
        case.mkdir(parents=True)
        write_case(
            case, spec, flavor=spec["flavor"], nprocs=ranks,
            end_time=args.iteration, write_interval=spec["write_interval"],
        )
        shutil.copytree(
            args.mesh / "constant/polyMesh", case / "constant/polyMesh",
            copy_function=shutil.copyfile, ignore=shutil.ignore_patterns("sets"),
        )
        shutil.copyfile(args.mesh / "system/blockMeshDict", case / "system/blockMeshDict")
        for filename in ("log.blockMesh", "log.checkMesh", "log.checkMesh.basic", "mesh.json"):
            shutil.copyfile(args.mesh / filename, case / filename)
        study.save_json(case / "run-status.json", status)
        for log in ("log.warmup", "log.solve"):
            (case / log).symlink_to(args.output / log)

    spec["iterations"] = args.iteration
    options = SimpleNamespace(resume=resume, mesh_only=False, launcher="mpirun", mpi_type=None)
    study.run_case(case, spec, options)


if __name__ == "__main__":
    main()

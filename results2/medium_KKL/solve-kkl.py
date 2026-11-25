import json
import os
from pathlib import Path
import re
import shutil
import sys
from types import SimpleNamespace

baseline, mesh, work, output = map(Path, sys.argv[1:5])
layout, step = sys.argv[5], int(sys.argv[6])
sys.path.insert(0, str(work / 'code'))
import study
from foam.dictionaries import write_case

case = work / 'cases/low_kkl_016'
if step == 1500:
    spec = json.loads((baseline / 'cases/low_sst_016/case.json').read_text())
    assert (spec['level'], spec['reynolds'], spec['alpha'], spec['tu'], spec['nprocs']) == (
        'medium', 130000, 10, 0.16, int(os.environ['SLURM_NTASKS']))
    spec.update(id='low_kkl_016', model='kkLOmega', iterations=step,
                mesh_backend='blockMesh' if layout == 'blockMesh' else 'gmsh',
                kkl_scalar_convection='bounded Gauss upwind', kkl_turbulence_relaxation=0.15)
    assert json.loads((mesh / 'mesh.json').read_text())['cells'] == 209_952
    status_name = 'mesh-status.json' if layout == 'blockMesh' else 'run-status.json'
    audit = json.loads((mesh / status_name).read_text())
    assert audit['mesh_checked'] is True
    case.mkdir(parents=True)
    write_case(case, spec, flavor=spec['flavor'], nprocs=spec['nprocs'],
               end_time=step, write_interval=spec['write_interval'])
    shutil.copytree(mesh / 'constant/polyMesh', case / 'constant/polyMesh',
                    copy_function=shutil.copyfile, ignore=shutil.ignore_patterns('sets'))
    for rel in ('mesh.json', 'log.checkMesh', 'log.checkMesh.basic',
                'log.blockMesh', 'log.mesh', 'system/blockMeshDict'):
        if (mesh / rel).is_file():
            shutil.copyfile(mesh / rel, case / rel)
    # The previous scalar-upwind diagnostic avoided an early model crash.
    # Apply identical convection and relaxation settings on both grids.
    schemes = case / 'system/fvSchemes.final'
    text, count = re.subn(r'(div\(phi,(?:kt|kl|omega)\) bounded Gauss )limitedLinear 1',
                          r'\g<1>upwind', schemes.read_text())
    assert count == 3
    schemes.write_text(text)
    solution = case / 'system/fvSolution'
    old = '"(U|kt|kl|omega)" 0.35;'
    assert solution.read_text().count(old) == 1
    solution.write_text(solution.read_text().replace(old, 'U 0.35;\n        kt 0.15;\n        kl 0.15;\n        omega 0.15;'))
    study.save_json(case / 'run-status.json', {
        'mesh_checked': True, 'mesh_geometry_audit': audit['mesh_geometry_audit'],
    })
    for log in ('log.warmup', 'log.solve'):
        (case / log).symlink_to(output / log)
else:
    spec = json.loads((case / 'case.json').read_text())
spec['iterations'] = step
options = SimpleNamespace(resume=step != 1500, mesh_only=False, launcher='mpirun', mpi_type=None)
study.run_case(case, spec, options)

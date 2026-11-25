"""Report wall-shear metrics from one reconstructed iteration."""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

root = Path(sys.argv[1]).resolve()
iteration = sys.argv[2]
sys.path.insert(0, str(root / 'code' / 'analysis'))
from surface import surface_metrics, surface_rows

case = root / 'cases' / 'low_sst_016'
config = json.loads((case / 'case.json').read_text())
with TemporaryDirectory() as scratch:
    work = Path(scratch)
    for name in ('constant', iteration):
        (work / name).symlink_to(case / name, target_is_directory=True)
    rows, notes = surface_rows(work, config)
    metrics = surface_metrics(rows, config['wall_treatment'])

for name in ('surface_time', 'x_sep_over_c', 'x_reatt_over_c', 'y_plus_max_upper',
             'upper_separation_intervals'):
    print(f'{name}: {metrics.get(name)}')
for note in notes:
    print(note)

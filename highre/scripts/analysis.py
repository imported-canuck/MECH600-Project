import csv
import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


CASES = ['sst_L0', 'sst_L1', 'sst_L2', 'sst_domain', 'sa_L2', 'sst_wf']
CELLS = {'sst_L0': 36000, 'sst_L1': 81000, 'sst_L2': 181980,
         'sst_domain': 81000, 'sa_L2': 181980, 'sst_wf': 181980}
WINDOW = 1000 
R = 1.5

CL_EXP = 1.067
CD_EXP = 0.0117


def read_coeffs(case):
    data = np.loadtxt(f'{case}/postProcessing/forceCoeffs/0/coefficient.dat', comments='#')
    return data[:, 0], data[:, 1], data[:, 4]


def window_stats(values):
    last = values[-WINDOW:]
    return last.mean(), 0.5 * (last.max() - last.min())


def gci(f1, f2, f3, r):
    e21 = f2 - f1
    e32 = f3 - f2
    p = math.log(abs(e32 / e21)) / math.log(r)
    f_ext = (r**p * f1 - f2) / (r**p - 1)
    gci_fine = 1.25 * abs((f1 - f2) / f1) / (r**p - 1)
    monotonic = (e32 / e21) > 0
    return p, f_ext, gci_fine, monotonic


if __name__ == '__main__':
    results = {}

    print(f"{'case':11s} {'cells':>7s} {'Cl avg':>8s} {'+/-':>7s} {'Cd avg':>8s} {'+/-':>8s} "
          f"{'Cl err%':>8s} {'Cd err%':>8s}")
    for case in CASES:
        it, cd, cl = read_coeffs(case)
        cl_avg, cl_amp = window_stats(cl)
        cd_avg, cd_amp = window_stats(cd)
        cl_err = 100 * (cl_avg - CL_EXP) / CL_EXP
        cd_err = 100 * (cd_avg - CD_EXP) / CD_EXP
        results[case] = dict(cl=cl_avg, cl_amp=cl_amp, cd=cd_avg, cd_amp=cd_amp,
                             cl_err=cl_err, cd_err=cd_err)
        print(f'{case:11s} {CELLS[case]:7d} {cl_avg:8.4f} {cl_amp:7.4f} {cd_avg:8.5f} {cd_amp:8.5f} '
              f'{cl_err:8.2f} {cd_err:8.2f}')
    print(f"{'Ladson':11s} {'':7s} {CL_EXP:8.4f} {'':7s} {CD_EXP:8.5f}")

    print('\nGrid study (SST, wall-resolved, averages)')
    for name in ['cl', 'cd']:
        f1 = results['sst_L2'][name]
        f2 = results['sst_L1'][name]
        f3 = results['sst_L0'][name]
        p, f_ext, g, mono = gci(f1, f2, f3, R)
        print(f'  {name}: coarse {f3:.5f}  medium {f2:.5f}  fine {f1:.5f}  '
              f'p = {p:.2f}  extrapolated {f_ext:.5f}  GCI_fine = {100 * g:.2f} %  monotonic {mono}')

    print('\nComparisons (difference in %)')
    pairs = [('sst_domain', 'sst_L1', 'domain 200c vs 100c'),
             ('sa_L2', 'sst_L2', 'SA vs SST'),
             ('sst_wf', 'sst_L2', 'wall functions vs resolved')]
    for a, b, label in pairs:
        dcl = 100 * (results[a]['cl'] - results[b]['cl']) / results[b]['cl']
        dcd = 100 * (results[a]['cd'] - results[b]['cd']) / results[b]['cd']
        print(f'  {label:28s} Cl {dcl:+6.2f} %   Cd {dcd:+6.2f} %')

    with open('results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case', 'cells', 'cl_avg', 'cl_amp', 'cd_avg', 'cd_amp', 'cl_err_pct', 'cd_err_pct'])
        for case in CASES:
            r = results[case]
            w.writerow([case, CELLS[case], r['cl'], r['cl_amp'], r['cd'], r['cd_amp'], r['cl_err'], r['cd_err']])

    fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    for case in CASES:
        it, cd, cl = read_coeffs(case)
        keep = it > it[-1] - 4000
        ax[0].plot(it[keep], cl[keep], lw=0.8, label=case)
        ax[1].plot(it[keep], cd[keep], lw=0.8, label=case)
    ax[0].axhline(CL_EXP, color='k', ls='--', lw=1, label='Ladson')
    ax[1].axhline(CD_EXP, color='k', ls='--', lw=1, label='Ladson')
    ax[0].set_ylabel('$C_L$')
    ax[1].set_ylabel('$C_D$')
    ax[1].set_xlabel('iteration')
    ax[0].legend(fontsize=7, ncol=4)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('force_history.png', dpi=200)
    print('\nwrote results.csv and force_history.png')
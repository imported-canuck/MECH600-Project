import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analysis import CASES, CELLS, WINDOW, R, CL_EXP, CD_EXP, read_coeffs, window_stats, gci


def mean_uncertainty(values, n_blocks=10):
    blocks = values[-WINDOW:].reshape(n_blocks, -1).mean(axis=1)
    return blocks.std(ddof=1) / math.sqrt(n_blocks)


def read_yplus(case):
    last = None
    for line in open(f'{case}/postProcessing/yPlus/0/yPlus.dat'):
        if line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 5 and parts[1] == 'airfoil':
            last = parts
    return float(last[2]), float(last[3]), float(last[4])


if __name__ == '__main__':
    res = {}
    for case in CASES:
        it, cd, cl = read_coeffs(case)
        cl_avg, cl_amp = window_stats(cl)
        cd_avg, cd_amp = window_stats(cd)
        res[case] = dict(cl=cl_avg, cl_amp=cl_amp, cd=cd_avg, cd_amp=cd_amp,
                         cl_unc=mean_uncertainty(cl), cd_unc=mean_uncertainty(cd))
    print(f"{'case':11s} {'Cl avg':>8s} {'+/- avg':>8s} {'Cd avg':>9s} {'+/- avg':>9s}")
    for case in CASES:
        r = res[case]
        print(f"{case:11s} {r['cl']:8.4f} {r['cl_unc']:8.5f} {r['cd']:9.6f} {r['cd_unc']:9.6f}")

    print(f"{'case':11s} {'y+ min':>8s} {'y+ max':>8s} {'y+ avg':>8s}")
    for case in CASES:
        ymin, ymax, yavg = read_yplus(case)
        print(f'{case:11s} {ymin:8.3f} {ymax:8.3f} {yavg:8.3f}')

    grid = ['sst_L0', 'sst_L1', 'sst_L2']
    h = [math.sqrt(CELLS['sst_L2'] / CELLS[c]) for c in grid]     # cell size relative to fine mesh
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for j, (name, label, exp) in enumerate([('cl', '$C_L$', CL_EXP), ('cd', '$C_D$', CD_EXP)]):
        vals = [res[c][name] for c in grid]
        uncs = [res[c][name + '_unc'] for c in grid]
        p, f_ext, g, mono = gci(vals[2], vals[1], vals[0], R)
        a = ax[j]
        a.errorbar(h, vals, yerr=uncs, fmt='o-', capsize=4, label='SST (1000-iteration average)')
        a.plot(0, f_ext, 's', color='C3', label=f'Richardson extrapolation (p = {p:.2f})')
        a.errorbar(0.93, vals[2], yerr=g * vals[2], fmt='none', ecolor='C3', capsize=8, lw=2,
                   label=f'GCI fine = {100 * g:.1f} %')
        a.axhline(exp, color='k', ls='--', lw=1, label='Ladson (1988)')
        a.set_xlabel('relative cell size $h/h_{fine}$')
        a.set_ylabel(label)
        a.set_xlim(-0.1, 2.4)
        a.grid(alpha=0.3)
        a.legend(fontsize=7, loc='center right' if name == 'cl' else 'upper left')
    fig.suptitle('Grid convergence study (L2 fine, L1 medium, L0 coarse; r = 1.5)', fontsize=11)
    fig.tight_layout()
    fig.savefig('grid_study.png', dpi=200)
    plt.close(fig)
    print('wrote grid_study.png')

    labels = ['SST L0', 'SST L1', 'SST L2', 'SST 200c', 'SA L2', 'SST wall fn']
    x = np.arange(len(CASES))
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for j, (name, label, exp) in enumerate([('cl', '$C_L$', CL_EXP), ('cd', '$C_D$', CD_EXP)]):
        vals = [res[c][name] for c in CASES]
        uncs = [res[c][name + '_unc'] for c in CASES]
        a = ax[j]
        for i in range(len(CASES)):
            a.errorbar(x[i], vals[i], yerr=uncs[i], fmt='o', ms=8, capsize=4, color=f'C{i}')
        a.axhline(exp, color='k', ls='--', lw=1.2, label='Ladson (1988)')
        a.set_xticks(x)
        a.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
        a.set_ylabel(label)
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    fig.suptitle('All cases vs experiment (1000-iteration averages)', fontsize=11)
    fig.tight_layout()
    fig.savefig('cases_vs_ladson.png', dpi=200)
    plt.close(fig)
    print('wrote cases_vs_ladson.png')
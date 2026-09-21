"""Finite-receptor proofreading and a hypothetical cytokine signaling network."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIGURE = 'fig15_receptor_proofreading_cytokine_balance.png'
CONFIG = {'proofreading_steps': 5, 'kp_h': 720., 'kon_L_h': [90., 90.],
          'koff_h': [36., 360.], 'antigen_decay_start_h': 24., 'antigen_decay_h': 0.08,
          'antagonist_start_h': 12., 'antagonist_half_life_h': 48.,
          'antagonist_initial_C_over_KD': [0., 1., 10., 100.], 'horizon_h': 96.,
          'kinetics_provenance': 'all numerical kinetic values are scenario assumptions',
          'cytokine_unit': 'normalized model concentration; no patient calibration',
          'clinical_recommendations': None}

def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')

def write_csv(path, names, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(names); w.writerows(rows)

def steady_proofreading(kon=(90., 90.), koff=(36., 360.), kp=720., steps=5):
    """Closed ligand reservoirs; only the finite receptor pool is conserved."""
    if kp <= 0 or steps < 1 or any(x <= 0 for x in koff) or any(x < 0 for x in kon):
        raise ValueError('positive rates and a nonzero cascade are required')
    free = 1. / (1. + sum(a / b for a, b in zip(kon, koff)))
    blocks = []
    for on, off in zip(kon, koff):
        flux = on * free; ratio = kp / (kp + off)
        blocks.append([flux / (kp + off) * ratio ** j for j in range(steps)] +
                      [flux / off * ratio ** steps])
    return free, np.asarray(blocks)

def receptor_rhs(receptor, kon, koff=(36., 360.), kp=720.):
    d = np.zeros(13); r = receptor[0]
    for i in range(2):
        s = 1 + 6 * i; c = receptor[s:s+6]; on = kon[i] * r
        d[0] += -on + koff[i] * c.sum()
        d[s] = on - (kp + koff[i]) * c[0]
        d[s+1:s+5] = kp * c[:4] - (kp + koff[i]) * c[1:5]
        d[s+5] = kp * c[4] - koff[i] * c[5]
    return d

def blockade(t, ratio):
    x = ratio * np.exp(-np.log(2) * max(t - 12., 0.) / 48.) if t >= 12. else 0.
    return x / (1. + x)

def derivative(t, y, ratio):
    antigen = np.exp(-0.08 * max(t - 24., 0.))
    d = np.zeros(19)
    d[:13] = receptor_rhs(y[:13], (90. * antigen, 90.))
    tcell, macrophage, il6, tnf, ifng, il1 = y[13:]
    signal = il6 / (1. + il6) * (1. - blockade(t, ratio))
    trigger = (y[6] + y[12]) / (0.15 + y[6] + y[12])
    d[13] = 0.22 * trigger * (1. - tcell) - 0.06 * tcell
    d[14] = (0.35 * ifng / (1. + ifng) + 0.20 * signal) * (1. - macrophage) - 0.1 * macrophage
    # Blocking IL-6R reduces signaling, not the independent ligand production.
    d[15] = 0.18 * tcell + 1.1 * macrophage - (0.12 + 0.08 * (1. - blockade(t, ratio))) * il6
    d[16] = 0.55 * tcell + 0.6 * macrophage - 0.3 * tnf
    d[17] = 0.9 * tcell + 0.1 * macrophage - 0.25 * ifng
    d[18] = 0.08 * tcell + 0.8 * macrophage - 0.2 * il1
    return d

def simulate(ratio=0., rtol=2e-8):
    if ratio < 0: raise ValueError('antagonist ratio cannot be negative')
    initial = np.zeros(19); initial[0] = 1.
    tt = np.linspace(0., 96., 961)
    # Resolve the intervention discontinuity explicitly.
    first = solve_ivp(lambda t, y: derivative(t, y, 0.), (0., 12.), initial,
                      method='Radau', rtol=rtol, atol=rtol * 1e-2, dense_output=True)
    second = solve_ivp(lambda t, y: derivative(t, y, ratio), (12., 96.), first.y[:, -1],
                       method='Radau', rtol=rtol, atol=rtol * 1e-2, dense_output=True)
    if not first.success or not second.success: raise RuntimeError('cytokine integration failed')
    yy = np.empty((19, len(tt))); mask = tt < 12.
    yy[:, mask] = first.sol(tt[mask]); yy[:, ~mask] = second.sol(tt[~mask])
    sig = yy[15] / (1. + yy[15]) * np.array([1. - blockade(t, ratio) for t in tt])
    return tt, yy, sig

def run(out: Path) -> dict:
    out = Path(out)
    if out.exists(): raise FileExistsError(f'Output must be new: {out}')
    out.mkdir(parents=True); (out / 'figures').mkdir()
    write_json(out / 'config.json', CONFIG)
    records = []; runs = []; stats = []; masserr = 0.; minimum = 0.
    for ratio in CONFIG['antagonist_initial_C_over_KD']:
        t, y, sig = simulate(ratio); runs.append((t, y, sig))
        masserr = max(masserr, float(np.max(np.abs(y[:13].sum(axis=0) - 1.))))
        minimum = min(minimum, float(y.min()))
        for j, tj in enumerate(t): records.append([ratio, tj, *y[:, j], sig[j]])
        stats.append({'initial_C_over_KD': ratio, 'IL6_peak_normalized': float(y[15].max()),
                      'IL6R_signal_AUC_h': float(trapezoid(sig, t)),
                      'macrophage_peak_fraction': float(y[14].max())})
    names = ['antagonist_C_over_KD', 'time_h', 'free_receptor'] + [f'C_{a}_{i}' for a in ['target', 'self'] for i in range(6)]
    write_csv(out / 'cytokine_trajectories.csv', names + ['T_active', 'M_active', 'IL6', 'TNFa', 'IFNg', 'IL1b', 'IL6R_signal'], records)
    free, blocks = steady_proofreading()
    sensitivity = []
    for off in np.geomspace(3.6, 3600., 65):
        f, b = steady_proofreading(koff=(36., float(off)))
        sensitivity.append([off, b[0, -1], b[1, -1], b[0, -1] / b[1, -1]])
    write_csv(out / 'proofreading_sensitivity.csv', ['self_koff_h', 'target_C5', 'self_C5', 'signal_discrimination'], sensitivity)
    _, yr, sr = simulate(10., 2e-10)
    refine = float(np.max(np.abs(yr - runs[2][1])) / max(1., np.max(np.abs(yr))))
    verification = {'receptor_mass_error': masserr, 'minimum_state': minimum,
                    'refinement_max_scaled_error': refine,
                    'steady_receptor_mass_error': float(abs(free + blocks.sum() - 1.)),
                    'passed': bool(masserr < 1e-7 and minimum > -1e-8 and refine < 2e-5)}
    if not verification['passed']: raise RuntimeError(str(verification))
    summary = {'task': 15, 'evidence': 'uncalibrated mechanistic scenario', 'runs': stats,
               'target_to_self_C5_ratio': float(blocks[0, -1] / blocks[1, -1]),
               'clinical_recommendations': None, 'homeostasis_restoration_established': None,
               'limitations': ['No measured homeostatic set point or patient cytokine calibration.',
                               'Anti-IL6R exposure is dimensionless C/KD, not a tocilizumab dose.',
                               'Receptor mass conserved; antigen is an imposed reservoir.']}
    write_json(out / 'summary.json', summary); write_json(out / 'verification.json', verification)
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)
    ax[0,0].semilogy(range(6), blocks[0], 'o-', label='target'); ax[0,0].semilogy(range(6), blocks[1], 'o-', label='self')
    ax[0,0].set(xlabel='Proofreading state', ylabel='Receptor fraction', title='A  Finite receptor competition'); ax[0,0].legend()
    for k, label in enumerate(['IL-6', 'TNF-alpha', 'IFN-gamma', 'IL-1beta']): ax[0,1].plot(runs[0][0], runs[0][1][15+k], label=label)
    ax[0,1].set(xlabel='Time (h)', ylabel='Normalized concentration', title='B  Untreated cytokine network'); ax[0,1].legend()
    for ratio, (t,y,sig) in zip(CONFIG['antagonist_initial_C_over_KD'], runs):
        ax[1,0].plot(t, sig, label=f'C/KD={ratio:g}'); ax[1,1].plot(t, y[15], label=f'C/KD={ratio:g}')
    ax[1,0].set(xlabel='Time (h)', ylabel='IL-6R signal (fraction)', title='C  Receptor signal inhibition at 12 h'); ax[1,0].legend()
    ax[1,1].set(xlabel='Time (h)', ylabel='IL-6 normalized concentration', title='D  Ligand need not decrease on blockade'); ax[1,1].legend()
    for a in ax.flat: a.grid(alpha=.2)
    fig.suptitle('Task 15 | Proofreading and cytokine signaling — hypothetical, normalized signaling', fontsize=14)
    fig.savefig(out / 'figures' / FIGURE, dpi=300); plt.close(fig)
    return summary

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--out', type=Path, default=Path(__file__).parent / 'outputs')
    print(json.dumps(run(p.parse_args().out), indent=2))

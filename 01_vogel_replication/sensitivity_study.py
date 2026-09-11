"""
Sensitivity study for thesis: shows the progression of test cases
needed to reproduce Vogel et al. Fig 10.

Each case adds or changes ONE physical effect to demonstrate its impact.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
import aurora

from profiles_fig4 import rho as rho_kin, ne_noRMP, Te_noRMP
from transport_fig11 import transport
from geometry import R0, a, kappa


# ---------- the cases to compare ----------
CASES = [
    dict(label='1. Base\nno CX, n0=0, src=1e20',
         te_scale=1.0, cxr=False, n0_edge=0.0,    lambda_n=0.05,
         recycling=False, src=1.0e20),
    dict(label='2. Te×0.5 diagnostic\nno CX, n0=0, src=1e20',
         te_scale=0.5, cxr=False, n0_edge=0.0,    lambda_n=0.05,
         recycling=False, src=1.0e20),
    dict(label='3. + CX enabled\nn0=1e18 m⁻³, λ=5cm, src=1e20',
         te_scale=1.0, cxr=True,  n0_edge=1.0e18, lambda_n=0.05,
         recycling=False, src=1.0e20),
    dict(label='4. Final tuned\nn0=1e17 m⁻³, λ=10cm, recycling, src=2.5e18',
         te_scale=1.0, cxr=True,  n0_edge=1.0e17, lambda_n=0.10,
         recycling=True,  src=2.5e18),
]


def build_namelist(cfg):
    nml = aurora.default_nml.load_default_namelist()
    nml['device']       = 'EAST'
    nml['imp']          = 'Fe'
    nml['main_element'] = 'D'
    nml['cxr_flag']     = cfg['cxr']
    nml['recycling_flag'] = cfg['recycling']
    nml['wall_recycling'] = 0.9 if cfg['recycling'] else 0.0

    nml['Raxis_cm']  = R0 * 100.0
    nml['rvol_lcfs'] = a * 100.0 * np.sqrt(kappa)

    nml['source_type'] = 'const'
    nml['source_rate'] = cfg['src']

    nml['timing'] = {
        'times':    [0.0, 2.0],
        'dt_start': [1.0e-5, 1.0e-5],
        'steps_per_cycle': [1, 1],
        'dt_increase': [1.005, 1.0],
    }

    ne_safe = np.maximum(ne_noRMP, 0.05)
    nml['kin_profs']['ne']['fun']  = 'interp'
    nml['kin_profs']['ne']['rhop'] = rho_kin
    nml['kin_profs']['ne']['vals'] = ne_safe * 1.0e19 * 1.0e-6

    Te_safe = np.maximum(Te_noRMP * cfg['te_scale'], 0.02)
    nml['kin_profs']['Te']['fun']  = 'interp'
    nml['kin_profs']['Te']['rhop'] = rho_kin
    nml['kin_profs']['Te']['vals'] = Te_safe * 1000.0

    if cfg['cxr']:
        n0 = cfg['n0_edge'] * np.exp(-(1.0 - rho_kin) * a / cfg['lambda_n'])
    else:
        n0 = np.full_like(rho_kin, 1.0)   # negligible
    nml['kin_profs']['n0']['fun']  = 'interpa'
    nml['kin_profs']['n0']['rhop'] = rho_kin
    nml['kin_profs']['n0']['vals'] = n0 * 1.0e-6

    return nml


def run(cfg):
    nml = build_namelist(cfg)
    asim = aurora.aurora_sim(nml, geqdsk=None)

    rho_DV  = transport["No RMP"]["rho"]
    D = np.interp(asim.rhop_grid, rho_DV, transport["No RMP"]["D"]) * 1.0e4
    V = np.interp(asim.rhop_grid, rho_DV, transport["No RMP"]["V"]) * 1.0e2

    out = asim.run_aurora(D, V)
    return asim.rhop_grid, out['nz'][:, :, -1] * 1.0e6 / 1.0e15


# ---------- run all cases ----------
results = []
for cfg in CASES:
    print(f"Running: {cfg['label']}")
    rho, nz = run(cfg)
    results.append((cfg['label'], rho, nz))


# ---------- plot: one column per case, rows = Fe20+/Fe21+/Fe22+ ----------
n_cases = len(CASES)
fig, axes = plt.subplots(3, n_cases, figsize=(4 * n_cases, 8), sharey='row')

charge_states = [20, 21, 22]

for col, (label, rho, nz) in enumerate(results):
    for row, Z in enumerate(charge_states):
        ax = axes[row, col]
        ax.plot(rho, nz[:, Z], 'k-', lw=2)
        ax.grid(alpha=0.3)
        ax.set_xlim(0, 1)
        if col == 0:
            ax.set_ylabel(f'Fe$^{{{Z}+}}$\n(10$^{{15}}$ m$^{{-3}}$)')
        if row == 0:
            ax.set_title(label, fontsize=9)
        if row == 2:
            ax.set_xlabel(r'$\rho$')

plt.suptitle('Sensitivity study: parameter exploration toward Vogel et al. Fig 10', y=1.00)
plt.tight_layout()
plt.savefig(HERE / 'sensitivity_study.png', dpi=120, bbox_inches='tight')
plt.show()


# ---------- second plot: overlay all cases for direct shape comparison ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
colors = ['C0', 'C1', 'C2', 'C3']

for row, Z in enumerate(charge_states):
    ax = axes[row]
    for (label, rho, nz), color in zip(results, colors):
        ax.plot(rho, nz[:, Z], color=color, lw=2, label=label)
    ax.set_title(f'Fe$^{{{Z}+}}$')
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$n_{Fe}$ (10$^{15}$ m$^{-3}$)')
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 1)
    if row == 0:
        ax.legend(fontsize=7, loc='best')

plt.suptitle('All sensitivity cases overlaid')
plt.tight_layout()
plt.savefig(HERE / 'sensitivity_overlay.png', dpi=120, bbox_inches='tight')
plt.show()

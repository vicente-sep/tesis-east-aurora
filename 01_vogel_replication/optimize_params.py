"""
Automatic parameter optimization to match Vogel 2021 Fig 10 profiles.

Free parameters (4):
  log10(N0_EDGE)        in [15, 19]
  lambda_n        [m]   in [0.005, 0.30]
  log10(source_rate)    in [16, 21]
  wall_recycling        in [0.0, 1.0]

Cost: RMS error of (sim - paper) summed over Fe20+/Fe21+/Fe22+ and both
cases (No RMP + With RMP), normalized by each charge state's peak.

Uses scipy.optimize.differential_evolution (global, no derivatives needed).
"""

import numpy as np
import matplotlib.pyplot as plt
import aurora
from scipy.optimize import differential_evolution
from pathlib import Path
import time

from profiles_fig4 import (
    rho as rho_kin,
    ne_noRMP, ne_withRMP,
    Te_noRMP, Te_withRMP,
)
from transport_fig11 import transport
from geometry import R0, a, kappa
from paper_targets import rho_target, targets_noRMP, targets_withRMP

HERE = Path(__file__).parent

CHARGE_STATES = [20, 21, 22]


# ---------- Aurora helpers ----------
def build_namelist(ne_profile, Te_profile, N0_EDGE, LAMBDA_N, SOURCE_RATE, WALL_RECYCLING):
    nml = aurora.default_nml.load_default_namelist()
    nml['device']       = 'EAST'
    nml['imp']          = 'Fe'
    nml['main_element'] = 'D'
    nml['cxr_flag']     = True
    nml['recycling_flag'] = True
    nml['wall_recycling'] = WALL_RECYCLING

    nml['Raxis_cm']  = R0 * 100.0
    nml['rvol_lcfs'] = a * 100.0 * np.sqrt(kappa)

    nml['source_type'] = 'const'
    nml['source_rate'] = SOURCE_RATE

    nml['timing'] = {
        'times':    [0.0, 2.0],
        'dt_start': [1.0e-5, 1.0e-5],
        'steps_per_cycle': [1, 1],
        'dt_increase': [1.005, 1.0],
    }

    nml['kin_profs']['ne']['fun']  = 'interp'
    nml['kin_profs']['ne']['rhop'] = rho_kin
    nml['kin_profs']['ne']['vals'] = np.maximum(ne_profile, 0.05) * 1.0e19 * 1.0e-6

    nml['kin_profs']['Te']['fun']  = 'interp'
    nml['kin_profs']['Te']['rhop'] = rho_kin
    nml['kin_profs']['Te']['vals'] = np.maximum(Te_profile, 0.02) * 1000.0

    n0 = N0_EDGE * np.exp(-(1.0 - rho_kin) * a / LAMBDA_N)
    nml['kin_profs']['n0']['fun']  = 'interpa'
    nml['kin_profs']['n0']['rhop'] = rho_kin
    nml['kin_profs']['n0']['vals'] = n0 * 1.0e-6

    return nml


def run_case(case_name, params):
    if case_name == "No RMP":
        ne_p, Te_p = ne_noRMP, Te_noRMP
    else:
        ne_p, Te_p = ne_withRMP, Te_withRMP

    rho_DV  = transport[case_name]["rho"]
    D_knots = transport[case_name]["D"]
    V_knots = transport[case_name]["V"]

    nml = build_namelist(ne_p, Te_p, **params)
    asim = aurora.aurora_sim(nml, geqdsk=None)

    D = np.interp(asim.rhop_grid, rho_DV, D_knots) * 1.0e4
    V = np.interp(asim.rhop_grid, rho_DV, V_knots) * 1.0e2

    out = asim.run_aurora(D, V)
    return asim.rhop_grid, out['nz'][:, :, -1] * 1.0e6 / 1.0e15   # 10^15 m^-3


# ---------- cost function ----------
eval_counter = {'n': 0, 'best': np.inf}


def cost(x):
    """x = [log10_N0, lambda_n, log10_src, wall_recycling]"""
    eval_counter['n'] += 1
    params = dict(
        N0_EDGE        = 10.0 ** x[0],
        LAMBDA_N       = x[1],
        SOURCE_RATE    = 10.0 ** x[2],
        WALL_RECYCLING = x[3],
    )
    try:
        rho_n, nz_n = run_case("No RMP",   params)
        rho_r, nz_r = run_case("With RMP", params)
    except Exception as e:
        print(f"  [eval {eval_counter['n']}] FAILED: {e}")
        return 1e6

    err2 = 0.0
    count = 0
    for Z in CHARGE_STATES:
        # No RMP
        sim = np.interp(rho_target, rho_n, nz_n[:, Z])
        ref = targets_noRMP[Z]
        err2 += np.mean((sim - ref) ** 2) / (ref.max() ** 2 + 1e-30)
        count += 1
        # With RMP
        sim = np.interp(rho_target, rho_r, nz_r[:, Z])
        ref = targets_withRMP[Z]
        err2 += np.mean((sim - ref) ** 2) / (ref.max() ** 2 + 1e-30)
        count += 1

    c = np.sqrt(err2 / count)
    if c < eval_counter['best']:
        eval_counter['best'] = c
        print(f"  [eval {eval_counter['n']:4d}] new best: cost={c:.4f}  "
              f"N0={10**x[0]:.2e}  λ={x[1]:.3f}  src={10**x[2]:.2e}  rec={x[3]:.2f}")
    return c


# ---------- optimization ----------
BOUNDS = [
    (15.0, 19.0),    # log10(N0_EDGE)
    (0.005, 0.30),   # lambda_n [m]
    (16.0, 21.0),    # log10(source_rate)
    (0.0, 1.0),      # wall_recycling
]


print("=" * 60)
print("Running differential evolution optimization...")
print("Free parameters: log10(N0), lambda, log10(src), wall_recycling")
print(f"Bounds: {BOUNDS}")
print("=" * 60)

t0 = time.time()
result = differential_evolution(
    cost,
    bounds=BOUNDS,
    maxiter=30,         # generations
    popsize=10,         # population per generation
    tol=1e-3,
    mutation=(0.5, 1.5),
    recombination=0.7,
    seed=42,
    polish=True,        # local refinement at the end
    workers=1,          # Aurora not thread-safe
    updating='deferred',
    disp=True,
)
elapsed = time.time() - t0

best = result.x
print("\n" + "=" * 60)
print(f"Done in {elapsed:.1f} s after {eval_counter['n']} Aurora evaluations")
print(f"Best cost: {result.fun:.4f}")
print(f"  N0_EDGE        = {10**best[0]:.3e} m^-3")
print(f"  lambda_n       = {best[1]:.3f} m  ({best[1]*100:.1f} cm)")
print(f"  source_rate    = {10**best[2]:.3e} atoms/s")
print(f"  wall_recycling = {best[3]:.3f}")
print("=" * 60)


# ---------- final plot with best params ----------
best_params = dict(
    N0_EDGE        = 10.0 ** best[0],
    LAMBDA_N       = best[1],
    SOURCE_RATE    = 10.0 ** best[2],
    WALL_RECYCLING = best[3],
)
rho_n, nz_n = run_case("No RMP",   best_params)
rho_r, nz_r = run_case("With RMP", best_params)

fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
for row, Z in enumerate(CHARGE_STATES):
    for col, (case, sim_rho, sim_nz, target, color) in enumerate([
        ("No RMP",   rho_n, nz_n, targets_noRMP[Z],   'k'),
        ("With RMP", rho_r, nz_r, targets_withRMP[Z], 'r'),
    ]):
        ax = axes[row, col]
        ax.plot(sim_rho, sim_nz[:, Z], color=color, lw=2, label='Aurora (best fit)')
        ax.plot(rho_target, target, color=color, ls='', marker='o', ms=6, alpha=0.7,
                label='Vogel Fig 10')
        ax.set_ylabel(f'Fe$^{{{Z}+}}$ (10$^{{15}}$ m$^{{-3}}$)')
        ax.set_xlim(0, 1); ax.grid(alpha=0.3)
        if row == 0:
            ax.set_title(case)
            ax.legend(fontsize=8)
axes[-1, 0].set_xlabel(r'$\rho$')
axes[-1, 1].set_xlabel(r'$\rho$')

plt.suptitle(f'Best fit  |  N0={10**best[0]:.1e}  λ={best[1]*100:.1f}cm  '
             f'src={10**best[2]:.1e}  rec={best[3]:.2f}  |  cost={result.fun:.3f}')
plt.tight_layout()
plt.show()

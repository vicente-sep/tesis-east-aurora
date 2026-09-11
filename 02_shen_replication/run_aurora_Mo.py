"""
Aurora driver to reproduce Shen et al. 2019 Figs 8 & 9.
Runs Mo transport before and during on-axis ECRH (shot 55339).

Background parameters benchmarked against the Vogel 2021 Fe replication:
  - EAST geometry (Miller, R0=1.85m, a=0.45m, kappa=1.80, delta=0.70)
  - CX with n0=1e17 m^-3 at LCFS, lambda=10cm
  - Wall recycling = 0.9
  - Source rate tuned for Mo magnitudes (~10^16 m^-3, 10x higher than Fe)
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Force interactive matplotlib popup window
import matplotlib.pyplot as plt
plt.ioff()  # Disable interactive mode to ensure plt.show() blocks
import aurora
from pathlib import Path

HERE = Path(__file__).parent

from profiles_fig4 import rho as rho_kin, ne_before, ne_ecrh, Te_before, Te_ecrh
from transport_fig10_11 import transport
from geometry import R0, a, kappa


# ---------- knobs (same as final Vogel benchmark) ----------
N0_EDGE      = 1.0e17   # m^-3, neutral H/D at LCFS for CX
LAMBDA_N     = 0.10     # m, neutral decay length
SOURCE_RATE  = 2.5e19   # Mo atoms/s at edge (start higher than Fe -> Mo magnitudes are 10x)


def build_namelist(ne_profile, Te_profile):
    nml = aurora.default_nml.load_default_namelist()
    nml['device']       = 'EAST'
    nml['imp']          = 'Mo'
    nml['main_element'] = 'D'
    nml['cxr_flag']     = True
    nml['recycling_flag'] = True
    nml['wall_recycling'] = 0.9

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

    # ne: floor zeros, 10^19 m^-3 -> cm^-3
    ne_safe = np.maximum(ne_profile, 0.05)
    nml['kin_profs']['ne']['fun']  = 'interp'
    nml['kin_profs']['ne']['rhop'] = rho_kin
    nml['kin_profs']['ne']['vals'] = ne_safe * 1.0e19 * 1.0e-6

    # Te: floor, keV -> eV
    Te_safe = np.maximum(Te_profile, 0.02)
    nml['kin_profs']['Te']['fun']  = 'interp'
    nml['kin_profs']['Te']['rhop'] = rho_kin
    nml['kin_profs']['Te']['vals'] = Te_safe * 1000.0

    # neutral H/D
    n0 = N0_EDGE * np.exp(-(1.0 - rho_kin) * a / LAMBDA_N)
    nml['kin_profs']['n0']['fun']  = 'interpa'
    nml['kin_profs']['n0']['rhop'] = rho_kin
    nml['kin_profs']['n0']['vals'] = n0 * 1.0e-6

    return nml


def run_case(case_name):
    if case_name == "Before ECRH":
        ne_p, Te_p = ne_before, Te_before
    else:
        ne_p, Te_p = ne_ecrh, Te_ecrh

    rho_DV  = transport[case_name]["rho"]
    D_knots = transport[case_name]["D"]
    V_knots = transport[case_name]["V"]

    nml = build_namelist(ne_p, Te_p)
    asim = aurora.aurora_sim(nml, geqdsk=None)

    D = np.interp(asim.rhop_grid, rho_DV, D_knots) * 1.0e4
    V = np.interp(asim.rhop_grid, rho_DV, V_knots) * 1.0e2

    out = asim.run_aurora(D, V)
    return asim.rhop_grid, out['nz'][:, :, -1] * 1.0e6 / 1.0e16  # 10^16 m^-3 (Mo paper units)


# ---------- run both cases ----------
print("Running Before ECRH...")
rho_b, nz_b = run_case("Before ECRH")
print("Running On-axis ECRH...")
rho_e, nz_e = run_case("On-axis ECRH")


# ---------- plot like Figs 8 & 9: one figure per case ----------
charge_states = [28, 29, 30, 31]


def plot_case(rho, nz, color, title, fname):
    fig, axes = plt.subplots(4, 1, figsize=(6, 11), sharex=True)
    for ax, Z in zip(axes, charge_states):
        ax.plot(rho, nz[:, Z], color=color, lw=2)
        ax.set_ylabel(f'Mo$^{{{Z}+}}$ (10$^{{16}}$ m$^{{-3}}$)')
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel(r'$\rho$')
    axes[-1].set_xlim(0, 0.8)
    plt.suptitle(f'{title}  (n0_edge={N0_EDGE:.0e}, src={SOURCE_RATE:.0e})')
    plt.tight_layout()
    plt.savefig(HERE / fname, dpi=120)


plot_case(rho_b, nz_b, 'b', 'Aurora reproduction of Shen 2019 — Before ECRH (Fig 8)',
          'shen_reproduction_before.png')
plot_case(rho_e, nz_e, 'r', 'Aurora reproduction of Shen 2019 — On-axis ECRH (Fig 9)',
          'shen_reproduction_ecrh.png')


# ---------- diagnostic: where does each charge state peak? ----------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (label, rho, nz) in zip(axes, [('Before ECRH', rho_b, nz_b), ('On-axis ECRH', rho_e, nz_e)]):
    nz_norm = nz / nz.max(axis=0, keepdims=True).clip(min=1e-30)
    for Z in [15, 20, 24, 28, 30, 32, 36, 42]:
        ax.plot(rho, nz_norm[:, Z], lw=1.5, label=f'Mo$^{{{Z}+}}$')
    ax.set_xlabel(r'$\rho$'); ax.set_ylabel('normalised density')
    ax.set_title(f'{label}: where does each charge state peak?')
    ax.set_xlim(0, 0.8); ax.legend(fontsize=8); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(HERE / 'shen_diagnostic.png', dpi=120)


# ---------- diagnostic: total Mo (sum of all charge states) ----------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(rho_b, nz_b.sum(axis=1), 'b-', lw=2, label='Before ECRH')
ax.plot(rho_e, nz_e.sum(axis=1), 'r-', lw=2, label='On-axis ECRH')
ax.set_xlabel(r'$\rho$'); ax.set_ylabel('total Mo (10$^{16}$ m$^{-3}$)')
ax.set_title('Total Mo density (all charge states)')
ax.set_xlim(0, 0.8); ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(HERE / 'shen_total_Mo.png', dpi=120)

plt.show(block=True)

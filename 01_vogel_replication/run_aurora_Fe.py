"""
Aurora driver to reproduce Vogel et al. Figure 10.
Runs both No RMP and With RMP cases, overlays Fe20+/Fe21+/Fe22+ for direct comparison.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
OUT  = HERE.parent.parent / "Thesis" / "figures" / "fig10_reproduction.png"
import aurora

from profiles_fig4 import (
    rho as rho_kin,
    ne_noRMP, ne_withRMP,
    Te_noRMP, Te_withRMP,
)
from transport_fig11 import transport
from geometry import R0, a, kappa


# ---------- knobs ----------
TE_SCALE       = 1.0      # 1.0 = paper value
N0_EDGE        = 1.0e17   # m^-3, neutral H/D at LCFS for CX
LAMBDA_N       = 0.10     # m, neutral decay length
SOURCE_RATE    = 2.5e18   # Fe atoms/s at edge (tuned to match paper magnitudes)
WALL_RECYCLING = 0.9      # fraction recycled at wall


def build_namelist(ne_profile, Te_profile):
    nml = aurora.default_nml.load_default_namelist()

    nml['device']       = 'EAST'
    nml['imp']          = 'Fe'
    nml['main_element'] = 'D'
    nml['cxr_flag']     = True

    # impurity recycling at wall: lets Fe come back into the plasma instead of
    # being lost. Should let it accumulate more centrally over time.
    nml['recycling_flag'] = True
    nml['wall_recycling'] = WALL_RECYCLING   # fraction recycled (0=no recycling, 1=full)

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

    # ne: floor zeros, convert 10^19 m^-3 -> cm^-3
    ne_safe = np.maximum(ne_profile, 0.05)
    nml['kin_profs']['ne']['fun']  = 'interp'
    nml['kin_profs']['ne']['rhop'] = rho_kin
    nml['kin_profs']['ne']['vals'] = ne_safe * 1.0e19 * 1.0e-6

    # Te: floor, scale, convert keV -> eV
    Te_safe = np.maximum(Te_profile * TE_SCALE, 0.02)
    nml['kin_profs']['Te']['fun']  = 'interp'
    nml['kin_profs']['Te']['rhop'] = rho_kin
    nml['kin_profs']['Te']['vals'] = Te_safe * 1000.0

    # neutral H/D, edge-peaked exponential
    n0 = N0_EDGE * np.exp(-(1.0 - rho_kin) * a / LAMBDA_N)
    nml['kin_profs']['n0']['fun']  = 'interpa'
    nml['kin_profs']['n0']['rhop'] = rho_kin
    nml['kin_profs']['n0']['vals'] = n0 * 1.0e-6

    return nml


def run_case(case_name):
    ne_p = ne_noRMP if case_name == "No RMP" else ne_withRMP
    Te_p = Te_noRMP if case_name == "No RMP" else Te_withRMP

    rho_DV  = transport[case_name]["rho"]
    D_knots = transport[case_name]["D"]
    V_knots = transport[case_name]["V"]

    nml = build_namelist(ne_p, Te_p)
    asim = aurora.aurora_sim(nml, geqdsk=None)

    D = np.interp(asim.rhop_grid, rho_DV, D_knots) * 1.0e4   # m^2/s -> cm^2/s
    V = np.interp(asim.rhop_grid, rho_DV, V_knots) * 1.0e2   # m/s -> cm/s

    out = asim.run_aurora(D, V)
    return asim.rhop_grid, out['nz'][:, :, -1] * 1.0e6 / 1.0e15  # 10^15 m^-3


# ---------- run both cases ----------
print("Running No RMP...")
rho_noRMP, nz_noRMP = run_case("No RMP")
print("Running With RMP...")
rho_RMP,   nz_RMP   = run_case("With RMP")


# ---------- plot like Fig 10 ----------
fig, axes = plt.subplots(3, 1, figsize=(6, 9), sharex=True)
charge_states = [20, 21, 22]

for ax, Z in zip(axes, charge_states):
    ax.plot(rho_noRMP, nz_noRMP[:, Z], 'k-',  lw=2, label='no RMP')
    ax.plot(rho_RMP,   nz_RMP[:, Z],   'r-',  lw=2, label='with RMP')
    ax.set_ylabel(f'Fe$^{{{Z}+}}$ (10$^{{15}}$ m$^{{-3}}$)')
    ax.grid(alpha=0.3)
    if Z == 20:
        ax.legend()

axes[-1].set_xlabel(r'$\rho$')
axes[-1].set_xlim(0, 1)
plt.tight_layout()

OUT.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT, dpi=200, bbox_inches='tight')
print(f"Saved: {OUT}")
plt.show()

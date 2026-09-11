"""
plot_total_Mo.py
================
Genera la figura 'shen_total_Mo.png' para la tesis:
densidad total de Mo (sumada sobre todos los estados de carga)
antes vs durante ECRH, usando los inputs de Shen 2019.

Corre esto en tu computador (necesita aurora instalado):
    python plot_total_Mo.py

Guarda el PNG directamente en Thesis/figures/.
"""

import numpy as np
import matplotlib.pyplot as plt
import aurora
from pathlib import Path

from profiles_fig4 import rho as rho_kin, ne_before, ne_ecrh, Te_before, Te_ecrh
from transport_fig10_11 import transport
from geometry import R0, a, kappa

HERE = Path(__file__).parent
OUT  = HERE.parent.parent / "Thesis" / "figures" / "shen_total_Mo.png"

# Mismos parametros que run_aurora.py
N0_EDGE     = 1.0e17
LAMBDA_N    = 0.10
SOURCE_RATE = 2.5e19


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
    ne_safe = np.maximum(ne_profile, 0.05)
    nml['kin_profs']['ne']['fun']  = 'interp'
    nml['kin_profs']['ne']['rhop'] = rho_kin
    nml['kin_profs']['ne']['vals'] = ne_safe * 1.0e19 * 1.0e-6
    Te_safe = np.maximum(Te_profile, 0.02)
    nml['kin_profs']['Te']['fun']  = 'interp'
    nml['kin_profs']['Te']['rhop'] = rho_kin
    nml['kin_profs']['Te']['vals'] = Te_safe * 1000.0
    n0 = N0_EDGE * np.exp(-(1.0 - rho_kin) * a / LAMBDA_N)
    nml['kin_profs']['n0']['fun']  = 'interpa'
    nml['kin_profs']['n0']['rhop'] = rho_kin
    nml['kin_profs']['n0']['vals'] = n0 * 1.0e-6
    return nml


def run(case):
    ne_p, Te_p = (ne_before, Te_before) if case == "Before ECRH" else (ne_ecrh, Te_ecrh)
    rho_DV  = transport[case]["rho"]
    D_k     = transport[case]["D"]
    V_k     = transport[case]["V"]
    nml  = build_namelist(ne_p, Te_p)
    asim = aurora.aurora_sim(nml, geqdsk=None)
    D = np.interp(asim.rhop_grid, rho_DV, D_k) * 1.0e4
    V = np.interp(asim.rhop_grid, rho_DV, V_k) * 1.0e2
    out = asim.run_aurora(D, V)
    return asim.rhop_grid, out['nz'][:, :, -1] * 1.0e6 / 1.0e16   # 10^16 m^-3


print("Running Before ECRH...")
rho_b, nz_b = run("Before ECRH")
print("Running On-axis ECRH...")
rho_e, nz_e = run("On-axis ECRH")

total_b = nz_b.sum(axis=1)
total_e = nz_e.sum(axis=1)

dC = 100 * (total_e[0] - total_b[0]) / total_b[0]

plt.rcParams.update({'font.size': 15, 'axes.labelsize': 17,
                     'legend.fontsize': 13,
                     'xtick.labelsize': 14, 'ytick.labelsize': 14,
                     'axes.linewidth': 1.4})

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(rho_b, total_b, 'b-', lw=2.8, label='Before ECRH')
ax.plot(rho_e, total_e, 'r-', lw=2.8, label='On-axis ECRH')
ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'Total Mo density (10$^{16}$ m$^{-3}$)')

ax.set_xlim(0, 0.8)
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()

OUT.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT, dpi=200, bbox_inches='tight')
print(f"\nSaved: {OUT}")
print(f"Center suppression: {dC:.1f}%")

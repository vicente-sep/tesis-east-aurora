"""
plot_inputs_figure.py
=====================
Genera 'shen_inputs.png' para el capitulo 5.2 de la tesis:
perfiles ne, Te (Fig 4) y D, v (Figs 10, 11) digitalizados de Shen 2019.

Uso:
    python plot_inputs_figure.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from pathlib import Path

from profiles_fig4 import rho as rho_k, ne_before, ne_ecrh, Te_before, Te_ecrh
from transport_fig10_11 import transport

HERE = Path(__file__).parent
OUT  = HERE.parent.parent / "Thesis" / "figures" / "shen_inputs.png"

D_bef, V_bef = transport["Before ECRH"]["D"],   transport["Before ECRH"]["V"]
D_ecr, V_ecr = transport["On-axis ECRH"]["D"],  transport["On-axis ECRH"]["V"]
rho_t = transport["Before ECRH"]["rho"]


def pchip(x_k, y_k, N=200):
    x_f = np.linspace(x_k.min(), x_k.max(), N)
    return x_f, PchipInterpolator(x_k, y_k)(x_f)


plt.rcParams.update({'font.size': 15, 'axes.labelsize': 17,
                     'axes.titlesize': 16, 'legend.fontsize': 13,
                     'xtick.labelsize': 14, 'ytick.labelsize': 14,
                     'axes.linewidth': 1.4})

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# (a) ne
ax = axes[0]
for prof, col, lbl in [(ne_before, 'b', 'before ECRH'), (ne_ecrh, 'r', 'on-axis ECRH')]:
    x, y = pchip(rho_k, prof); ax.plot(x, y, color=col, lw=2, label=lbl)
    ax.plot(rho_k, prof, 'o', color=col, ms=4, alpha=0.7)
ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$n_e$ (10$^{19}$ m$^{-3}$)')
ax.set_title('(a) Electron density — Fig. 4 of paper')
ax.set_xlim(0, 0.8); ax.set_ylim(0, 5); ax.legend(); ax.grid(alpha=0.3)

# (b) Te
ax = axes[1]
for prof, col, lbl in [(Te_before, 'b', 'before ECRH'), (Te_ecrh, 'r', 'on-axis ECRH')]:
    x, y = pchip(rho_k, prof); ax.plot(x, y, color=col, lw=2, label=lbl)
    ax.plot(rho_k, prof, 'o', color=col, ms=4, alpha=0.7)
ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$T_e$ (keV)')
ax.set_title('(b) Electron temperature — Fig. 4')
ax.set_xlim(0, 0.8); ax.set_ylim(0, 2); ax.legend(); ax.grid(alpha=0.3)

# (c) D, v twin
ax = axes[2]; ax2 = ax.twinx()
for D, V, col, lbl in [(D_bef, V_bef, 'b', 'before ECRH'),
                       (D_ecr, V_ecr, 'r', 'on-axis ECRH')]:
    x, y = pchip(rho_t, D);  ax.plot(x, y, color=col, ls='-', lw=2, label=f'D {lbl}')
    x, y = pchip(rho_t, V); ax2.plot(x, y, color=col, ls='--', lw=1.5)
ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$D_\mathrm{Mo}$ (m$^2$/s) — solid')
ax2.set_ylabel(r'$v_\mathrm{Mo}$ (m/s) — dashed')
ax.set_title('(c) Transport coefficients — Figs. 10 & 11')
ax.set_xlim(0, 0.8); ax.grid(alpha=0.3); ax.legend(fontsize=12, loc='lower right')
ax2.axhline(0, color='gray', lw=0.5, ls=':')

plt.tight_layout()

OUT.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT, dpi=200, bbox_inches='tight')
print(f"Saved: {OUT}")

"""
plot_time_traces.py
===================
Genera 'fig_time_traces.png' para el capitulo 6 de la tesis:
time traces del shot 160869 (Ip, ECRH, WMHD, Dalpha).

Datos requeridos:
    shot_160869.npz  - generado con mdsplus_time_traces.py

Uso:
    python plot_time_traces.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pathlib import Path

HERE     = Path(__file__).parent
NPZ_PATH = HERE / "shot_160869.npz"
OUT_PATH = HERE.parent.parent / "Thesis" / "figures" / "fig_time_traces.png"

# ---- config ----
T_A   = (2.0, 5.0)
T_B   = (6.0, 9.0)
T_OFF = 5.58
XLIM  = (-0.5, 10.0)

# ---- load data ----
d = np.load(NPZ_PATH, allow_pickle=True)
t = d['t']
mask = (t >= XLIM[0]) & (t <= XLIM[1])

# ---- style ----
plt.rcParams.update({
    'font.size':        16, 'axes.labelsize': 18, 'axes.titlesize': 18,
    'legend.fontsize':  14, 'xtick.labelsize': 15, 'ytick.labelsize': 15,
    'axes.linewidth':   1.5, 'lines.linewidth': 2.4,
})

fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
fig.subplots_adjust(hspace=0.15)

# (1) Ip
axes[0].plot(t[mask], d['ip'][mask], color='#1f4e79')
axes[0].set_ylabel(r'$I_p$ (kA)')
axes[0].set_ylim(0, 350)

# (2) ECRH
axes[1].plot(t[mask], d['ecrh'][mask]/1000, color='#c00000')
axes[1].fill_between(t[mask], 0, d['ecrh'][mask]/1000, color='#c00000', alpha=0.15)
axes[1].set_ylabel(r'$P_{\mathrm{ECRH}}$ (MW)')
axes[1].set_ylim(0, 1.3)

# (3) WMHD
axes[2].plot(t[mask], d['wmhd'][mask]/1000, color='#2e7d32')
axes[2].set_ylabel(r'$W_{\mathrm{MHD}}$ (kJ)')
axes[2].set_ylim(0, 65)

# (4) Dalpha
axes[3].plot(t[mask], d['dalpha_mid'][mask], color='#7b1fa2', lw=1.3)
axes[3].set_ylabel(r'$D_\alpha$ midplane (a.u.)')
axes[3].set_xlabel('Time (s)')
axes[3].set_xlim(*XLIM)
axes[3].set_ylim(0, 0.25)

for ax in axes:
    ax.axvspan(*T_A, alpha=0.11, color='#c00000', zorder=0)
    ax.axvspan(*T_B, alpha=0.11, color='#1f4e79', zorder=0)
    ax.axvline(T_OFF, color='k', ls='--', lw=1.5, alpha=0.7)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.tick_params(axis='both', direction='in', length=5)

handles = [
    Patch(facecolor='#c00000', alpha=0.11, label='Phase A: ECRH ON'),
    Patch(facecolor='#1f4e79', alpha=0.11, label='Phase B: ECRH OFF'),
    Line2D([0], [0], color='k', ls='--', lw=1.5, label=f'ECRH switch-off (t = {T_OFF} s)'),
]
axes[0].legend(handles=handles, loc='upper right', fontsize=13, framealpha=0.95)

plt.tight_layout()

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_PATH, dpi=220, bbox_inches='tight')
print(f"Saved: {OUT_PATH}")

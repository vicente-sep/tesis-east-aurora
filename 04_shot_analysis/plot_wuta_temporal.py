"""
plot_wuta_temporal.py
=====================
Genera 'fig_wuta_temporal.png' para el capitulo 6 de la tesis:
evolucion temporal de la emision W-UTA en la cuerda midplane del shot 160869,
con la traza de ECRH para referencia.

Datos requeridos (local):
    160869.sif          - datos crudos del espectrometro XEUV_Long
    shot_160869.npz     - time traces MDSplus (generado con mdsplus_time_traces.py)

Uso:
    python plot_wuta_temporal.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import sif_parser
from pathlib import Path
from z_channels import Z_CM, CH_MIDPLANE

# ------------ paths ------------
HERE = Path(__file__).parent
SIF_PATH = HERE / "160869.sif"
NPZ_PATH = HERE / "shot_160869.npz"
OUT_PATH = HERE.parent.parent / "Thesis" / "figures" / "fig_wuta_temporal.png"

# ------------ config ------------
# CH_MIDPLANE viene de z_channels
WUTA_NM     = (4.80, 5.30)   # rango del W-UTA en nm
T_A         = (2.0, 5.0)     # fase ECRH ON
T_B         = (6.0, 9.0)     # fase ECRH OFF
T_OFF       = 5.58           # instante del switch-off
T_OFFSET    = 0.12           # offset de tiempo entre el .sif y el reloj MDSplus

# ------------ calibracion wavelength (biaodingbochang.m) ------------
_COEFS = np.polyfit([248, 376, 841], [33.73, 40.268, 33.73*2], 2)

def px_from_nm(nm):
    a, b, c = _COEFS
    disc = b**2 - 4*a*(c - nm*10)
    return int((-b + np.sqrt(disc)) / (2*a))


# ------------ 1) cargar .sif y extraer W-UTA(t) en midplane ------------
data, info = sif_parser.np_open(str(SIF_PATH))     # shape (n_t, n_sp, n_wl)
n_t, n_sp, n_wl = data.shape
cycle = info['CycleTime']

# background subtraction (ultimo frame, post-plasma)
data_bg = data - data[-1][None, :, :]

# eje temporal de los frames alineado con MDSplus
t_frame = np.arange(n_t) * cycle + T_OFFSET

# rango del W-UTA en pixeles
p_lo, p_hi = px_from_nm(WUTA_NM[0]), px_from_nm(WUTA_NM[1])
print(f"W-UTA pixel range: {p_lo}--{p_hi}")

# integrar sobre el rango de wavelength -> (n_t, n_sp)
wuta_signal = data_bg[:, :, p_lo:p_hi].sum(axis=2)
# seleccionar la cuerda midplane
wuta_mid = wuta_signal[:, CH_MIDPLANE]


# ------------ 2) cargar time traces MDSplus (para plot de ECRH) ------------
d = np.load(NPZ_PATH, allow_pickle=True)
t_mds = d['t']
ecrh_MW = d['ecrh'] / 1000.0


# ------------ 3) figura ------------
plt.rcParams.update({
    'font.size':        16,
    'axes.labelsize':   18,
    'axes.titlesize':   18,
    'legend.fontsize':  14,
    'xtick.labelsize':  15,
    'ytick.labelsize':  15,
    'axes.linewidth':   1.5,
    'lines.linewidth':  2.4,
})

fig, axes = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True)
fig.subplots_adjust(hspace=0.15)

# Panel 1: W-UTA temporal
ax = axes[0]
ax.plot(t_frame, wuta_mid, 'o-', color='#7b1fa2', ms=5, lw=1.8,
        label='W-UTA emission (midplane chord)')
ax.set_ylabel('W-UTA intensity (counts)')

# Panel 2: ECRH
ax = axes[1]
ax.plot(t_mds, ecrh_MW, color='#c00000', lw=2)
ax.fill_between(t_mds, 0, ecrh_MW, color='#c00000', alpha=0.2)
ax.set_ylabel(r'$P_{\mathrm{ECRH}}$ (MW)')
ax.set_xlabel('Time (s)')
ax.set_xlim(-0.5, 10)
ax.set_ylim(0, 1.3)

for ax in axes:
    ax.axvspan(*T_A, alpha=0.11, color='#c00000', zorder=0)
    ax.axvspan(*T_B, alpha=0.11, color='#1f4e79', zorder=0)
    ax.axvline(T_OFF, color='k', ls='--', lw=1.5, alpha=0.7)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.tick_params(axis='both', which='major', direction='in', length=5)

handles = [
    Patch(facecolor='#c00000', alpha=0.11, label='Phase A: ECRH ON'),
    Patch(facecolor='#1f4e79', alpha=0.11, label='Phase B: ECRH OFF'),
    Line2D([0], [0], color='k', ls='--', lw=1.5, label='ECRH switch-off (t = 5.58 s)'),
]
axes[0].legend(handles=handles, loc='upper right', fontsize=13, framealpha=0.95)

plt.tight_layout()

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_PATH, dpi=220, bbox_inches='tight')
print(f"\nSaved: {OUT_PATH}")


# ------------ 4) numeros clave ------------
mask_A = (t_frame >= T_A[0]) & (t_frame <= T_A[1])
mask_B = (t_frame >= T_B[0]) & (t_frame <= T_B[1])
mean_A = wuta_mid[mask_A].mean()
mean_B = wuta_mid[mask_B].mean()

print("\n" + "="*50)
print("KEY NUMBERS")
print("="*50)
print(f"W-UTA at midplane, mean over phase:")
print(f"  Phase A (ECRH ON):   {mean_A:.2e} counts")
print(f"  Phase B (ECRH OFF):  {mean_B:.2e} counts")
print(f"  Ratio B/A:           {mean_B/mean_A:.3f}")
print(f"  Reduction:           {100*(1-mean_B/mean_A):.1f}%")

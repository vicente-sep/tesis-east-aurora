"""
check_fe_lines.py
=================
Verifica si las lineas de Fe (12.88, 13.29, 13.51 nm) del paper de Vogel
son visibles sobre el background en el shot 160869.

Cada linea:
  - Fe XXI  = Fe20+ = 12.88 nm  (multiplet)
  - Fe XXIII = Fe22+ = 13.29 nm
  - Fe XXII  = Fe21+ = 13.51 nm  (posible blending con Li III a 13.50 nm)

Uso:
    python check_fe_lines.py

Salida:
  fig_fe_check.png  - zoom del espectro en la region 12-15 nm
  reporte por linea: peak/background ratio
"""

import numpy as np
import matplotlib.pyplot as plt
import sif_parser
from pathlib import Path
from z_channels import Z_CM, CH_MIDPLANE

HERE     = Path(__file__).parent
SIF_PATH = HERE / "160869.sif"
OUT_FIG  = HERE.parent.parent / "Thesis" / "figures" / "fig_fe_check.png"

# ---- config ----
T_A      = (2.0, 5.0)
T_B      = (6.0, 9.0)
T_OFFSET = 0.12
CH_MID   = 9

# Lineas Fe segun Vogel 2021
FE_LINES = {
    'Fe XXI (Fe20+)':   12.88,
    'Fe XXIII (Fe22+)': 13.29,
    'Fe XXII (Fe21+)':  13.51,   # posible blending con Li III
}
LINE_HALF_WIDTH_NM = 0.06        # ventana total ~0.12 nm para integrar (~10 pixels)

# ---- calibracion ----
_COEFS = np.polyfit([248, 376, 841], [33.73, 40.268, 33.73*2], 2)

def pix_to_wl_nm(p):
    return np.polyval(_COEFS, p) / 10.0

def wl_nm_to_pix(nm):
    a_, b_, c_ = _COEFS
    disc = b_**2 - 4*a_*(c_ - nm*10)
    return (-b_ + np.sqrt(disc)) / (2*a_)


# ---- load .sif ----
data, info = sif_parser.np_open(str(SIF_PATH))
n_t, n_sp, n_wl = data.shape
cycle = info['CycleTime']
t_frame = np.arange(n_t) * cycle + T_OFFSET
data_bg = data - data[-1][None, :, :]

wl_nm = pix_to_wl_nm(np.arange(n_wl))

# Espectros midplane, promediados por fase
mask_A = (t_frame >= T_A[0]) & (t_frame <= T_A[1])
mask_B = (t_frame >= T_B[0]) & (t_frame <= T_B[1])
spec_A = data_bg[mask_A, CH_MID, :].mean(axis=0)
spec_B = data_bg[mask_B, CH_MID, :].mean(axis=0)


# ---- analisis por linea ----
def line_stats(spec, wl_nm, target_nm, half_width):
    """Devuelve (peak, background, snr) para una linea."""
    mask_line = np.abs(wl_nm - target_nm) < half_width
    mask_bg   = (np.abs(wl_nm - target_nm) < 3*half_width) & ~mask_line
    peak = float(spec[mask_line].max())
    bg   = float(np.median(spec[mask_bg]))
    bg_std = float(np.std(spec[mask_bg]))
    snr  = (peak - bg) / max(bg_std, 1e-9)
    integ = float(spec[mask_line].sum())
    return peak, bg, bg_std, snr, integ


print("=" * 70)
print(f"Fe line diagnostics for shot 160869 (midplane channel {CH_MID})")
print("=" * 70)
print(f"{'Line':<22} {'nm':>6} | {'Phase':<8} {'Peak':>10} {'BG':>10} {'BG std':>10} {'SNR':>8}")
print("-" * 90)

fits = {}
for name, wl in FE_LINES.items():
    for phase, spec in [('A_ECRHon', spec_A), ('B_ECRHoff', spec_B)]:
        p, b, bs, snr, integ = line_stats(spec, wl_nm, wl, LINE_HALF_WIDTH_NM)
        note = ''
        if snr > 10:    note = '  <- STRONG'
        elif snr > 3:   note = '  <- weak but usable'
        elif snr > 0:   note = '  <- marginal'
        else:           note = '  <- not visible'
        print(f"{name:<22} {wl:>6.2f} | {phase:<8} {p:>10.1f} {b:>10.1f} {bs:>10.1f} {snr:>8.1f}{note}")
    print()


# ---- figura: zoom en 12.5 - 14 nm ----
plt.rcParams.update({
    'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 13,
    'legend.fontsize': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'axes.linewidth': 1.2,
})

fig, ax = plt.subplots(figsize=(12, 5.5))
ax.plot(wl_nm, spec_A, color='#c00000', lw=1.4, label='Phase A (ECRH ON)')
ax.plot(wl_nm, spec_B, color='#1f4e79', lw=1.4, label='Phase B (ECRH OFF)')

for name, wl in FE_LINES.items():
    ax.axvline(wl, color='#2e7d32', ls=':', lw=1.5, alpha=0.7)
    ax.text(wl, ax.get_ylim()[1] if False else 1.02, name,
            transform=ax.get_xaxis_transform(),
            fontsize=9, ha='center', va='bottom', color='#2e7d32',
            fontweight='bold', rotation=0)

# Li III para referencia
ax.axvline(13.50, color='gray', ls=':', lw=1, alpha=0.5)
ax.text(13.50, 1.02, 'Li III', transform=ax.get_xaxis_transform(),
        fontsize=9, ha='center', va='bottom', color='gray', style='italic')

ax.set_xlim(12.4, 14.2)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Intensity (counts, midplane)')
ax.set_title('EAST #160869 - Zoom on Fe line region (midplane chord)')
ax.legend(loc='upper right', framealpha=0.95)
ax.grid(True, alpha=0.3, linestyle=':')
ax.tick_params(direction='in', length=5)

plt.tight_layout()
OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_FIG, dpi=140, bbox_inches='tight')
print(f"Saved: {OUT_FIG}")

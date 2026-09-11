"""
plot_FeXXIII_raw_zprofile.py
============================
Visualiza el perfil ESPACIAL crudo (intensity vs z, sin folding a rho_tan)
del Fe XXIII para las dos fases. Ayuda a diagnosticar si el flat-top del
fit v/D viene de la fisica o de un z_axis mal elegido.

Salida:
    fig_FeXXIII_raw_zprofile.png
    + reporte de peak/centroid por fase

Uso:
    python plot_FeXXIII_raw_zprofile.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sif_parser
from pathlib import Path
from z_channels import Z_CM

HERE     = Path(__file__).parent
SIF_PATH = HERE / "160869.sif"
OUT_FIG  = HERE.parent.parent / "Thesis" / "figures" / "fig_FeXXIII_raw_zprofile.png"

# ---- config ----
FE_LINE_NM  = 13.29
LINE_HW_NM  = 0.06
T_A, T_B    = (2.0, 5.0), (6.0, 9.0)
T_OFFSET    = 0.12

# ---- calibracion ----
_COEFS = np.polyfit([248, 376, 841], [33.73, 40.268, 33.73*2], 2)
def pix_to_wl_nm(p): return np.polyval(_COEFS, p) / 10.0


# ---- load & extract Fe XXIII brightness ----
data, info = sif_parser.np_open(str(SIF_PATH))
n_t, n_sp, n_wl = data.shape
cycle = info['CycleTime']
t_frame = np.arange(n_t) * cycle + T_OFFSET
data_bg = data - data[-1][None, :, :]
wl_nm   = pix_to_wl_nm(np.arange(n_wl))

mask_line = np.abs(wl_nm - FE_LINE_NM) < LINE_HW_NM
mask_bg   = ((np.abs(wl_nm - FE_LINE_NM) > 2*LINE_HW_NM) &
             (np.abs(wl_nm - FE_LINE_NM) < 5*LINE_HW_NM))
signal    = data_bg[:, :, mask_line].sum(axis=2)
bg        = data_bg[:, :, mask_bg].mean(axis=2) * mask_line.sum()
fe = signal - bg

mask_A = (t_frame >= T_A[0]) & (t_frame <= T_A[1])
mask_B = (t_frame >= T_B[0]) & (t_frame <= T_B[1])
prof_A = fe[mask_A].mean(axis=0)
prof_B = fe[mask_B].mean(axis=0)


# ---- diagnostics: peak & centroid ----
z_peak_A = Z_CM[np.argmax(prof_A)]
z_peak_B = Z_CM[np.argmax(prof_B)]
z_cent_A = np.sum(Z_CM * np.maximum(prof_A, 0)) / np.sum(np.maximum(prof_A, 0))
z_cent_B = np.sum(Z_CM * np.maximum(prof_B, 0)) / np.sum(np.maximum(prof_B, 0))

print("=" * 60)
print("Fe XXIII raw z-profile diagnostics")
print("=" * 60)
print(f"Phase A (ECRH ON):  peak @ z = {z_peak_A:+6.1f} cm  |  centroid = {z_cent_A:+6.1f} cm")
print(f"Phase B (ECRH OFF): peak @ z = {z_peak_B:+6.1f} cm  |  centroid = {z_cent_B:+6.1f} cm")
print()
print(f"Peak intensities:  A = {prof_A.max():.1f}  |  B = {prof_B.max():.1f}   (ratio A/B = {prof_A.max()/prof_B.max():.2f})")
print()
print("Interpretation guide:")
print("  * If peak & centroid ~ same (within few cm)  ->  emission profile is symmetric")
print("  * If peak sharp but centroid different       ->  asymmetry or long tail")
print("  * z_peak ~ z_axis (the vertical shift of the magnetic axis)")
print("  * A flat top around z_peak means eps(rho) is hollow / plateau in the core")


# ---- figure ----
plt.rcParams.update({
    'font.size': 16, 'axes.labelsize': 18, 'axes.titlesize': 17,
    'legend.fontsize': 13, 'xtick.labelsize': 15, 'ytick.labelsize': 15,
    'axes.linewidth': 1.5,
})

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

for ax, prof, tag, color, z_p, z_c in [
    (axes[0], prof_A, 'A (ECRH ON)',  '#c00000', z_peak_A, z_cent_A),
    (axes[1], prof_B, 'B (ECRH OFF)', '#1f4e79', z_peak_B, z_cent_B),
]:
    ax.plot(Z_CM, prof, 'o-', color=color, ms=10, lw=2.0,
            label='Fe XXIII raw $B(z)$')
    ax.axvline(0.0, color='gray', ls=':',  lw=1.2, label='$z = 0$')
    ax.axvline(z_p, color='k',    ls='--', lw=1.8, label=f'peak @ {z_p:+.1f} cm')
    ax.axvline(z_c, color='green',ls='--', lw=1.8, label=f'centroid @ {z_c:+.1f} cm')
    ax.set_xlabel('z (cm)')
    ax.set_ylabel('Fe XXIII brightness (counts)')
    ax.set_title(f'Phase {tag}')
    ax.grid(alpha=0.3, ls=':')
    ax.legend(loc='best', framealpha=0.95)
    ax.tick_params(direction='in', length=5)

plt.tight_layout()

OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_FIG, dpi=220, bbox_inches='tight')
print(f"\nSaved: {OUT_FIG}")

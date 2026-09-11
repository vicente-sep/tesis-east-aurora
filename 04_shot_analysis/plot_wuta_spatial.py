"""
plot_wuta_spatial.py
====================
Genera 'fig_wuta_spatial.png' para el capitulo 6 de la tesis:
perfil espacial del W-UTA en coordenada rho tangent-point,
promediado en cada fase de ECRH.

El eje X es rho_tan = |z - z_axis|/(kappa*a), el rho minimo (tangent point)
que atraviesa cada chord vertical. Como z y (-z) dan el mismo rho_tan,
los datos aparecen "doblados" (upper y lower channels marcados distinto).

Datos requeridos:
    160869.sif  - datos crudos del espectrometro XEUV_Long

Uso:
    python plot_wuta_spatial.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sif_parser
from pathlib import Path
from z_channels import Z_CM, CH_MIDPLANE

HERE     = Path(__file__).parent
SIF_PATH = HERE / "160869.sif"
OUT_PATH = HERE.parent.parent / "Thesis" / "figures" / "fig_wuta_spatial.png"

# ---- config ----
WUTA_NM  = (4.80, 5.30)
T_A      = (2.0, 5.0)
T_B      = (6.0, 9.0)
T_OFFSET = 0.12

# EAST geometry
R0, a, kappa = 1.85, 0.45, 1.80
Z_AXIS_CM = 10.7   # obtenido del fit de Fe XXIII (Cap 6, seccion fit v/D)

# 20 canales espaciales cubren aproximadamente z = -27 a +33 cm
# Z_MIN, Z_MAX ya no se usan; ver Z_CM en z_channels.py

# ---- wavelength calibration ----
_COEFS = np.polyfit([248, 376, 841], [33.73, 40.268, 33.73*2], 2)

def px_from_nm(nm):
    a, b, c = _COEFS
    disc = b**2 - 4*a*(c - nm*10)
    return int((-b + np.sqrt(disc)) / (2*a))


# ---- load .sif ----
data, info = sif_parser.np_open(str(SIF_PATH))
n_t, n_sp, n_wl = data.shape
cycle = info['CycleTime']
t_frame = np.arange(n_t) * cycle + T_OFFSET

# background subtraction
data_bg = data - data[-1][None, :, :]

# rango W-UTA en pixeles
p_lo, p_hi = px_from_nm(WUTA_NM[0]), px_from_nm(WUTA_NM[1])
print(f"W-UTA pixel range: {p_lo}--{p_hi}")

# integrar sobre wavelength -> (n_t, n_sp)
wuta = data_bg[:, :, p_lo:p_hi].sum(axis=2)

# perfiles espaciales por fase (promedio temporal, normalizado por # de frames)
mask_A = (t_frame >= T_A[0]) & (t_frame <= T_A[1])
mask_B = (t_frame >= T_B[0]) & (t_frame <= T_B[1])
prof_A = wuta[mask_A].mean(axis=0)
prof_B = wuta[mask_B].mean(axis=0)

# eje espacial en cm
z_cm = Z_CM
ch   = np.arange(n_sp)

# ---- mapeo a rho tangent-point ----
rho_tan = np.abs(z_cm - Z_AXIS_CM) / (kappa * a * 100.0)  # dimensionless
rho_tan = np.clip(rho_tan, 0, 1)

# separar upper (z > z_axis) y lower (z < z_axis) para plotear con marcadores distintos
upper = z_cm >= Z_AXIS_CM
lower = z_cm <  Z_AXIS_CM


# ---- style ----
plt.rcParams.update({
    'font.size':        16, 'axes.labelsize': 18, 'axes.titlesize': 17,
    'legend.fontsize':  13, 'xtick.labelsize': 15, 'ytick.labelsize': 15,
    'axes.linewidth':   1.5, 'lines.linewidth': 2.4,
})

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Peak positions per phase (for legend)
iA_peak = int(np.argmax(prof_A))
iB_peak = int(np.argmax(prof_B))

# Panel 1: perfil absoluto vs z (cm)
ax = axes[0]
ax.plot(z_cm, prof_A, 'o-', color='#c00000', ms=8, lw=2.0,
        label='Phase A: ECRH ON')
ax.plot(z_cm, prof_B, 's-', color='#1f4e79', ms=8, lw=2.0,
        label='Phase B: ECRH OFF')
ax.axvline(0.0, color='gray', ls=':', lw=1.2)

ax.set_xlabel(r'Vertical position $z$ (cm)')
ax.set_ylabel('W-UTA intensity per frame (counts)')
ax.set_title('(a) Absolute spatial profile')
ax.legend(loc='upper left', framealpha=0.95)
ax.grid(True, alpha=0.3, linestyle=':')
ax.tick_params(direction='in', length=5)

# Panel 2: perfiles normalizados
ax = axes[1]
pA_n = prof_A / prof_A.max()
pB_n = prof_B / prof_B.max()
ax.plot(z_cm, pA_n, 'o-', color='#c00000', ms=8, lw=2.0,
        label=f'Phase A  (peak @ $z = {z_cm[iA_peak]:+.1f}$ cm)')
ax.plot(z_cm, pB_n, 's-', color='#1f4e79', ms=8, lw=2.0,
        label=f'Phase B  (peak @ $z = {z_cm[iB_peak]:+.1f}$ cm)')
ax.axvline(0.0, color='gray', ls=':', lw=1.2)
ax.set_xlabel(r'Vertical position $z$ (cm)')
ax.set_ylabel('W-UTA intensity (normalised to peak)')
ax.set_title('(b) Normalised profile (shape)')
ax.set_ylim(0, 1.15)
ax.legend(loc='lower left', framealpha=0.95)
ax.grid(True, alpha=0.3, linestyle=':')
ax.tick_params(direction='in', length=5)

plt.tight_layout()

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_PATH, dpi=220, bbox_inches='tight')
print(f"\nSaved: {OUT_PATH}")

# ---- key numbers ----
iA, iB = int(np.argmax(prof_A)), int(np.argmax(prof_B))
print("\n" + "="*50)
print("KEY NUMBERS")
print("="*50)
print(f"z_axis assumed: {Z_AXIS_CM:+.1f} cm  (from Fe XXIII fit)")
print(f"\nPeak of W-UTA profile:")
print(f"  Phase A (ECRH ON):   z = {z_cm[iA]:+6.1f} cm  ->  rho_tan = {rho_tan[iA]:.3f}  |  I = {prof_A[iA]:.2e}")
print(f"  Phase B (ECRH OFF):  z = {z_cm[iB]:+6.1f} cm  ->  rho_tan = {rho_tan[iB]:.3f}  |  I = {prof_B[iB]:.2e}")
print(f"\nIntegrated (sum over all channels):")
print(f"  Phase A: {prof_A.sum():.2e}")
print(f"  Phase B: {prof_B.sum():.2e}")
print(f"  Ratio B/A: {prof_B.sum()/prof_A.sum():.3f}")

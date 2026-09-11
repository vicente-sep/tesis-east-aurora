"""
verify_forward_projection_CVI.py
================================
Test de la metodologia de forward projection usando la linea C VI (3.373 nm).

LOGICA:
  C es hidrogeneno (C^6+) practicamente en TODO el plasma cuando Te > 500 eV
  (que es toda la zona de emision). Entonces:
     - la fraccion de carbono en C6+ es ~100% en todos lados
     - la densidad de C se puede aproximar como constante en el core
     - el shape de la emisividad C VI depende SOLO de n_e(rho) y PEC(Te)
     - NO depende del transporte de C

Si el forward projection funciona correctamente, entonces:
     B_sim(z)  computado con  eps(rho) proportional a n_e(rho)
     deberia matchear
     B_meas(z) de la linea C VI

Si matchean -> nuestra chord integration + geometria (Miller, Z channels,
z_axis) esta bien -> podemos confiar en los fits de Fe XXIII y W-UTA.

Si NO matchean -> hay un error sistematico que estamos arrastrando (chord
width finita, z_axis mal, calibracion espacial, etc.).

Uso:
    python verify_forward_projection_CVI.py
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import sif_parser

HERE     = Path(__file__).parent
SIF_PATH = HERE / "160869.sif"
OUT_FIG  = HERE.parent.parent / "Thesis" / "figures" / "fig_CVI_forward_projection_check.png"

SHEN_DIR = HERE.parent / "02_shen_replication"
sys.path.insert(0, str(SHEN_DIR))
from profiles_fig4 import rho as rho_shen, ne_before, ne_ecrh
from z_channels import Z_CM

# ---- geometria ----
R0, a, kappa = 1.85, 0.45, 1.80

# ---- config ----
CVI_NM      = 3.373
LINE_HW_NM  = 0.03      # C VI es un pico agudo, ventana estrecha
T_A, T_B    = (2.0, 5.0), (6.0, 9.0)
T_OFFSET    = 0.12

# ---- calibracion ----
_COEFS = np.polyfit([248, 376, 841], [33.73, 40.268, 33.73*2], 2)
def pix_to_wl_nm(p): return np.polyval(_COEFS, p) / 10.0


# =============================================================================
# 1) medir B(z) para C VI
# =============================================================================

def measure_CVI_brightness():
    data, info = sif_parser.np_open(str(SIF_PATH))
    n_t, n_sp, n_wl = data.shape
    cycle = info['CycleTime']
    t_frame = np.arange(n_t) * cycle + T_OFFSET
    data_bg = data - data[-1][None, :, :]
    wl_nm   = pix_to_wl_nm(np.arange(n_wl))

    mask_line = np.abs(wl_nm - CVI_NM) < LINE_HW_NM
    mask_bg   = ((np.abs(wl_nm - CVI_NM) > 2*LINE_HW_NM) &
                 (np.abs(wl_nm - CVI_NM) < 5*LINE_HW_NM))

    signal    = data_bg[:, :, mask_line].sum(axis=2)
    bg_per_px = data_bg[:, :, mask_bg].mean(axis=2)
    cvi       = signal - bg_per_px * mask_line.sum()

    mask_A = (t_frame >= T_A[0]) & (t_frame <= T_A[1])
    mask_B = (t_frame >= T_B[0]) & (t_frame <= T_B[1])
    B_A = cvi[mask_A].mean(axis=0)
    B_B = cvi[mask_B].mean(axis=0)
    return B_A, B_B


# =============================================================================
# 2) emisividad esperada  eps(rho) ~ n_e(rho)  (asumiendo n_C, PEC constantes)
# =============================================================================

def eps_CVI_prediction(ne_prof):
    """
    eps(rho) proporcional a n_e(rho) x n_C x PEC(Te)
                ~ n_e(rho)                (constantes: n_C fully stripped,
                                            PEC casi cte para Te >> Ei_C6)
    Devuelve rho (0..1 fine grid), eps.
    """
    rho = np.linspace(0, 1, 200)
    ne_at_rho = np.interp(rho, rho_shen, np.maximum(ne_prof, 0.05))   # 10^19 m^-3
    eps = ne_at_rho.copy()      # unidades arbitrarias, importa solo el shape
    return rho, eps


# =============================================================================
# 3) chord integration (mismo codigo que Fe XXIII)
# =============================================================================

def chord_integrate(rho_grid, eps, z_cm_channels, z_axis_cm=0.0, n_u=200):
    B = np.zeros_like(z_cm_channels, dtype=float)
    for i, z in enumerate(z_cm_channels):
        rho_min = abs((z - z_axis_cm)/100.0) / (kappa * a)
        if rho_min >= 1.0:
            B[i] = 0.0; continue
        u_max = np.sqrt(1.0 - rho_min**2)
        u   = np.linspace(0.0, u_max, n_u)
        rho = np.sqrt(u**2 + rho_min**2)
        eps_at_u = np.interp(rho, rho_grid, eps, left=0, right=0)
        B[i] = 2.0 * a * np.trapz(eps_at_u, u)
    return B


# =============================================================================
# 4) fit: solo el z_axis (dejamos que el fit encuentre el shift optimo)
# =============================================================================

def fit_zaxis(B_meas, ne_prof, tag):
    print(f"\nFitting z_axis for phase {tag}...")
    rho_g, eps_g = eps_CVI_prediction(ne_prof)
    B_meas_n = B_meas / B_meas.max()

    def cost(params):
        z_axis = params[0]
        B_sim = chord_integrate(rho_g, eps_g, Z_CM, z_axis_cm=z_axis)
        if B_sim.max() <= 0:
            return 1.0e6
        return float(np.mean((B_sim/B_sim.max() - B_meas_n)**2))

    res = differential_evolution(cost, [(-15.0, 15.0)], maxiter=40, popsize=15,
                                 seed=42, tol=1e-4)
    z_axis = res.x[0]
    B_sim  = chord_integrate(rho_g, eps_g, Z_CM, z_axis_cm=z_axis)
    print(f"  z_axis fitted = {z_axis:+.2f} cm")
    print(f"  best cost = {res.fun:.4e}")
    return z_axis, B_sim, rho_g, eps_g


# =============================================================================
# 5) main
# =============================================================================

print("Loading measured C VI brightness B(z)...")
B_A, B_B = measure_CVI_brightness()
print(f"  Phase A peak: {B_A.max():.1f}")
print(f"  Phase B peak: {B_B.max():.1f}")

zA, B_A_sim, rho_g, eps_A = fit_zaxis(B_A, ne_ecrh,   'A (ECRH ON)')
zB, B_B_sim, _,     eps_B = fit_zaxis(B_B, ne_before, 'B (ECRH OFF)')


# =============================================================================
# 6) figura
# =============================================================================

plt.rcParams.update({
    'font.size': 16, 'axes.labelsize': 18, 'axes.titlesize': 17,
    'legend.fontsize': 13, 'xtick.labelsize': 15, 'ytick.labelsize': 15,
    'axes.linewidth': 1.5,
})

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

ax = axes[0]
ax.plot(Z_CM, B_A/B_A.max(), 'o', color='#c00000', ms=10, label='Measured C VI (norm.)')
ax.plot(Z_CM, B_A_sim/B_A_sim.max(), '-', color='k', lw=2.5,
        label=f'Forward-projected  ($z_{{ax}}={zA:+.1f}$ cm)')
ax.set_xlabel('z (cm)'); ax.set_ylabel('Brightness (norm.)')
ax.set_title('(a) Phase A: ECRH ON')
ax.set_ylim(0, 1.15); ax.legend(); ax.grid(alpha=0.3, ls=':')
ax.tick_params(direction='in', length=5)

ax = axes[1]
ax.plot(Z_CM, B_B/B_B.max(), 's', color='#1f4e79', ms=10, label='Measured C VI (norm.)')
ax.plot(Z_CM, B_B_sim/B_B_sim.max(), '-', color='k', lw=2.5,
        label=f'Forward-projected  ($z_{{ax}}={zB:+.1f}$ cm)')
ax.set_xlabel('z (cm)'); ax.set_ylabel('Brightness (norm.)')
ax.set_title('(b) Phase B: ECRH OFF')
ax.set_ylim(0, 1.15); ax.legend(); ax.grid(alpha=0.3, ls=':')
ax.tick_params(direction='in', length=5)

plt.tight_layout()

OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_FIG, dpi=220, bbox_inches='tight')
print(f"\nSaved: {OUT_FIG}")

print("\n" + "="*60)
print("VALIDATION SUMMARY")
print("="*60)
print(f"z_axis converge a:")
print(f"  Phase A (ECRH ON):  {zA:+.2f} cm")
print(f"  Phase B (ECRH OFF): {zB:+.2f} cm")
print(f"\nComparar con el z_axis obtenido con Fe XXIII (~ +10.7 cm)")
print(f"Si estos valores son similares -> la geometria/methodo esta validada.")
print(f"Si son distintos -> algo esta mal en el analisis de una de las lineas.")

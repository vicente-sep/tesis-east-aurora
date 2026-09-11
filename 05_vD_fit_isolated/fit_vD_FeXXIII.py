"""
fit_FeXXIII_shen_method.py
==========================
Ajuste v/D del transporte de Fe usando la linea Fe XXIII (13.29 nm),
siguiendo el procedimiento de Shen 2019:

  1) D = constante = 1 m^2/s
  2) Escanear DIFERENTES perfiles v(rho) (parametrizados por knots)
     Para cada uno, correr Aurora, computar nz(Fe22+), comparar con la
     forma espacial medida de la emision Fe XXIII.
  3) Elegir el mejor v(rho)  ==> obtenemos v/D
  4) Fijar el v/D optimo. Escanear D constante (0.1, 0.5, 1, 2, 5 m^2/s)
     para caracterizar sensibilidad.

Con una sola linea NO podemos separar D y v independientemente (necesitariamos
>= 3 estados de carga consecutivos), pero podemos:
  - determinar v/D
  - reportar un rango razonable de D

Datos requeridos:
    160869.sif
    Modulos importables desde ../02_shen_replication/

Uso:
    python fit_FeXXIII_shen_method.py
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sif_parser
import aurora
from scipy.optimize import differential_evolution
from scipy.interpolate import PchipInterpolator

# ---- paths ----
HERE     = Path(__file__).parent
SIF_PATH = HERE / "160869.sif"
OUT_FIG  = HERE.parent.parent / "Thesis" / "figures" / "fig_vD_result.png"

# All dependencies are in the same directory
sys.path.insert(0, str(HERE))
from profiles_fig4 import rho as rho_shen, ne_before, ne_ecrh, Te_before, Te_ecrh
from z_channels import Z_CM, CH_MIDPLANE

# ---- geometria EAST (hardcodeada) ----
R0, a, kappa = 1.85, 0.45, 1.80

# ---- config ----
FE_LINE_NM  = 13.29
LINE_HW_NM  = 0.06
T_A, T_B    = (2.0, 5.0), (6.0, 9.0)
T_OFFSET    = 0.12
# Z_MIN, Z_MAX ya no se usan; ver Z_CM en z_channels.py

# z_axis por fase, tomado del fit de C VI (transport-independent)
Z_AXIS_A_CM = 8.4    # ECRH ON
Z_AXIS_B_CM = 5.5    # ECRH OFF

D_CONST_M2S = 1.0                     # Shen initial guess
V_KNOT_RHOS = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
V_BOUNDS    = [(-3.0, 3.0)] * len(V_KNOT_RHOS)   # m/s por knot
D_SCAN      = np.array([0.1, 0.3, 0.5, 1.0, 2.0, 5.0])

# baseline Aurora
N0_EDGE, LAMBDA_N, SOURCE_RATE = 1.0e17, 0.10, 2.5e18


# ---- calibracion ----
_COEFS = np.polyfit([248, 376, 841], [33.73, 40.268, 33.73*2], 2)
def pix_to_wl_nm(p): return np.polyval(_COEFS, p) / 10.0


# =============================================================================
# 1) medir Fe XXIII(z) del .sif
# =============================================================================

def measure_FeXXIII_profile():
    data, info = sif_parser.np_open(str(SIF_PATH))
    n_t, n_sp, n_wl = data.shape
    cycle = info['CycleTime']
    t_frame = np.arange(n_t) * cycle + T_OFFSET
    data_bg = data - data[-1][None, :, :]
    wl_nm   = pix_to_wl_nm(np.arange(n_wl))

    mask_line = np.abs(wl_nm - FE_LINE_NM) < LINE_HW_NM
    # background continuo: promediar los 20 pixeles a cada lado sin la linea
    mask_bg = ((np.abs(wl_nm - FE_LINE_NM) > 2*LINE_HW_NM) &
               (np.abs(wl_nm - FE_LINE_NM) < 5*LINE_HW_NM))

    # emisividad Fe XXIII por (frame, canal) = suma en la linea - background local
    signal = data_bg[:, :, mask_line].sum(axis=2)              # (n_t, n_sp)
    bg     = data_bg[:, :, mask_bg].mean(axis=2) * mask_line.sum()
    fe_xxiii = signal - bg                                     # (n_t, n_sp)

    mask_A = (t_frame >= T_A[0]) & (t_frame <= T_A[1])
    mask_B = (t_frame >= T_B[0]) & (t_frame <= T_B[1])
    prof_A = fe_xxiii[mask_A].mean(axis=0)
    prof_B = fe_xxiii[mask_B].mean(axis=0)

    z_cm = Z_CM
    return z_cm, prof_A, prof_B


def z_to_rho(z_cm, z_ax_cm=0.0):
    """Fold z-> rho_tan usando el z_axis correcto para cada fase."""
    return np.clip(np.abs((z_cm - z_ax_cm)/100.0) / (kappa * a), 0.0, 0.99)


# =============================================================================
# 2) Aurora forward: dado v(rho) knots -> Fe22+ (rho)
# =============================================================================

def build_namelist(ne_prof, Te_prof):
    nml = aurora.default_nml.load_default_namelist()
    nml['device']         = 'EAST'
    nml['imp']            = 'Fe'
    nml['main_element']   = 'D'
    nml['cxr_flag']       = True
    nml['recycling_flag'] = True
    nml['wall_recycling'] = 0.9
    nml['Raxis_cm']  = R0 * 100.0
    nml['rvol_lcfs'] = a  * 100.0 * np.sqrt(kappa)
    nml['source_type']    = 'const'
    nml['source_rate']    = SOURCE_RATE
    nml['timing'] = dict(times=[0.0, 2.0], dt_start=[1e-5, 1e-5],
                         steps_per_cycle=[1, 1], dt_increase=[1.005, 1.0])
    nml['kin_profs']['ne']['fun']  = 'interp'
    nml['kin_profs']['ne']['rhop'] = rho_shen
    nml['kin_profs']['ne']['vals'] = np.maximum(ne_prof, 0.05) * 1.0e19 * 1.0e-6
    nml['kin_profs']['Te']['fun']  = 'interp'
    nml['kin_profs']['Te']['rhop'] = rho_shen
    nml['kin_profs']['Te']['vals'] = np.maximum(Te_prof, 0.02) * 1000.0
    n0 = N0_EDGE * np.exp(-(1.0 - rho_shen) * a / LAMBDA_N)
    nml['kin_profs']['n0']['fun']  = 'interpa'
    nml['kin_profs']['n0']['rhop'] = rho_shen
    nml['kin_profs']['n0']['vals'] = n0 * 1.0e-6
    return nml


def fe22_profile(ne_prof, Te_prof, v_knots, D_const_m2s):
    """Corre Aurora y devuelve rho, n_Fe22+(rho) [cm^-3]."""
    nml  = build_namelist(ne_prof, Te_prof)
    asim = aurora.aurora_sim(nml, geqdsk=None)

    # v(rho) via PCHIP en la grilla de Aurora
    spline = PchipInterpolator(V_KNOT_RHOS, v_knots)
    v_rho  = spline(asim.rhop_grid)                # m/s
    V = v_rho * 1.0e2                              # -> cm/s
    D = np.full_like(asim.rhop_grid, D_const_m2s * 1.0e4)   # cm^2/s

    out = asim.run_aurora(D, V)
    nz  = out['nz'][:, :, -1]                      # (rho, Z+1)
    return asim.rhop_grid, nz[:, 22]               # Fe22+


# =============================================================================
# 3) fit: encontrar v(rho) que reproduce la forma medida
# =============================================================================

def make_cost(rho_meas, y_meas, ne_prof, Te_prof, D_const=D_CONST_M2S):
    y_meas_n = y_meas / y_meas.max()

    def cost(v_knots):
        try:
            rho_sim, n22 = fe22_profile(ne_prof, Te_prof, v_knots, D_const)
        except Exception:
            return 1.0e6
        y_sim = np.interp(rho_meas, rho_sim, n22)
        y_sim_n = y_sim / max(y_sim.max(), 1e-30)
        return float(np.mean((y_sim_n - y_meas_n)**2))
    return cost


def fit_phase(rho_meas, y_meas, ne_prof, Te_prof, tag):
    print(f"\nFitting phase {tag}...  (D fixed = {D_CONST_M2S} m^2/s, "
          f"scanning v({V_KNOT_RHOS}) with bounds {V_BOUNDS[0]})")
    cost = make_cost(rho_meas, y_meas, ne_prof, Te_prof)
    res = differential_evolution(cost, V_BOUNDS, maxiter=25, popsize=8,
                                 seed=42, workers=1, tol=1e-3, updating='deferred')
    v_best = res.x
    print(f"  best cost = {res.fun:.4e}")
    print(f"  v_best (m/s) at rho {V_KNOT_RHOS}: {v_best}")

    # Perfil final con v_best
    rho_sim, n22 = fe22_profile(ne_prof, Te_prof, v_best, D_CONST_M2S)

    return {'v_best': v_best, 'cost': res.fun,
            'rho_sim': rho_sim, 'n22': n22,
            'v_smooth_rho': np.linspace(0, 1, 100),
            'v_smooth': PchipInterpolator(V_KNOT_RHOS, v_best)(np.linspace(0, 1, 100)),
            'ne_prof': ne_prof, 'Te_prof': Te_prof}


# =============================================================================
# 4) sensibilidad a D
# =============================================================================

def D_scan_sensitivity(fit_result, rho_meas, y_meas):
    """Con v/D fijo (i.e. v_best) escaneamos D constante."""
    y_meas_n = y_meas / y_meas.max()
    v_best   = fit_result['v_best']
    ne, Te   = fit_result['ne_prof'], fit_result['Te_prof']

    results = []
    for D in D_SCAN:
        # cuando escalas D, v/D se mantiene si escalas v tambien: v_new = D * (v_best/D0)
        v_scaled = v_best * (D / D_CONST_M2S)
        try:
            rho_sim, n22 = fe22_profile(ne, Te, v_scaled, D)
        except Exception:
            results.append({'D': D, 'cost': np.inf, 'rho': None, 'n22': None})
            continue
        y_sim = np.interp(rho_meas, rho_sim, n22)
        y_sim_n = y_sim / max(y_sim.max(), 1e-30)
        cost = float(np.mean((y_sim_n - y_meas_n)**2))
        results.append({'D': D, 'cost': cost, 'rho': rho_sim, 'n22': n22})
    return results


# =============================================================================
# 5) main
# =============================================================================

print("Loading measured Fe XXIII profile...")
z_cm, prof_A, prof_B = measure_FeXXIII_profile()

# folding con z_axis distinto por fase (de C VI)
rho_meas_A = z_to_rho(z_cm, z_ax_cm=Z_AXIS_A_CM)
rho_meas_B = z_to_rho(z_cm, z_ax_cm=Z_AXIS_B_CM)
print(f"  Phase A (z_ax={Z_AXIS_A_CM:+.1f}): rho range {rho_meas_A.min():.2f} to {rho_meas_A.max():.2f}")
print(f"  Phase B (z_ax={Z_AXIS_B_CM:+.1f}): rho range {rho_meas_B.min():.2f} to {rho_meas_B.max():.2f}")
print(f"  Phase A peak: {prof_A.max():.1f}  |  Phase B peak: {prof_B.max():.1f}")

fit_A = fit_phase(rho_meas_A, prof_A, ne_ecrh,   Te_ecrh,   'A (ECRH ON)')
fit_B = fit_phase(rho_meas_B, prof_B, ne_before, Te_before, 'B (ECRH OFF)')

print("\nD scan sensitivity (Phase A):")
scan_A = D_scan_sensitivity(fit_A, rho_meas_A, prof_A)
for r in scan_A:
    print(f"  D = {r['D']:>4.1f} m2/s  ->  cost = {r['cost']:.4e}")

print("\nD scan sensitivity (Phase B):")
scan_B = D_scan_sensitivity(fit_B, rho_meas_B, prof_B)
for r in scan_B:
    print(f"  D = {r['D']:>4.1f} m2/s  ->  cost = {r['cost']:.4e}")


# =============================================================================
# 6) figura
# =============================================================================

plt.rcParams.update({
    'font.size': 16, 'axes.labelsize': 18, 'axes.titlesize': 17,
    'legend.fontsize': 13, 'xtick.labelsize': 15, 'ytick.labelsize': 15,
    'axes.linewidth': 1.5,
})

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# mostrar solo la rama lower (z < z_axis) por fase — es la que tiene mas chords
# y llega mas afuera en rho, entonces cubre mejor el perfil
lower_A = z_cm < Z_AXIS_A_CM
lower_B = z_cm < Z_AXIS_B_CM

# (a) fit shape phase A
ax = axes[0]
prof_A_n = prof_A/prof_A.max()
ax.plot(rho_meas_A[lower_A], prof_A_n[lower_A], 'o', color='#c00000', ms=10,
        label='Measured Fe XXIII (norm.)')
ax.plot(fit_A['rho_sim'], fit_A['n22']/fit_A['n22'].max(), '-', color='k', lw=2,
        label='Aurora best fit')
ax.set_xlabel(r'$\rho$'); ax.set_ylabel('Fe XXIII intensity (norm.)')
ax.set_title('(a) Phase A: ECRH ON')
ax.set_xlim(0, 0.6); ax.set_ylim(0, 1.15)
ax.legend(loc='upper right', fontsize=12); ax.grid(alpha=0.3, ls=':')
ax.tick_params(direction='in', length=5)

# (b) fit shape phase B
ax = axes[1]
prof_B_n = prof_B/prof_B.max()
ax.plot(rho_meas_B[lower_B], prof_B_n[lower_B], 's', color='#1f4e79', ms=10,
        label='Measured Fe XXIII (norm.)')
ax.plot(fit_B['rho_sim'], fit_B['n22']/fit_B['n22'].max(), '-', color='k', lw=2,
        label='Aurora best fit')
ax.set_xlabel(r'$\rho$'); ax.set_ylabel('Fe XXIII intensity (norm.)')
ax.set_title('(b) Phase B: ECRH OFF')
ax.set_xlim(0, 0.6); ax.set_ylim(0, 1.15)
ax.legend(loc='upper right', fontsize=12); ax.grid(alpha=0.3, ls=':')
ax.tick_params(direction='in', length=5)

# (c) v/D profile obtenido
ax = axes[2]
vD_A = fit_A['v_smooth'] / D_CONST_M2S
vD_B = fit_B['v_smooth'] / D_CONST_M2S
ax.plot(fit_A['v_smooth_rho'], vD_A, '-', color='#c00000', lw=2.5,
        label='Phase A (ECRH ON)')
ax.plot(fit_B['v_smooth_rho'], vD_B, '-', color='#1f4e79', lw=2.5,
        label='Phase B (ECRH OFF)')
# knots
ax.plot(V_KNOT_RHOS, fit_A['v_best']/D_CONST_M2S, 'o', color='#c00000', ms=10)
ax.plot(V_KNOT_RHOS, fit_B['v_best']/D_CONST_M2S, 's', color='#1f4e79', ms=10)
ax.axhline(0, color='k', lw=0.5, ls='--')
ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$v/D$ (1/m)')
ax.set_title(f'(c) Extracted v/D profiles')
ax.set_xlim(0, 1); ax.legend(); ax.grid(alpha=0.3, ls=':')
ax.tick_params(direction='in', length=5)

plt.tight_layout()

OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_FIG, dpi=220, bbox_inches='tight')
print(f"\nSaved: {OUT_FIG}")


# =============================================================================
# 7) numeros clave
# =============================================================================

print("\n" + "="*60)
print("CAPTION-READY NUMBERS")
print("="*60)
for tag, fit in [('Phase A (ECRH ON)', fit_A), ('Phase B (ECRH OFF)', fit_B)]:
    print(f"\n{tag}:  (with D fixed = {D_CONST_M2S} m2/s)")
    print(f"  v(rho) knots (m/s): {dict(zip(V_KNOT_RHOS, np.round(fit['v_best'], 2)))}")
    print(f"  v/D at rho=0.2: {fit['v_smooth'][20]/D_CONST_M2S:+.2f} 1/m")
    print(f"  v/D at rho=0.5: {fit['v_smooth'][50]/D_CONST_M2S:+.2f} 1/m")
    print(f"  cost = {fit['cost']:.4e}")

best_D_A = D_SCAN[np.argmin([r['cost'] for r in scan_A])]
best_D_B = D_SCAN[np.argmin([r['cost'] for r in scan_B])]
print(f"\nD scan best values:")
print(f"  Phase A: D = {best_D_A} m2/s")
print(f"  Phase B: D = {best_D_B} m2/s")

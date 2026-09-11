"""
Pump-out verification for tungsten using Shen 2019 D, V coefficients.
Computes total W density (sum over all charge states) before and during ECRH.

If pump-out is real, we should see:
  1. Central total W decrease with ECRH
  2. Profile flatten (peaked -> flat)
  3. Volume-integrated W content decrease
"""

import matplotlib
matplotlib.use('TkAgg')   # force interactive backend
import numpy as np
import matplotlib.pyplot as plt
import aurora

from profiles_fig4 import rho as rho_kin, ne_before, ne_ecrh, Te_before, Te_ecrh
from transport_fig10_11 import transport
from geometry import R0, a, kappa


# ---------- knobs (same as run_aurora_W.py) ----------
N0_EDGE      = 1.0e17
LAMBDA_N     = 0.10
SOURCE_RATE  = 2.5e18


def build_namelist(ne_profile, Te_profile):
    nml = aurora.default_nml.load_default_namelist()
    nml['device']       = 'EAST'
    nml['imp']          = 'W'
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

    nml['kin_profs']['ne']['fun']  = 'interp'
    nml['kin_profs']['ne']['rhop'] = rho_kin
    nml['kin_profs']['ne']['vals'] = np.maximum(ne_profile, 0.05) * 1.0e19 * 1.0e-6

    nml['kin_profs']['Te']['fun']  = 'interp'
    nml['kin_profs']['Te']['rhop'] = rho_kin
    nml['kin_profs']['Te']['vals'] = np.maximum(Te_profile, 0.02) * 1000.0

    n0 = N0_EDGE * np.exp(-(1.0 - rho_kin) * a / LAMBDA_N)
    nml['kin_profs']['n0']['fun']  = 'interpa'
    nml['kin_profs']['n0']['rhop'] = rho_kin
    nml['kin_profs']['n0']['vals'] = n0 * 1.0e-6

    return nml


def run_case(case_name):
    if case_name == "Before ECRH":
        ne_p, Te_p = ne_before, Te_before
    else:
        ne_p, Te_p = ne_ecrh, Te_ecrh

    rho_DV  = transport[case_name]["rho"]
    D_knots = transport[case_name]["D"]
    V_knots = transport[case_name]["V"]

    nml = build_namelist(ne_p, Te_p)
    asim = aurora.aurora_sim(nml, geqdsk=None)

    D = np.interp(asim.rhop_grid, rho_DV, D_knots) * 1.0e4
    V = np.interp(asim.rhop_grid, rho_DV, V_knots) * 1.0e2

    out = asim.run_aurora(D, V)
    return asim.rhop_grid, out['nz'][:, :, -1] * 1.0e6 / 1.0e16   # 10^16 m^-3


# ---------- run both cases ----------
print("Running W transport - Before ECRH...")
rho_b, nz_b = run_case("Before ECRH")
print("Running W transport - On-axis ECRH...")
rho_e, nz_e = run_case("On-axis ECRH")


# ---------- pump-out metrics ----------
total_W_before = nz_b.sum(axis=1)
total_W_ecrh   = nz_e.sum(axis=1)

central_b    = total_W_before[0]
central_e    = total_W_ecrh[0]
integrated_b = np.trapz(total_W_before * rho_b, rho_b)   # proxy for volume-integrated content
integrated_e = np.trapz(total_W_ecrh   * rho_e, rho_e)

dC = 100 * (central_e - central_b) / central_b
dI = 100 * (integrated_e - integrated_b) / integrated_b

print("\n========== W PUMP-OUT CHECK ==========")
print(f"Central total W (rho=0):")
print(f"  Before ECRH:  {central_b:.4f} x 10^16 m^-3")
print(f"  On-axis ECRH: {central_e:.4f} x 10^16 m^-3")
print(f"  Change at center: {dC:+.1f}%")
print(f"\nVolume-weighted integrated W (proxy, int(n*rho dr)):")
print(f"  Before ECRH:  {integrated_b:.4f}")
print(f"  On-axis ECRH: {integrated_e:.4f}")
print(f"  Change: {dI:+.1f}%")
print("======================================\n")
if dC < -5:
    print("==> Pump-out CONFIRMED at center for total W.")
elif dC > 5:
    print("==> Center ENHANCEMENT - no pump-out for total W.")
else:
    print("==> Roughly unchanged at center - weak / no pump-out.")


# ---------- plots ----------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

axes[0].plot(rho_b, total_W_before, 'b-', lw=2, label='Before ECRH')
axes[0].plot(rho_e, total_W_ecrh,   'r-', lw=2, label='On-axis ECRH')
axes[0].set_xlabel(r'$\rho$'); axes[0].set_ylabel('total W (10$^{16}$ m$^{-3}$)')
axes[0].set_title(f'Total W density   (Δ at center: {dC:+.0f}%)')
axes[0].set_xlim(0, 0.8); axes[0].legend(); axes[0].grid(alpha=0.3)

ratio = total_W_ecrh / total_W_before.clip(min=1e-30)
axes[1].plot(rho_b, ratio, 'g-', lw=2)
axes[1].axhline(1.0, color='gray', ls='--', lw=0.8, label='no change')
axes[1].set_xlabel(r'$\rho$'); axes[1].set_ylabel('n_W(ECRH) / n_W(Before)')
axes[1].set_title('Ratio: <1 = pump-out,  >1 = enhancement')
axes[1].set_xlim(0, 0.8); axes[1].legend(); axes[1].grid(alpha=0.3)

plt.suptitle('W pump-out verification using Shen 2019 D, V')
plt.tight_layout()
plt.show(block=True)

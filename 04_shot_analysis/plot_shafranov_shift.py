"""
plot_shafranov_shift.py
=======================
Figura 2D de la seccion poloidal de EAST con superficies magneticas anidadas,
mostrando el corrimiento vertical del eje magnetico en el shot 160869.

Todo esta ploteado en CM (para que set_aspect('equal') tenga sentido).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from pathlib import Path
from z_channels import Z_CM

HERE     = Path(__file__).parent
OUT_FIG  = HERE.parent.parent / "Thesis" / "figures" / "fig_shafranov_shift.png"

# ---- EAST Miller geometry (en cm) ----
R0_cm    = 185.0
a_cm     = 45.0
kappa    = 1.80
delta    = 0.70

# ---- observed magnetic-axis positions (cm) ----
Z_AXIS_A = 8.4
Z_AXIS_B = 5.5
Z_PEAK_FE_A = 6.2
Z_PEAK_FE_B = 9.2


def miller_surface(rho, theta, z_shift=0.0):
    """R,Z in cm."""
    alpha = theta + np.arcsin(delta) * np.sin(theta)
    R = R0_cm + a_cm * rho * np.cos(alpha)
    Z = kappa * a_cm * rho * np.sin(theta) + z_shift
    return R, Z


theta = np.linspace(0, 2*np.pi, 400)
rhos  = np.arange(0.1, 1.01, 0.1)

plt.rcParams.update({
    'font.size':        16, 'axes.labelsize':  18, 'axes.titlesize': 18,
    'legend.fontsize':  13, 'xtick.labelsize': 14, 'ytick.labelsize': 14,
    'axes.linewidth':   1.6,
    'font.family':      'DejaVu Sans',
})

# figsize: subplot ~135 cm ancho x 200 cm alto → ratio 1:1.48
fig, axes = plt.subplots(1, 2, figsize=(16, 11))

for ax, z_shift, tag, color, z_peak in [
    (axes[0], Z_AXIS_A, 'A (ECRH ON)',  '#c00000', Z_PEAK_FE_A),
    (axes[1], Z_AXIS_B, 'B (ECRH OFF)', '#1f4e79', Z_PEAK_FE_B),
]:
    # -- nominal Miller surfaces (Z_axis = 0) --
    for rho in rhos:
        R_n, Z_n = miller_surface(rho, theta, z_shift=0.0)
        lw = 1.5 if rho == 1.0 else 0.8
        ls = '-'  if rho == 1.0 else '--'
        c  = 'dimgray' if rho == 1.0 else 'lightgrey'
        ax.plot(R_n, Z_n, color=c, lw=lw, ls=ls, zorder=1)

    # -- shifted flux surfaces --
    for rho in rhos:
        R_s, Z_s = miller_surface(rho, theta, z_shift=z_shift)
        lw = 2.2 if rho == 1.0 else 1.1
        ax.plot(R_s, Z_s, color=color, lw=lw, ls='-', alpha=0.85, zorder=2)

    # -- magnetic axes markers --
    ax.plot(R0_cm, 0,       'x', color='dimgray', ms=16, mew=3.0, zorder=6)
    ax.plot(R0_cm, z_shift, 'o', color=color,     ms=16, mec='k',  mew=1.8, zorder=7)

    # -- shift arrow: placed to the LEFT of the plasma so it doesn't overlap the surfaces --
    R_arrow = 135.0                                # x-position of the arrow (left of plasma)
    # thin horizontal reference lines at z=0 and z=z_shift, connecting arrow to axis
    ax.plot([R_arrow, R0_cm], [0, 0],
            color='dimgray', ls=':', lw=1.2, alpha=0.7, zorder=3)
    ax.plot([R_arrow, R0_cm], [z_shift, z_shift],
            color=color, ls=':', lw=1.2, alpha=0.7, zorder=3)
    # vertical arrow between the two horizontal reference lines
    ax.annotate('', xy=(R_arrow, z_shift), xytext=(R_arrow, 0),
                arrowprops=dict(arrowstyle='<->', color=color, lw=2.8, mutation_scale=22),
                zorder=5)
    ax.text(R_arrow - 3, z_shift/2, f'$\\Delta z = {z_shift:+.1f}$ cm',
            color=color, fontweight='bold', fontsize=15,
            ha='right', va='center', rotation=90)

    # -- Fe XXIII peak: short marker on the right side + label, without a full horizontal line --
    ax.plot([236, 246], [z_peak, z_peak], color='k', lw=2.0, alpha=0.75, zorder=4)
    ax.text(248, z_peak, f'Fe XXIII peak\n$z = {z_peak:+.1f}$ cm',
            fontsize=12, ha='left', va='center', color='k', fontweight='bold')

    # -- labels & title --
    ax.set_xlabel('R (cm)')
    ax.set_title(f'Phase {tag} — magnetic axis at $z = {z_shift:+.1f}$ cm',
                 color=color, fontweight='bold')

    legend_handles = [
        Line2D([0], [0], color='lightgrey', ls='--', lw=1.6,
               label='Nominal surfaces ($z_{\\rm ax}=0$)'),
        Line2D([0], [0], color=color, lw=2.6,
               label='Shifted surfaces (observed)'),
        Line2D([0], [0], marker='x', color='dimgray', mew=3.0, ms=13,
               lw=0, label='Geometric axis'),
        Line2D([0], [0], marker='o', color=color, mec='k', mew=1.8, ms=13,
               lw=0, label='Magnetic axis (C VI fit)'),
    ]
    ax.legend(handles=legend_handles, loc='lower left',
              framealpha=0.95, fontsize=13)

    # NOW both axes are in cm -> aspect equal works properly
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(115, 285)
    ax.set_ylim(-95, 105)

    ax.grid(alpha=0.25, ls=':')
    ax.tick_params(direction='in', length=5, top=True, right=True)

axes[0].set_ylabel('Z (cm)')


plt.tight_layout()
OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_FIG, dpi=220, bbox_inches='tight')
print(f"Saved: {OUT_FIG}")

print()
print("=" * 60)
print("Vertical axis displacement (from C VI forward-projection fit)")
print("=" * 60)
print(f"Phase A (ECRH ON):  z_axis = {Z_AXIS_A:+.1f} cm")
print(f"Phase B (ECRH OFF): z_axis = {Z_AXIS_B:+.1f} cm")
print(f"Delta (A - B)      = {Z_AXIS_A - Z_AXIS_B:+.1f} cm")

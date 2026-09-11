"""
Shen et al. 2019 (Phys. Plasmas 26, 032507) - Figure 4 profiles
Mo impurity suppression with on-axis ECRH at EAST. Shot 55339.

ne and Te before ECRH (ohmic) and during on-axis ECRH.
Values digitized by eye from Fig 4 - adjust as needed.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from pathlib import Path

HERE = Path(__file__).parent


# Radial coordinate
rho = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

# ne [10^19 m^-3]
# Before ECRH: peaked, higher in core
# During ECRH: flattened, slightly lower in core ("particle pump-out")
ne_before = np.array([4.50, 4.50, 4.40, 4.10, 3.65, 3.20, 2.85, 2.55, 2.30])
ne_ecrh   = np.array([3.80,3.80, 3.80, 3.75, 3.65, 3.40, 3, 2.55, 2.30])

# Te [keV]
# Before ECRH: lower, less peaked
# During ECRH: higher and more peaked in core
Te_before = np.array([0.80, 0.75, 0.70, 0.62, 0.55, 0.50, 0.45, 0.40, 0.30])
Te_ecrh   = np.array([1.70, 1.55, 1.30, 1.05, 0.85, 0.70, 0.55, 0.45, 0.35])


def interpolate(rho_knots, y_knots, n=200):
    rho_fine = np.linspace(rho_knots.min(), rho_knots.max(), n)
    spline = PchipInterpolator(rho_knots, y_knots)
    return rho_fine, spline(rho_fine)


def _draw(ax, rho_k, y_k, color, label):
    rho_f, y_f = interpolate(rho_k, y_k)
    ax.plot(rho_f, y_f, '-', color=color, lw=2, label=label)
    ax.plot(rho_k, y_k, 'o', color=color, ms=4)


def plot_profiles():
    fig, axes = plt.subplots(2, 1, figsize=(6, 7), sharex=True)

    _draw(axes[0], rho, ne_before, 'r', 'before ECRH')
    _draw(axes[0], rho, ne_ecrh,   'g', 'on-axis ECRH')
    axes[0].set_ylabel(r'$n_e$ (10$^{19}$ m$^{-3}$)')
    axes[0].set_ylim(2, 5)
    axes[0].legend(); axes[0].grid(alpha=0.3)

    _draw(axes[1], rho, Te_before, 'r', 'before ECRH')
    _draw(axes[1], rho, Te_ecrh,   'g', 'on-axis ECRH')
    axes[1].set_ylabel(r'$T_e$ (keV)')
    axes[1].set_xlabel(r'$\rho$')
    axes[1].set_ylim(0, 2)
    axes[1].set_xlim(0, 0.8)
    axes[1].legend(); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(HERE / 'profiles_fig4_check.png', dpi=120)
    plt.show()


if __name__ == '__main__':
    plot_profiles()

"""
Vogel et al. - Figure 4 profiles digitized
ne, Te, Ti vs rho, with and without RMP.
Values read by eye from the figure - adjust as needed.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
from scipy.interpolate import PchipInterpolator

# Radial coordinate
rho = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0])

# ne [10^19 m^-3]
ne_noRMP   = np.array([3.05, 3.00, 2.85, 2.70, 2.55, 2.40, 2.25, 2.05, 1.90, 1.55, 1.00, 0, 0])
ne_withRMP = np.array([3.50, 3.45, 3.30, 3.10, 2.90, 2.70, 2.55, 2.30, 2.10, 1.70, 1.10, 0, 0])

# Te [keV]
Te_noRMP   = np.array([3.15, 3, 2.20, 1.85, 1.55, 1.30, 1.10, 0.90, 0.70, 0.60, 0.50, 0.35, 0.15])
Te_withRMP = np.array([3.50, 3.54, 2.60, 1.9, 1.80, 1.55, 1.30, 1.10, 0.95, 0.85, 0.75, 0.55, 0.40])

# Ti [keV] - profiles overlap in figure, so same array used; truncate at rho ~0.8 (no data beyond)
rho_Ti = np.array([0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80])
Ti_noRMP   = np.array([1.10, 1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75])
Ti_withRMP = np.array([1.10, 1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75])


def interpolate(rho_knots, y_knots, n=200):
    """Smooth PCHIP interpolation onto a fine grid."""
    rho_fine = np.linspace(rho_knots.min(), rho_knots.max(), n)
    spline = PchipInterpolator(rho_knots, y_knots)
    return rho_fine, spline(rho_fine)


def _draw(ax, rho_k, y_k, color, label):
    rho_f, y_f = interpolate(rho_k, y_k)
    ax.plot(rho_f, y_f, '-', color=color, lw=2, label=label)
    ax.plot(rho_k, y_k, 'o', color=color, ms=4)


def plot_profiles():
    fig, axes = plt.subplots(3, 1, figsize=(6, 9), sharex=True)

    _draw(axes[0], rho, ne_noRMP,   'k', 'no RMP')
    _draw(axes[0], rho, ne_withRMP, 'c', 'with RMP')
    axes[0].set_ylabel(r'$n_e$ (10$^{19}$ m$^{-3}$)')
    axes[0].set_ylim(0, 4)
    axes[0].legend(); axes[0].grid(alpha=0.3)

    _draw(axes[1], rho, Te_noRMP,   'k', 'no RMP')
    _draw(axes[1], rho, Te_withRMP, 'c', 'with RMP')
    axes[1].set_ylabel(r'$T_e$ (keV)')
    axes[1].set_ylim(0, 4)
    axes[1].legend(); axes[1].grid(alpha=0.3)

    _draw(axes[2], rho_Ti, Ti_noRMP,   'k', 'no RMP')
    _draw(axes[2], rho_Ti, Ti_withRMP, 'c', 'with RMP')
    axes[2].set_ylabel(r'$T_i$ (keV)')
    axes[2].set_xlabel(r'$\rho$')
    axes[2].set_ylim(0.4, 1.2)
    axes[2].set_xlim(0, 1)
    axes[2].legend(); axes[2].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(HERE / 'profiles_fig4_check.png', dpi=120)
    plt.show()


if __name__ == '__main__':
    plot_profiles()

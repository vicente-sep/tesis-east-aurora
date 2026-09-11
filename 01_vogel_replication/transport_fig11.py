"""
Vogel et al. - Figure 11 transport coefficients digitized
D_Fe and v_Fe vs rho, with and without RMP. Calculated with STRAHL.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
from scipy.interpolate import PchipInterpolator

rho = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

transport = {
    "No RMP": {
        "rho": rho,
        # m^2/s
        "D":   np.array([0.55, 0.60, 0.60, 0.75, 0.85, 1.05, 1.30, 1.60, 1.90, 1.95, 2.00]),
        # m/s
        "V":   np.array([0.00, -0.10, -0.25, -0.60, -1.00, -0.40, 0.20, 0.30, 0.35, 0.30, 0.20]),
    },
    "With RMP": {
        "rho": rho,
        "D":   np.array([0.65, 0.75, 0.90, 1.00, 1.10, 1.50, 1.90, 2.50, 3.00, 3.15, 3.30]),
        "V":   np.array([0.00, -0.05, -0.10, -0.30, -0.50, 0.70, 0.75, 1.20, 0.7, 0.55, 0.50]),
    },
}


def interpolate(rho_knots, y_knots, n=200):
    """Smooth PCHIP interpolation onto a fine grid."""
    rho_fine = np.linspace(rho_knots.min(), rho_knots.max(), n)
    spline = PchipInterpolator(rho_knots, y_knots)
    return rho_fine, spline(rho_fine)


def plot_transport():
    fig, axes = plt.subplots(2, 1, figsize=(6, 7), sharex=True)

    colors = {"No RMP": "k", "With RMP": "r"}

    for label, data in transport.items():
        c = colors[label]
        # smooth curves
        rho_f, D_f = interpolate(data["rho"], data["D"])
        _,     V_f = interpolate(data["rho"], data["V"])
        axes[0].plot(rho_f, D_f, '-', color=c, lw=2, label=label)
        axes[1].plot(rho_f, V_f, '-', color=c, lw=2, label=label)
        # original knots
        axes[0].plot(data["rho"], data["D"], 'o', color=c, ms=5)
        axes[1].plot(data["rho"], data["V"], 'o', color=c, ms=5)

    axes[0].set_ylabel(r'$D_{Fe}$ (m$^2$/s)')
    axes[0].set_ylim(0, 3.5)
    axes[0].legend(); axes[0].grid(alpha=0.3)

    axes[1].axhline(0, color='gray', ls='--', lw=0.8)
    axes[1].set_ylabel(r'$v_{Fe}$ (m/s)')
    axes[1].set_xlabel(r'$\rho$')
    axes[1].set_ylim(-1.2, 1.3)
    axes[1].set_xlim(0, 1)
    axes[1].legend(); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(HERE / 'transport_fig11_check.png', dpi=120)
    plt.show()


if __name__ == '__main__':
    plot_transport()

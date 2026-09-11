"""
Shen et al. 2019 - Figures 10 (D) and 11 (V) transport coefficients
for Mo impurity, before and during on-axis ECRH. Shot 55339.

Values read by eye from the figures - adjust as needed.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from pathlib import Path

HERE = Path(__file__).parent


rho = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

transport = {
    "Before ECRH": {
        "rho": rho,
        # D [m^2/s] - Fig 10, blue
        "D":   np.array([0.30, 0.34, 0.43, 0.50, 0.95, 1.00, 1.11, 1.22, 1.40]),
        # V [m/s] - Fig 11, blue (inward pinch everywhere)
        "V":   np.array([0.00, -0.50, -0.85, -0.60, -0.34, -0.25, -0.18, -0.14, -0.10]),
    },
    "On-axis ECRH": {
        "rho": rho,
        # D [m^2/s] - Fig 10, red (higher in core than before ECRH)
        "D":   np.array([0.50, 0.56, 0.70, 0.91, 1.10, 1.10, 1.21, 1.30, 1.40]),
        # V [m/s] - Fig 11, red (outward in core, inward outside)
        "V":   np.array([0.00, 0.12, 0.25, -0.19, -0.20, -0.18, -0.10, -0.11, -0.10]),
    },
}


def interpolate(rho_knots, y_knots, n=200):
    rho_fine = np.linspace(rho_knots.min(), rho_knots.max(), n)
    spline = PchipInterpolator(rho_knots, y_knots)
    return rho_fine, spline(rho_fine)


def plot_transport():
    fig, axes = plt.subplots(2, 1, figsize=(6, 7), sharex=True)

    colors = {"Before ECRH": "b", "On-axis ECRH": "r"}

    for label, data in transport.items():
        c = colors[label]
        rho_f, D_f = interpolate(data["rho"], data["D"])
        _,     V_f = interpolate(data["rho"], data["V"])
        axes[0].plot(rho_f, D_f, '-', color=c, lw=2, label=label)
        axes[1].plot(rho_f, V_f, '-', color=c, lw=2, label=label)
        axes[0].plot(data["rho"], data["D"], 'o', color=c, ms=5)
        axes[1].plot(data["rho"], data["V"], 'o', color=c, ms=5)

    axes[0].set_ylabel(r'$D_{Mo}$ (m$^2$/s)')
    axes[0].set_ylim(0, 2)
    axes[0].legend(); axes[0].grid(alpha=0.3)

    axes[1].axhline(0, color='gray', ls='--', lw=0.8)
    axes[1].set_ylabel(r'$v_{Mo}$ (m/s)')
    axes[1].set_xlabel(r'$\rho$')
    axes[1].set_ylim(-1, 0.5)
    axes[1].set_xlim(0, 0.8)
    axes[1].legend(); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(HERE / 'transport_fig10_11_check.png', dpi=120)
    plt.show()


if __name__ == '__main__':
    plot_transport()

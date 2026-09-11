"""
Digitized Mo charge state profiles from Shen et al. 2019, Figs 8 & 9.
These are the targets the parameter sweep tries to match.

Units: 10^16 m^-3
Edit values to match the paper figures more precisely.
"""

import numpy as np
import matplotlib.pyplot as plt

# Common rho grid for digitized targets
rho_target = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

# ---- Fig 8: Before ECRH (peaked, higher central density) ----
targets_before = {
    28: np.array([1.500, 1.450, 1.300, 0.850, 0.250, 0.100, 0.070, 0.060, 0.050]),
    29: np.array([1.250, 1.220, 1.100, 0.700, 0.250, 0.120, 0.080, 0.060, 0.040]),
    30: np.array([0.900, 0.870, 0.750, 0.400, 0.120, 0.060, 0.040, 0.030, 0.020]),
    31: np.array([0.600, 0.580, 0.480, 0.250, 0.100, 0.060, 0.040, 0.030, 0.030]),
}

# ---- Fig 9: On-axis ECRH (flatter, lower central density) ----
targets_ecrh = {
    28: np.array([0.075, 0.080, 0.076, 0.066, 0.049, 0.030, 0.023, 0.018, 0.013]),
    29: np.array([0.100, 0.105, 0.100, 0.075, 0.040, 0.025, 0.020, 0.015, 0.012]),
    30: np.array([0.209, 0.210, 0.190, 0.150, 0.090, 0.060, 0.040, 0.030, 0.020]),
    31: np.array([0.180, 0.180, 0.160, 0.117, 0.075, 0.051, 0.041, 0.036, 0.020]),
}


def plot_targets():
    fig, axes = plt.subplots(4, 2, figsize=(11, 11), sharex=True)
    charge_states = [28, 29, 30, 31]
    for row, Z in enumerate(charge_states):
        axes[row, 0].plot(rho_target, targets_before[Z], 'b-o', lw=2, ms=5)
        axes[row, 0].set_ylabel(f'Mo$^{{{Z}+}}$ (10$^{{16}}$ m$^{{-3}}$)')
        axes[row, 0].grid(alpha=0.3)

        axes[row, 1].plot(rho_target, targets_ecrh[Z], 'r-o', lw=2, ms=5)
        axes[row, 1].grid(alpha=0.3)

    axes[0, 0].set_title('Before ECRH (Fig 8 target)')
    axes[0, 1].set_title('On-axis ECRH (Fig 9 target)')
    axes[-1, 0].set_xlabel(r'$\rho$')
    axes[-1, 1].set_xlabel(r'$\rho$')
    axes[-1, 0].set_xlim(0, 0.8)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plot_targets()

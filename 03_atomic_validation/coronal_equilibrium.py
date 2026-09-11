"""
Coronal equilibrium fractional abundance vs Te.
Used to validate Aurora's atomic data against published curves
(e.g. Hu et al. 2018 Figs 3-5 for Ar and Xe).

Also can include charge exchange recombination effect.
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Force interactive matplotlib popup window
import matplotlib.pyplot as plt
plt.ioff()  # Disable interactive mode to ensure plt.show() blocks
import aurora
from pathlib import Path

HERE = Path(__file__).parent


# ---------- impurities and charge states of interest ----------
IMPURITIES = {
    'Ar': {
        'Z':           18,
        'Te_range':    (0.1e3, 10e3),
        'highlight':   [14, 15, 16, 17, 18],
        'reference':   'Hu et al. 2018, Fig 3',
    },
    'Mo': {
        'Z':           42,
        'Te_range':    (0.1e3, 10e3),
        'highlight':   [28, 29, 30, 31, 32],
        'reference':   'no paper reference - using IE check',
        # NIST ionization energies (eV) for Mo Z+ -> (Z+1)+
        'IE': {28: 1290, 29: 1480, 30: 1726, 31: 1815, 32: 4280},
    },
}


def coronal_fractions(imp, Te_eV, ne_cm3=1e14, n0_ne_ratio=0.0):
    """
    Compute coronal equilibrium fractional abundances for impurity `imp`.

    Args:
        imp: element symbol ('Ar', 'Fe', 'Mo', 'Xe', ...)
        Te_eV: array of electron temperatures [eV]
        ne_cm3: electron density [cm^-3] (scalar)
        n0_ne_ratio: ratio of neutral density to electron density,
                     used for charge exchange. 0 disables CX.

    Returns:
        Te_eV: input Te (for convenience)
        frac:  array shape (Z+1, len(Te)) - fractional abundance per charge state
    """
    atom_data = aurora.atomic.get_atom_data(imp, files=['scd', 'acd', 'ccd'])

    ne_arr = np.full_like(Te_eV, ne_cm3)
    n0_arr = ne_arr * n0_ne_ratio if n0_ne_ratio > 0 else None

    # get_frac_abundances signature: (atom_data, ne_cm3, Te_eV, n0_by_ne=None, ...)
    Te_out, frac = aurora.atomic.get_frac_abundances(
        atom_data, ne_arr, Te_eV,
        n0_by_ne=n0_arr / ne_arr if n0_arr is not None else None,
        plot=False,
    )
    return Te_out, frac.T   # transpose so shape is (charge_state, Te)


def plot_for_impurity(imp_name, with_cx=False, n0_by_ne=1e-4):
    cfg = IMPURITIES[imp_name]
    Te = np.logspace(np.log10(cfg['Te_range'][0]),
                     np.log10(cfg['Te_range'][1]),
                     200)

    Te_out, frac_noCX = coronal_fractions(imp_name, Te, n0_ne_ratio=0.0)
    if with_cx:
        _, frac_CX = coronal_fractions(imp_name, Te, n0_ne_ratio=n0_by_ne)

    fig, ax = plt.subplots(figsize=(8, 5))

    # plot all charge states faintly
    for Z in range(cfg['Z'] + 1):
        ax.plot(Te / 1000, frac_noCX[Z], color='lightgray', lw=0.5)

    # highlight key charge states
    cmap = plt.cm.tab10
    for i, Z in enumerate(cfg['highlight']):
        c = cmap(i)
        ax.plot(Te / 1000, frac_noCX[Z], color=c, lw=2,
                label=f'{imp_name}$^{{{Z}+}}$')
        if with_cx:
            ax.plot(Te / 1000, frac_CX[Z], color=c, lw=1.5, ls='--',
                    alpha=0.7)

    # overlay ionization-energy markers if available
    if 'IE' in cfg:
        for Z, Ei in cfg['IE'].items():
            if Z in cfg['highlight']:
                ax.axvline(Ei / 1000, ls=':', color='gray', alpha=0.5)
                ax.text(Ei / 1000, 1.02, f'IE({Z}+)', rotation=90,
                        fontsize=7, ha='right', va='bottom', color='gray')

    ax.set_xlabel('Electron temperature (keV)')
    ax.set_xlim(0, 10)
    ax.set_ylabel('Fractional abundance')
    ax.set_title(f'{imp_name} coronal equilibrium  '
                 f'(ref: {cfg["reference"]})'
                 + (f'\nDashed: with CX, n0/ne={n0_by_ne}' if with_cx else ''))
    ax.set_ylim(0, 1.0)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig_path = HERE / f'coronal_{imp_name}{"_withCX" if with_cx else ""}.png'
    plt.savefig(fig_path, dpi=120)
    return fig


if __name__ == '__main__':
    # Run for all impurities, with and without CX
    for imp in IMPURITIES:
        plot_for_impurity(imp, with_cx=False)
        plot_for_impurity(imp, with_cx=True, n0_by_ne=1e-4)
    plt.show(block=True)


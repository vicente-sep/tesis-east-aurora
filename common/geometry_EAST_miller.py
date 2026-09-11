"""
Miller analytic geometry for EAST.
Builds flux surfaces R(rho, theta), Z(rho, theta) and computes
volume V(rho) and surface area A(rho) for Aurora.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent

# ---- EAST Miller parameters ----
R0    = 1.85    # major radius [m]
a     = 0.45    # minor radius [m]
kappa = 1.80    # elongation (EAST typical 1.6-2)
delta = 0.70    # triangularity (EAST typical 0.6-0.8)


def flux_surface(rho, theta):
    """R(rho,theta), Z(rho,theta) for Miller geometry."""
    alpha = theta + np.arcsin(delta) * np.sin(theta)
    R = R0 + a * rho * np.cos(alpha)
    Z = kappa * a * rho * np.sin(theta)
    return R, Z


def _derivatives(rho, theta):
    """Partial derivatives needed for volume/area integrals."""
    s = np.arcsin(delta)
    alpha = theta + s * np.sin(theta)
    dalpha_dtheta = 1.0 + s * np.cos(theta)

    dR_drho   = a * np.cos(alpha)
    dR_dtheta = -a * rho * np.sin(alpha) * dalpha_dtheta
    dZ_drho   = kappa * a * np.sin(theta)
    dZ_dtheta = kappa * a * rho * np.cos(theta)
    return dR_drho, dR_dtheta, dZ_drho, dZ_dtheta


def volume_and_area(rho_grid, n_theta=400):
    """V(rho), A(rho) by integration over poloidal angle."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta)
    dtheta = theta[1] - theta[0]

    V = np.zeros_like(rho_grid)
    A = np.zeros_like(rho_grid)

    for i, rho in enumerate(rho_grid):
        R, _ = flux_surface(rho, theta)
        dRdr, dRdt, dZdr, dZdt = _derivatives(rho, theta)

        # Jacobian for cross-section element
        J = np.abs(dRdr * dZdt - dRdt * dZdr)
        # cross-sectional area enclosed (for V) - we need ∫∫ 2πR·J drho dtheta
        # but since V is enclosed up to rho, integrate over rho' from 0 to rho:
        # easier: use ∫(2πR · dA_cross) where dA_cross at surface = (1/2)|R dZ - Z dR| dtheta -> Pappus-like
        # Direct way: compute V(rho) = ∫₀^rho ∫ 2πR J dtheta drho' -> we accumulate below

        # Surface arc length element ds (poloidal)
        ds = np.sqrt(dRdt**2 + dZdt**2)

        # Toroidal surface area of the flux surface
        A[i] = np.sum(2.0 * np.pi * R * ds) * dtheta

    # Volume: integrate cross-section area enclosed * 2π·R_centroid...
    # cleanest: cumulative integral of dV/drho = ∫ 2πR J dtheta
    dV_drho = np.zeros_like(rho_grid)
    for i, rho in enumerate(rho_grid):
        R, _ = flux_surface(rho, theta)
        dRdr, dRdt, dZdr, dZdt = _derivatives(rho, theta)
        J = np.abs(dRdr * dZdt - dRdt * dZdr)
        dV_drho[i] = np.sum(2.0 * np.pi * R * J) * dtheta

    V = np.concatenate(([0.0], np.cumsum(0.5 * (dV_drho[:-1] + dV_drho[1:]) * np.diff(rho_grid))))
    return V, A


def plot_geometry():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # cross-section: nested flux surfaces
    theta = np.linspace(0, 2 * np.pi, 200)
    for rho in np.linspace(0.1, 1.0, 6):
        R, Z = flux_surface(rho, theta)
        axes[0].plot(R, Z, lw=1)
    axes[0].set_aspect('equal')
    axes[0].set_xlabel('R [m]'); axes[0].set_ylabel('Z [m]')
    axes[0].set_title(f'EAST: R0={R0}, a={a}, κ={kappa}, δ={delta}')
    axes[0].grid(alpha=0.3)

    # V(rho), A(rho)
    rho = np.linspace(0.001, 1.0, 50)
    V, A = volume_and_area(rho)
    axes[1].plot(rho, V, 'b-', lw=2)
    axes[1].set_xlabel(r'$\rho$'); axes[1].set_ylabel(r'$V$ [m$^3$]')
    axes[1].set_title(f'Total V = {V[-1]:.2f} m³')
    axes[1].grid(alpha=0.3)

    axes[2].plot(rho, A, 'r-', lw=2)
    axes[2].set_xlabel(r'$\rho$'); axes[2].set_ylabel(r'$A$ [m$^2$]')
    axes[2].set_title(f'LCFS A = {A[-1]:.2f} m²')
    axes[2].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(HERE / 'geometry_check.png', dpi=120)
    plt.show()


if __name__ == '__main__':
    plot_geometry()

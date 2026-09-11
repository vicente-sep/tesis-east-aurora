"""
Wavelength calibration for the EAST XEUV_Long spectrometer.
From biaodingbochang.m (3 anchor points, quadratic fit).
"""

import numpy as np

# Anchor points (from biaodingbochang.m)
_PIX_ANCHOR    = np.array([248, 376, 841])
_WL_A_ANCHOR   = np.array([33.73, 40.268, 33.73 * 2])   # Ångström (33.73*2 = 2nd order C VI)

# Polynomial coefficients: wavelength_Å = a*pix^2 + b*pix + c
_COEFS = np.polyfit(_PIX_ANCHOR, _WL_A_ANCHOR, 2)


def pix_to_wl_A(pix):
    """Convert pixel number to wavelength in Ångström."""
    return np.polyval(_COEFS, np.asarray(pix))


def pix_to_wl_nm(pix):
    """Convert pixel number to wavelength in nm."""
    return pix_to_wl_A(pix) / 10.0


def wl_nm_to_pix(wl_nm):
    """Inverse: from wavelength in nm to (fractional) pixel."""
    # invert quadratic: c + b*pix + a*pix^2 = wl_A  →  a*pix^2 + b*pix + (c - wl_A) = 0
    wl_A = np.asarray(wl_nm) * 10.0
    a, b, c = _COEFS
    disc = b**2 - 4*a*(c - wl_A)
    return (-b + np.sqrt(disc)) / (2*a)


# Standard pixel grid for a 2048-pixel detector
PIXELS_2048 = np.arange(2048)
WL_NM_2048  = pix_to_wl_nm(PIXELS_2048)


if __name__ == '__main__':
    print("EAST XEUV_Long calibration")
    print(f"Polynomial coefficients (Å = a*p^2 + b*p + c):")
    print(f"  a = {_COEFS[0]:.6e}")
    print(f"  b = {_COEFS[1]:.6e}")
    print(f"  c = {_COEFS[2]:.6e}")
    print(f"\nWavelength range: {WL_NM_2048.min():.2f} - {WL_NM_2048.max():.2f} nm")
    print(f"Resolution at midrange: {WL_NM_2048[1024] - WL_NM_2048[1023]:.4f} nm/pixel")

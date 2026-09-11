"""
Digitized Fe charge state profiles from Vogel 2021, Fig 10.
Units: 10^15 m^-3
"""

import numpy as np

rho_target = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

# Fig 10: black (no RMP) curves - centrally peaked
targets_noRMP = {
    20: np.array([3.30, 3.20, 3.00, 2.70, 2.30, 1.90, 1.50, 1.10, 0.80, 0.60, 0.50]),
    21: np.array([3.00, 2.90, 2.70, 2.40, 2.00, 1.60, 1.25, 0.90, 0.60, 0.45, 0.35]),
    22: np.array([2.50, 2.40, 2.10, 1.70, 1.25, 0.85, 0.55, 0.35, 0.20, 0.13, 0.10]),
}

# Fig 10: red (with RMP) curves - lower and slightly flatter
targets_withRMP = {
    20: np.array([2.00, 1.95, 1.80, 1.65, 1.40, 1.15, 0.90, 0.70, 0.55, 0.45, 0.40]),
    21: np.array([1.80, 1.75, 1.60, 1.45, 1.20, 0.95, 0.75, 0.55, 0.40, 0.30, 0.25]),
    22: np.array([1.50, 1.45, 1.30, 1.05, 0.80, 0.55, 0.35, 0.22, 0.13, 0.08, 0.06]),
}

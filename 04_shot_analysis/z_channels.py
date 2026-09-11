"""
z_channels.py
=============
Posiciones verticales Z de los 20 canales espaciales del espectrometro XEUV_Long
de EAST, tomadas del array Z_m del codigo canonico de Zichao.

IMPORTANTE: el array NO es uniforme. Hay un gap entre los canales 2 y 3
(probablemente por un canal muerto o un binning diferente del CCD).

El canal mas cercano al midplane (Z = 0) es el canal 8 (Z = +0.09 cm).
"""

import numpy as np

# Z de cada canal en METROS (de Zichao)
Z_M = np.array([
    -0.2721, -0.2418, -0.2114, -0.1508, -0.1205,
    -0.0901, -0.0598, -0.0295,  0.0009,  0.0312,
     0.0616,  0.0918,  0.1222,  0.1525,  0.1828,
     0.2132,  0.2435,  0.2738,  0.3041,  0.3344,
])

# En centimetros (mas conveniente para graficar)
Z_CM = Z_M * 100.0

# Numero de canales
N_CH = len(Z_M)   # 20

# Canal mas cercano al midplane
CH_MIDPLANE = int(np.argmin(np.abs(Z_M)))   # = 8, en z = +0.09 cm

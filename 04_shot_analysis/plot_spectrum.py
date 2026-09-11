"""
plot_spectrum.py
================
Genera 'fig_spectrum.png' para el capitulo 6 de la tesis:
espectro EUV del shot 160869 en la cuerda midplane, comparando
ECRH ON vs ECRH OFF, con lineas espectrales conocidas etiquetadas.

Datos requeridos:
    160869.sif  - datos crudos del espectrometro XEUV_Long

Uso:
    python plot_spectrum.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sif_parser
from pathlib import Path

HERE     = Path(__file__).parent
SIF_PATH = HERE / "160869.sif"
OUT_PATH = HERE.parent.parent / "Thesis" / "figures" / "fig_spectrum.png"

# ---- config ----
CH_MIDPLANE = 9              # canal ~z=0
T_A         = (2.0, 5.0)
T_B         = (6.0, 9.0)
T_OFFSET    = 0.12           # offset entre .sif y MDSplus

# ---- wavelength calibration (biaodingbochang.m) ----
_COEFS = np.polyfit([248, 376, 841], [33.73, 40.268, 33.73*2], 2)

def pix_to_wl_nm(p):
    return np.polyval(_COEFS, p) / 10.0

# ---- load .sif ----
data, info = sif_parser.np_open(str(SIF_PATH))
n_t, n_sp, n_wl = data.shape
cycle = info['CycleTime']
t_frame = np.arange(n_t) * cycle + T_OFFSET

# background subtraction (ultimo frame)
data_bg = data - data[-1][None, :, :]

# ejes
wl_nm  = pix_to_wl_nm(np.arange(n_wl))
mask_A = (t_frame >= T_A[0]) & (t_frame <= T_A[1])
mask_B = (t_frame >= T_B[0]) & (t_frame <= T_B[1])

# espectros midplane, promediados en cada fase
spec_A = data_bg[mask_A, CH_MIDPLANE, :].mean(axis=0)
spec_B = data_bg[mask_B, CH_MIDPLANE, :].mean(axis=0)

# ---- lineas conocidas: (nombre, nm, dx_text_nm, dy_text_offset_data) ----
# dx_text_nm: horizontal offset of the text label in wavelength units
# dy_text_offset_data: vertical offset above the local peak in COUNTS
LINE_ENTRIES = [
    ('C VI',      3.373, +0.80, 800),
    ('Fe XXIII', 13.29,  -0.55, 800),
    ('Li III',   13.50,  +0.55, 800),
]
WUTA = (4.80, 5.30)

# ---- style ----
plt.rcParams.update({
    'font.size':        13, 'axes.labelsize': 14, 'axes.titlesize': 14,
    'legend.fontsize':  12, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'axes.linewidth':   1.2,
})

fig, ax = plt.subplots(figsize=(12, 6))
fig.subplots_adjust(top=0.80)

# W-UTA band
ax.axvspan(*WUTA, alpha=0.18, color='#7b1fa2', zorder=0)

# spectra
ax.plot(wl_nm, spec_A, color='#c00000', lw=1.1, label='Phase A: ECRH ON  (2--5 s)')
ax.plot(wl_nm, spec_B, color='#1f4e79', lw=1.1, label='Phase B: ECRH OFF (6--9 s)')

ax.set_xlim(2, 16)
ax.set_ylim(0, max(spec_A.max(), spec_B.max()) * 1.05)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Intensity (counts, midplane chord)')

# ---- label placement: keep labels inside the plot area ----
y_max_plot = max(spec_A.max(), spec_B.max()) * 1.05   # matches set_ylim
y_cap      = y_max_plot * 0.95                        # max height for a label

# line markers with annotate: label sits above the peak with a small connector,
# capped so it never leaves the plot area
for name, wl, dx_nm, dy_counts in LINE_ENTRIES:
    mask_win = (wl_nm >= wl - 0.05) & (wl_nm <= wl + 0.05)
    peak_y = spec_A[mask_win].max() if mask_win.any() else 0.0
    text_y = min(peak_y + dy_counts, y_cap)
    text_x = wl + dx_nm
    ax.annotate(name,
                xy=(wl, peak_y),
                xytext=(text_x, text_y),
                fontsize=11, ha='center', va='bottom',
                color='#2e7d32', fontweight='bold',
                arrowprops=dict(arrowstyle='-', color='#2e7d32',
                                lw=0.9, alpha=0.7,
                                shrinkA=0, shrinkB=3))

# W-UTA label: connector arrow from the band peak to a label offset to the right
wuta_mask    = (wl_nm >= WUTA[0]) & (wl_nm <= WUTA[1])
wuta_peak    = spec_A[wuta_mask].max()
wuta_peak_wl = wl_nm[wuta_mask][np.argmax(spec_A[wuta_mask])]
wuta_text_y  = min(wuta_peak + 800, y_cap)
ax.annotate('W-UTA',
            xy=(wuta_peak_wl, wuta_peak),
            xytext=(wuta_peak_wl + 1.2, wuta_text_y),
            fontsize=11, ha='center', va='bottom',
            color='#7b1fa2', fontweight='bold',
            arrowprops=dict(arrowstyle='-', color='#7b1fa2',
                            lw=0.9, alpha=0.7,
                            shrinkA=0, shrinkB=3))

ax.legend(loc='upper right', framealpha=0.95, fontsize=12,
          bbox_to_anchor=(0.99, 0.98))
ax.grid(True, alpha=0.3, linestyle=':', which='both')
ax.tick_params(direction='in', length=5, which='both')

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_PATH, dpi=200, bbox_inches='tight')
print(f"Saved: {OUT_PATH}")

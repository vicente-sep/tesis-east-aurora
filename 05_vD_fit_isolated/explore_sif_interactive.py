"""
explore_sif_interactive.py
==========================
Visor interactivo del cubo de datos del espectrómetro XEUV_Long.

Abre una ventana con:
  - Panel superior izquierdo: espectro para el (tiempo, canal) seleccionado
  - Panel superior derecho: mapa 2D (t, z) integrando una banda espectral
  - Panel inferior izquierdo: evolución temporal en el canal seleccionado
  - Panel inferior derecho: perfil vertical (chord) en la ventana espectral
  - Dos sliders: tiempo (frame) y canal espacial (z)

Uso:
    python explore_sif_interactive.py

Requiere: numpy, matplotlib, sif_parser
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, RadioButtons
from pathlib import Path
import sif_parser

# ---- config ----
HERE = Path(__file__).parent
SIF_PATH = HERE / "160869.sif"

T_OFFSET = 0.12   # segundos: offset entre el reloj del .sif y MDSplus
Z_M = np.array([
    -0.2721, -0.2418, -0.2114, -0.1508, -0.1205,
    -0.0901, -0.0598, -0.0295,  0.0009,  0.0312,
     0.0616,  0.0918,  0.1222,  0.1525,  0.1828,
     0.2132,  0.2435,  0.2738,  0.3041,  0.3344,
]) * 100.0   # cm

# Presets de bandas espectrales conocidas (nm)
LINE_PRESETS = {
    'W-UTA':    (4.8, 5.3),
    'C VI':     (3.33, 3.42),
    'Fe XXIII': (13.24, 13.34),
    'Full':     (2.2, 16.3),
}

# ---- calibración wavelength ----
_COEFS = np.polyfit([248, 376, 841], [33.73, 40.268, 33.73*2], 2)
def pix_to_wl_nm(p):
    return np.polyval(_COEFS, np.asarray(p)) / 10.0


# ---- load .sif ----
print(f"Loading {SIF_PATH}...")
data, info = sif_parser.np_open(str(SIF_PATH))
n_t, n_sp, n_wl = data.shape
cycle = info['CycleTime']
print(f"  Shape: (t={n_t}, z={n_sp}, wl={n_wl}), cycle={cycle:.3f} s")

# background: último frame (post-plasma)
data_bg = data - data[-1][None, :, :]

# ejes
t_frame = np.arange(n_t) * cycle + T_OFFSET
wl_nm   = pix_to_wl_nm(np.arange(n_wl))
z_cm    = Z_M


# ---- helpers ----
def band_indices(band_nm):
    lo, hi = band_nm
    return np.where((wl_nm >= lo) & (wl_nm <= hi))[0]


def band_integrated_2d(band_nm):
    """Devuelve mapa 2D (t, z) integrando la banda espectral seleccionada."""
    idx = band_indices(band_nm)
    if len(idx) == 0:
        return np.zeros((n_t, n_sp))
    return data_bg[:, :, idx].sum(axis=2)


# ---- figure ----
plt.rcParams.update({
    'font.size': 11,
    'axes.linewidth': 1.2,
})

fig = plt.figure(figsize=(15, 9))
gs = gridspec.GridSpec(3, 3, height_ratios=[1, 1, 0.08],
                       hspace=0.35, wspace=0.30,
                       left=0.07, right=0.98, top=0.94, bottom=0.05)

ax_spec  = fig.add_subplot(gs[0, 0:2])   # espectro (arriba izq, ancho)
ax_map   = fig.add_subplot(gs[0, 2])     # mapa 2D t-z (arriba der)
ax_time  = fig.add_subplot(gs[1, 0:2])   # evolución temporal (abajo izq, ancho)
ax_prof  = fig.add_subplot(gs[1, 2])     # perfil vertical (abajo der)

# sliders al fondo
ax_t     = fig.add_subplot(gs[2, 0])
ax_ch    = fig.add_subplot(gs[2, 1])
ax_radio = fig.add_subplot(gs[2, 2])

fig.suptitle(f'XEUV\\_Long interactive viewer  |  {SIF_PATH.name}',
             fontsize=13, fontweight='bold', y=0.99)


# ---- estado inicial ----
current_band = 'W-UTA'
current_t = n_t // 2
current_ch = 8   # midplane approx

def get_spectrum(t_idx, ch_idx):
    return data_bg[t_idx, ch_idx, :]


# ---- initial plots ----

# Panel 1: spectrum
line_spec, = ax_spec.plot(wl_nm, get_spectrum(current_t, current_ch),
                          color='#c00000', lw=1.2)
band_lo, band_hi = LINE_PRESETS[current_band]
band_span = ax_spec.axvspan(band_lo, band_hi, alpha=0.20, color='#7b1fa2',
                            label=current_band)
ax_spec.set_xlim(wl_nm.min(), wl_nm.max())
ax_spec.set_xlabel('Wavelength (nm)')
ax_spec.set_ylabel('Intensity (counts, bg-subtracted)')
ax_spec.set_title(
    f'Spectrum  |  t = {t_frame[current_t]:.2f} s  |  ch = {current_ch}  '
    f'(z = {z_cm[current_ch]:+.1f} cm)'
)
ax_spec.grid(alpha=0.3, ls=':')
ax_spec.legend(loc='upper right')

# Panel 2: mapa 2D t-z de la banda seleccionada
band_map = band_integrated_2d(LINE_PRESETS[current_band])
im = ax_map.imshow(band_map.T,
                   aspect='auto', origin='lower',
                   extent=[t_frame[0], t_frame[-1],
                           z_cm.min(), z_cm.max()],
                   cmap='viridis',
                   interpolation='nearest')
line_t_map = ax_map.axvline(t_frame[current_t], color='w', lw=1.2, alpha=0.8)
line_z_map = ax_map.axhline(z_cm[current_ch], color='w', lw=1.2, alpha=0.8)
ax_map.set_xlabel('Time (s)')
ax_map.set_ylabel('z (cm)')
ax_map.set_title(f'{current_band} integrated  (t vs z)')
cbar = plt.colorbar(im, ax=ax_map, pad=0.02)
cbar.set_label('counts')

# Panel 3: evolución temporal en el canal seleccionado, banda integrada
temporal = band_map[:, current_ch]
line_time, = ax_time.plot(t_frame, temporal, 'o-', color='#1f4e79', ms=3, lw=1.2)
marker_t, = ax_time.plot([t_frame[current_t]], [temporal[current_t]],
                          'o', color='#c00000', ms=10, zorder=5)
ax_time.set_xlabel('Time (s)')
ax_time.set_ylabel(f'{current_band} intensity (counts)')
ax_time.set_title(f'Time evolution at ch = {current_ch} (z = {z_cm[current_ch]:+.1f} cm)')
ax_time.grid(alpha=0.3, ls=':')

# Panel 4: perfil vertical (todos los canales, tiempo seleccionado)
profile = band_map[current_t, :]
line_prof, = ax_prof.plot(profile, z_cm, 'o-', color='#2e7d32', ms=6, lw=1.5)
marker_ch, = ax_prof.plot([profile[current_ch]], [z_cm[current_ch]],
                          'o', color='#c00000', ms=10, zorder=5)
ax_prof.axhline(0.0, color='gray', ls=':', lw=0.8)
ax_prof.set_xlabel(f'{current_band} intensity (counts)')
ax_prof.set_ylabel('z (cm)')
ax_prof.set_title(f'Vertical profile at t = {t_frame[current_t]:.2f} s')
ax_prof.grid(alpha=0.3, ls=':')


# ---- sliders ----
slider_t = Slider(ax_t, 'Frame',
                  0, n_t - 1, valinit=current_t, valstep=1,
                  color='#c00000')
slider_t.label.set_size(10)

slider_ch = Slider(ax_ch, 'Chord',
                   0, n_sp - 1, valinit=current_ch, valstep=1,
                   color='#1f4e79')
slider_ch.label.set_size(10)

# Radio buttons para seleccionar la banda espectral
radio = RadioButtons(ax_radio, list(LINE_PRESETS.keys()),
                     active=list(LINE_PRESETS.keys()).index(current_band))


def redraw_spectrum(t_idx, ch_idx):
    line_spec.set_ydata(get_spectrum(t_idx, ch_idx))
    ax_spec.set_ylim(auto=True)
    ax_spec.relim()
    ax_spec.autoscale_view(scalex=False, scaley=True)
    ax_spec.set_title(
        f'Spectrum  |  t = {t_frame[t_idx]:.2f} s  |  ch = {ch_idx}  '
        f'(z = {z_cm[ch_idx]:+.1f} cm)'
    )


def redraw_map_lines(t_idx, ch_idx):
    line_t_map.set_xdata([t_frame[t_idx], t_frame[t_idx]])
    line_z_map.set_ydata([z_cm[ch_idx], z_cm[ch_idx]])


def redraw_temporal(band_name, ch_idx):
    band_map_local = band_integrated_2d(LINE_PRESETS[band_name])
    temporal_local = band_map_local[:, ch_idx]
    line_time.set_ydata(temporal_local)
    marker_t.set_ydata([temporal_local[int(slider_t.val)]])
    ax_time.set_ylim(auto=True)
    ax_time.relim()
    ax_time.autoscale_view(scalex=False, scaley=True)
    ax_time.set_ylabel(f'{band_name} intensity (counts)')
    ax_time.set_title(
        f'Time evolution at ch = {ch_idx} (z = {z_cm[ch_idx]:+.1f} cm)'
    )
    return band_map_local


def redraw_profile(band_map_local, t_idx):
    profile_local = band_map_local[t_idx, :]
    line_prof.set_xdata(profile_local)
    marker_ch.set_xdata([profile_local[int(slider_ch.val)]])
    ax_prof.set_xlim(auto=True)
    ax_prof.relim()
    ax_prof.autoscale_view(scalex=True, scaley=False)
    ax_prof.set_title(f'Vertical profile at t = {t_frame[t_idx]:.2f} s')


def redraw_map(band_name):
    band_map_local = band_integrated_2d(LINE_PRESETS[band_name])
    im.set_data(band_map_local.T)
    im.set_clim(vmin=band_map_local.min(), vmax=band_map_local.max())
    ax_map.set_title(f'{band_name} integrated  (t vs z)')
    return band_map_local


def redraw_band_span(band_name):
    global band_span
    band_span.remove()
    lo, hi = LINE_PRESETS[band_name]
    band_span = ax_spec.axvspan(lo, hi, alpha=0.20, color='#7b1fa2',
                                label=band_name)
    ax_spec.legend(loc='upper right')


# ---- callbacks ----
def on_slider(val):
    t_idx = int(slider_t.val)
    ch_idx = int(slider_ch.val)
    band_name = current_band

    redraw_spectrum(t_idx, ch_idx)
    redraw_map_lines(t_idx, ch_idx)

    band_map_local = band_integrated_2d(LINE_PRESETS[band_name])
    temporal_local = band_map_local[:, ch_idx]
    line_time.set_ydata(temporal_local)
    marker_t.set_xdata([t_frame[t_idx]])
    marker_t.set_ydata([temporal_local[t_idx]])
    ax_time.set_ylabel(f'{band_name} intensity (counts)')
    ax_time.set_title(
        f'Time evolution at ch = {ch_idx} (z = {z_cm[ch_idx]:+.1f} cm)'
    )
    ax_time.relim()
    ax_time.autoscale_view(scalex=False, scaley=True)

    profile_local = band_map_local[t_idx, :]
    line_prof.set_xdata(profile_local)
    marker_ch.set_xdata([profile_local[ch_idx]])
    marker_ch.set_ydata([z_cm[ch_idx]])
    ax_prof.set_title(f'Vertical profile at t = {t_frame[t_idx]:.2f} s')
    ax_prof.relim()
    ax_prof.autoscale_view(scalex=True, scaley=False)

    fig.canvas.draw_idle()


def on_radio(label):
    global current_band
    current_band = label
    redraw_band_span(label)
    band_map_local = redraw_map(label)

    t_idx = int(slider_t.val)
    ch_idx = int(slider_ch.val)

    temporal_local = band_map_local[:, ch_idx]
    line_time.set_ydata(temporal_local)
    marker_t.set_ydata([temporal_local[t_idx]])
    ax_time.set_ylabel(f'{label} intensity (counts)')
    ax_time.set_title(
        f'Time evolution at ch = {ch_idx} (z = {z_cm[ch_idx]:+.1f} cm)'
    )
    ax_time.relim()
    ax_time.autoscale_view(scalex=False, scaley=True)

    profile_local = band_map_local[t_idx, :]
    line_prof.set_xdata(profile_local)
    marker_ch.set_xdata([profile_local[ch_idx]])
    ax_prof.set_xlabel(f'{label} intensity (counts)')
    ax_prof.relim()
    ax_prof.autoscale_view(scalex=True, scaley=False)

    fig.canvas.draw_idle()


slider_t.on_changed(on_slider)
slider_ch.on_changed(on_slider)
radio.on_clicked(on_radio)

print("\nInteractive viewer open. Close the window to exit.\n")
plt.show()

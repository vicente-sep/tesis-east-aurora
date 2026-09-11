"""
plot_shots.py
=============
Grafica y guarda las señales temporales clave de un shot EAST.
Se conecta directamente al servidor MDSplus.

Escribe un numero de shot y presiona Enter para graficar.
Presiona "Save NPZ" para exportar toda la data a un archivo .npz.

Basado en el node reference verificado (Manual_Datos_EAST.md, Vicente Sepulveda).
"""

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox, Button
import mdsthin
warnings.filterwarnings("ignore")

# =============================================================================
#  CONFIGURACION
# =============================================================================

SERVER   = "202.127.204.12"
T_FLAT   = (2.0, 6.0)
SAVE_DIR = os.path.expanduser("~/Desktop")

# =============================================================================

def fetch(shot):
    conn = mdsthin.Connection(SERVER)

    def get(tree, node):
        try:
            conn.openTree(tree, shot)
            t = np.array(conn.get(f"dim_of({node})").data(), float).ravel()
            y = np.array(conn.get(node).data(), float).ravel()
            return t, y
        except Exception:
            return None, None

    # Ip (base temporal) - probar varios nodos comunes
    IP_CANDIDATES = [
        ("pcs_east", "\\pcrl01"),
        ("pcs_east", "\\pcrl02"),
        ("pcs_east", "\\pcrlm"),
        ("pcs_east", "\\ip"),
        ("east",     "\\ip"),
        ("east",     "\\ipm"),
        ("east_1",   "\\ip"),
    ]
    t_ip, ip = None, None
    for tree, node in IP_CANDIDATES:
        t_ip, ip = get(tree, node)
        if t_ip is not None:
            print(f"[shot {shot}] Ip: {tree}::{node}")
            break

    # Fallback: si Ip no aparece, usar cualquier otra senal para tener base temporal
    if t_ip is None:
        print(f"[shot {shot}] Ip no disponible, buscando base temporal alternativa...")
        FALLBACK_TIME = [
            ("efit_east", "\\WMHD"),
            ("pcs_east",  "\\dfsdev"),
            ("east",      "\\Dam1"),
            ("east",      "\\pxuv32"),
        ]
        for tree, node in FALLBACK_TIME:
            t_ip, _ = get(tree, node)
            if t_ip is not None:
                print(f"[shot {shot}] base temporal: {tree}::{node}")
                break
        if t_ip is None:
            raise ValueError(f"Shot #{shot}: no se pudo leer NINGUNA senal - shot inexistente?")
        ip = np.full_like(t_ip, np.nan)   # Ip queda como NaN
    else:
        ip = ip / 1000.0

    def interp(tree, node, default=np.nan):
        t, y = get(tree, node)
        if t is None:
            return np.full_like(t_ip, default)
        return np.interp(t_ip, t, y, left=default, right=default)

    # ECRH total (suma de 4 girotrones)
    ecrh = np.zeros_like(t_ip)
    for i in range(1, 5):
        t_ec, y_ec = get("ecrh_east", f"\\PECRH{i}I")
        if t_ec is not None:
            ecrh += np.interp(t_ip, t_ec, y_ec, left=0, right=0)

    # LHW neta: (PLHI1 - PLHR1) + PLHI2
    lhw_inj = interp("east", "\\PLHI1", default=0.0)
    lhw_ref = interp("east", "\\PLHR1", default=0.0)
    lhw2    = interp("east", "\\PLHI2", default=0.0)
    lhw     = np.nan_to_num(lhw_inj) - np.nan_to_num(lhw_ref) + np.nan_to_num(lhw2)

    # Plasma parameters
    ne   = interp("pcs_east",  "\\dfsdev")
    te   = interp("hrs_east",  "\\te0_hrs")
    ti   = interp("txcs_east", "\\ti0_txcs")
    wmhd = interp("efit_east", "\\WMHD")

    # Dalpha (dos nodos: midplane y divertor)
    dalpha_mid = interp("east", "\\Dam1")
    dalpha_div = interp("east", "\\Dal1")

    # Filterscopes y radiación
    w_filter = interp("east", "\\wu1")
    c_filter = interp("east", "\\ciiil1")
    zeff     = interp("east", "\\vbm1")
    prad     = interp("east", "\\pxuv32")

    # aminor para Greenwald
    a_i = interp("efit_east", "\\aminor")
    a_i = np.where(np.isnan(a_i) | (a_i < 0.1), 0.45, a_i)
    ip_ma = ip / 1000.0
    nG    = ip_ma / (np.pi * a_i**2) * 10
    ratio = ne / np.where(nG > 0, nG, np.nan)

    return dict(
        shot=shot, t=t_ip,
        ip=ip, ecrh=ecrh, lhw=lhw,
        ne=ne, nG=nG, ratio=ratio,
        te=te, ti=ti, wmhd=wmhd,
        dalpha_mid=dalpha_mid, dalpha_div=dalpha_div,
        w_filter=w_filter, c_filter=c_filter,
        zeff=zeff, prad=prad,
    )


# ── State ─────────────────────────────────────────────────────────────────────
last_data = {"d": None}

# ── Layout ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(7, 1, figsize=(11, 14), sharex=True)
fig.subplots_adjust(left=0.11, right=0.96, top=0.93, bottom=0.11, hspace=0.10)

ax_box  = fig.add_axes([0.20, 0.03, 0.14, 0.030])
textbox = TextBox(ax_box, "Shot: ", initial="", color="lightyellow", hovercolor="lightyellow")
ax_save = fig.add_axes([0.40, 0.03, 0.12, 0.032])
btn_save = Button(ax_save, "Save NPZ", color="lightgreen", hovercolor="green")
status  = fig.text(0.56, 0.038, "", fontsize=8, color="gray", va="center")
msg     = fig.text(0.56, 0.018, "", fontsize=8, color="red",  va="center")

YLABELS = [
    "Ip (kA)",
    "Heating (kW)",
    "ne / nG",
    "ne (10¹⁹ m⁻³)",
    "Te₀ / Ti₀ (keV)",
    "WMHD (kJ)",
    "Dα (a.u.)",
]

for ax, yl in zip(axes, YLABELS):
    ax.set_ylabel(yl, fontsize=10)
    ax.grid(True, alpha=0.25); ax.tick_params(labelsize=9)
axes[-1].set_xlabel("Tiempo (s)", fontsize=11)
fig.suptitle("Escribe un shot y presiona Enter", fontsize=12, color="gray")


def draw(shot_str):
    shot_str = shot_str.strip()
    if not shot_str:
        return
    try:
        shot = int(shot_str)
    except ValueError:
        msg.set_text(f'"{shot_str}" no es numero valido'); fig.canvas.draw_idle(); return

    for ax, yl in zip(axes, YLABELS):
        ax.cla(); ax.set_ylabel(yl, fontsize=10); ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("Tiempo (s)", fontsize=11)
    msg.set_text(""); status.set_text(f"Descargando shot #{shot}...")
    fig.suptitle(f"Cargando shot #{shot}...", fontsize=12, color="gray")
    fig.canvas.draw_idle(); plt.pause(0.05)

    try:
        d = fetch(shot)
    except Exception as e:
        msg.set_text(str(e)[:90]); status.set_text("")
        fig.suptitle(f"Shot #{shot} - error", fontsize=12, color="red")
        fig.canvas.draw_idle(); return

    last_data["d"] = d
    t = d["t"]

    axes[0].plot(t, d["ip"], color="steelblue", lw=1.4)

    axes[1].plot(t, d["ecrh"], color="firebrick",  lw=1.4, label="ECRH")
    axes[1].plot(t, d["lhw"],  color="royalblue",  lw=1.4, label="LHW net")
    axes[1].legend(fontsize=8, loc="upper right")

    axes[2].plot(t, d["ratio"], color="darkorange", lw=1.4)
    axes[2].axhline(1.0, color="k", lw=1.5, alpha=0.7)
    axes[2].axhline(1.1, color="r", lw=1.0, ls="--", alpha=0.5)

    axes[3].plot(t, d["ne"], color="seagreen", lw=1.4, label="ne")
    axes[3].plot(t, d["nG"], color="k", lw=1.0, ls="--", alpha=0.6, label="nG")
    axes[3].legend(fontsize=8, loc="upper right")

    axes[4].plot(t, d["te"], color="mediumpurple", lw=1.4, label="Te₀")
    if not np.all(np.isnan(d["ti"])):
        # Ti se guarda en eV (raw/1000 = keV según el manual)
        axes[4].plot(t, d["ti"] / 1000.0, color="orchid", lw=1.4, label="Ti₀")
    axes[4].legend(fontsize=8, loc="upper right")

    axes[5].plot(t, d["wmhd"] / 1000.0, color="darkslategray", lw=1.4)

    axes[6].plot(t, d["dalpha_mid"], color="crimson",   lw=0.9, label="midplane")
    axes[6].plot(t, d["dalpha_div"], color="darkred",   lw=0.9, alpha=0.6, label="divertor")
    axes[6].legend(fontsize=8, loc="upper right")

    for ax in axes:
        ax.axvspan(*T_FLAT, alpha=0.07, color="gold", zorder=0)
        ax.grid(True, alpha=0.25); ax.tick_params(labelsize=9)
    axes[-1].set_xlabel("Tiempo (s)", fontsize=11)

    mask = (t >= T_FLAT[0]) & (t <= T_FLAT[1])
    def fm(a): v = a[mask]; return np.nanmean(v) if np.any(~np.isnan(v)) else np.nan
    def fs(v, fmt): return f"{v:{fmt}}" if not np.isnan(v) else "---"

    fig.suptitle(
        f"Shot #{shot}   |   "
        f"Ip={fs(fm(d['ip']),'.0f')} kA   "
        f"ECRH={fs(fm(d['ecrh']),'.0f')} kW   "
        f"LHW={fs(fm(d['lhw']),'.0f')} kW   "
        f"ne={fs(fm(d['ne']),'.2f')}   "
        f"Te₀={fs(fm(d['te']),'.2f')} keV   "
        f"W={fs(fm(d['wmhd'])/1000,'.1f')} kJ",
        fontsize=9, fontweight="bold", color="black"
    )
    status.set_text("Listo. 'Save NPZ' para exportar.")
    fig.canvas.draw_idle()


def save_data(event):
    d = last_data["d"]
    if d is None:
        msg.set_text("Primero descarga un shot"); fig.canvas.draw_idle(); return
    os.makedirs(SAVE_DIR, exist_ok=True)
    path = os.path.join(SAVE_DIR, f"shot_{d['shot']}.npz")
    np.savez_compressed(path, **d)
    status.set_text(f"Guardado: {path}"); msg.set_text("")
    fig.canvas.draw_idle()
    print(f"Saved: {path}")


textbox.on_submit(draw)
btn_save.on_clicked(save_data)
plt.show()

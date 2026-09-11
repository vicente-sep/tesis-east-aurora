"""
fetch_ts_profiles.py
====================
Descarga el perfil radial de temperatura electrónica medido por ECE
(Electron Cyclotron Emission, heterodino 32-ch) para los 7 shots de
Liu et al. (Sci. Adv. 12, eadz3040, 2026).

Nodos confirmados en servidor 202.127.204.12:
  hrs_east / \\te_hrs   → shape (11000, 12)  Te(t, canal) en keV
                           t = -0.5 → 10.5 s, resolución ~1 ms
  hrs_east / \\r_hrs    → shape (12,)         R [m] de cada canal ECE
                           R = 1.925 → 2.288 m (región de núcleo)

Nota sobre unidades:
  \te0_hrs está confirmado en keV (sin dividir).
  \te_hrs también viene en keV (misma calibración).
  Valores negativos a t < 0 son ruido pre-plasma — se enmascaran.

Outputs en Output/ECE/:
  shot_<N>_ece_profile.csv    → Time(s), R_ch01..R_ch12 (m), Te_ch01..Te_ch12 (keV)
  shot_<N>_ece_map.png        → mapa Te(R, t)
  shot_<N>_ece_profiles.png   → perfiles Te(R) cada 0.5 s
  summary_ece_Te.png          → todos los shots @ t ~ 2 s
"""

import numpy as np
import pandas as pd
import mdsthin
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

# ─── Configuración ────────────────────────────────────────────────────────────

SERVER  = "202.127.204.12"
SHOTS   = [143064, 143069, 143073, 143074, 143075, 143077, 143079]
T_START = 0.0    # descartar pre-plasma (ruido antes del breakdown)
T_END   = 10.0

OUTPUT_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ECE"))


# ─── Descarga de un shot ──────────────────────────────────────────────────────

def fetch_shot_ece(shot):
    print(f"[#{shot}] Conectando...")
    try:
        conn = mdsthin.Connection(SERVER)
    except Exception as e:
        print(f"[#{shot}] ERROR conexión: {e}")
        return shot, None, None

    # ── Perfil Te(t, canal) ───────────────────────────────────────────────────
    try:
        conn.openTree("hrs_east", shot)
        raw  = np.array(conn.get("\\te_hrs").data(),  dtype=float)  # (11000, 12)
        time = np.array(conn.get("dim_of(\\te_hrs)").data(), dtype=float).ravel()
        r    = np.array(conn.get("\\r_hrs").data(),   dtype=float).ravel()
    except Exception as e:
        print(f"[#{shot}] ERROR descargando te_hrs: {e}")
        return shot, None, None

    # ── Normalizar shape a (n_t, n_ch) ───────────────────────────────────────
    if raw.ndim == 1:
        raw = raw.reshape(-1, 1)
    # Si viene transpuesto (n_ch, n_t):
    if raw.shape[0] == len(r) and raw.shape[1] == len(time):
        raw = raw.T   # → (n_t, n_ch)

    # ── Filtrar ventana temporal ──────────────────────────────────────────────
    mask = (time >= T_START) & (time <= T_END)
    time = time[mask]
    raw  = raw[mask, :]

    # ── Limpiar valores no físicos ────────────────────────────────────────────
    # Valores negativos = ruido/calibración; Te > 30 keV = artefacto
    raw = np.where((raw < 0) | (raw > 30), np.nan, raw)

    n_ok = np.sum(~np.isnan(raw))
    print(f"[#{shot}] ✓  te_hrs  shape={raw.shape}  R=[{r.min():.3f},{r.max():.3f}] m  "
          f"Te=[{np.nanmin(raw):.2f},{np.nanmax(raw):.2f}] keV  datos_ok={n_ok}")

    return shot, time, raw, r


# ─── Guardado y visualización ─────────────────────────────────────────────────

def save_csv(shot, time, raw, r):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    n_ch = raw.shape[1]
    df = pd.DataFrame({"Time(s)": time})
    for i in range(n_ch):
        df[f"R_ch{i+1:02d}(m)"]   = r[i] if i < len(r) else np.nan
    for i in range(n_ch):
        df[f"Te_ch{i+1:02d}(keV)"] = raw[:, i]
    path = os.path.join(OUTPUT_DIR, f"shot_{shot}_ece_profile.csv")
    df.to_csv(path, index=False)
    return path


def plot_ece_map(shot, time, raw, r):
    """Mapa de color Te(R, t)."""
    fig, ax = plt.subplots(figsize=(11, 5))
    vmax = min(np.nanpercentile(raw, 98), 15)
    vmin = 0
    pcm = ax.pcolormesh(time, r, raw.T,
                        shading="auto", cmap="inferno",
                        norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label("Te (keV)", fontsize=11)
    ax.set_xlabel("Tiempo (s)", fontsize=11)
    ax.set_ylabel("R (m)", fontsize=11)
    ax.set_title(f"Shot #{shot} — ECE Te(R, t)  [hrs_east/\\te_hrs, 12 canales]", fontsize=12)
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"shot_{shot}_ece_map.png")
    plt.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_ece_profiles(shot, time, raw, r):
    """Perfiles radiales en varios instantes de tiempo."""
    t_slices = np.arange(1.0, min(time[-1], 9.0), 0.5)
    colors = plt.cm.viridis(np.linspace(0, 1, len(t_slices)))

    fig, ax = plt.subplots(figsize=(8, 5))
    for color, t_ref in zip(colors, t_slices):
        idx = np.argmin(np.abs(time - t_ref))
        profile = raw[idx, :]
        valid = ~np.isnan(profile)
        if valid.sum() >= 2:
            ax.plot(r[valid], profile[valid], "o-", color=color,
                    label=f"t={t_ref:.1f}s", ms=4, lw=1.5)

    ax.set_xlabel("R (m)", fontsize=11)
    ax.set_ylabel("Te (keV)", fontsize=11)
    ax.set_title(f"Shot #{shot} — Perfiles ECE Te(R)", fontsize=12)
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"shot_{shot}_ece_profiles.png")
    plt.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_summary(all_data):
    """Perfiles Te(R) de todos los shots superpuestos a t ~ 2 s."""
    T_REP   = 2.0
    colors  = plt.cm.tab10(np.linspace(0, 0.9, len(SHOTS)))
    fig, ax = plt.subplots(figsize=(9, 5))

    for color, shot in zip(colors, SHOTS):
        if shot not in all_data or all_data[shot] is None:
            continue
        time, raw, r = all_data[shot]
        idx = np.argmin(np.abs(time - T_REP))
        profile = raw[idx, :]
        valid = ~np.isnan(profile)
        if valid.sum() >= 2:
            ax.plot(r[valid], profile[valid], "o-", color=color,
                    label=f"#{shot}", ms=5, lw=2)

    ax.set_xlabel("R (m)", fontsize=12)
    ax.set_ylabel("Te (keV)", fontsize=12)
    ax.set_title(f"Perfil ECE Te(R) @ t ≈ {T_REP} s — Liu et al. 2026", fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

    # Marcar eje magnético de EAST (~R = 1.85 m) con línea de referencia
    ax.axvline(x=1.85, color="gray", linestyle="--", lw=1, alpha=0.6, label="Eje magnético")

    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, "summary_ece_Te.png")
    plt.savefig(path, dpi=150)
    plt.close(fig)
    print(f"\n  Resumen → {path}")
    return path


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    print(f"\n{'='*60}")
    print(f"  ECE Profile Fetch — Liu et al. 2026")
    print(f"  Nodo: hrs_east/\\te_hrs  (12 canales, ~1 ms)")
    print(f"  Shots   : {SHOTS}")
    print(f"  Output  : {OUTPUT_DIR}")
    print(f"{'='*60}\n")

    all_data  = {}
    all_saved = {}
    failed    = []

    with ThreadPoolExecutor(max_workers=len(SHOTS)) as pool:
        futures = {pool.submit(fetch_shot_ece, s): s for s in SHOTS}
        for future in as_completed(futures):
            shot = futures[future]
            try:
                result = future.result()
                sid = result[0]
                if len(result) == 4 and result[1] is not None:
                    _, time, raw, r = result
                    all_data[sid] = (time, raw, r)

                    saved = []
                    saved.append(save_csv(sid, time, raw, r))
                    saved.append(plot_ece_map(sid, time, raw, r))
                    saved.append(plot_ece_profiles(sid, time, raw, r))
                    all_saved[sid] = saved
                    print(f"[#{sid}] Guardado ✓")
                else:
                    failed.append(sid)
            except Exception as exc:
                print(f"[#{shot}] ERROR: {exc}")
                failed.append(shot)

    if all_data:
        plot_summary(all_data)

    print(f"\n{'='*60}")
    print(f"  OK      ({len(all_data)}/{len(SHOTS)}): {sorted(all_data.keys())}")
    if failed:
        print(f"  Fallidos({len(failed)}/{len(SHOTS)}): {sorted(failed)}")
    print(f"  Archivos en: {OUTPUT_DIR}")
    for shot, paths in sorted(all_saved.items()):
        for p in paths:
            print(f"    #{shot}: {os.path.basename(p)}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

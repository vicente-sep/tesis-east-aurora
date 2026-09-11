"""
Análisis completo de los 7 shots Liu et al. 2026
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os, warnings
warnings.filterwarnings("ignore")

SHOTS   = [143064, 143069, 143073, 143074, 143075, 143077, 143079]
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "CSV_Shots"))
OUTDIR  = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "Output"))
os.makedirs(OUTDIR, exist_ok=True)

COLORS  = plt.cm.tab10(np.linspace(0, 0.9, len(SHOTS)))
CMAP    = dict(zip(SHOTS, COLORS))

def load(shot):
    return pd.read_csv(f"{DATADIR}/shot_{shot}_data.csv")

def col(df, name, default=np.nan):
    return df[name].values if name in df.columns else np.full(len(df), default)

def ft(df, t0=2.0, t1=6.0):
    return df[(df["Time (s)"] >= t0) & (df["Time (s)"] <= t1)]

dfs = {s: load(s) for s in SHOTS}

# ─── nG: Límite de Greenwald ───────────────────────────────────────────────
# nG [10¹⁹ m⁻³] = Ip[MA] / (π × a[m]²) × 10
def greenwald(df):
    ip = col(df, "Ip (kA)") / 1000          # → MA
    a  = col(df, "aminor (m)")
    a  = np.where(a < 0.1, 0.45, a)         # fallback si NaN
    nG = ip / (np.pi * a**2)                 # × 10²⁰ m⁻³ → × 10 en 10¹⁹ m⁻³
    return nG * 10

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 1 — Series temporales principales (6 paneles)
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(6, 1, figsize=(13, 16), sharex=True)
fig.suptitle("Liu et al. 2026 — Series temporales principales", fontsize=14, y=1.001)

panels = [
    ("Ip (kA)",         "Ip (kA)",        None),
    ("ne (1e19 m-3)",   "ne (10¹⁹ m⁻³)", None),
    ("Te0 (keV)",       "Te₀ (keV)",      None),
    ("PEC (kW)",        "ECRH (kW)",      None),
    ("PLHW_net (kW)",   "LHW net (kW)",   None),
    ("Da1 (a.u.)",      "Dα (a.u.)",      None),
]
fallbacks = ["Da1 (a.u.)", "Da (a.u.)"]

for ax, (col_name, ylabel, _) in zip(axes, panels):
    for shot in SHOTS:
        df = dfs[shot]
        t  = df["Time (s)"].values
        if col_name == "Da1 (a.u.)":
            y = col(df, "Da1 (a.u.)") if "Da1 (a.u.)" in df.columns else col(df, "Da (a.u.)")
        else:
            y = col(df, col_name)
        ax.plot(t, y, lw=1.2, color=CMAP[shot], label=f"#{shot}", alpha=0.85)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.grid(True, alpha=0.25)
    if col_name == "ne (1e19 m-3)":
        # Greenwald medio
        nG_vals = [greenwald(dfs[s]) for s in SHOTS]
        nG_mean = np.nanmean([np.nanmean(n[(dfs[s]["Time (s)"]>=2)&(dfs[s]["Time (s)"]<=6)])
                              for n, s in zip(nG_vals, SHOTS)])
        ax.axhline(nG_mean, color="k", ls="--", lw=1.2, alpha=0.6, label=f"nG≈{nG_mean:.1f}")

axes[-1].set_xlabel("Tiempo (s)", fontsize=11)
axes[0].legend(fontsize=8, loc="upper right", ncol=4)
axes[1].legend(fontsize=8, loc="upper right", ncol=4)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig1_timeseries.png", dpi=150, bbox_inches="tight")
plt.close(); print("fig1 OK")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 2 — EFIT: equilibrio MHD
# ══════════════════════════════════════════════════════════════════════════════
efit_panels = [
    ("q0",        "q₀",          "Factor de seguridad central"),
    ("q95",       "q₉₅",         "Factor de seguridad edge"),
    ("WMHD (kJ)", "Wmhd (kJ)",   "Energía almacenada"),
    ("betaN",     "βN",           "Beta normalizada"),
    ("kappa",     "κ",            "Elongación"),
    ("li",        "li",           "Inductancia interna"),
]
fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
axes = axes.flatten()
fig.suptitle("Liu et al. 2026 — Equilibrio MHD (EFIT)", fontsize=13)
for ax, (cn, yl, title) in zip(axes, efit_panels):
    for shot in SHOTS:
        df = dfs[shot]; t = df["Time (s)"].values; y = col(df, cn)
        ax.plot(t, y, lw=1.5, color=CMAP[shot], label=f"#{shot}", alpha=0.85)
    ax.set_ylabel(yl, fontsize=10); ax.set_title(title, fontsize=9)
    ax.grid(True, alpha=0.25)
for ax in axes[-2:]: ax.set_xlabel("Tiempo (s)", fontsize=10)
axes[0].legend(fontsize=8, ncol=4)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig2_efit.png", dpi=150, bbox_inches="tight")
plt.close(); print("fig2 OK")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 3 — Cinética: Te0, Ti0, Te/Ti, Vt0
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(4, 1, figsize=(13, 12), sharex=True)
fig.suptitle("Liu et al. 2026 — Parámetros cinéticos", fontsize=13)
for shot in SHOTS:
    df = dfs[shot]; t = df["Time (s)"].values
    te = col(df, "Te0 (keV)"); ti = col(df, "Ti0 (keV)")
    vt = col(df, "Vt0 (km/s)")
    ratio = np.where((ti > 0) & (te > 0), te/ti, np.nan)
    axes[0].plot(t, te, lw=1.3, color=CMAP[shot], label=f"#{shot}")
    axes[1].plot(t, ti, lw=1.3, color=CMAP[shot], marker="o", ms=3, ls="none" if np.sum(~np.isnan(ti))<50 else "-")
    axes[2].plot(t, ratio, lw=1.3, color=CMAP[shot], marker="o", ms=3, ls="none" if np.sum(~np.isnan(ratio))<50 else "-")
    axes[3].plot(t, vt, lw=1.3, color=CMAP[shot], marker="o", ms=3, ls="none" if np.sum(~np.isnan(vt))<50 else "-")

for ax, yl in zip(axes, ["Te₀ (keV)", "Ti₀ (keV)", "Te/Ti", "Vt₀ (km/s)"]):
    ax.set_ylabel(yl, fontsize=10); ax.grid(True, alpha=0.25)
axes[2].axhline(1, color="k", ls="--", lw=1, alpha=0.5, label="Te=Ti")
axes[3].axhline(0, color="k", ls="--", lw=1, alpha=0.4)
axes[-1].set_xlabel("Tiempo (s)", fontsize=11)
axes[0].legend(fontsize=8, ncol=4, loc="upper right")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig3_kinetics.png", dpi=150, bbox_inches="tight")
plt.close(); print("fig3 OK")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 4 — ne/nG: régimen density-free
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
fig.suptitle("Liu et al. 2026 — Régimen density-free (ne / nG)", fontsize=13)
for shot in SHOTS:
    df = dfs[shot]; t = df["Time (s)"].values
    ne = col(df, "ne (1e19 m-3)")
    nG = greenwald(df)
    ratio = ne / np.where(nG > 0, nG, np.nan)
    axes[0].plot(t, ne,    lw=1.3, color=CMAP[shot], label=f"#{shot}")
    axes[0].plot(t, nG,    lw=1.0, color=CMAP[shot], ls="--", alpha=0.5)
    axes[1].plot(t, ratio, lw=1.5, color=CMAP[shot], label=f"#{shot}")

axes[0].set_ylabel("ne, nG (10¹⁹ m⁻³)\n(sólido=ne, guión=nG)", fontsize=10)
axes[1].set_ylabel("ne / nG", fontsize=10)
axes[1].axhline(1.0, color="k", ls="-",  lw=2,   alpha=0.8, label="nG límite")
axes[1].axhline(1.1, color="r", ls="--", lw=1.2, alpha=0.6, label="110% nG")
axes[-1].set_xlabel("Tiempo (s)", fontsize=11)
for ax in axes: ax.grid(True, alpha=0.25)
axes[0].legend(fontsize=8, ncol=4)
axes[1].legend(fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig4_nG_ratio.png", dpi=150, bbox_inches="tight")
plt.close(); print("fig4 OK")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 5 — Impurezas y radiación
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(5, 1, figsize=(13, 14), sharex=True)
fig.suptitle("Liu et al. 2026 — Impurezas y radiación", fontsize=13)
imp_panels = [
    ("Da1 (a.u.)",      "Da (a.u.)",       "Dα divertor inf."),
    ("W_filter (a.u.)", "W (a.u.)",        "Filterscope W"),
    ("C_filter (a.u.)", "CIII (a.u.)",     "Filterscope CIII"),
    ("Prad_axuv (a.u.)","Prad AXUV (a.u.)","Prad AXUV"),
    ("Zeff_brem (a.u.)","Zeff (a.u.)",     "VBM Bremsstrahlung"),
]
for ax, (cn, yl, title) in zip(axes, imp_panels):
    for shot in SHOTS:
        df = dfs[shot]; t = df["Time (s)"].values
        if cn == "Da1 (a.u.)":
            y = col(df, "Da1 (a.u.)") if "Da1 (a.u.)" in df.columns else col(df, "Da (a.u.)")
        else:
            y = col(df, cn)
        ax.plot(t, y, lw=1.2, color=CMAP[shot], label=f"#{shot}", alpha=0.85)
    ax.set_ylabel(yl, fontsize=9); ax.set_title(title, fontsize=9, pad=2)
    ax.grid(True, alpha=0.25)
axes[-1].set_xlabel("Tiempo (s)", fontsize=11)
axes[0].legend(fontsize=8, ncol=4, loc="upper right")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig5_impurities.png", dpi=150, bbox_inches="tight")
plt.close(); print("fig5 OK")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 6 — Divertor: triple probe + Tt calculada
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(4, 1, figsize=(13, 11), sharex=True)
fig.suptitle("Liu et al. 2026 — Divertor (triple probe)", fontsize=13)
for shot in SHOTS:
    df = dfs[shot]; t = df["Time (s)"].values
    vp = col(df, "Tt_Vprobe (V)")
    vf = col(df, "Tt_Vfloat (V)")
    ii = col(df, "Tt_Isat (A)")
    tt = np.clip((vp - vf) / np.log(2), 0, 150)
    axes[0].plot(t, vp, lw=1.0, color=CMAP[shot], label=f"#{shot}", alpha=0.8)
    axes[1].plot(t, vf, lw=1.0, color=CMAP[shot], alpha=0.8)
    axes[2].plot(t, ii, lw=1.0, color=CMAP[shot], alpha=0.8)
    axes[3].plot(t, tt, lw=1.0, color=CMAP[shot], alpha=0.8)
for ax, yl in zip(axes, ["Vprobe (V)", "Vfloat (V)", "Isat (A)", "Tt = (Vp-Vf)/ln2 (eV)"]):
    ax.set_ylabel(yl, fontsize=9); ax.grid(True, alpha=0.25)
axes[-1].set_xlabel("Tiempo (s)", fontsize=11)
axes[0].legend(fontsize=8, ncol=4)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig6_divertor.png", dpi=150, bbox_inches="tight")
plt.close(); print("fig6 OK")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 7 — Calentamiento total
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
fig.suptitle("Liu et al. 2026 — Potencias de calentamiento", fontsize=13)
for shot in SHOTS:
    df = dfs[shot]; t = df["Time (s)"].values
    pec  = col(df, "PEC (kW)",      0)
    plhw = col(df, "PLHW_net (kW)", 0)
    ptot = np.where(np.isnan(pec), 0, pec) + np.where(np.isnan(plhw), 0, plhw)
    axes[0].plot(t, pec,  lw=1.3, color=CMAP[shot], label=f"#{shot}")
    axes[1].plot(t, plhw, lw=1.3, color=CMAP[shot])
    axes[2].plot(t, ptot, lw=1.3, color=CMAP[shot])
for ax, yl in zip(axes, ["ECRH (kW)", "LHW net (kW)", "Ptot = ECRH+LHW (kW)"]):
    ax.set_ylabel(yl, fontsize=10); ax.grid(True, alpha=0.25)
axes[-1].set_xlabel("Tiempo (s)", fontsize=11)
axes[0].legend(fontsize=8, ncol=4)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig7_heating.png", dpi=150, bbox_inches="tight")
plt.close(); print("fig7 OK")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 8 — Tabla resumen flat-top (t=2–6 s)
# ══════════════════════════════════════════════════════════════════════════════
rows = []
for shot in SHOTS:
    df = dfs[shot]
    f  = ft(df)
    nG = greenwald(f)
    row = {
        "Shot":          shot,
        "Ip (kA)":       f["Ip (kA)"].mean() if "Ip (kA)" in f else np.nan,
        "ne (10¹⁹)":     f["ne (1e19 m-3)"].mean() if "ne (1e19 m-3)" in f else np.nan,
        "nG (10¹⁹)":     np.nanmean(nG),
        "ne/nG":          np.nanmean(f["ne (1e19 m-3)"].values / np.where(nG>0,nG,np.nan)) if "ne (1e19 m-3)" in f else np.nan,
        "Te0 (keV)":     f["Te0 (keV)"].mean() if "Te0 (keV)" in f else np.nan,
        "Ti0 (keV)":     f["Ti0 (keV)"].mean() if "Ti0 (keV)" in f else np.nan,
        "Te/Ti":         (f["Te0 (keV)"].mean()/f["Ti0 (keV)"].mean()) if ("Te0 (keV)" in f and "Ti0 (keV)" in f and f["Ti0 (keV)"].mean()>0) else np.nan,
        "Vt0 (km/s)":    f["Vt0 (km/s)"].mean() if "Vt0 (km/s)" in f else np.nan,
        "q0":            f["q0"].mean() if "q0" in f else np.nan,
        "q95":           f["q95"].mean() if "q95" in f else np.nan,
        "WMHD (kJ)":     f["WMHD (kJ)"].mean() if "WMHD (kJ)" in f else np.nan,
        "βN":            f["betaN"].mean() if "betaN" in f else np.nan,
        "ECRH (kW)":     f["PEC (kW)"].mean() if "PEC (kW)" in f else np.nan,
        "LHW (kW)":      f["PLHW_net (kW)"].mean() if "PLHW_net (kW)" in f else np.nan,
    }
    rows.append(row)

summary = pd.DataFrame(rows).set_index("Shot").round(3)
print("\n" + summary.to_string())
summary.to_csv(f"{OUTDIR}/summary_flatop.csv")

# Plot tabla
fig, ax = plt.subplots(figsize=(16, 3.5))
ax.axis("off")
tbl = ax.table(
    cellText=summary.reset_index().values.round(2),
    colLabels=["Shot"] + list(summary.columns),
    loc="center", cellLoc="center"
)
tbl.auto_set_font_size(False); tbl.set_fontsize(9)
tbl.scale(1, 1.6)
for j in range(len(summary.columns)+1):
    tbl[0, j].set_facecolor("#2c7bb6"); tbl[0, j].set_text_props(color="white", fontweight="bold")
for i, shot in enumerate(SHOTS):
    c = CMAP[shot]
    tbl[i+1, 0].set_facecolor((*c[:3], 0.3))
plt.title("Flat-top estadístico (t = 2–6 s) — Liu et al. 2026", fontsize=12, pad=10)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig8_summary_table.png", dpi=150, bbox_inches="tight")
plt.close(); print("fig8 OK")

print("\n✓ Todas las figuras generadas en", OUTDIR)

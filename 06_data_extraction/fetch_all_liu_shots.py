"""
fetch_all_liu_shots.py
======================
Descarga todas las señales verificadas en el servidor 202.127.204.12
para los 7 disparos de Liu et al. (Sci. Adv. 12, eadz3040, 2026):
    143064, 143069, 143073, 143074, 143075, 143077, 143079

Inventario completo verificado con test_all_nodes.py (2026-04-19):
─────────────────────────────────────────────────────────────────────
PARÁMETROS GLOBALES (pcs_east)
  ✓ Ip           \\pcrl01          → /1000 = kA
  ✓ ne central   \\dfsdev          → 10¹⁹ m⁻³  (cuerda HCN central)
  ✓ ne2          \\dfsdev2         → 10¹⁹ m⁻³  (cuerda HCN 2)
  ✓ ne3          \\dfsdev3         → 10¹⁹ m⁻³  (cuerda HCN 3)
  ✓ Vloop        \\pcvloop         → V

EQUILIBRIO MHD (efit_east, ~106 puntos temporales)
  ✓ q0           \\q0              → factor de seguridad central
  ✓ q95          \\q95             → factor de seguridad edge
  ✓ li           \\li              → inductancia interna normalizada
  ✓ WMHD         \\WMHD            → J  (energía almacenada MHD)
  ✓ Rmaxis       \\rmaxis          → m  (posición eje magnético)
  ✓ aminor       \\aminor          → m  (radio menor)
  ✓ kappa        \\kappa           → elongación
  ✓ betaN        \\betan           → beta normalizada
  ✓ betap        \\betap           → beta poloidal
  ✓ drsep        \\drsep           → m  (separación d_rsep)

TEMPERATURA ELECTRÓNICA — ECE (hrs_east)
  ✓ Te0          \\te0_hrs         → keV (raw ya en keV, NO dividir)

TEMPERATURA IÓNICA Y ROTACIÓN — XCS (txcs_east, ~20 pts t=[2.6,6.6]s)
  ✓ Ti0          \\ti0_txcs        → eV raw /1000 = keV
  ✓ Ti0err       \\ti0_txcserr     → eV raw /1000 = keV
  ✓ Vt0          \\vt0_txcs        → km/s (negativo = counter-current)
  ✓ Vt0err       \\vt0_txcserr     → km/s

RADIACIÓN E IMPUREZAS (east)
  ✓ Dα           \\Dal1/2/3, \\Dam1 → a.u. (divertor inf. y plano medio)
  ✓ Prad_axuv    \\pxuv32          → a.u. (canal central AXUV)
  ✓ Prad_bolo    \\bolo1           → a.u. (bolómetro)
  ✓ Zeff_brem    \\vbm1            → a.u. (bremsstrahlung)
  ✓ W_filter     \\wu1             → a.u. (Filterscope W)
  ✓ C_filter     \\ciiil1          → a.u. (Filterscope CIII)

DIVERTOR (east)
  ✓ Tt_Vprobe    \\lovp01          → V   (sonda triple Vprobe)
  ✓ Tt_Isat      \\lois01          → A   (sonda triple Isat)
  ✓ Tt_Vfloat    \\lovf01          → V   (sonda triple Vfloat)
  ✓ ne_horiz     \\point_n3        → a.u. (interferómetro POINT cuerda 3)

GAS E INYECCIÓN (east)
  ✓ Vgas_max     \\vodpev1-3 + \\vhdpev1 + \\voupev1-3 + \\vcdpev1  → V (máximo)
  ✓ SMBI2        \\smbi2           → V   (Super Molecular Beam Injection)

CALENTAMIENTO (ecrh_east, east)
  ✓ PEC          \\PECRH1I-4I      → kW  (ECRH total, suma 4 girotrones)
  ✓ PLHI1        \\PLHI1           → kW  (LHW 2.45GHz inyectada)
  ✓ PLHR1        \\PLHR1           → kW  (LHW 2.45GHz reflejada)
  ✓ PLHW_net     PLHI1 - PLHR1    → kW  (LHW potencia neta ≥ 0)

NO DISPONIBLES REMOTAMENTE:
  ✗ CXRS (Ti/Vt perfiles)  — E-FOPENR (árbol no existe en servidor público)
  ✗ XEUV01/02 (EUV espectrómetro) — NNF (datos en servidor interno ASIPP)
  ✗ MHD/Mirnovs — NNF
  ✗ CXRS — NNF
─────────────────────────────────────────────────────────────────────
"""

import numpy as np
import pandas as pd
import mdsthin
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

# ─── Configuración ────────────────────────────────────────────────────────────

SERVER  = "202.127.204.12"
T_START = 0.0
T_END   = 10.0
DT      = 0.001   # resolución 1 ms → 10 000 puntos por shot

SHOTS = [143064, 143069, 143073, 143074, 143075, 143077, 143079]

OUTPUT_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "CSV_Shots"))


# ─── Señales escalares 1D confirmadas ────────────────────────────────────────
# (col_csv, árbol, nodo, escala)

SCALAR_SIGNALS = [
    # ── Parámetros globales ───────────────────────────────────────────────
    ("Ip (kA)",          "pcs_east",   "\\pcrl01",      1/1000),
    ("Vloop (V)",        "pcs_east",   "\\pcvloop",     1     ),
    ("ne (1e19 m-3)",    "pcs_east",   "\\dfsdev",      1     ),
    ("ne2 (1e19 m-3)",   "pcs_east",   "\\dfsdev2",     1     ),
    ("ne3 (1e19 m-3)",   "pcs_east",   "\\dfsdev3",     1     ),
    # ── ECE temperatura electrónica ───────────────────────────────────────
    # te0_hrs: raw en keV — NO dividir por 1000
    ("Te0 (keV)",        "hrs_east",   "\\te0_hrs",     1     ),
    # ── XCS temperatura iónica y rotación (~20 pts, t=[2.6,6.6]s) ────────
    ("Ti0 (keV)",        "txcs_east",  "\\ti0_txcs",    1/1000),
    ("Ti0err (keV)",     "txcs_east",  "\\ti0_txcserr", 1/1000),
    ("Vt0 (km/s)",       "txcs_east",  "\\vt0_txcs",    1     ),
    ("Vt0err (km/s)",    "txcs_east",  "\\vt0_txcserr", 1     ),
    # ── Equilibrio MHD EFIT (~106 pts) ───────────────────────────────────
    ("q0",               "efit_east",  "\\q0",          1     ),
    ("q95",              "efit_east",  "\\q95",         1     ),
    ("li",               "efit_east",  "\\li",          1     ),
    ("WMHD (kJ)",        "efit_east",  "\\WMHD",        1/1000),
    ("Rmaxis (m)",       "efit_east",  "\\rmaxis",      1     ),
    ("aminor (m)",       "efit_east",  "\\aminor",      1     ),
    ("kappa",            "efit_east",  "\\kappa",       1     ),
    ("betaN",            "efit_east",  "\\betan",       1     ),
    ("betap",            "efit_east",  "\\betap",       1     ),
    ("drsep (m)",        "efit_east",  "\\drsep",       1     ),
    # ── Dα emisión ───────────────────────────────────────────────────────
    ("Da1 (a.u.)",       "east",       "\\Dal1",        1     ),
    ("Da2 (a.u.)",       "east",       "\\Dal2",        1     ),
    ("Da3 (a.u.)",       "east",       "\\Dal3",        1     ),
    ("Dam1 (a.u.)",      "east",       "\\Dam1",        1     ),
    # ── SMBI ─────────────────────────────────────────────────────────────
    ("SMBI2 (V)",        "east",       "\\smbi2",       1     ),
]

# ─── Señales con fallback (primer nodo que responda) ─────────────────────────
# (col_csv, escala, [(árbol, nodo), ...])

FALLBACK_SIGNALS = [
    ("Prad_bolo (a.u.)", 1, [("east", "\\bolo1"),  ("east", "\\bolo01")]),
    ("Prad_axuv (a.u.)", 1, [("east", "\\pxuv32"), ("east", "\\pxuv1")]),
    ("Zeff_brem (a.u.)", 1, [("east", "\\vbm1"),   ("east", "\\VBM1")]),
    ("W_filter (a.u.)",  1, [("east", "\\wu1"),    ("east", "\\WU1")]),
    ("C_filter (a.u.)",  1, [("east", "\\ciiil1"), ("east", "\\CIIIL1")]),
    ("Tt_Vprobe (V)",    1, [("east", "\\lovp01"), ("east", "\\lovp1")]),
    ("Tt_Isat (A)",      1, [("east", "\\lois01"), ("east", "\\lois1")]),
    ("Tt_Vfloat (V)",    1, [("east", "\\lovf01"), ("east", "\\lovf1")]),
    ("ne_horiz (a.u.)",  1, [("east", "\\point_n3"), ("pcs_east", "\\point_n3")]),
]

# ─── Válvulas de gas (máximo sobre todos los canales activos) ─────────────────
VGAS_NODES = [
    ("east", "\\vodpev1"), ("east", "\\vodpev2"), ("east", "\\vodpev3"),
    ("east", "\\vhdpev1"),
    ("east", "\\voupev1"), ("east", "\\voupev2"), ("east", "\\voupev3"),
    ("east", "\\vcdpev1"),
]

# ─── ECRH: suma de 4 girotrones [kW raw] ─────────────────────────────────────
ECRH_NODES = [
    ("ecrh_east", "\\PECRH1I"), ("ecrh_east", "\\PECRH2I"),
    ("ecrh_east", "\\PECRH3I"), ("ecrh_east", "\\PECRH4I"),
]

# ─── LHW 2.45 GHz: potencia neta = inyectada − reflejada ─────────────────────
LHW_NODES = {
    "inj": ("east", "\\PLHI1"),
    "ref": ("east", "\\PLHR1"),
}


# ─── Funciones ────────────────────────────────────────────────────────────────

def fetch_node(conn, tree, node, shot):
    try:
        conn.openTree(tree, shot)
        t = np.array(conn.get(f"dim_of({node})").data(), dtype=float)
        y = np.array(conn.get(node).data(),              dtype=float)
        if t.ndim > 1: t = t.ravel()
        if y.ndim > 1: y = y.ravel()
        if len(t) == 0 or len(y) == 0:
            return None, None
        return t, y
    except Exception as e:
        print(f"  [#{shot}] WARN {tree}/{node}: {str(e)[:60]}")
        return None, None


def interp_to_common(t_raw, y_raw, t_common, fill=np.nan):
    if t_raw is None or len(t_raw) == 0:
        return np.full(len(t_common), fill)
    idx = np.argsort(t_raw)
    return np.interp(t_common, t_raw[idx], y_raw[idx], left=fill, right=fill)


def fetch_shot(shot):
    print(f"[#{shot}] Conectando...")
    try:
        conn = mdsthin.Connection(SERVER)
    except Exception as e:
        print(f"[#{shot}] ERROR conexión: {e}")
        return shot, None

    t_common = np.arange(T_START, T_END, DT)
    df = pd.DataFrame({"Time (s)": t_common})

    # ── Escalares directos ───────────────────────────────────────────────────
    for col, tree, node, scale in SCALAR_SIGNALS:
        t, y = fetch_node(conn, tree, node, shot)
        if t is not None:
            df[col] = interp_to_common(t, y * scale, t_common)
            print(f"  [#{shot}] {col} ✓")
        else:
            df[col] = np.nan

    # ── Fallback ─────────────────────────────────────────────────────────────
    for col, scale, node_list in FALLBACK_SIGNALS:
        found = False
        for tree, node in node_list:
            t, y = fetch_node(conn, tree, node, shot)
            if t is not None:
                df[col] = interp_to_common(t, y * scale, t_common)
                print(f"  [#{shot}] {col} ← {tree}/{node} ✓")
                found = True
                break
        if not found:
            df[col] = np.nan

    # ── Vgas: máximo sobre válvulas activas ──────────────────────────────────
    vgas_accum = None
    for tree, node in VGAS_NODES:
        t, y = fetch_node(conn, tree, node, shot)
        if t is not None:
            yi = interp_to_common(t, y, t_common, fill=0.0)
            vgas_accum = yi if vgas_accum is None else np.maximum(vgas_accum, yi)
    df["Vgas_max (V)"] = vgas_accum if vgas_accum is not None else np.nan

    # ── ECRH total ───────────────────────────────────────────────────────────
    ecrh_total = np.zeros(len(t_common))
    has_ecrh = False
    for tree, node in ECRH_NODES:
        t, y = fetch_node(conn, tree, node, shot)
        if t is not None:
            has_ecrh = True
            ecrh_total += interp_to_common(t, y, t_common, fill=0.0)
    df["PEC (kW)"] = ecrh_total if has_ecrh else np.nan
    if has_ecrh:
        print(f"  [#{shot}] PEC ✓")

    # ── LHW 2.45 GHz ────────────────────────────────────────────────────────
    t_inj, y_inj = fetch_node(conn, *LHW_NODES["inj"], shot)
    t_ref, y_ref = fetch_node(conn, *LHW_NODES["ref"], shot)
    if t_inj is not None:
        plhi = interp_to_common(t_inj, y_inj, t_common, fill=0.0)
        plhr = interp_to_common(t_ref, y_ref, t_common, fill=0.0) if t_ref is not None else 0.0
        df["PLHI1 (kW)"]    = plhi
        df["PLHR1 (kW)"]    = plhr
        df["PLHW_net (kW)"] = np.maximum(plhi - plhr, 0.0)
        print(f"  [#{shot}] LHW ✓  max_net={df['PLHW_net (kW)'].max():.0f} kW")
    else:
        df["PLHI1 (kW)"] = df["PLHR1 (kW)"] = df["PLHW_net (kW)"] = np.nan

    return shot, df


def save_shot(shot, df):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    path = os.path.join(OUTPUT_DIR, f"shot_{shot}_data.csv")
    df.to_csv(path, index=False)
    cols_ok = df.notna().any().sum() - 1   # excluye Time
    print(f"[#{shot}] ✓  {cols_ok} señales  →  {path}")
    return path


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    n_signals = len(SCALAR_SIGNALS) + len(FALLBACK_SIGNALS) + 3  # +Vgas,PEC,LHW
    print(f"\n{'='*60}")
    print(f"  Liu et al. 2026 — Descarga completa EAST")
    print(f"  Shots    : {SHOTS}")
    print(f"  Servidor : {SERVER}")
    print(f"  Señales  : {n_signals} columnas por shot")
    print(f"  Resolución: {DT*1000:.0f} ms  ({int((T_END-T_START)/DT)} puntos)")
    print(f"{'='*60}\n")

    results, failed = {}, []

    with ThreadPoolExecutor(max_workers=len(SHOTS)) as pool:
        futures = {pool.submit(fetch_shot, s): s for s in SHOTS}
        for future in as_completed(futures):
            shot = futures[future]
            try:
                shot_id, df = future.result()
                if df is not None:
                    results[shot_id] = save_shot(shot_id, df)
                else:
                    failed.append(shot_id)
            except Exception as exc:
                print(f"[#{shot}] ERROR: {exc}")
                failed.append(shot)

    print(f"\n{'='*60}")
    print(f"  OK      ({len(results)}/{len(SHOTS)}): {sorted(results.keys())}")
    if failed:
        print(f"  Fallidos({len(failed)}/{len(SHOTS)}): {sorted(failed)}")
    print(f"  Archivos en: {OUTPUT_DIR}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

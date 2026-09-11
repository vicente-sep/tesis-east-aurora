"""
test_ts_nodes.py  (v3)
======================
Resultados anteriores:
  ✓ east/\\TS3 y \\TS5 existen  →  shape (4250000,), 250 kHz crudo  ← señal raw, no perfil
  ✗ hrs_east/\\HRS01H–HRS32H   →  NNF

Ahora buscamos:
  1. Perfil ECE procesado      →  \\te_hrs (pista del MATLAB Mydata_V24.m)
  2. Nodos alternativos ECE    →  distintas convenciones de nombre
  3. Perfiles TS procesados    →  \\TS3_Te, \\TS3_NE, \\te_ts3, etc.
  4. Eje radial del ECE/TS     →  \\R, \\r_hrs, \\freq_hrs, etc.
"""

import mdsthin
import numpy as np

SERVER = "202.127.204.12"
SHOT   = 143073

print(f"\n{'='*65}")
print(f"  Test ECE + TS perfiles procesados — shot #{SHOT}")
print(f"{'='*65}\n")

try:
    conn = mdsthin.Connection(SERVER)
    print("  Conexión OK\n")
except Exception as e:
    print(f"  ERROR: {e}")
    exit(1)

nodes = [
    # ── ECE perfil completo (pista del MATLAB: \te_hrs) ──────────────────
    ("hrs_east", "\\te_hrs",       "ECE perfil Te (Mydata_V24 hint) → /1000 = keV"),
    ("hrs_east", "\\te_hrs_pro",   "ECE perfil Te procesado"),
    ("hrs_east", "\\te0_hrs",      "ECE Te0 central [ya confirmado]"),
    # Eje espacial / frecuencias ECE
    ("hrs_east", "\\freq_hrs",     "Frecuencias ECE (mapean a R)"),
    ("hrs_east", "\\r_hrs",        "Posiciones R de canales ECE"),
    ("hrs_east", "\\rho_hrs",      "rho normalizado canales ECE"),
    # Canales individuales con otra convención
    ("hrs_east", "\\HRS1H",        "ECE canal 1 (sin cero)"),
    ("hrs_east", "\\HRS01",        "ECE canal 01 (sin H)"),
    ("hrs_east", "\\te_01",        "ECE te_01"),
    ("hrs_east", "\\te01",         "ECE te01"),
    ("hrs_east", "\\ch01",         "ECE ch01"),
    # ── TS perfiles procesados ────────────────────────────────────────────
    ("east",     "\\TS3_TE",       "TS3 Te procesada"),
    ("east",     "\\TS3_NE",       "TS3 ne procesada"),
    ("east",     "\\TS5_TE",       "TS5 Te procesada"),
    ("east",     "\\TS5_NE",       "TS5 ne procesada"),
    ("east",     "\\te_ts3",       "Te TS core"),
    ("east",     "\\ne_ts3",       "ne TS core"),
    ("east",     "\\te_ts5",       "Te TS edge"),
    ("east",     "\\ne_ts5",       "ne TS edge"),
    # Eje radial de TS
    ("east",     "\\R_TS3",        "R posiciones TS3"),
    ("east",     "\\R_TS5",        "R posiciones TS5"),
    ("east",     "\\rho_ts3",      "rho TS3"),
    # ── ECE en árbol east (alternativo) ───────────────────────────────────
    ("east",     "\\te_hrs",       "east/te_hrs"),
    ("east",     "\\HRS01H",       "east/HRS01H"),
]

found = []
for tree, node, desc in nodes:
    try:
        conn.openTree(tree, SHOT)
        raw  = conn.get(node).data()
        time = conn.get(f"dim_of({node})").data()
        d = np.array(raw,  dtype=float)
        t = np.array(time, dtype=float).ravel()
        print(f"  ✓  {tree}/{node}")
        print(f"       [{desc}]")
        print(f"       shape = {d.shape}")
        if len(t) > 0:
            print(f"       t = {t[0]:.4f} → {t[-1]:.4f} s   n_t = {len(t)}")
        print(f"       min={d.min():.4g}  max={d.max():.4g}  mean={d.mean():.4g}\n")
        found.append((tree, node, desc, d.shape))
    except Exception as e:
        print(f"  ✗  {tree:12s}/{node:18s}  → {e}")

print(f"\n{'='*65}")
print(f"  RESULTADO: {len(found)}/{len(nodes)} nodos con datos")
print(f"{'='*65}")
for tree, node, desc, shape in found:
    print(f"  {tree}/{node}  shape={shape}")
print()

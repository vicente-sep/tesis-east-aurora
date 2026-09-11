"""
test_xcs_nodes.py
=================
Prueba rápida de todos los nodos XCS disponibles en el servidor EAST.
Corre este script para verificar qué nodos están accesibles antes de
hacer la descarga masiva.

Diagnóstico XCS en EAST:
  - 切向弯晶谱仪 (Tangential Bent Crystal Spectrometer)
  - Ventana C, extremo del tubo de vacío
  - Mide: Ti (temperatura iónica), Vt (rotación toroidal)
  - Árbol MDSplus: txcs_east
  - Señales crudas: XCS01, XCS02, XEUV01, XEUV02 (1 kS/s)
"""

import mdsthin
import numpy as np

SERVER = "202.127.204.12"
SHOT   = 143073   # Shot de referencia

print(f"\n{'='*65}")
print(f"  Test nodos XCS — shot #{SHOT} — servidor {SERVER}")
print(f"{'='*65}\n")

try:
    conn = mdsthin.Connection(SERVER)
except Exception as e:
    print(f"ERROR conexión: {e}")
    exit(1)

# ─── Nodos a probar ─────────────────────────────────────────────────────────
nodes = [
    # ── Árbol txcs_east (confirmado en Mydata_V24.m) ─────────────────────
    # Ti0 central: raw en eV → /1000 = keV
    ("txcs_east", "\\ti0_txcs",      "Ti0 central (eV raw → /1000 keV)"),
    ("txcs_east", "\\ti0_txcserr",   "Ti0 error"),
    # Velocidad de rotación toroidal
    ("txcs_east", "\\vt0_txcs",      "Vt0 rotación toroidal"),
    ("txcs_east", "\\vt0_txcserr",   "Vt0 error"),
    # Perfiles radiales (si están procesados)
    ("txcs_east", "\\ti_txcs",       "Ti perfil radial"),
    ("txcs_east", "\\vt_txcs",       "Vt perfil radial"),
    ("txcs_east", "\\ti_txcs_t",     "Ti perfil (variante _t)"),
    ("txcs_east", "\\vt_txcs_t",     "Vt perfil (variante _t)"),
    # Señales crudas (del documento oficial de diagnósticos)
    ("txcs_east", "\\xcs01",         "XCS raw canal 1"),
    ("txcs_east", "\\xcs02",         "XCS raw canal 2"),
    ("txcs_east", "\\xeuv01",        "XEUV canal 1"),
    ("txcs_east", "\\xeuv02",        "XEUV canal 2"),
    # Posiciones radiales de la medición
    ("txcs_east", "\\r_txcs",        "Posiciones R de la medición"),
    ("txcs_east", "\\rho_txcs",      "Posiciones rho normalizadas"),
    # ── Árbol east (alternativo) ──────────────────────────────────────────
    ("east",      "\\xcs01",         "east/XCS01"),
    ("east",      "\\xcs02",         "east/XCS02"),
    ("east",      "\\ti0_txcs",      "east/ti0_txcs"),
    ("east",      "\\vt0_txcs",      "east/vt0_txcs"),
]

found = []
for tree, node, desc in nodes:
    try:
        conn.openTree(tree, SHOT)
        data = conn.get(node).data()
        time = conn.get(f"dim_of({node})").data()
        y = np.array(data, dtype=float).ravel()
        t = np.array(time, dtype=float).ravel()
        status = f"✓  {tree}/{node}"
        print(f"  {status}")
        print(f"     [{desc}]")
        print(f"     t = {t[0]:.3f} → {t[-1]:.3f} s   n = {len(t)}")
        print(f"     min = {y.min():.4g}   max = {y.max():.4g}   mean = {y.mean():.4g}\n")
        found.append((tree, node, desc, t, y))
    except Exception as e:
        print(f"  ✗  {tree}/{node:25s}  → {e}")

print(f"\n{'='*65}")
print(f"  RESULTADO: {len(found)}/{len(nodes)} nodos con datos")
print(f"{'='*65}\n")

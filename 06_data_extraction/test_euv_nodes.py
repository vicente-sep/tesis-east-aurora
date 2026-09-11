"""
test_euv_nodes.py
=================
Explora los nodos MDSplus del espectrómetro EUV de EAST.

Según la literatura, EAST tiene 4 espectrómetros EUV:
  EUV_Short   →  5–45 Å   (W muy ionizado, core caliente)
  EUV_Long_a  →  40–180 Å (W43+-W45+, líneas principales de W)
  EUV_Long_b  →  270–480 Å (W5+ en 382/394 Å → flujo de W desde pared)
  EUV_Long_c  →  130–330 Å

Lo que buscamos:
  - Espectros resueltos en longitud de onda: intensidad(λ, t)
  - Eje de longitud de onda: λ en Å
  - Si hay resolución espacial (perfiles radiales integrados)

Cómo correr:
  1. Conéctate a EAST por VSCode remoto
  2. python test_euv_nodes.py
  3. Los nodos que aparezcan con ✓ son los que hay que descargar

Referencia:
  Xiao et al. 2016 (RSI) — EUV_Long en EAST
  Zhang et al. 2021 (NIM) — EUV_Short 10-130 Å en EAST
  arXiv:2410.02669        — W5+ a 382/394 Å en EAST
"""

import mdsthin
import numpy as np

SERVER = "202.127.204.12"
SHOT   = 143073   # Shot de referencia (alta ECRH, más señal de W)

print(f"\n{'='*65}")
print(f"  Test nodos EUV — shot #{SHOT} — servidor {SERVER}")
print(f"  Buscando espectrómetros EUV_Short / EUV_Long_a/b/c")
print(f"{'='*65}\n")

try:
    conn = mdsthin.Connection(SERVER)
    print("  Conexión OK\n")
except Exception as e:
    print(f"  ERROR conexión: {e}")
    exit(1)

nodes = [
    # ── euv_east confirmado en Mydata_V24.m ──────────────────────────────────
    # \wuta = W-UTA (45-50 Å), traza temporal de intensidad — NUNCA extraído
    # porque Tungsten_plot=0 estaba apagado en el MATLAB del grupo chino.
    ("euv_east",    "\\wuta",           "W-UTA intensidad (Mydata_V24 confirmado)"),
    # Explorar otras señales del árbol euv_east
    ("euv_east",    "\\wuta2",          "W-UTA canal 2"),
    ("euv_east",    "\\euv_short",      "EUV_Short espectro (5-45 Å)"),
    ("euv_east",    "\\euv_long_a",     "EUV_Long_a (40-180 Å) W43+-W45+"),
    ("euv_east",    "\\euv_long_b",     "EUV_Long_b (270-480 Å) W5+ 382/394 Å"),
    ("euv_east",    "\\euv_long_c",     "EUV_Long_c (130-330 Å)"),
    ("euv_east",    "\\euv_s",          "EUV_Short abreviado"),
    ("euv_east",    "\\euv_la",         "EUV_Long_a abreviado"),
    ("euv_east",    "\\euv_lb",         "EUV_Long_b abreviado"),
    ("euv_east",    "\\euv_lc",         "EUV_Long_c abreviado"),
    ("euv_east",    "\\spec_short",     "spec_short"),
    ("euv_east",    "\\spec_long",      "spec_long"),
    ("euv_east",    "\\spec1",          "spec1"),
    ("euv_east",    "\\spec2",          "spec2"),
    ("euv_east",    "\\wavelength",     "eje de λ"),
    ("euv_east",    "\\intensity",      "intensidad espectral"),
    ("euv_east",    "\\signal",         "señal genérica"),
    ("euv_east",    "\\w_filter",       "W filter en euv_east"),
    ("euv_east",    "\\c_filter",       "C filter en euv_east"),
    # Árbol 'east' con nombres alternativos
    ("east",        "\\wuta",           "east/wuta"),
    ("east",        "\\XEUV01",         "east/XEUV01"),
    ("east",        "\\XEUV02",         "east/XEUV02"),
    ("txcs_east",   "\\xeuv01",         "txcs_east/xeuv01"),
    ("txcs_east",   "\\xeuv02",         "txcs_east/xeuv02"),
    # ── Thomson Scattering ne(r) — bonus ─────────────────────────────────────
    ("east",        "\\ne_ts",          "TS ne (árbol east)"),
    ("east",        "\\te_ts",          "TS Te (árbol east)"),
    ("ts_east",     "\\ne",             "ts_east/ne"),
    ("ts_east",     "\\te",             "ts_east/te"),
    ("ts_east",     "\\r_ts",           "ts_east/r_ts (eje radial)"),
    ("ts_east",     "\\ne_core",        "ts_east/ne_core"),
    ("ts_east",     "\\te_core",        "ts_east/te_core"),
    ("east",        "\\TS3_TE",         "TS3 Te (árbol east)"),
    ("east",        "\\TS3_NE",         "TS3 ne (árbol east)"),
    ("east",        "\\TS5_TE",         "TS5 Te (árbol east)"),
    ("east",        "\\TS5_NE",         "TS5 ne (árbol east)"),
]

found   = []
missing = []

for tree, node, desc in nodes:
    try:
        conn.openTree(tree, SHOT)
        data = conn.get(node).data()
        dims = conn.get(f"dim_of({node})").data()
        d = np.array(data, dtype=float)
        t = np.array(dims, dtype=float).ravel()
        print(f"  ✓  {tree}/{node}")
        print(f"       [{desc}]")
        print(f"       shape = {d.shape}")
        if len(t) > 0:
            print(f"       dim   = {t[0]:.4g} → {t[-1]:.4g}   n = {len(t)}")
        print(f"       rango = {d.min():.4g} … {d.max():.4g}\n")
        found.append((tree, node, desc, d.shape))
    except Exception as e:
        err = str(e)[:60]
        print(f"  ✗  {tree:14s}/{node:22s}  → {err}")
        missing.append((tree, node))

# ─── Resumen ──────────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"  RESULTADO: {len(found)}/{len(nodes)} nodos encontrados")
print(f"{'='*65}")

if found:
    print("\n  NODOS CON DATOS:")
    for tree, node, desc, shape in found:
        print(f"    ✓  {tree}/{node}  shape={shape}")
    print()
    print("  → Copia los nodos encontrados y avisa para escribir")
    print("    el script de descarga completo (fetch_euv.py).")
else:
    print("\n  No se encontró ningún nodo EUV.")
    print("  Opciones:")
    print("    1. Preguntarle al grupo de EAST cuál es el árbol MDSplus")
    print("       del espectrómetro EUV (EUV_Short / EUV_Long_a/b/c).")
    print("    2. Buscar en el archivo Mydata_V24.m qué nodos usa el")
    print("       grupo chino para leer los datos EUV en MATLAB.")
    print("    3. Revisar si el espectrómetro EUV estaba activo durante")
    print("       los shots 143064–143079 (no siempre está habilitado).")
print()

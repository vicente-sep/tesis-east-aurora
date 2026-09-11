"""
test_all_nodes_ipp.py
=====================
Igual que test_all_nodes.py pero apunta al servidor alternativo mds.ipp.ac.cn.
Corre esto para ver si EUV/XCS están disponibles en este servidor.

Cambia SHOT al número que quieras probar.
"""

import mdsthin
import numpy as np
import sys
from datetime import datetime

SERVER = "mds.ipp.ac.cn"   # ← servidor alternativo ASIPP
SHOT   = 161075            # ← CAMBIA AQUÍ

ALL_NODES = [
    # Parámetros globales
    ("pcs_east", "\\pcrl01",    "Ip corriente de plasma [A]"),
    ("pcs_east", "\\ipm",       "Ip (nodo ipm — sintaxis pmds)"),
    ("pcs_east", "\\dfsdev",    "ne HCN cuerda 1"),
    ("pcs_east", "\\dfsdev2",   "ne HCN cuerda 2"),
    ("pcs_east", "\\pcvloop",   "Vloop [V]"),
    ("pcs_east", "\\bcentr",    "Bt [T]"),
    # EFIT
    ("efit_east", "\\q95",      "q95"),
    ("efit_east", "\\WMHD",     "Wmhd [J]"),
    ("efit_east", "\\betan",    "betaN"),
    # ECE
    ("hrs_east",  "\\te0_hrs",  "Te0 ECE [keV]"),
    ("hrs_east",  "\\te_hrs",   "Te perfil ECE"),
    # XCS / Ti
    ("txcs_east", "\\ti0_txcs",    "Ti0 XCS"),
    ("txcs_east", "\\vt0_txcs",    "Vt0 XCS"),
    ("txcs_east", "\\xcs01",       "XCS raw ch01"),
    ("txcs_east", "\\xcs02",       "XCS raw ch02"),
    ("txcs_east", "\\xeuv01",      "XEUV ch01 (txcs_east)"),
    ("txcs_east", "\\xeuv02",      "XEUV ch02 (txcs_east)"),
    # EUV
    ("east",      "\\XEUV",        "EUV total"),
    ("east",      "\\XEUV01",      "EUV ch01"),
    ("east",      "\\XEUV02",      "EUV ch02"),
    ("euv_east",  "\\xeuv01",      "EUV ch01 (euv_east)"),
    ("euv_east",  "\\wuta",        "W línea"),
    ("euv_east",  "\\cuta",        "C línea"),
    # CXRS
    ("cxrs_east", "\\ti0_cxrs",    "Ti0 CXRS"),
    ("cxrs_east", "\\vt0_cxrs",    "Vt0 CXRS"),
    ("cxrs_east", "\\int0_cxrs",   "Intensidad central CXRS"),
    # Filterscopes
    ("east",  "\\wu1",       "Filterscope W ch1"),
    ("east",  "\\wu7",       "Filterscope W ch7"),
    ("east",  "\\ciiil1",    "Filterscope CIII ch1"),
    # SXR
    ("east",  "\\SXR15V",    "SXR vertical ch15"),
    ("east",  "\\SXR23U",    "SXR horizontal ch23"),
    # Radiación
    ("prad_east", "\\pradtot",      "Prad total"),
    ("prad_east", "\\pradtot_axuv", "Prad AXUV"),
    # ECRH
    ("ecrh_east", "\\PECRH1I", "ECRH girotron 1"),
    ("ecrh_east", "\\PECRH2I", "ECRH girotron 2"),
    # Thomson
    ("east",  "\\te_ts3",    "Te Core TS procesado"),
    ("east",  "\\ne_ts3",    "ne Core TS procesado"),
]


print(f"\n{'='*65}")
print(f"  Servidor  : {SERVER}")
print(f"  Shot      : #{SHOT}")
print(f"  Timestamp : {datetime.now():%Y-%m-%d %H:%M}")
print(f"{'='*65}\n")

try:
    conn = mdsthin.Connection(SERVER)
    print("  Conexión OK\n")
except Exception as e:
    print(f"  ERROR conexión: {e}")
    sys.exit(1)

found, missing = [], []

for tree, node, desc in ALL_NODES:
    try:
        conn.openTree(tree, SHOT)
        d = np.array(conn.get(node).data(), dtype=float).ravel()
        t = np.array(conn.get(f"dim_of({node})").data(), dtype=float).ravel()
        n = len(t)
        tag = "PROC" if n < 50_000 else ("MED " if n < 500_000 else "RAW ")
        t_str = f"t=[{t[0]:.2f},{t[-1]:.2f}]s" if len(t) > 1 else "sin_t"
        print(f"  [{tag}] {tree}/{node}")
        print(f"         {desc}")
        print(f"         shape={d.shape}  {t_str}  min={d.min():.3g}  max={d.max():.3g}\n")
        found.append((tag, tree, node, desc))
    except Exception as e:
        err = str(e).split(",")[0].replace("%TREE-", "").strip()
        missing.append((tree, node, err))
        print(f"  [----] {tree}/{node:25s} → {err}")

print(f"\n{'='*65}")
print(f"  RESULTADO: {len(found)}/{len(ALL_NODES)} nodos con datos")
print(f"{'='*65}\n")
if found:
    print("  DISPONIBLES:")
    for tag, tree, node, desc in found:
        print(f"    [{tag}] {tree}/{node}  [{desc}]")

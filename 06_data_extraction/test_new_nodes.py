"""
test_new_nodes.py
=================
Verifica los nodos identificados en el documento "Mapeo de Diagnósticos
EAST (MDSplus)" que aún no habían sido testeados.

Prioridad alta para tesis de impurezas/transporte:
  - HKSO1H–HKSO32H  : Edge Thomson Scattering 32 canales (Te, ne borde)
  - SXR1V–SXR30V    : Soft X-ray array vertical (perfil impurezas)
  - SXR1U–SXR46U    : Soft X-ray array horizontal
  - XEUV, XEUV01/02 : EUV spectrometer (espectrómetro del usuario)
  - JJHG1, OUG1, ODFG1 : Gas puffing (Vgas real)
  - LIVP01–LIVP15, UOVP, LOVF : Sondas Langmuir divertor
"""

import mdsthin
import numpy as np

SERVER = "202.127.204.12"
SHOT   = 143073

print(f"\n{'='*65}")
print(f"  Test nodos nuevos — shot #{SHOT} — servidor {SERVER}")
print(f"{'='*65}\n")

conn = mdsthin.Connection(SERVER)
print("  Conexión OK\n")

nodes = [
    # ── Edge Thomson Scattering (HKSO) ────────────────────────────────────
    ("east",     "\\HKSO1H",    "Edge TS canal 1"),
    ("east",     "\\HKSO8H",    "Edge TS canal 8"),
    ("east",     "\\HKSO16H",   "Edge TS canal 16"),
    ("east",     "\\HKSO32H",   "Edge TS canal 32"),
    ("east",     "\\hkso1h",    "Edge TS canal 1 (minúsculas)"),
    ("ts_east",  "\\HKSO1H",    "ts_east/HKSO1H"),
    # ── Soft X-ray vertical (SXR1V–SXR30V) ───────────────────────────────
    ("east",     "\\SXR1V",     "SXR vertical canal 1"),
    ("east",     "\\SXR15V",    "SXR vertical canal 15"),
    ("east",     "\\SXR30V",    "SXR vertical canal 30"),
    ("east",     "\\sxr1v",     "SXR vertical 1 (minúsculas)"),
    # ── Soft X-ray horizontal (SXR1U–SXR46U) ─────────────────────────────
    ("east",     "\\SXR1U",     "SXR horizontal canal 1"),
    ("east",     "\\SXR23U",    "SXR horizontal canal 23"),
    ("east",     "\\SXR46U",    "SXR horizontal canal 46"),
    # ── EUV spectrometer ──────────────────────────────────────────────────
    ("east",     "\\XEUV",      "EUV espectrómetro (señal total)"),
    ("east",     "\\XEUV01",    "EUV canal 01"),
    ("east",     "\\XEUV02",    "EUV canal 02"),
    ("east",     "\\xeuv",      "EUV (minúsculas)"),
    ("east",     "\\xeuv01",    "EUV canal 01 (minúsculas)"),
    ("txcs_east","\\xeuv01",    "txcs_east/xeuv01"),
    ("txcs_east","\\XEUV01",    "txcs_east/XEUV01"),
    # ── Gas puffing ───────────────────────────────────────────────────────
    ("east",     "\\JJHG1",     "Gas puffing JJHG1"),
    ("east",     "\\OUG1",      "Gas puffing OUG1"),
    ("east",     "\\ODFG1",     "Gas puffing ODFG1"),
    ("east",     "\\jjhg1",     "Gas puffing (minúsculas)"),
    ("east",     "\\oug1",      "Gas puffing oug1"),
    # ── Sondas Langmuir divertor ──────────────────────────────────────────
    ("east",     "\\LIVP01",    "Sonda Langmuir LIVP01"),
    ("east",     "\\LIVP08",    "Sonda Langmuir LIVP08"),
    ("east",     "\\LIVP15",    "Sonda Langmuir LIVP15"),
    ("east",     "\\UOVP",      "Sonda Upper Outer Vprobe"),
    ("east",     "\\LOVF",      "Sonda Lower Outer Vfloat"),
    ("east",     "\\livp01",    "LIVP01 (minúsculas)"),
    ("east",     "\\uovp",      "UOVP (minúsculas)"),
    ("east",     "\\lovf",      "LOVF (minúsculas)"),
    # ── LHW (Lower Hybrid Heating) ────────────────────────────────────────
    ("east",     "\\PLHI1",     "LHW 2.45GHz potencia inyectada"),
    ("east",     "\\PLHR1",     "LHW 2.45GHz potencia reflejada"),
    ("east",     "\\PLHI2",     "LHW 4.6GHz potencia inyectada"),
    ("east",     "\\PLHR2",     "LHW 4.6GHz potencia reflejada"),
]

found = []
for tree, node, desc in nodes:
    try:
        conn.openTree(tree, SHOT)
        d = np.array(conn.get(node).data(), dtype=float)
        t = np.array(conn.get(f"dim_of({node})").data(), dtype=float).ravel()
        print(f"  ✓  {tree}/{node}")
        print(f"       [{desc}]  shape={d.shape}")
        if len(t) > 0:
            print(f"       t=[{t[0]:.3f},{t[-1]:.3f}] s  n_t={len(t)}")
        print(f"       min={d.min():.4g}  max={d.max():.4g}  mean={d.mean():.4g}\n")
        found.append((tree, node, desc, d.shape, len(t)))
    except Exception as e:
        print(f"  ✗  {tree:12s}/{node:14s}  → {e}")

print(f"\n{'='*65}")
print(f"  RESULTADO: {len(found)}/{len(nodes)} nodos con datos")
print(f"{'='*65}")
for tree, node, desc, shape, nt in found:
    processed = "PROCESADO" if nt < 100000 else "CRUDO"
    print(f"  [{processed}] {tree}/{node}  shape={shape}  [{desc}]")
print()

"""
test_ne_profile_nodes.py
========================
Verifica qué nodos de perfil de densidad están disponibles en el servidor.

Candidatos:
  POINT interferómetro  →  POINT_N1–POINT_N5  (5 cuerdas horizontales ne)
  HCN interferómetro    →  HCN1, HCN2, HCN3   (3 cuerdas verticales ne)
  efit_east             →  perfiles reconstruidos ne, Te, etc.
  analysis              →  perfiles procesados
"""

import mdsthin
import numpy as np

SERVER = "202.127.204.12"
SHOT   = 143073

print(f"\n{'='*65}")
print(f"  Test perfiles ne — shot #{SHOT} — servidor {SERVER}")
print(f"{'='*65}\n")

conn = mdsthin.Connection(SERVER)
print("  Conexión OK\n")

nodes = [
    # ── POINT: 5 cuerdas horizontales (pcs_east o east) ─────────────────
    ("pcs_east", "\\point_n1",   "POINT cuerda 1"),
    ("pcs_east", "\\point_n2",   "POINT cuerda 2"),
    ("pcs_east", "\\point_n3",   "POINT cuerda 3 (central)"),
    ("pcs_east", "\\point_n4",   "POINT cuerda 4"),
    ("pcs_east", "\\point_n5",   "POINT cuerda 5"),
    ("east",     "\\point_n1",   "east/point_n1"),
    ("east",     "\\point_n2",   "east/point_n2"),
    ("east",     "\\point_n3",   "east/point_n3"),
    ("east",     "\\point_n4",   "east/point_n4"),
    ("east",     "\\point_n5",   "east/point_n5"),
    # Ángulo de Faraday
    ("east",     "\\point_f1",   "POINT Faraday F1"),
    ("east",     "\\point_f3",   "POINT Faraday F3"),
    # ── HCN: 3 cuerdas verticales ────────────────────────────────────────
    ("pcs_east", "\\hcn1",       "HCN cuerda 1"),
    ("pcs_east", "\\hcn2",       "HCN cuerda 2 (central)"),
    ("pcs_east", "\\hcn3",       "HCN cuerda 3"),
    ("east",     "\\hcn1",       "east/hcn1"),
    ("east",     "\\hcn2",       "east/hcn2"),
    ("east",     "\\hcn3",       "east/hcn3"),
    # ── dfsdev: ya confirmado ────────────────────────────────────────────
    ("pcs_east", "\\dfsdev",     "dfsdev central [confirmado]"),
    ("pcs_east", "\\dfsdev2",    "dfsdev cuerda 2"),
    ("pcs_east", "\\dfsdev3",    "dfsdev cuerda 3"),
    # ── efit_east: perfiles reconstruidos ────────────────────────────────
    ("efit_east","\\ne_pro",     "efit ne perfil"),
    ("efit_east","\\te_pro",     "efit Te perfil"),
    ("efit_east","\\pressure",   "efit presión perfil"),
    ("efit_east","\\ne",         "efit ne"),
    # ── analysis: perfiles procesados ────────────────────────────────────
    ("analysis", "\\ne_pro",     "analysis ne perfil"),
    ("analysis", "\\ne_ts",      "analysis ne TS"),
    ("analysis", "\\ne0",        "analysis ne0"),
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
        found.append((tree, node, desc, d.shape))
    except Exception as e:
        print(f"  ✗  {tree:12s}/{node:16s}  → {e}")

print(f"\n{'='*65}")
print(f"  RESULTADO: {len(found)}/{len(nodes)} nodos disponibles")
for t, n, d, s in found:
    print(f"    {t}/{n}  shape={s}  [{d}]")
print()

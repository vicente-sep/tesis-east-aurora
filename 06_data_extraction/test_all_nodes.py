"""
test_all_nodes.py
=================
Escaneo exhaustivo de todos los diagnósticos EAST disponibles en el
servidor 202.127.204.12 para el shot de referencia #143073.

Cubre todos los árboles y nodos identificados en:
  - Catálogo oficial de diagnósticos EAST (documento interno ASIPP)
  - Documento "Mapeo de Diagnósticos EAST (MDSplus)"
  - MATLAB Mydata_V24.m (Mydata script del laboratorio)
  - Sesiones previas de verificación

Clasifica cada nodo como:
  [PROC]  < 50,000 muestras  → procesado/usable directamente
  [MED]   50k–500k muestras  → muestreo medio, probablemente usable
  [RAW]   > 500k muestras    → señal cruda de alta velocidad

Output: imprime tabla completa y guarda resumen en test_all_nodes_result.txt
"""

import mdsthin
import numpy as np
import sys
from datetime import datetime

SERVER = "202.127.204.12"
SHOT   = 143073

# ─── Definición de todos los nodos a probar ───────────────────────────────────
# (árbol, nodo, descripción)
ALL_NODES = [

    # ════════════════════════════════════════════════════════════════════════
    # 1. PARÁMETROS GLOBALES DE PLASMA (pcs_east)
    # ════════════════════════════════════════════════════════════════════════
    ("pcs_east", "\\pcrl01",    "Ip corriente de plasma [A]"),
    ("pcs_east", "\\dfsdev",    "ne central HCN cuerda 1 [10¹⁹ m⁻³]"),
    ("pcs_east", "\\dfsdev2",   "ne HCN cuerda 2"),
    ("pcs_east", "\\dfsdev3",   "ne HCN cuerda 3"),
    ("pcs_east", "\\pcvloop",   "Voltaje de lazo Vloop [V]"),
    ("pcs_east", "\\bcentr",    "Campo toroidal central Bt [T]"),
    ("pcs_east", "\\point_n1",  "POINT interferómetro cuerda 1"),
    ("pcs_east", "\\point_n2",  "POINT interferómetro cuerda 2"),
    ("pcs_east", "\\point_n3",  "POINT interferómetro cuerda 3"),
    ("pcs_east", "\\point_n4",  "POINT interferómetro cuerda 4"),
    ("pcs_east", "\\point_n5",  "POINT interferómetro cuerda 5"),

    # ════════════════════════════════════════════════════════════════════════
    # 2. EQUILIBRIO MHD (efit_east)
    # ════════════════════════════════════════════════════════════════════════
    ("efit_east", "\\q0",       "Factor de seguridad central q0"),
    ("efit_east", "\\q95",      "Factor de seguridad q95"),
    ("efit_east", "\\li",       "Inductancia interna normalizada li"),
    ("efit_east", "\\WMHD",     "Energía almacenada Wmhd [J]"),
    ("efit_east", "\\rmaxis",   "Posición eje magnético R [m]"),
    ("efit_east", "\\zmaxis",   "Posición eje magnético Z [m]"),
    ("efit_east", "\\aminor",   "Radio menor a [m]"),
    ("efit_east", "\\kappa",    "Elongación κ"),
    ("efit_east", "\\delta",    "Triangularidad δ"),
    ("efit_east", "\\betan",    "Beta normalizada βN"),
    ("efit_east", "\\betap",    "Beta poloidal βp"),
    ("efit_east", "\\fpol",     "Flujo poloidal en borde"),
    ("efit_east", "\\drsep",    "Separación d_rsep [m]"),
    ("efit_east", "\\ne_pro",   "Perfil ne reconstruido EFIT"),
    ("efit_east", "\\te_pro",   "Perfil Te reconstruido EFIT"),
    ("efit_east", "\\pressure", "Perfil de presión EFIT"),
    ("efit_east", "\\pprime",   "dP/dΨ perfil"),
    ("efit_east", "\\ffprime",  "FF' perfil"),

    # ════════════════════════════════════════════════════════════════════════
    # 3. TEMPERATURA ELECTRÓNICA — ECE (hrs_east)
    # ════════════════════════════════════════════════════════════════════════
    ("hrs_east",  "\\te0_hrs",  "Te0 central ECE [keV]"),
    ("hrs_east",  "\\te_hrs",   "Te perfil ECE 12 canales [keV]"),
    ("hrs_east",  "\\r_hrs",    "Posiciones R canales ECE [m]"),

    # ════════════════════════════════════════════════════════════════════════
    # 4. TEMPERATURA IÓNICA Y ROTACIÓN — XCS (txcs_east)
    # ════════════════════════════════════════════════════════════════════════
    ("txcs_east", "\\ti0_txcs",    "Ti0 central XCS [eV raw → /1000 keV]"),
    ("txcs_east", "\\ti0_txcserr", "Ti0 error XCS"),
    ("txcs_east", "\\vt0_txcs",    "Vt0 rotación toroidal XCS [km/s]"),
    ("txcs_east", "\\vt0_txcserr", "Vt0 error XCS"),
    ("txcs_east", "\\ti_txcs",     "Ti perfil radial XCS"),
    ("txcs_east", "\\vt_txcs",     "Vt perfil radial XCS"),

    # ════════════════════════════════════════════════════════════════════════
    # 5. CXRS — Charge Exchange Recombination Spectroscopy (cxrs_east)
    # ════════════════════════════════════════════════════════════════════════
    ("cxrs_east", "\\ti0_cxrs",    "Ti0 central CXRS"),
    ("cxrs_east", "\\vt0_cxrs",    "Vt0 central CXRS"),
    ("cxrs_east", "\\ti_cxrs",     "Ti perfil radial CXRS"),
    ("cxrs_east", "\\vt_cxrs",     "Vt perfil radial CXRS"),
    ("cxrs_east", "\\vp_cxrs",     "Vp velocidad poloidal CXRS"),
    ("cxrs_east", "\\ti0",         "Ti0 (nodo corto)"),
    ("cxrs_east", "\\vt0",         "Vt0 (nodo corto)"),
    ("cxrs_east", "\\ti",          "Ti perfil (nodo corto)"),
    ("cxrs_east", "\\vt",          "Vt perfil (nodo corto)"),
    ("cxrs_east", "\\r_cxrs",      "R posiciones CXRS"),
    ("cxrs_east", "\\rho_cxrs",    "rho CXRS"),
    ("cxrs_east", "\\int_cxrs",    "Intensidad línea CXRS"),
    ("cxrs_east", "\\int0_cxrs",   "Intensidad central CXRS"),

    # ════════════════════════════════════════════════════════════════════════
    # 6. ANÁLISIS / PERFILES PROCESADOS (analysis)
    # ════════════════════════════════════════════════════════════════════════
    ("analysis",  "\\eng",      "Energía almacenada análisis [J]"),
    ("analysis",  "\\ne_pro",   "ne perfil análisis"),
    ("analysis",  "\\te_pro",   "Te perfil análisis"),
    ("analysis",  "\\ti_pro",   "Ti perfil análisis"),
    ("analysis",  "\\vt_pro",   "Vt perfil análisis"),
    ("analysis",  "\\q_pro",    "q perfil análisis"),
    ("analysis",  "\\ne0",      "ne0 análisis"),
    ("analysis",  "\\te0",      "Te0 análisis"),

    # ════════════════════════════════════════════════════════════════════════
    # 7. RADIACIÓN E IMPUREZAS (east)
    # ════════════════════════════════════════════════════════════════════════
    ("east",  "\\Dal1",      "Dα divertor inferior canal 1"),
    ("east",  "\\Dal2",      "Dα canal 2"),
    ("east",  "\\Dal3",      "Dα canal 3"),
    ("east",  "\\Dam1",      "Dα plano medio canal 1"),
    # AXUV array
    ("east",  "\\pxuv1",     "AXUV canal 1"),
    ("east",  "\\pxuv16",    "AXUV canal 16"),
    ("east",  "\\pxuv32",    "AXUV canal 32 (central)"),
    ("east",  "\\pxuv48",    "AXUV canal 48"),
    ("east",  "\\pxuv64",    "AXUV canal 64"),
    ("east",  "\\cxuv1v",    "AXUV-C vertical canal 1"),
    ("east",  "\\cxuv15v",   "AXUV-C vertical canal 15"),
    ("east",  "\\cxuv30v",   "AXUV-C vertical canal 30"),
    # Bolómetros
    ("east",  "\\bolo1",     "Bolómetro canal 1"),
    ("east",  "\\bolo24",    "Bolómetro canal 24"),
    ("east",  "\\bolo48",    "Bolómetro canal 48"),
    # Prad total árbol prad_east
    ("prad_east", "\\pradtot",      "Prad total"),
    ("prad_east", "\\pradtot_axuv", "Prad AXUV total"),
    # VBM: Bremsstrahlung (Zeff)
    ("east",  "\\vbm1",      "VBM Zeff bremsstrahlung canal 1"),
    ("east",  "\\vbm7",      "VBM canal 7"),
    ("east",  "\\vbm13",     "VBM canal 13"),
    # Filterscopes
    ("east",  "\\wu1",       "Filterscope W divertor sup canal 1"),
    ("east",  "\\wu7",       "Filterscope W canal 7"),
    ("east",  "\\ciiil1",    "Filterscope CIII divertor inf canal 1"),
    ("east",  "\\ciiil7",    "Filterscope CIII canal 7"),
    # Soft X-ray
    ("east",  "\\SXR1V",     "SXR vertical canal 1"),
    ("east",  "\\SXR15V",    "SXR vertical canal 15"),
    ("east",  "\\SXR30V",    "SXR vertical canal 30"),
    ("east",  "\\SXR1U",     "SXR horizontal canal 1"),
    ("east",  "\\SXR23U",    "SXR horizontal canal 23"),
    ("east",  "\\SXR46U",    "SXR horizontal canal 46"),
    # EUV
    ("east",      "\\XEUV",    "EUV espectrómetro total"),
    ("east",      "\\XEUV01",  "EUV canal 01"),
    ("east",      "\\XEUV02",  "EUV canal 02"),
    ("txcs_east", "\\xeuv01",  "txcs_east/XEUV01"),
    ("euv_east",  "\\xeuv01",  "euv_east/XEUV01"),
    ("euv_east",  "\\wuta",    "euv_east/wuta W línea"),
    ("euv_east",  "\\cuta",    "euv_east/cuta C línea"),

    # ════════════════════════════════════════════════════════════════════════
    # 8. DIVERTOR Y PARED (east)
    # ════════════════════════════════════════════════════════════════════════
    ("east",  "\\lovp01",    "Triple probe lower outer Vprobe"),
    ("east",  "\\lois01",    "Triple probe Isat"),
    ("east",  "\\lovf01",    "Triple probe Vfloat"),
    ("east",  "\\LIVP01",    "Langmuir LIVP01"),
    ("east",  "\\LIVP02",    "Langmuir LIVP02"),
    ("east",  "\\LIVP03",    "Langmuir LIVP03"),
    ("east",  "\\LIVP04",    "Langmuir LIVP04"),
    ("east",  "\\LIVP05",    "Langmuir LIVP05"),
    ("east",  "\\LIVP06",    "Langmuir LIVP06"),
    ("east",  "\\LIVP07",    "Langmuir LIVP07"),
    ("east",  "\\LIVP08",    "Langmuir LIVP08"),
    ("east",  "\\LIVP09",    "Langmuir LIVP09"),
    ("east",  "\\LIVP10",    "Langmuir LIVP10"),
    ("east",  "\\LIVP11",    "Langmuir LIVP11"),
    ("east",  "\\LIVP12",    "Langmuir LIVP12"),
    ("east",  "\\LIVP13",    "Langmuir LIVP13"),
    ("east",  "\\LIVP14",    "Langmuir LIVP14"),
    ("east",  "\\LIVP15",    "Langmuir LIVP15"),
    ("east",  "\\uovp",      "Upper outer Vprobe"),
    ("east",  "\\lovf",      "Lower outer Vfloat"),

    # ════════════════════════════════════════════════════════════════════════
    # 9. GAS PUFFING Y VACÍO (east)
    # ════════════════════════════════════════════════════════════════════════
    ("east",  "\\vodpev1",   "Gas valve outer D lower 1"),
    ("east",  "\\vodpev2",   "Gas valve outer D lower 2"),
    ("east",  "\\vodpev3",   "Gas valve outer D lower 3"),
    ("east",  "\\vhdpev1",   "Gas valve outer D upper 1"),
    ("east",  "\\voupev1",   "Gas valve outer upper 1"),
    ("east",  "\\voupev2",   "Gas valve outer upper 2"),
    ("east",  "\\voupev3",   "Gas valve outer upper 3"),
    ("east",  "\\vcdpev1",   "Gas valve center D 1"),
    ("east",  "\\JJHG1",     "Gas puffing JJHG1"),
    ("east",  "\\OUG1",      "Gas/vacuum outer upper 1"),
    ("east",  "\\ODFG1",     "Gas puffing ODFG1"),
    ("east",  "\\smbi1",     "SMBI inyector 1"),
    ("east",  "\\smbi2",     "SMBI inyector 2"),
    ("east",  "\\smbi3",     "SMBI inyector 3"),

    # ════════════════════════════════════════════════════════════════════════
    # 10. CALENTAMIENTO RF (east + árboles dedicados)
    # ════════════════════════════════════════════════════════════════════════
    # LHW
    ("east",  "\\PLHI1",     "LHW 2.45GHz inyectada [kW]"),
    ("east",  "\\PLHR1",     "LHW 2.45GHz reflejada [kW]"),
    ("east",  "\\PLHW1",     "LHW 2.45GHz potencia neta"),
    ("east",  "\\PLHI2",     "LHW 4.6GHz inyectada [kW]"),
    ("east",  "\\PLHR2",     "LHW 4.6GHz reflejada [kW]"),
    ("east",  "\\PLHW2",     "LHW 4.6GHz potencia neta"),
    # ECRH
    ("ecrh_east", "\\PECRH1I", "ECRH girotron 1 [kW]"),
    ("ecrh_east", "\\PECRH2I", "ECRH girotron 2 [kW]"),
    ("ecrh_east", "\\PECRH3I", "ECRH girotron 3 [kW]"),
    ("ecrh_east", "\\PECRH4I", "ECRH girotron 4 [kW]"),
    # ICRF
    ("icrf_east", "\\PICRF1",  "ICRF potencia 1 [kW]"),
    ("icrf_east", "\\picrf1",  "ICRF potencia 1 (minúsc)"),
    ("east",      "\\PICRF1",  "east/PICRF1"),
    ("east",      "\\icrf1",   "east/icrf1"),
    # NBI
    ("nbi_east",  "\\PNBI1",   "NBI potencia 1 [kW]"),
    ("nbi_east",  "\\pnbi1",   "NBI potencia 1 (minúsc)"),
    ("nbi_east",  "\\PNBI_TOT","NBI potencia total"),
    ("east",      "\\PNBI1",   "east/PNBI1"),

    # ════════════════════════════════════════════════════════════════════════
    # 11. THOMSON SCATTERING (east)
    # ════════════════════════════════════════════════════════════════════════
    ("east",  "\\TS3",       "Core TS raw signal"),
    ("east",  "\\TS5",       "Edge TS raw signal"),
    # Procesados (si existen)
    ("east",  "\\te_ts3",    "Te Core TS procesado"),
    ("east",  "\\ne_ts3",    "ne Core TS procesado"),
    ("east",  "\\te_ts5",    "Te Edge TS procesado"),
    ("east",  "\\ne_ts5",    "ne Edge TS procesado"),

    # ════════════════════════════════════════════════════════════════════════
    # 12. MHD — MIRNOVS Y OTROS (east)
    # ════════════════════════════════════════════════════════════════════════
    ("east",  "\\mhd01",     "Bobina Mirnov 01"),
    ("east",  "\\mhd08",     "Bobina Mirnov 08"),
    ("east",  "\\mhd16",     "Bobina Mirnov 16"),
    ("east",  "\\bpmb01",    "Sonda magnética poloidal 01"),
    ("east",  "\\bpmt01",    "Sonda magnética toroidal 01"),
    ("east",  "\\vloop",     "Voltaje de lazo east"),

    # ════════════════════════════════════════════════════════════════════════
    # 13. DIAGNÓSTICOS ADICIONALES
    # ════════════════════════════════════════════════════════════════════════
    # Reflectómetro (perfil ne de borde)
    ("east",  "\\refl01",    "Reflectómetro canal 01"),
    ("east",  "\\ne_refl",   "ne borde reflectómetro"),
    # Neutrones (performance nuclear)
    ("east",  "\\neut",      "Tasa de neutrones"),
    ("east",  "\\neutron",   "Tasa de neutrones (alt)"),
    # Espejo IR (temperatura de pared)
    ("east",  "\\ir01",      "Cámara IR canal 01"),
    # CXRS en árbol east
    ("east",  "\\ti0_cxrs",  "east/Ti0 CXRS"),
    ("east",  "\\vt0_cxrs",  "east/Vt0 CXRS"),
]


# ─── Test ─────────────────────────────────────────────────────────────────────

def classify(n_samples):
    if n_samples < 50_000:
        return "PROC"
    elif n_samples < 500_000:
        return "MED "
    else:
        return "RAW "

print(f"\n{'='*70}")
print(f"  Escaneo completo — shot #{SHOT} — {datetime.now():%Y-%m-%d %H:%M}")
print(f"  Servidor: {SERVER}")
print(f"  Nodos a probar: {len(ALL_NODES)}")
print(f"{'='*70}\n")

try:
    conn = mdsthin.Connection(SERVER)
    print("  Conexión OK\n")
except Exception as e:
    print(f"  ERROR: {e}")
    sys.exit(1)

found   = []
missing = []

for tree, node, desc in ALL_NODES:
    try:
        conn.openTree(tree, SHOT)
        d = np.array(conn.get(node).data(), dtype=float)
        t = np.array(conn.get(f"dim_of({node})").data(), dtype=float).ravel()
        cls = classify(d.size)
        t_info = f"t=[{t[0]:.2f},{t[-1]:.2f}]s" if len(t) > 1 else "sin_t"
        print(f"  [{cls}] {tree}/{node}")
        print(f"          {desc}")
        print(f"          shape={d.shape}  {t_info}  "
              f"min={d.min():.3g}  max={d.max():.3g}\n")
        found.append((cls, tree, node, desc, d.shape, t_info,
                      d.min(), d.max()))
    except Exception as e:
        err = str(e).split(",")[0].replace("%TREE-", "")
        missing.append((tree, node, err))
        print(f"  [----] {tree}/{node:20s} → {err}")

# ─── Resumen ──────────────────────────────────────────────────────────────────

print(f"\n{'='*70}")
print(f"  RESUMEN: {len(found)}/{len(ALL_NODES)} nodos con datos")
print(f"  PROC (<50k): {sum(1 for r in found if r[0]=='PROC')}  "
      f"MED (50k-500k): {sum(1 for r in found if r[0]=='MED ')}  "
      f"RAW (>500k): {sum(1 for r in found if r[0]=='RAW ')}")
print(f"{'='*70}\n")

print("  NODOS CON DATOS:")
for cls, tree, node, desc, shape, t_info, vmin, vmax in found:
    print(f"  [{cls}] {tree}/{node:22s} shape={str(shape):12s} {t_info}  [{desc}]")

# ─── Guardar resultado ────────────────────────────────────────────────────────

import os
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "test_all_nodes_result.txt")
with open(out_path, "w", encoding="utf-8") as f:
    f.write(f"Escaneo completo EAST — shot #{SHOT} — {datetime.now():%Y-%m-%d %H:%M}\n")
    f.write(f"Servidor: {SERVER}\n")
    f.write(f"Total: {len(found)}/{len(ALL_NODES)} nodos con datos\n\n")
    f.write("NODOS DISPONIBLES:\n")
    for cls, tree, node, desc, shape, t_info, vmin, vmax in found:
        f.write(f"[{cls}] {tree}/{node}  shape={shape}  {t_info}  "
                f"min={vmin:.3g}  max={vmax:.3g}  [{desc}]\n")
    f.write("\nNODOS NO DISPONIBLES:\n")
    for tree, node, err in missing:
        f.write(f"  {tree}/{node}  → {err}\n")

print(f"\n  Resultado guardado en: test_all_nodes_result.txt\n")

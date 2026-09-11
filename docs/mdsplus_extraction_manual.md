# Manual de Uso — Datos EAST 

---

## 1. Acceso remoto al computador de EAST

Para acceder al computador del laboratorio EAST desde fuera del campus se usa **RustDesk** (cliente de escritorio remoto). El código es: (Ver referencia 1 en Discord)

En Settings -> Network -> ID/Relay Server pegar esto: (Ver referencia 2 en Discord)

OJO: Tener copiado el código pero apretar el botón el portapapeles de la pestaña, no hacer Ctrl+v

Una vez conectado al computador EAST, abrir el script `Extract_PUC` 

Lo único que importa realmente es poder correr códigos desde el computador allá, solo que el archivo mencionado está en un entonrno donde hice pip install de las librerías necesarias, no tuve necesariamente el cuidado de ponerlo en el mismo entorno donde tengo todo listo para usar.

---
Los nodos fueron extraídos del Mdsplus que se abre en Chrome, qué es una especie de interfaz donde se puede explorar gráficos de datos pero no lo encontré muy útil y además no puedo descargar los datos de esa forma.
 Para llegar a estos entrar a Chrome al link y credenciales en (Ver referencia 3 en Discord)

---
## 2. Referencia de nodos MDSplus

Esta tabla consolida el mapeo de diagnósticos del tokamak EAST, verificado experimentalmente contra los shots 143064–143079. Se indica el árbol MDSplus, el nombre del nodo, las unidades y si fue confirmado (✓) o no disponible remotamente (✗).

NOTA: EL ESTADO FUE INDICADO POR CLAUDE QUE HIZO EL ANÁLISIS, HE ENCONTRADO QUE NO TODOS LOS NODOS ESTUVIERON DISPONIBLES PARA USAR DESPUÉS. 

Los nombres de los nodos fueron sacaados de MDsplus

### 2.1 Parámetros globales del plasma

| Señal | Árbol | Nodo | Unidades | Estado |
|-------|-------|------|----------|--------|
| Corriente de plasma Ip | `pcs_east` | `\pcrl01` | kA (raw/1000) | ✓ |
| Voltaje de lazo Vloop | `pcs_east` | `\pcvloop` | V | ✓ |
| Densidad nₑ — cuerda central | `pcs_east` | `\dfsdev` | 10¹⁹ m⁻³ | ✓ |
| Densidad nₑ — cuerda 2 | `pcs_east` | `\dfsdev2` | 10¹⁹ m⁻³ | ✓ |
| Densidad nₑ — cuerda 3 | `pcs_east` | `\dfsdev3` | 10¹⁹ m⁻³ | ✓ |
| Densidad horizontal POINT | `east` | `\point_n3` | a.u. | ✓ |

> Las tres cuerdas de densidad (`dfsdev`, `dfsdev2`, `dfsdev3`) corresponden al interferómetro HCN en el puerto K vertical. No son perfiles radiales completos, sino tres cuerdas de línea integrada a distintas posiciones.

### 2.2 Equilibrio MHD (EFIT)

| Señal | Árbol | Nodo | Unidades | Estado |
|-------|-------|------|----------|--------|
| Factor de seguridad central q₀ | `efit_east` | `\q0` | — | ✓ |
| Factor de seguridad q₉₅ | `efit_east` | `\q95` | — | ✓ |
| Inductancia interna normalizada lᵢ | `efit_east` | `\li` | — | ✓ |
| Energía almacenada W_MHD | `efit_east` | `\WMHD` | J | ✓ |
| Posición eje magnético R_axis | `efit_east` | `\Rmaxis` | m | ✓ |
| Radio menor a | `efit_east` | `\aminor` | m | ✓ |
| Elongación κ | `efit_east` | `\kappa` | — | ✓ |
| Beta normalizada βₙ | `efit_east` | `\betan` | — | ✓ |
| Beta poloidal βₚ | `efit_east` | `\betap` | — | ✓ |
| Separación d_rsep | `efit_east` | `\drsep` | m | ✓ |

> EFIT tiene resolución temporal de ~106 puntos por shot (resolución ~50 ms).

### 2.3 Temperatura electrónica — ECE

| Señal | Árbol | Nodo | Unidades | Estado |
|-------|-------|------|----------|--------|
| Te central (escalar) | `hrs_east` | `\te0_hrs` | keV (NO dividir) | ✓ |
| Te perfil radial — 12 canales | `hrs_east` | `\te_hrs` | keV, shape (11000, 12) | ✓ |
| Posiciones R de canales ECE | `hrs_east` | `\r_hrs` | m, shape (12,) | ✓ |

Los 12 canales ECE cubren R = 1.925 → 2.288 m (lado de bajo campo), con resolución temporal de ~1 ms. El canal 1 (R ≈ 1.925 m) está a ~7.5 cm del eje magnético (R_axis ≈ 1.85 m). Para convertir R a ρ_p se requieren los datos de EFIT (`Rmaxis`, `aminor`).

### 2.4 Temperatura iónica y rotación — XCS

| Señal | Árbol | Nodo | Unidades | Estado |
|-------|-------|------|----------|--------|
| Ti central | `txcs_east` | `\ti0_txcs` | eV (raw/1000 = keV) | ✓ |
| Error Ti central | `txcs_east` | `\ti0_txcserr` | eV | ✓ |
| Rotación toroidal Vt | `txcs_east` | `\vt0_txcs` | km/s | ✓ |
| Error Vt | `txcs_east` | `\vt0_txcserr` | km/s | ✓ |

> Los perfiles radiales de Ti y Vt del sistema CXRS (`cxrs_east`) **no están disponibles**

### 2.5 Radiación e impurezas

| Señal | Árbol | Nodo | Unidades | Estado |
|-------|-------|------|----------|--------|
| Radiación AXUV (canal central) | `east` | `\pxuv32` | a.u. | ✓ |
| Radiación bolómetro | `east` | `\bolo1` | a.u. | ✓ |
| Bremsstrahlung (Zeff proxy) | `east` | `\vbm1` | a.u. | ✓ |
| Filterscope W (W-UTA broadband) | `east` | `\wu1` | a.u. | ✓ |
| Filterscope C (CIII broadband) | `east` | `\ciiil1` | a.u. | ✓ |
| Dα divertor inferior | `east` | `\Dal1` | a.u. | ✓ |
| Dα plano medio | `east` | `\Dam1` | a.u. | ✓ |

> **Importante:** `W_filter` (`\wu1`) y `C_filter` (`\ciiil1`) son señales de fotodiodos con filtros de película delgada — miden una banda ancha de longitudes de onda, **no** una línea espectral específica. No permiten identificar estados de carga individuales de W.

**Sobre el espectrómetro EUV:**  
El árbol `euv_east` existe en MDSplus y contiene el nodo `\wuta` (W-UTA integrado). Sin embargo, **no es accesible**. Los datos del espectrómetro EUV (EUV_Short 5–45 Å, EUV_Long_a 40–180 Å, EUV_Long_b 270–480 Å) Para acceder a ellos es necesario contactar directamente al grupo de diagnósticos de EAST.

### 2.6 Divertor

| Señal | Árbol | Nodo | Unidades | Estado |
|-------|-------|------|----------|--------|
| Sonda triple — Vprobe | `east` | `\lovp01` | V | ✓ |
| Sonda triple — Isat | `east` | `\lois01` | A | ✓ |
| Sonda triple — Vfloat | `east` | `\lovf01` | V | ✓ |

### 2.7 Gas e inyección

| Señal | Árbol | Nodo | Unidades | Estado |
|-------|-------|------|----------|--------|
| SMBI canal 2 | `east` | `\smbi2` | V | ✓ |
| Válvulas gas (máx sobre canales) | `east` | `\vodpev1–3`, `\vhdpev1`, `\voupev1–3`, `\vcdpev1` | V | ✓ |

### 2.8 Calentamiento

| Señal | Árbol | Nodo | Unidades | Estado |
|-------|-------|------|----------|--------|
| ECRH total (4 girotrones) | `ecrh_east` | `\PECRH1I` + `\PECRH2I` + `\PECRH3I` + `\PECRH4I` | kW | ✓ |
| LHW 2.45 GHz — inyectada | `east` | `\PLHI1` | kW | ✓ |
| LHW 2.45 GHz — reflejada | `east` | `\PLHR1` | kW | ✓ |
| LHW 4.6 GHz — inyectada | `east` | `\PLHI2` | kW | ✓ |

> LHW neta = `PLHI1 - PLHR1` (se calcula en el script, no existe como nodo propio).

---

## 3. Scripts de extracción de datos

Todos los scripts de extracción viven en la carpeta `Extracted_EAST/Scripts/`. Están diseñados para correr localmente, estos fueron hechos por Claude a medida que iba probando códigos para ver cómo extraer datos desde su computador. Una vez que pude correr códigos en su computador, estos se transiferen al computador propio mediante la misma aplicación de RustDesk

### `fetch_all_liu_shots.py` — Descarga principal

Descarga todas las señales verificadas para los 7 shots y las guarda como CSV.

```bash
cd Extracted_EAST/Scripts/
python fetch_all_liu_shots.py
```

**Output:** `Extracted_EAST/CSV_Shots/shot_<N>_data.csv` — un CSV por shot con todas las señales temporales alineadas a una grilla común de 1 ms.

**Columnas generadas:**

| Columna | Descripción |
|---------|-------------|
| `Time (s)` | Tiempo |
| `Ip (kA)` | Corriente de plasma |
| `Vloop (V)` | Voltaje de lazo |
| `ne (1e19 m-3)` | Densidad cuerda central |
| `ne2`, `ne3` | Densidad cuerdas 2 y 3 |
| `Te0 (keV)` | Temperatura electrónica central |
| `Ti0 (keV)` | Temperatura iónica central |
| `q0`, `q95`, `li` | Parámetros de equilibrio |
| `WMHD (kJ)` | Energía almacenada |
| `Rmaxis (m)`, `aminor (m)`, `kappa` | Geometría |
| `Da1/2/3 (a.u.)` | Dα divertor |
| `Prad_bolo`, `Prad_axuv` | Radiación (a.u.) |
| `Zeff_brem (a.u.)` | Proxy de Zeff |
| `W_filter`, `C_filter` | Filterscopes W y C |
| `Tt_Vprobe/Isat/Vfloat` | Sonda triple divertor |
| `PEC (kW)` | Potencia ECRH total |
| `PLHI1`, `PLHR1`, `PLHW_net (kW)` | Potencia LHW |

---

### `fetch_ts_profiles.py` — Perfiles radiales de Te (ECE)

Descarga los perfiles radiales de temperatura electrónica del sistema ECE de 12 canales para todos los shots.

```bash
python fetch_ts_profiles.py
```

**Output:** en `Extracted_EAST/ECE/`

| Archivo | Descripción |
|---------|-------------|
| `shot_<N>_ece_profile.csv` | Perfil Te(R,t): columnas Time(s), R_ch01–R_ch12(m), Te_ch01–Te_ch12(keV) |
| `shot_<N>_ece_map.png` | Mapa de color Te(R,t) |
| `shot_<N>_ece_profiles.png` | Perfiles Te(R) en snapshots cada 0.5 s |
| `summary_ece_Te.png` | Todos los shots superpuestos a t ≈ 2 s |

**Cobertura radial:** R = 1.925 → 2.288 m (12 canales uniformemente distribuidos en el lado de bajo campo).

---

### `analyze_liu_shots.py` — Análisis y visualización

Genera todas las figuras de análisis a partir de los CSVs descargados.

```bash
python analyze_liu_shots.py
```

**Requiere:** que `fetch_all_liu_shots.py` haya corrido primero.

**Figuras generadas:**

| Figura | Contenido |
|--------|-----------|
| `fig1_timeseries.png` | Trazas temporales de Ip, nₑ, Te0, ECRH, LHW, Prad |
| `fig2_impurities.png` | W_filter, C_filter, Zeff_brem vs tiempo |
| `fig3_ne_nG.png` | nₑ/n_Greenwald por shot |
| `fig4_divertor.png` | Señales de divertor (sonda triple, Dα) |

---

### Scripts de exploración de nodos (`test_*.py`)

Estos scripts fueron escritos durante la exploración inicial para identificar qué nodos están disponibles en el servidor. **No son necesarios para el análisis rutinario**, pero son útiles si se quiere agregar nuevas señales.

| Script | Propósito |
|--------|-----------|
| `test_all_nodes.py` | Inventario general de todos los nodos probados |
| `test_ts_nodes.py` | Exploración de nodos ECE y Thomson scattering |
| `test_xcs_nodes.py` | Exploración de nodos XCS (Ti, Vt) y XEUV |
| `test_ne_profile_nodes.py` | Exploración de perfiles radiales de nₑ |
| `test_euv_nodes.py` | Exploración de nodos del espectrómetro EUV |
| `test_new_nodes.py` | Exploración de nodos adicionales |

---

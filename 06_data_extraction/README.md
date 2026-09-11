# Data extraction scripts (MDSplus)

Python scripts to extract time traces and diagnostic signals from the EAST
MDSplus server. All scripts are designed to be **executed on the ASIPP
workstation** (Hefei), which is the only machine with direct access to the
MDSplus tree. Once the CSV / NPZ outputs are produced, they are transferred
to the local analysis machine.

See [`docs/mdsplus_extraction_manual.md`](../docs/mdsplus_extraction_manual.md)
for the full workflow and node reference.

## Main extraction scripts

| Script | Role | Output |
|---|---|---|
| `fetch_all_liu_shots.py`   | Downloads global plasma parameters + heating + radiation + divertor for a list of shots. | `shot_<N>_data.csv` |
| `fetch_ts_profiles.py`     | Downloads 12-channel ECE radial profiles of `Te(R, t)`. | `shot_<N>_ece_profile.csv` + plots |
| `analyze_liu_shots.py`     | Generates the standard analysis figures from the CSVs produced by `fetch_all_liu_shots.py`. | 4 PNGs (time series, impurities, ne/nG, divertor) |

## ECE-specific plotting scripts

| Script | Role |
|---|---|
| `plot_all_ece_profiles.py` | Overlay of `Te(R)` at fixed times for all shots. |
| `plot_te0_ece.py`          | Central `Te0` time trace for a single shot. |
| `plot_te_3d.py`            | 3D visualisation of `Te(R, t)`. |

## Node-exploration scripts (`test_*.py`)

These were used to identify which MDSplus nodes are actually available
remotely. They are not needed for routine analysis, but are useful when
adding new signals or verifying node availability in a new campaign.

| Script | Purpose |
|---|---|
| `test_all_nodes.py`         | General inventory across trees. |
| `test_all_nodes_ipp.py`     | Variant intended to run inside the IPP network. |
| `test_ts_nodes.py`          | ECE and Thomson scattering nodes. |
| `test_xcs_nodes.py`         | XCS (Ti, Vt) and XEUV nodes. |
| `test_ne_profile_nodes.py`  | Radial density profiles. |
| `test_euv_nodes.py`         | EUV spectrometer nodes. |
| `test_new_nodes.py`         | Miscellaneous exploration of additional signals. |

## Usage

```bash
# On the ASIPP workstation:
cd 06_data_extraction/
python fetch_all_liu_shots.py     # writes CSVs to ../CSV_Shots/
python fetch_ts_profiles.py        # writes ECE profiles to ../ECE/
python analyze_liu_shots.py        # generates analysis figures

# Then transfer the CSV / PNG outputs to the local analysis machine
# and continue with the Aurora pipeline in the other folders.
```

## Requirements

- `mdsthin`: Python client for MDSplus (`pip install mdsthin`)
- Direct network access to the MDSplus server (`202.127.204.12`), only available
  from within the ASIPP local network in Hefei

## Verified shot range

All node names were verified experimentally against shots
`#143064`–`#143079` (the density-free discharges of Liu et al. 2026).
Node availability may vary in later campaigns; use the `test_*.py` scripts
to re-verify if working with a different shot range.

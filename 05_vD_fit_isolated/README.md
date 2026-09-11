# v/D fit — Isolated pipeline (Fe XXIII, shot #160869)

Self-contained version of the Fe XXIII single-line inversion (Shen 2019 method)
applied to EAST shot #160869. All dependencies are inside this folder — nothing
else in the repository is required to run it.

**Base for future work** on multi-line inversions, alternative parametrisations
of `v(ρ)`, or application to other shots.

## Contents

| File | Role |
|---|---|
| `fit_vD_FeXXIII.py`         | Main script. Runs the fit and produces `fig_vD_result.png`. |
| `profiles_fig4.py`          | Kinetic profiles `ne(ρ), Te(ρ)` (digitised from Shen 2019 Fig. 4, used as proxy). |
| `z_channels.py`             | XEUV chord vertical positions `Z_M` (canonical array). |
| `explore_sif_interactive.py`| Interactive viewer for the `.sif` cube (matplotlib sliders). |

## Data not included

The raw spectrometer file `160869.sif` (17 MB) is not distributed. Place your own
copy in this folder before running:

```bash
cp /path/to/160869.sif .
```

## Pipeline overview

```
160869.sif        profiles_fig4.py       z_channels.py
    │                    │                    │
    ▼                    ▼                    ▼
    Fe XXIII line  ─── Aurora ───→  n_Fe22+(ρ) ─── × ne · PEC ─→  ε(ρ)
                                                                    │
                                                                    ▼
                                              Chord integration (Miller κ, a, z_ax)
                                                                    │
                                                                    ▼
                                                              B_sim(z)
                                                                    │
                                                                    ▼
                                                      Differential evolution fit
                                                       (min ‖B_sim − B_meas‖²)
                                                                    │
                                                                    ▼
                                                   v(ρ) best fit (6 PCHIP knots)
                                                        + D sensitivity scan
```

## Fit configuration (defaults)

| Parameter        | Value                                        |
|------------------|----------------------------------------------|
| `D_CONST_M2S`    | 1.0 m²/s (fixed)                             |
| `V_KNOT_RHOS`    | `{0, 0.2, 0.4, 0.6, 0.8, 1.0}`               |
| `V_BOUNDS`       | `[−3, +3]` m/s per knot                      |
| `D_SCAN`         | `{0.1, 0.3, 0.5, 1.0, 2.0, 5.0}` m²/s        |
| `Z_AXIS_A_CM`    | +8.4 (from C VI fit)                          |
| `Z_AXIS_B_CM`    | +5.5 (from C VI fit)                          |
| Optimiser        | `scipy.optimize.differential_evolution`      |

## Usage

```bash
python fit_vD_FeXXIII.py
```

Runtime: ~5–10 minutes.

## Known result

The single-line fit does not converge to a physical `v/D` (residual cost ~0.14).
This is not a bug — it is the empirical demonstration of the single-line
degeneracy: transport and `T_e`-driven ionisation shift are indistinguishable
when only Fe XXIII is available. The only meaningful contrast between phases
appears at `ρ ≈ 0.6`, qualitatively consistent with the ECRH pump-out signature.

## Interactive exploration

```bash
python explore_sif_interactive.py
```

Opens a matplotlib window with sliders for time and spatial channel, plus
selectable spectral bands (W-UTA, C VI, Fe XXIII, Full). Useful for
understanding the data cube before running any fit.

## Suggested extensions

- **Multi-line inversion** on discharges with Fe XXI, XXII and XXIII all above
  the noise floor → breaks the single-line degeneracy following Shen 2019.
- **Shot-specific kinetic profiles** from integrated data analysis
  (Liu 2024) → replaces the Shen proxy.
- **Alternative parametrisations of v(ρ)**: more knots, physics-motivated shapes
  (neoclassical + turbulent), or `z_axis` as a joint free parameter.

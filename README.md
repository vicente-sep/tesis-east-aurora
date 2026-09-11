# Impurity transport analysis in EAST via EUV spectroscopy

Python codes accompanying the undergraduate thesis by **Vicente Sepúlveda Bustos** (2026),
supervised by German Vogel, at Pontificia Universidad Católica de Chile and Institute of
Plasma Physics, Chinese Academy of Sciences (ASIPP).

The repository implements an [Aurora](https://github.com/fsciortino/Aurora)-based
framework for high-Z impurity transport analysis, validated against
Vogel~2021 (Fe with RMPs), Shen~2019 (Mo with ECRH) and Hu~2018 (Ar coronal equilibrium),
and applied to EAST shot~#160869 (ECRH switch-off event).

## Repository structure

```
.
├── common/                     Shared Miller geometry module
├── 01_vogel_replication/       Fe transport under RMPs — validation
├── 02_shen_replication/        Mo transport under ECRH — validation
├── 03_atomic_validation/       Aurora atomic data vs Hu 2018 (Ar)
├── 04_shot_analysis/           Analysis of EAST shot #160869
├── 05_vD_fit_isolated/         Self-contained v/D fit pipeline
├── 06_data_extraction/         MDSplus extraction scripts (ASIPP-side)
├── docs/                       Practical guides (Aurora, MDSplus, chord integration)
└── examples/                   Sample outputs (PNGs of reference figures)
```

## Getting started

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

Note that [Aurora](https://github.com/fsciortino/Aurora) requires a Fortran compiler
and downloads OpenADAS atomic data on first run (internet needed).

### 2. Choose a starting point

- **New to the framework?** Start by running the replication of Vogel~2021:
  ```bash
  cd 01_vogel_replication/
  python plot_inputs_figure.py     # visualise the digitised inputs
  python run_aurora_Fe.py           # reproduce Fig. 10 of the paper
  ```

- **Working on a new shot?** Use `05_vD_fit_isolated/` as the base for a
  self-contained forward-projection + fit pipeline. Everything is in one folder,
  documented in its own `README.md`.

- **Extending the theory?** The Aurora practical guide and the chord-integration
  derivation are in `docs/`.

## Data not included

The raw `.sif` spectrometer files and MDSplus `.npz` archives belong to EAST/ASIPP
and are not distributed with this repository. Contact the ASIPP impurity spectroscopy
group to request access:

- MDSplus server: `202.127.204.12` (ASIPP local network)
- Space-resolved EUV spectrometer (XEUV_Long): described in
  [Vogel et al., IEEE Trans. Plasma Sci. 46, 1350 (2018)](https://doi.org/10.1109/TPS.2018.2814586)
  and [Cheng et al., Plasma Phys. Control. Fusion 68, 045039 (2026)](https://doi.org/10.1088/1361-6587/ae5d6b)

## Key references

- Sciortino et al., *Nuclear Fusion* **60**, 126014 (2020) — Aurora
- Vogel et al., *IEEE Trans. Plasma Sci.* **46**, 1350 (2018) — XEUV_Long spectrometer
- Vogel et al., *J. Plasma Phys.* **87**, 905870213 (2021) — Fe/RMP transport
- Shen et al., *Phys. Plasmas* **26**, 032507 (2019) — Mo/ECRH transport
- Hu et al., *Rev. Sci. Instrum.* **89**, 10F110 (2018) — Ar coronal equilibrium

## Citation

If you use this code in your work, please cite:

> Sepúlveda Bustos, V. (2026). *Impurity transport analysis in the EAST tokamak
> via EUV spectroscopy and Aurora forward modeling*. Undergraduate thesis,
> Pontificia Universidad Católica de Chile.

## License

MIT — see [LICENSE](LICENSE).

## Contact

Vicente Sepúlveda Bustos — vicentesepulvedabustos@gmail.com

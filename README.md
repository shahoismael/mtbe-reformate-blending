# MTBE-blended Al-Daura catalytic reformate (version 2.1)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24427383.svg)](https://doi.org/10.5281/zenodo.24427383)

Data, open-source Python models and Aspen HYSYS cases for:

**Suitability of MTBE-blended Al-Daura catalytic reformate as a finished gasoline: a specification-based simulation study**

Serwan Ibrahem Abdulkhader, Shaho Ismael Hassen
Chemical and Petrochemical Engineering Department, College of Engineering, Salahaddin University-Erbil, Kurdistan Region, Iraq
serwan.abdulkhader@su.edu.krd, shaho.hassen@su.edu.krd

Zenodo (all versions): https://doi.org/10.5281/zenodo.24427383
Version 2.0: https://doi.org/10.5281/zenodo.24444962 · Version 1.0 (HYSYS workflow): https://doi.org/10.5281/zenodo.24427384

## What is new in version 2.1

- One fitted surrogate composition (`surrogate_<A|B>_rayleigh_fit.json`) is used for both vapour pressure and distillation.
- Reflux-ratio bound in the surrogate fit widened from 20 to 100 (surrogate A no longer at the bound; R = 37.2).
- Calibration and blending data checked against Hamid and Ali (1995), Mustafa (1989), Kopp et al. (1993) and Ershov et al. (2020).
- Figures renumbered: `fig3_validation`, `fig4_windows`.

## What is new in version 2.0

| Part | Script | Result |
|---|---|---|
| 1 Surrogate fit | `src/surrogate_fit.py` | Surrogates A and B fitted to measured D86, RVP, density, PIONA |
| 2 Vapour pressure | `src/blend_vapor_pressure.py` | RVP and TVP of blends: ideal, UNIFAC, NRTL fitted to binary VLE |
| 3 D86 model | `src/surrogate_fit_rayleigh.py`, `src/blend_distillation.py` | Batch distillation with one equilibrium stage, delta-corrected |
| 4 Reflux calibration | `src/calibrate_reflux.py` | Reflux fitted to measured MTBE blends (T50, E70) |
| 5 Octane uncertainty | `src/octane_uncertainty.py` | Monte Carlo grade bounds |
| 6 Windows | `src/windows.py` | MTBE blending windows |
| 7 Figures, tables | `src/figures.py`, `src/tables.py` | `figures/`, `results/tables/` |

## Folders

| Folder | Content |
|---|---|
| `data/processed` | Base fuel, specifications, surrogates, component properties (JSON) |
| `data/raw` | Version 1.0 results, parameters and compositions |
| `src` | Python code |
| `tests` | Code checks (pytest) |
| `results` | Outputs of all parts |
| `figures` | Figures (PNG and TIFF, 600 dpi) |
| `cases/hysys` | Version 1.0 Aspen HYSYS V11 case files |

## Run

```
pip install -r requirements.txt
python src/component_props.py
python src/surrogate_fit.py A
python src/surrogate_fit.py B
python src/blend_vapor_pressure.py A_fitted
python src/blend_vapor_pressure.py B_fitted
python src/blend_vapor_pressure.py A_original
python src/blend_vapor_pressure.py B_original
python src/surrogate_fit_rayleigh.py A
python src/surrogate_fit_rayleigh.py B
python src/calibrate_reflux.py A
python src/calibrate_reflux.py B
python src/blend_distillation.py A
python src/blend_distillation.py B
python src/octane_uncertainty.py
python src/windows.py
python src/figures.py
python src/tables.py
python -m pytest tests -q
```

Run all commands from the repository root.

## Citation

Abdulkhader, S. I., & Hassen, S. I. (2026). Suitability of MTBE-blended Al-Daura catalytic reformate as a finished gasoline: data, open-source models and Aspen HYSYS cases (Version 2.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.24427383

## License

CC BY 4.0

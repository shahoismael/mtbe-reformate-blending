# MTBE-blended Al-Daura catalytic reformate

Data, Aspen HYSYS cases and Python scripts for:

**Suitability of MTBE-blended Al-Daura catalytic reformate as a finished gasoline: a specification-based simulation study**

Serwan Ibrahem Abdulkhader, Shaho Ismael Hassen
Chemical and Petrochemical Engineering Department, College of Engineering, Salahaddin University-Erbil, Kurdistan Region, Iraq
serwan.abdulkhader@su.edu.krd, shaho.hassen@su.edu.krd

Zenodo: https://doi.org/10.5281/zenodo.23438232

## Contents

| Path | Content |
|---|---|
| `HYSYS_cases/` | Aspen HYSYS V11 case files (routes 1–7, surrogates A and B) |
| `Route1*` – `Route7*` | Results by distillation route (D86, assay, EFV, batch) |
| `Route8_rectification/` | Batch cascade reproduced in Python (Table S16) |
| `Route9_RVP_fitted/` | RVP with fitted MTBE parameters |
| `MTBE_binary_parameters_*` | Fitted NRTL and Peng–Robinson MTBE–hydrocarbon parameters |
| `Surrogate_*_composition.xlsx` | Surrogate compositions |
| `Octane_and_E70_rates.xlsx` | Octane sensitivity and measured-blend comparison |
| `scripts/` | Python scripts |
| `Supporting_Information_Tables_S1-S16.pdf` | Supplementary tables |

Full description: `README.txt`. File list: `FILE_LIST.txt`.

## Requirements

Python 3 with `numpy`, `scipy`, `openpyxl`; Route 8 also needs `thermo` and `chemicals`.

```
pip install numpy scipy openpyxl thermo chemicals
python scripts/<name>.py [arguments]
```

Run scripts from the repository root.

## Citation

Abdulkhader, S. I., & Hassen, S. I. (2026). Suitability of MTBE-blended Al-Daura catalytic reformate as a finished gasoline: data and simulation files [Data set]. Zenodo. https://doi.org/10.5281/zenodo.23438232

## License

CC BY 4.0

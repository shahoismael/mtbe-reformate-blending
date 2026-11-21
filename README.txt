Suitability of MTBE-blended Al-Daura catalytic reformate as a finished gasoline:
a specification-based simulation study
Data, calculation files and Aspen HYSYS cases (version 1.0)

Serwan Ibrahem Abdulkhader, Shaho Ismael Hassen
Chemical and Petrochemical Engineering Department, College of Engineering,
Salahaddin University-Erbil, Kurdistan Region, Iraq (serwan.abdulkhader@su.edu.krd, shaho.hassen@su.edu.krd)
License: CC BY 4.0

HOW TO CITE
  Abdulkhader, S. I., & Hassen, S. I. (2026). Suitability of MTBE-blended Al-Daura catalytic reformate as a finished gasoline:
  data and simulation files [Data set]. Zenodo. https://doi.org/10.5281/zenodo.23438232

All section, table and figure numbers below refer to the article. Supplementary tables (S1-S16) are in
Supporting_Information_Tables_S1-S16.pdf. Some files cover analyses that the article reports only in text
or in the supplementary tables; this is noted where it applies.


1. OVERVIEW
MTBE (0-20 vol%) was blended in Aspen HYSYS V11 into two real-component surrogates of one measured sample of
Al-Daura catalytic reformate, and every simulated response was used only as a change from the measured base fuel
(delta method, article Section 2.3, Eqs. 1-2).
  Surrogate A: 22 components (Table S1); used for all routes.
  Surrogate B: 30 components, built independently (Table S10); used for the vapour-pressure and batch routes.
Distillation routes: (1) surrogate + D86 utility, (2) Oil Manager assay of the measured D86 curve + D86 utility,
(3) surrogate + equilibrium flash vaporization (EFV), (4) assay + EFV, and (5-7) surrogate + 23-stage batch
(Rayleigh) cascade with four MTBE parameter sets: NRTL with UNIFAC-estimated parameters, Peng-Robinson kij = 0,
NRTL fitted to binary VLE data, Peng-Robinson with fitted kij.

Base fuel: measured properties from Abdulwahhab et al. (2025), J. Pet. Res. Stud. 15(2) 74-95, Table 1
(doi:10.52716/jprs.v15i2.935), CC BY 4.0. Base MON (82) and LHV (43 MJ/kg) are assumed (article, Section 2.1).


2. ASPEN HYSYS CASE FILES  (folder HYSYS_cases/, list in HSC_FILE_LIST.txt)
  01-05  Routes 1-4: surrogate A D86 (NRTL and PR), assay D86, surrogate A EFV, assay EFV.
  06-09  Batch route, surrogate A: NRTL-UNIFAC, NRTL fitted, PR kij = 0, PR fitted.
  10-13  Batch route, surrogate B: NRTL-UNIFAC, NRTL fitted, PR kij = 0, PR fitted.
  superseded/  Two earlier batch cases, kept for transparency; not used for any reported result.
  HYSYS_case_check_8_files.xlsx  Identifies the eight batch cases (TVP at 0 vol% gives the property package,
      Delta-TVP and RVP at 20 vol% give the parameter set) and lists the MTBE row and column of the NRTL and
      PR matrices in each file's component order.


3. RESULTS BY ROUTE
Route1_Surrogate_D86_NRTL.xlsx
    Surrogate A, D86 utility, NRTL (route S). Article: Tables 3-4, Figs 1-3; Tables S1-S3.
    Sheets: HYSYS_results (raw output), D86_curves, Blend_calc (delta correction and derived properties),
    Window, Surrogate composition, Base validation, Inputs. In Blend_calc, T10-T90 (blue) are the article values
    of Table S3 (evaporated-volume correction, Eq. 2); E70, E100, E150, RVP, density, FBP, DI, VLI, oxygen and
    LHV (base 43 MJ/kg, assumed) are computed by the workbook formulas and agree with Tables 3 and S3.
Route1b_Surrogate_D86_PengRobinson.xlsx
    Same with Peng-Robinson (kij = 0 for MTBE pairs). Article: Sections 3.1 and 3.4; Table S5.
Route2_Assay_OilManager_D86.xlsx, Route2_Assay_OilManager_D86_results.xlsx
    Oil Manager assay of the measured D86 curve and density (route OM); delta-corrected results (Table S6),
    computed by scripts/route2_assay_delta.py.
Route3_Surrogate_EFV_NRTL.xlsx
    Surrogate A flashed at 101.325 kPa (molar vapour fractions 0.02-0.98), route EFV. Evaporated volume = standard
    liquid volume of vapour / feed. Article: Section 3.4; Table S7.
Route4_Assay_OilManager_EFV_NRTL.xlsx
    Assay flashed at equilibrium (NRTL). Table S8.

Batch (Rayleigh) route (route B, article Sections 2.3 and 3.4), 23 heater-separator stages at 101.325 kPa (40-90 C in 5 C steps, 100-200 C in
10 C steps and 125 C), 21 states of 0-20 vol% MTBE. Columns: vapour standard liquid volume of each stage,
blend flow, TVP at 25/35/45 C. -32767 marks a non-converged stage. Delta-correct each file as supplied with
  python scripts/route5_rayleigh_delta.py <file>.csv
(states with a non-converged stage are skipped). For Route6_SurrogateA (Table S11 note) use
  python scripts/route5_rayleigh_delta.py Route6_SurrogateA_Batch_TVP_PR_kij0.csv --fill --base Route7b_SurrogateA_Batch_TVP_PR_fitted.csv
  Route5_SurrogateA_Batch_TVP_NRTL_UNIFAC.csv       Table S9
  Route5b_SurrogateB_Batch_TVP_NRTL_UNIFAC.csv      Table S10
  Route6_SurrogateA_Batch_TVP_PR_kij0.csv           Table S11 (stages above 100 C not converged in 9 of 21 states)
  Route6b_SurrogateB_Batch_TVP_PR_kij0.csv          Table S11
  Route7_SurrogateA_Batch_TVP_NRTL_fitted.csv       Table S15
  Route7b_SurrogateA_Batch_TVP_PR_fitted.csv        Table S15
  Route7c_SurrogateB_Batch_TVP_NRTL_fitted.csv      Table S15
  Route7d_SurrogateB_Batch_TVP_PR_fitted.csv        Table S15
  Base check at 0 vol% MTBE (TVP at 25/35/45 C): surrogate A 32.46/45.65/62.82 kPa (NRTL) and
  33.24/46.30/63.06 kPa (PR); surrogate B 33.55/46.96/64.35 kPa (NRTL) and 34.34/47.52/64.36 kPa (PR).
  Route5_StepCheck_2.5C_SurrogateA_NRTL.csv  Earlier step-size check on a superseded surrogate A version
      (2,2-dimethylbutane in place of 2-methyl-2-butene); replaced by the check in Route8_rectification/.

Route8_rectification/  (Table S16; Table S9 note)
    Batch cascade of surrogate A reproduced outside HYSYS (Python thermo 0.6.1 / chemicals 1.5.2, ideal vapour,
    UNIFAC or ideal liquid). reflux.json and reflux_eq1.json: one partial-condenser stage returning 30% or 50% of
    each stage vapour, evaluated on the evaporated-volume (Eq. 2) and temperature (Eq. 1) axes.
    py_uni_g5.csv and py_uni_g25.csv: step-size check (5 C against 2.5 C steps); py_uni_loss.csv: n-butane
    vapour counted as D86 loss. Same CSV layout as the HYSYS batch files.

Route9_RVP_fitted/  (Section 3.4; Table S14 note)
    RVP at 0-20 vol% MTBE read with the HYSYS stream calculator "Reid VP at 37.8 C" on the blended stream,
    because the Cold Properties utility returned no value with the fitted parameters. Surrogate A: NRTL-UNIFAC,
    NRTL fitted, PR kij = 0, PR fitted; surrogate B: NRTL-UNIFAC, NRTL fitted, PR fitted.
    Delta-RVP at 20 vol% MTBE: surrogate A +3.9, +2.5, +1.2 and +1.5 kPa; surrogate B +3.6, +2.2 and +1.2 kPa.
    For the base-case parameters the calculator gives +3.9 and +1.2 kPa against +3.7 and +1.2 kPa from the
    Cold Properties utility (ASTM D323) used for the main results.


4. MTBE-HYDROCARBON PARAMETERS  (Tables S14-S15)
MTBE_binary_parameters_fitted.csv       NRTL A12/A21 (cal/mol, alpha 0.3, Bij = 0) and Peng-Robinson kij by
                                        hydrocarbon class (as entered, with the unrounded regression values),
                                        with fitted and UNIFAC gamma-infinity and pressure deviations.
MTBE_binary_parameters_by_component.csv The same values assigned to every surrogate component.
MTBE_binary_parameters_HYSYS_entry.xlsx The values in the component order of the surrogate A (23 components incl.
                                        MTBE) and surrogate B (31) case files, ready to paste; Class_values sheet.
Surrogate_A_composition.xlsx            Initial and final composition of surrogate A and its fit to the base fuel.
Surrogate_B_composition.xlsx            Initial and final composition of surrogate B and its fit to the base fuel.


5. OCTANE AND MEASURED-BLEND COMPARISON
Octane_and_E70_rates.xlsx  Table S4 (MON sensitivity of the grade bounds), Badra et al. (2022) model
    extrapolation, Monte Carlo summary, E70 rates of measured blends and back-calculated blending RON,
    MON and RVP of MTBE for the reformate-rich bases of Table 5 (including Gasoline A of Hamid and Ali, 1995).
    The Kopp et al. (1993) pair of Table 5 is computed in scripts/octane_crosschecks.py.


6. SCRIPTS  (folder scripts/; Python 3 with numpy, scipy, openpyxl; Route8 also needs thermo and chemicals)
  route2_assay_delta.py           Delta correction of the assay D86 curves (Table S6).
  route5_rayleigh_delta.py        Delta correction (Eq. 2) of the batch curves (Tables S9-S11, S15).
  second_sample_routes.py         Route changes applied to the second Al-Daura sample (Table S12).
  e70_rates.py                    Linear, PCHIP and cubic interpolation of measured D86 data (Table 5).
  monte_carlo_octane.py           Monte Carlo of octane bounds (Section 3.2, Table S13), numpy seed 1.
  badra_extrapolation.py          Badra et al. (2022) models applied to the reformate (Table S13).
  nrtl_fit_bennett1993.py         NRTL fit to MTBE + 1-heptene at 311 K (Bennett et al., 1993).
  nrtl_fit_mtbe_hydrocarbons.py   NRTL regression on the binary VLE data and UNIFAC check (Table S14).
  pr_kij_fit_chamorro2004.py      Peng-Robinson kij regression for C6+ paraffins and aromatics (Table S14).
  octane_crosschecks.py           Ghosh et al. (2006) octane model on surrogate A, Stratiev et al. (2026)
                                  RON-MON and RON-aromatics correlations, and the Kopp et al. (1993) Base 2/5B
                                  pair (Sections 3.2, 3.4 and 3.5; Table 5). No extra packages needed.
  m2_cascade_core.py, m2_rectification_cascade.py, m2_stepsize_lightend_check.py   Route 8.
  Run all scripts from the package folder: python scripts/<name>.py [arguments]; outputs are written there.
  m2_stepsize_lightend_check.py takes one argument: g5, g25 or loss.


7. OTHER FILES
Supporting_Information_Tables_S1-S16.pdf  Supporting Information of the article, included for convenience.
FILE_LIST.txt                             List of all files in this package.

NOTE: the "README" sheet inside the Route 1, 1b and 2 workbooks is an earlier working template. Where it
differs, the article and this file apply.

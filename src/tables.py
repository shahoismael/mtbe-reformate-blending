"""Part 7: tables for the revised manuscript (CSV).

  table3_blend_properties.csv   properties at 0, 5, 10, 15, 20 vol% MTBE (main case, surrogate A)
  table4_windows.csv            MTBE windows (copy of results/windows/windows_summary.csv)
  table5_validation.csv         model responses against published measurements

Run from the project folder (a few seconds):
    python src/tables.py
Outputs in results/tables/.
"""
import csv, json, os, shutil
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "tables")
LEVELS = (0, 5, 10, 15, 20)


def rows_of(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def pick(rows, key, val, col):
    return {float(r["MTBE vol%"]): float(r[col]) for r in rows if r[key] == val}


def main():
    os.makedirs(OUT, exist_ok=True)
    base = json.load(open(os.path.join(ROOT, "data/processed/base_fuel.json")))
    vp = [r for r in rows_of("results/vapor_pressure/vp_A_fitted.csv")]
    d86 = rows_of("results/distillation/d86_A.csv")
    rb, rm = base["measured"]["density_kg_m3"]["value"], base["MTBE"]["density_kg_m3"]
    wo, lhv_b, lhv_m = base["MTBE"]["oxygen_mass_fraction"], base["assumed"]["LHV_MJ_kg"], base["MTBE"]["LHV_MJ_kg"]
    ron0, mon0 = base["measured"]["RON"]["value"], base["assumed"]["MON"]
    bron, bmon = base["MTBE"]["blending_RON"], base["MTBE"]["blending_MON"]
    aro = base["measured"]["composition_vol_pct"]["aromatics"]

    lines = [["Property"] + [f"{x} vol%" for x in LEVELS]]
    def add(name, f, nd=1):
        lines.append([name] + [round(f(x), nd) for x in LEVELS])
    for mdl, lab in (("nrtl", "RVP, NRTL fitted (kPa)"), ("unifac", "RVP, UNIFAC (kPa)"), ("ideal", "RVP, ideal (kPa)")):
        d = pick(vp, "model", mdl, "RVP delta-corrected (kPa)"); add(lab, lambda x: d[x])
    rho = lambda x: (1 - x / 100) * rb + x / 100 * rm
    add("Density at 15 °C, ideal mixing (kg/m3)", rho)
    add("Oxygen (wt%)", lambda x: 100 * wo * (x / 100 * rm) / rho(x), 2)
    for col, lab, nd in (("T10", "D86 T10 (°C) a", 1), ("T30", "D86 T30 (°C)", 1), ("T50", "D86 T50 (°C)", 1),
                         ("T70", "D86 T70 (°C)", 1), ("T90", "D86 T90 (°C)", 1), ("E70", "E70 (vol%)", 1),
                         ("E100", "E100 (vol%)", 1), ("E150", "E150 (vol%)", 1), ("DI (F)", "DI (°F)", 0)):
        d = pick(d86, "case", "nrtl_R_blend", col); add(lab, lambda x: d[x], nd)
    add("Aromatics (vol%)", lambda x: aro * (1 - x / 100), 2)
    add("RON, linear blending b", lambda x: (1 - x / 100) * ron0 + x / 100 * bron)
    add("ΔMON from base, linear blending b", lambda x: x / 100 * (bmon - mon0))
    def dlhv(x):
        mb, mm = (1 - x / 100) * rb, x / 100 * rm
        return 100 * ((mb * lhv_b + mm * lhv_m) / (rb * lhv_b) - 1)
    add("Δ volumetric LHV (%)", dlhv)
    lines.append(["a T10 response overstated by the model (Section on validation); b reference case, base MON assumed 82, "
                  "blending RON 119, MON 103; octane uncertainty in table4_windows.csv"])
    with open(os.path.join(OUT, "table3_blend_properties.csv"), "w", newline="", encoding="utf-8-sig") as f:
        csv.writer(f).writerows(lines)

    shutil.copy(os.path.join(ROOT, "results/windows/windows_summary.csv"), os.path.join(OUT, "table4_windows.csv"))

    from blend_distillation import MEASURED
    val = [["Quantity", "Measured (range)", "Model, calibrated reflux (A / B)", "Model, no reflux (A / B)",
            "Model, base-fit reflux (A / B)", "Fitted in calibration"]]
    cases = {s: json.load(open(os.path.join(ROOT, f"results/distillation/d86_{s}_summary.json")))["cases"] for s in ("A", "B")}
    for k_model, k_meas, name, fitted in (("dT10_10vol_C", "dT10_10vol_C", "ΔT10 at 10 vol% (°C)", "no"),
                                          ("dT50_10vol_C", "dT50_10vol_C", "ΔT50 at 10 vol% (°C)", "yes"),
                                          ("dE70_per_vol_pct_0to15", "dE70_per_vol_pct", "ΔE70 per vol% MTBE", "yes")):
        m = MEASURED[k_meas]
        f = lambda case: " / ".join(str(cases[s][case][k_model]) for s in ("A", "B"))
        val.append([name, f"{min(m)} to {max(m)}", f("nrtl_R_blend"), f("nrtl_R_0"), f("nrtl_R_calibrated"), fitted])
    for s in ("A", "B"):
        d = json.load(open(os.path.join(ROOT, f"results/vapor_pressure/vp_{s}_fitted_summary.json")))
        val.append([f"Blending RVP of MTBE (kPa), surrogate {s}", "54 to 76 (Ershov, Mustafa, Kopp)",
                    d["nrtl"]["blending_RVP_MTBE_20vol_kPa"], f"ideal {d['ideal']['blending_RVP_MTBE_20vol_kPa']}",
                    f"UNIFAC {d['unifac']['blending_RVP_MTBE_20vol_kPa']}", "no"])
    with open(os.path.join(OUT, "table5_validation.csv"), "w", newline="", encoding="utf-8-sig") as f:
        csv.writer(f).writerows(val)
    for p in ("table3_blend_properties.csv", "table4_windows.csv", "table5_validation.csv"):
        print("saved results/tables/" + p)


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()

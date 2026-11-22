"""Part 6: MTBE blending windows from the results of Parts 2-5.

Reads:
  results/vapor_pressure/vp_<A|B>_fitted.csv          (Part 2)
  results/distillation/d86_<A|B>_summary.json          (Parts 3-4; main case nrtl_R_blend)
  results/octane/octane_summary.json                   (Part 5)
  data/processed/base_fuel.json, specifications.json

Upper bound: oxygen limit (3.7 wt%), checked against RVP, density, aromatics and olefins.
Compliance window (Kurdistan specification) per grade: lower bound from octane (median and
5th-95th percentile) and from aromatics, counting the 1.28 vol% unknown as aromatic.
Performance window (EN 228 E100, ASTM D4814 DI): lower bound from the main distillation case,
with the range over all distillation cases and both surrogates. E70 and T10 are reported
separately (E70 not reached in the main case; T10 response not reproduced by the model).

Run from the project folder (a few seconds):
    python src/windows.py
Outputs in results/windows/.
"""
import csv, json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "windows")
J = lambda p: json.load(open(os.path.join(ROOT, p)))


def oxygen_limit(base, spec):
    rb, rm, wo = base["measured"]["density_kg_m3"]["value"], base["MTBE"]["density_kg_m3"], base["MTBE"]["oxygen_mass_fraction"]
    w = spec["common"]["oxygen_wt_pct_max"] / 100 / wo           # MTBE mass fraction at the limit
    return 100 * w / rm / (w / rm + (1 - w) / rb)               # vol% MTBE


def rvp_check(spec):
    out = {}
    for s in ("A", "B"):
        with open(os.path.join(ROOT, f"results/vapor_pressure/vp_{s}_fitted.csv")) as f:
            for r in csv.DictReader(f):
                if float(r["MTBE vol%"]) == 20:
                    out[f"{s}_{r['model']}"] = float(r["RVP delta-corrected (kPa)"])
    return out


def main():
    base, spec = J("data/processed/base_fuel.json"), J("data/processed/specifications.json")
    o2 = oxygen_limit(base, spec)
    comp = base["measured"]["composition_vol_pct"]
    aro_lim, ole_lim = spec["common"]["aromatics_vol_pct_max"], spec["common"]["olefins_vol_pct_max"]
    aro_known, aro_unknown = comp["aromatics"], comp["aromatics"] + comp["unknown"]
    aro_bound = max(0.0, 100 * (1 - aro_lim / aro_unknown))      # MTBE dilution needed if unknown is aromatic
    rvp = rvp_check(spec); rvp_max = max(rvp.values())
    rho_b, rho_m = base["measured"]["density_kg_m3"]["value"], base["MTBE"]["density_kg_m3"]
    rho_o2 = (1 - o2 / 100) * rho_b + o2 / 100 * rho_m

    d86 = {s: J(f"results/distillation/d86_{s}_summary.json")["cases"] for s in ("A", "B")}
    main_case = {s: d86[s]["nrtl_R_blend"] for s in ("A", "B")}
    allv = lambda key: [c[key] for s in d86 for c in d86[s].values()]
    rng = lambda key: (min(x for x in allv(key) if x is not None), max(x for x in allv(key) if x is not None))
    perf_lo = max(max(main_case[s]["first_vol_E100_ge_46"], main_case[s]["first_vol_DI_le_1250"]) for s in main_case)
    oc = J("results/octane/octane_summary.json")["bounds"]

    rows = [["Window", "Lower bound (vol% MTBE)", "Set by", "Range / uncertainty", "Upper bound (vol% MTBE)", "Set by"]]
    for g in ("regular", "premium", "super"):
        b = oc[g]
        lo = max(b["median"], aro_bound)
        setby = "octane" if b["median"] >= aro_bound else "aromatics (unknown counted as aromatic)"
        note = f"octane 5-95%: {b['p5']:.1f}-{b['p95']:.1f}; P(> O2 limit) {b['P_above_O2_limit']:.2f}"
        rows.append([f"Compliance, {g}", round(lo, 1), setby, note, round(o2, 1), "oxygen 3.7 wt%"])
    e100, di = rng("first_vol_E100_ge_46"), rng("first_vol_DI_le_1250")
    rows.append(["Performance (E100, DI)", perf_lo, "E100 >= 46 vol% / DI <= 1250 F (main case)",
                 f"all cases: E100 {e100[0]:.0f}-{e100[1]:.0f}, DI {di[0]:.0f}-{di[1]:.0f}", round(o2, 1), "oxygen 3.7 wt%"])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "windows_summary.csv"), "w", newline="") as f:
        csv.writer(f).writerows(rows)

    checks = {
        "oxygen_limit_vol_pct": round(o2, 2),
        "RVP_at_20vol_kPa": rvp, "RVP_max_kPa": round(rvp_max, 2),
        "RVP_summer_max_kPa": spec["common"]["RVP_kPa"]["summer"][1], "RVP_margin_kPa": round(spec["common"]["RVP_kPa"]["summer"][1] - rvp_max, 2),
        "density_at_O2_limit_kg_m3_ideal_mixing": round(rho_o2, 1),
        "aromatics_base_known_unknown": [aro_known, round(aro_unknown, 2)], "aromatics_limit": aro_lim,
        "MTBE_needed_for_aromatics_if_unknown_aromatic_vol_pct": round(aro_bound, 2),
        "aromatics_at_O2_limit_vol_pct": round(aro_unknown * (1 - o2 / 100), 2),
        "olefins_at_0vol": comp["olefins"], "olefins_limit": ole_lim,
        "E70_main_case_at_20vol": {s: main_case[s]["E70_at_20vol"] for s in main_case},
        "E70_first_vol_ge_20_all_cases": sorted({str(x) for x in allv("first_vol_E70_ge_20")}),
        "T10_note": "model overstates the measured T10 drop about 4x; T10 bound not used",
    }
    json.dump({"windows": rows, "checks": checks}, open(os.path.join(OUT, "windows_summary.json"), "w"), indent=2)
    for r in rows: print(" | ".join(str(x) for x in r))
    for k, v in checks.items(): print(f"  {k}: {v}")


if __name__ == "__main__":
    main()

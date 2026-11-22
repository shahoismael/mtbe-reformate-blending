"""Part 2: vapour pressure of reformate-MTBE blends, 0-20 vol% MTBE, three liquid models.

For each surrogate (original and fitted, A and B) and each model (ideal, unifac, nrtl):
  RVP (D323 basis, 37.8 degC, V/L = 4), TVP at 25, 35 and 45 degC,
  delta-corrected RVP = 49.0 + (RVP_blend - RVP_base), and the volumetric blending RVP of MTBE.

Check against HYSYS (article, Table 3 and Route9): with the original surrogate A, the HYSYS
NRTL case with UNIFAC-estimated parameters gave +3.8 kPa at 20 vol% MTBE, Peng-Robinson
(kij = 0) +1.2 kPa and NRTL with the fitted parameters +2.5 kPa.

Run from the project folder (one surrogate per run, about 1-2 min each):
    python src/blend_vapor_pressure.py A_fitted
    python src/blend_vapor_pressure.py B_fitted
    python src/blend_vapor_pressure.py A_original
    python src/blend_vapor_pressure.py B_original
Outputs in results/vapor_pressure/.
"""
import csv, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_models import Mixture, load_props, ROOT

OUT = os.path.join(ROOT, "results", "vapor_pressure")
MODELS = ("ideal", "unifac", "nrtl")
MTBE_PCT = np.arange(0, 21, 1.0)
RVP_BASE = 49.0


def surrogate(tag):
    which, kind = tag.split("_")
    sur = json.load(open(os.path.join(ROOT, "data/processed/surrogates.json")))["surrogate_" + which]
    names = [c["component"] for c in sur]; classes = [c["class"] for c in sur]
    if kind == "original":
        v = [c["final_vol_pct"] for c in sur]
    else:
        fit = json.load(open(os.path.join(ROOT, "results/surrogate_fit", f"surrogate_{which}_rayleigh_fit.json")))  # same composition as the D86 model
        v = [fit["fitted_vol_pct"][n] for n in names]
    return names, classes, np.array(v, float)


def run(tag):
    names, classes, v = surrogate(tag)
    mix = Mixture(names + ["mtbe"], classes + ["Ether"], load_props())
    rows, summary = [], {}
    for model in MODELS:
        res = []
        for x in MTBE_PCT:
            vol = np.append(v * (100 - x) / v.sum(), x)
            res.append([mix.rvp(vol, model)] + [mix.tvp(vol, t, model) for t in (25.0, 35.0, 45.0)])
        res = np.array(res); d = res - res[0]
        for x, r, dd in zip(MTBE_PCT, res, d):
            rows.append([tag, model, x, round(r[0], 2), round(RVP_BASE + dd[0], 2),
                         round(r[1], 2), round(r[2], 2), round(r[3], 2), round(dd[3], 2)])
        rv20 = RVP_BASE + d[-1, 0]
        summary[model] = {"RVP_base_model_kPa": round(res[0, 0], 2), "dRVP_20vol_kPa": round(d[-1, 0], 2),
                          "RVP_20vol_corrected_kPa": round(rv20, 2),
                          "blending_RVP_MTBE_20vol_kPa": round((rv20 - 0.8 * RVP_BASE) / 0.2, 1),
                          "dTVP45_20vol_kPa": round(d[-1, 3], 2)}
        print(f"{tag:11s} {model:6s} base RVP {res[0,0]:5.1f} | dRVP(20%) {d[-1,0]:+5.2f} kPa | "
              f"blending RVP {summary[model]['blending_RVP_MTBE_20vol_kPa']:5.1f} | dTVP45(20%) {d[-1,3]:+5.2f}")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"vp_{tag}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["surrogate", "model", "MTBE vol%", "RVP model (kPa)", "RVP delta-corrected (kPa)",
                    "TVP 25C (kPa)", "TVP 35C (kPa)", "TVP 45C (kPa)", "dTVP 45C (kPa)"])
        w.writerows(rows)
    json.dump(summary, open(os.path.join(OUT, f"vp_{tag}_summary.json"), "w"), indent=2)


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "A_fitted")

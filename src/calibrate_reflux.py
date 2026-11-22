"""Part 4: calibrate the reflux ratio R on measured MTBE blend responses.

The surrogate (surrogate_fit_rayleigh.py) fixes the base fuel; R controls how MTBE moves the
front of the D86 curve. Here R is chosen so that the delta-corrected model reproduces the mean
measured changes caused by 10-11 vol% MTBE in published blends (article Table 5):
    dT10 at 10 vol%      mean of [-1.0, -1.5, -0.6]  degC
    dT50 at 10 vol%      mean of [-12.0, -4.5, -7.2] degC
    dE70 per vol% MTBE   mean of [0.33, 0.55, 0.56, 0.40]
Each deviation is divided by the spread (standard deviation) of the measured values.
FIT_KEYS sets which responses are fitted. dT10 is excluded: no reflux ratio reproduces the small
measured T10 drop (equilibrium models overstate it; first run of this script), so dT10 is
reported as a check only and the T10 bound is treated as unresolved.
Liquid model: NRTL with fitted MTBE parameters.

Run from the project folder (about 1 min each):
    python src/calibrate_reflux.py A
    python src/calibrate_reflux.py B
Outputs: results/distillation/reflux_calibration_<A|B>.csv and .json
"""
import csv, json, os, sys
import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_models import Mixture, load_props, ROOT
from distillation import D86Model
from blend_distillation import corrected_curves, props_from, MEASURED, OUT, STEP

GRID = [0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
TARGET = {"dT10": float(np.mean(MEASURED["dT10_10vol_C"])), "dT50": float(np.mean(MEASURED["dT50_10vol_C"])),
          "dE70": float(np.mean(MEASURED["dE70_per_vol_pct"]))}
SPREAD = {"dT10": float(np.std(MEASURED["dT10_10vol_C"])), "dT50": float(np.std(MEASURED["dT50_10vol_C"])),
          "dE70": float(np.std(MEASURED["dE70_per_vol_pct"]))}


def setup(which):
    fit = json.load(open(os.path.join(ROOT, "results/surrogate_fit", f"surrogate_{which}_rayleigh_fit.json")))
    sur = json.load(open(os.path.join(ROOT, "data/processed/surrogates.json")))["surrogate_" + which]
    names = [c["component"] for c in sur]; classes = [c["class"] for c in sur]
    v = np.array([fit["fitted_vol_pct"][n] for n in names])
    P = load_props()
    return Mixture(names + ["mtbe"], classes + ["Ether"], P), P, v


def response(R, mix, P, v):
    m = D86Model(mix, P, model="nrtl", R=R)
    Es, Ts = [], []
    for x in (0.0, 10.0, 15.0):
        E, T = m.run(np.append(v * (100 - x) / v.sum(), x), step=STEP)
        Es.append(E); Ts.append(T)
    grid, Ec, _ = corrected_curves(Es, Ts)
    b, r10, r15 = (props_from(grid, E) for E in Ec)
    return {"dT10": r10["T10"] - b["T10"], "dT50": r10["T50"] - b["T50"],
            "dE70": (r15["E70"] - b["E70"]) / 15}


FIT_KEYS = ("dT50", "dE70")


def score(r):
    return float(sum(((r[k] - TARGET[k]) / SPREAD[k]) ** 2 for k in FIT_KEYS))


def main(which):
    mix, P, v = setup(which)
    rows = []
    for R in GRID:
        r = response(R, mix, P, v); s = score(r)
        rows.append([R, round(r["dT10"], 2), round(r["dT50"], 2), round(r["dE70"], 3), round(s, 2)])
        print(f"{which} R {R:5.2f}  dT10 {r['dT10']:+5.2f}  dT50 {r['dT50']:+6.2f}  dE70/vol {r['dE70']:.3f}  score {s:6.2f}")
    i = int(np.argmin([x[-1] for x in rows]))
    lo = GRID[max(i - 1, 0)]; hi = GRID[min(i + 1, len(GRID) - 1)]
    opt = minimize_scalar(lambda R: score(response(R, mix, P, v)), bounds=(lo, hi), method="bounded",
                          options={"xatol": 0.01})
    R_best = float(opt.x); r = response(R_best, mix, P, v)
    print(f"{which} best R {R_best:.3f}  dT10 {r['dT10']:+.2f}  dT50 {r['dT50']:+.2f}  dE70/vol {r['dE70']:.3f}  "
          f"score {opt.fun:.2f}")
    print("target", {k: round(x, 3) for k, x in TARGET.items()})
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"reflux_calibration_{which}.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["R", "dT10 10vol (C)", "dT50 10vol (C)", "dE70 per vol%", "score"]); w.writerows(rows)
    json.dump({"surrogate": which, "R_blend": round(R_best, 4), "fitted_responses": FIT_KEYS, "score": round(float(opt.fun), 3),
               "model_response": {k: round(x, 3) for k, x in r.items()}, "target": TARGET, "spread": SPREAD},
              open(os.path.join(OUT, f"reflux_calibration_{which}.json"), "w"), indent=2)


if __name__ == "__main__":
    main(sys.argv[1].upper() if len(sys.argv) > 1 else "A")

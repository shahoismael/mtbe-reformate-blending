"""Part 3, step 1: refit the surrogate with the same distillation model that is used for the blends.

Same targets, weights and constraints as surrogate_fit.py, but the D86 curve is computed by
the batch-distillation model of distillation.py (flask + one equilibrium stage with reflux
ratio R) instead of the API TBP->D86 conversion, and is matched at T10, T20, ..., T90.
The composition and R are fitted together (R between 0 and 20). The fitted R is the D86
apparatus calibration and is used unchanged for all MTBE blends (blend_distillation.py),
so the base fuel and the blends are computed with one consistent model.

Run from the project folder (about 1 min each):
    python src/surrogate_fit_rayleigh.py A
    python src/surrogate_fit_rayleigh.py B
Outputs in results/surrogate_fit/ (files *_rayleigh_*).
"""
import csv, json, os, sys
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import surrogate_fit as sf
from models import density_15C, rvp_d323, class_totals
from thermo_models import Mixture, load_props
from distillation import D86Model, T_at

PTS = [10, 20, 30, 40, 50, 60, 70, 80, 90]
STEP = 1.0  # vol% per distillation step


def make_model(c):
    names = c["names"] + ["mtbe"]
    mix = Mixture(names, c["classes"] + ["Ether"], load_props())
    return D86Model(mix, load_props(), model="ideal", R=0.0)


def d86(m, v):
    E, T = m.run(np.append(v, 0.0), step=STEP)
    return [T_at(E, T, p) for p in [0] + PTS]


def targets(base):
    m = base["measured"]; d = dict(zip(m["D86"]["evaporated_vol_pct"], m["D86"]["temperature_C"]))
    t = sf.targets(base); t["D86"] = [d[p] for p in PTS]
    return t


def residuals(v, c, t, rho_offset, m):
    p = d86(m, v)
    r = [sf.W["D86"] * (np.array(p[1:]) - t["D86"]),
         [sf.W["RVP"] * (rvp_d323(v, c["rho"], c["MW"], c["psat"]) - t["RVP"]),
          sf.W["density"] * (density_15C(v, c["rho"]) + rho_offset - t["density"])]]
    cl = class_totals(v, c["classes"]); cr = []
    for k, tv in t["classes"].items():
        got = cl.get(k, 0.0)
        if k == "Aromatic":
            lo, hi = t["aro_band"]; cr.append(sf.W["class"] * (min(0.0, got - lo) + max(0.0, got - hi)))
        else:
            cr.append(sf.W["class"] * (got - tv))
    r.append(cr); r.append(sf.W["reg"] * (v - c["v0"]))
    return np.concatenate([np.ravel(x) for x in r])


R_BOUNDS = (0.0, 20.0)
R_START = 1.0


def main(which):
    base, c = sf.load(which); t = targets(base); m = make_model(c)
    rho_offset = sf.HYSYS_ORIGINAL[which]["density"] - density_15C(c["v0"], c["rho"])

    def obj(z):
        m.R = z[-1]
        return float(np.sum(residuals(z[:-1], c, t, rho_offset, m) ** 2))

    bounds = [(x, x) if n.lower() in sf.FIXED else (0.0, 25.0) for n, x in zip(c["names"], c["v0"])] + [R_BOUNDS]
    res = minimize(obj, np.append(c["v0"], R_START), method="SLSQP", bounds=bounds,
                   constraints=[{"type": "eq", "fun": lambda z: np.sum(z[:-1]) - 100.0}],
                   options={"maxiter": 150, "ftol": 1e-8, "eps": 0.01})
    R = float(res.x[-1])
    v = np.clip(res.x[:-1], 0, None); v = 100 * v / v.sum()
    m.R = 0.0; p0 = d86(m, c["v0"])          # original surrogate, simple Rayleigh
    m.R = R; p1 = d86(m, v)
    rows = [["Property", "Measured", "Original surrogate (R = 0)", f"Fitted surrogate (R = {R:.2f})", "Fitted - measured"]]
    meas0 = base["measured"]["D86"]["temperature_C"][0]
    rows.append(["D86 IBP (degC, not fitted)", meas0, round(p0[0], 1), round(p1[0], 1), round(p1[0] - meas0, 1)])
    for i, pt in enumerate(PTS):
        rows.append([f"D86 T{pt} (degC)", t["D86"][i], round(p0[i + 1], 1), round(p1[i + 1], 1), round(p1[i + 1] - t["D86"][i], 1)])
    for name, f, tv in (("RVP (kPa)", lambda x: rvp_d323(x, c["rho"], c["MW"], c["psat"]), t["RVP"]),
                        ("Density 15 degC (kg/m3)", lambda x: density_15C(x, c["rho"]) + rho_offset, t["density"])):
        rows.append([name, tv, round(f(c["v0"]), 1), round(f(v), 1), round(f(v) - tv, 1)])
    cl0, cl1 = class_totals(c["v0"], c["classes"]), class_totals(v, c["classes"])
    for k, tv in t["classes"].items():
        rows.append([f"{k} (vol%)", tv, round(cl0.get(k, 0), 2), round(cl1.get(k, 0), 2), round(cl1.get(k, 0) - tv, 2)])
    os.makedirs(sf.OUT, exist_ok=True)
    with open(os.path.join(sf.OUT, f"surrogate_{which}_rayleigh_fit_residuals.csv"), "w", newline="") as f:
        csv.writer(f).writerows(rows)
    with open(os.path.join(sf.OUT, f"surrogate_{which}_rayleigh_fitted_composition.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["#", "Component", "Class", "NBP (degC)", "Original (vol%)", "Fitted (vol%)", "Change"])
        for i, n in enumerate(c["names"]):
            w.writerow([i + 1, n, c["classes"][i], c["nbp"][i], c["v0"][i], round(v[i], 2), round(v[i] - c["v0"][i], 2)])
    json.dump({"surrogate": which, "d86_model": "flask + 1 stage, ideal liquid", "reflux_ratio_R": round(R, 4),
               "R_at_bound": bool(abs(R - R_BOUNDS[1]) < 1e-3), "weights": sf.W,
               "density_offset_kg_m3": round(rho_offset, 2), "objective": res.fun, "optimizer_message": res.message,
               "fitted_vol_pct": dict(zip(c["names"], [round(x, 3) for x in v]))},
              open(os.path.join(sf.OUT, f"surrogate_{which}_rayleigh_fit.json"), "w"), indent=2)
    print(f"Surrogate {which} (D86 model): R = {R:.2f}, objective {res.fun:.2f} ({res.message})")
    for r in rows: print("  " + " | ".join(str(x) for x in r))


if __name__ == "__main__":
    main(sys.argv[1].upper() if len(sys.argv) > 1 else "A")

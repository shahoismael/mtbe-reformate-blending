"""Part 1, module 1: fit the surrogate composition to the measured base fuel.

Keeps the components of the existing surrogate (so the result can be entered directly in the
HYSYS case) and adjusts their liquid volume fractions by constrained least squares to:
  - D86 T10, T30, T50, T70, T90        (API 1994 TBP->D86 from component boiling points)
  - RVP 49.0 kPa                        (dry vapour pressure, 37.8 degC, V/L = 4)
  - density 762.8 kg/m3                 (ideal volume mixing, calibrated to HYSYS, see below)
  - PIONA class totals (ASTM D6729)     (aromatics allowed between 34.65 and 35.93 vol%,
                                         i.e. with or without the 1.28 vol% unknown)
A regularisation term keeps the composition close to the original surrogate, so that only the
changes the data require are made. Benzene (not measured) is held at its original value.
IBP and end point are reported but not fitted: the D86 IBP lies well above the thermodynamic
bubble point and the end point is distorted by the 5.7 vol% residue.

Density calibration: the open-source densities differ slightly from the HYSYS library, so one
constant offset is applied so that the original surrogate gives the HYSYS value (Base_fit sheet).

Run from the project folder:
    python src/surrogate_fit.py A
    python src/surrogate_fit.py B
Outputs in results/surrogate_fit/.
"""
import json, os, sys, csv
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import density_15C, rvp_d323, d86_from_composition, class_totals

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "surrogate_fit")

# HYSYS results for the original surrogates (data/raw/Surrogate_*_composition.xlsx, Base_fit)
HYSYS_ORIGINAL = {"A": {"density": 761.4, "RVP": 48.4, "D86": [17.3, 62.7, 94.4, 114.2, 137.0, 169.2]},
                  "B": {"density": 763.5, "RVP": 49.4, "D86": [None, 57.3, 91.8, 112.5, 136.7, 170.2]}}

# weights: residual units are degC, kPa, kg/m3, vol%
W = {"D86": 1.0, "RVP": 3.0, "density": 2.0, "class": 4.0, "reg": 0.3}
# benzene was not measured: kept at the original surrogate value, not fitted
FIXED = {"benzene"}
CLASS_KEY = {"n-Paraffin": "n_paraffins", "iso-Paraffin": "iso_paraffins", "Olefin": "olefins",
             "Naphthene": "naphthenes", "Aromatic": "aromatics"}
PTS = [10, 30, 50, 70, 90]


def load(which):
    base = json.load(open(os.path.join(ROOT, "data/processed/base_fuel.json")))
    sur = json.load(open(os.path.join(ROOT, "data/processed/surrogates.json")))["surrogate_" + which]
    P = json.load(open(os.path.join(ROOT, "data/processed/component_props.json")))
    names = [c["component"] for c in sur]
    p = [P[n.lower()] for n in names]
    comp = {"names": names, "classes": [c["class"] for c in sur],
            "v0": np.array([c["final_vol_pct"] for c in sur], float),
            "nbp": np.array([x["Tb_C"] for x in p]), "rho": np.array([x["rho15_kg_m3"] for x in p]),
            "MW": np.array([x["MW"] for x in p]), "psat": np.array([x["psat378_Pa"] for x in p])}
    return base, comp


def predict(v, c, rho_offset):
    d86 = d86_from_composition(v, c["nbp"])
    return {"D86": d86,  # 0,10,30,50,70,90,100
            "RVP": rvp_d323(v, c["rho"], c["MW"], c["psat"]),
            "density": density_15C(v, c["rho"]) + rho_offset,
            "classes": class_totals(v, c["classes"])}


def targets(base):
    m = base["measured"]; d = dict(zip(m["D86"]["evaporated_vol_pct"], m["D86"]["temperature_C"]))
    comp = m["composition_vol_pct"]
    cls = {k: comp[CLASS_KEY[k]] for k in CLASS_KEY}
    return {"D86": [d[p] for p in PTS], "IBP": d[0], "EP": m["D86"]["end_point_C"],
            "RVP": m["vapor_pressure_kPa"]["value"], "density": m["density_kg_m3"]["value"],
            "classes": cls, "aro_band": (comp["aromatics"], comp["aromatics"] + comp["unknown"])}


def residuals(v, c, t, rho_offset, v_ref):
    p = predict(v, c, rho_offset)
    r = [W["D86"] * (np.array(p["D86"][1:6]) - t["D86"])]
    r.append([W["RVP"] * (p["RVP"] - t["RVP"]), W["density"] * (p["density"] - t["density"])])
    cr = []
    for k, tv in t["classes"].items():
        got = p["classes"].get(k, 0.0)
        if k == "Aromatic":
            lo, hi = t["aro_band"]; dev = min(0.0, got - lo) + max(0.0, got - hi)
        else:
            dev = got - tv
        cr.append(W["class"] * dev)
    r.append(cr)
    r.append(W["reg"] * (v - v_ref))
    return np.concatenate([np.ravel(x) for x in r])


def fit(which):
    base, c = load(which)
    t = targets(base)
    p0 = predict(c["v0"], c, 0.0)
    rho_offset = HYSYS_ORIGINAL[which]["density"] - p0["density"]
    obj = lambda v: float(np.sum(residuals(v, c, t, rho_offset, c["v0"]) ** 2))
    cons = [{"type": "eq", "fun": lambda v: np.sum(v) - 100.0}]
    bounds = [(x, x) if n.lower() in FIXED else (0.0, 25.0) for n, x in zip(c["names"], c["v0"])]
    rng = np.random.default_rng(1)
    starts = [c["v0"]] + [c["v0"] * rng.uniform(0.7, 1.3, len(c["v0"])) for _ in range(2)]
    best = None
    for s in starts:
        s = np.array([lo if lo == hi else x for (lo, hi), x in zip(bounds, 100 * s / s.sum())])
        res = minimize(obj, s, method="SLSQP", bounds=bounds, constraints=cons,
                       options={"maxiter": 500, "ftol": 1e-9})
        if res.success and (best is None or res.fun < best.fun):
            best = res
    v = np.clip(best.x, 0, None); v = 100 * v / v.sum()
    return base, c, t, rho_offset, p0, predict(v, c, rho_offset), v, best


def report(which):
    base, c, t, rho_offset, p0, p1, v, res = fit(which)
    os.makedirs(OUT, exist_ok=True)
    h = HYSYS_ORIGINAL[which]
    rows = [["Property", "Measured", "Original (HYSYS)", "Original (this model)", "Fitted (this model)",
             "Fitted - measured"]]
    rows.append(["D86 IBP (degC)", t["IBP"], h["D86"][0], round(p0["D86"][0], 1), round(p1["D86"][0], 1),
                 round(p1["D86"][0] - t["IBP"], 1)])
    for i, pt in enumerate(PTS):
        rows.append([f"D86 T{pt} (degC)", t["D86"][i], h["D86"][i + 1], round(p0["D86"][i + 1], 1),
                     round(p1["D86"][i + 1], 1), round(p1["D86"][i + 1] - t["D86"][i], 1)])
    rows.append(["RVP (kPa)", t["RVP"], h["RVP"], round(p0["RVP"], 1), round(p1["RVP"], 1),
                 round(p1["RVP"] - t["RVP"], 1)])
    rows.append(["Density 15 degC (kg/m3)", t["density"], h["density"], round(p0["density"] + rho_offset, 1),
                 round(p1["density"], 1), round(p1["density"] - t["density"], 1)])
    for k, tv in t["classes"].items():
        rows.append([f"{k} (vol%)", tv, "", round(p0["classes"].get(k, 0), 2), round(p1["classes"].get(k, 0), 2),
                     round(p1["classes"].get(k, 0) - tv, 2)])
    with open(os.path.join(OUT, f"surrogate_{which}_fit_residuals.csv"), "w", newline="") as f:
        csv.writer(f).writerows(rows)
    with open(os.path.join(OUT, f"surrogate_{which}_fitted_composition.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["#", "Component", "Class", "NBP (degC)", "Original (vol%)", "Fitted (vol%)", "Change"])
        for i, n in enumerate(c["names"]):
            w.writerow([i + 1, n, c["classes"][i], c["nbp"][i], c["v0"][i], round(v[i], 2), round(v[i] - c["v0"][i], 2)])
        w.writerow(["", "Total", "", "", round(c["v0"].sum(), 2), round(v.sum(), 2), ""])
    json.dump({"surrogate": which, "weights": W, "density_offset_kg_m3": round(rho_offset, 2),
               "objective": res.fun, "optimizer_message": res.message,
               "fitted_vol_pct": dict(zip(c["names"], [round(x, 3) for x in v]))},
              open(os.path.join(OUT, f"surrogate_{which}_fit.json"), "w"), indent=2)
    print(f"Surrogate {which}: objective {res.fun:.2f}  ({res.message})")
    for r in rows: print("  " + " | ".join(str(x) for x in r))


if __name__ == "__main__":
    report(sys.argv[1].upper() if len(sys.argv) > 1 else "A")

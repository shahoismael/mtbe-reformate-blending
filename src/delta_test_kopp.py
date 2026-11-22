"""Independent test of the delta method on a fuel with measured base and blend data.

Kopp et al. (1993), SAE 930143: AQIRP Base 2 (Table 11) and fuel 5B = 89.0 vol% Base 2 +
11.0 vol% MTBE (Table 12; analysed 10.3 vol%, Table 13). Values checked against the paper.

Steps
  1. Fit surrogate A (same components) and its reflux ratio to Base 2: D86 T10-T90, RVP,
     density and class totals (aromatics, olefins, total saturates; Kopp gives no n/iso split).
  2. Calibrate the blend reflux ratio on the reformate WITHOUT the Kopp data
     (Hamid and Ali, Mustafa only), so the test does not use its own answer.
  3. Predict 5B with the delta method (Eq. 2 for RVP, Eq. 4 for D86) and compare with Table 13.

Run from the project folder (about 2-3 min):
    python src/delta_test_kopp.py
Outputs: results/validation/kopp_delta_test.csv and .json
"""
import csv, json, os, sys
import numpy as np
from scipy.optimize import minimize, minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import surrogate_fit as sf
from models import density_15C, rvp_d323, class_totals
from thermo_models import Mixture, load_props, ROOT
from distillation import D86Model, T_at
import calibrate_reflux as cr
from blend_distillation import STEP, DT

OUT = os.path.join(ROOT, "results", "validation")
F = lambda f: (f - 32) / 1.8
PSI = 6.89476
PTS_ALL = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
FIT_PTS = [10, 20, 30, 40, 50, 60, 70, 80, 90]

# Kopp et al. (1993) Table 11 (Base 2) and Table 13 (5B); D86 in degF
BASE2 = {"D86_F": [93, 115, 125, 135, 147, 165, 189, 216, 238, 253, 270], "RVP_psi": 8.3, "API": 60.9,
         "aromatics": 31.0, "olefins": 19.3, "saturates": 49.7, "benzene": 0.92}
FUEL5B = {"D86_F": [93, 116, 124, 134, 144, 157, 176, 203, 231, 251, 267], "RVP_psi": 8.6}
MTBE_RECIPE, MTBE_ANALYSED = 11.0, 10.3

# reformate calibration targets without Kopp: Hamid and Ali [8], Mustafa [23]
TARGET = {"dT50": float(np.mean([-12.0, -4.5])), "dE70": float(np.mean([0.33, 0.55, 0.56]))}
SPREAD = {"dT50": float(np.std([-12.0, -4.5])), "dE70": float(np.std([0.33, 0.55, 0.56]))}


def density_api(api):
    return 141.5 / (131.5 + api) * 999.016


def fit_base2():
    _, c = sf.load("A")
    v0 = c["v0"].copy()
    ib = [n.lower() for n in c["names"]].index("benzene")
    v0[ib] = BASE2["benzene"]
    others = [i for i in range(len(v0)) if i != ib]
    v0[others] *= (100 - BASE2["benzene"]) / v0[others].sum()
    off = json.load(open(os.path.join(ROOT, "results/surrogate_fit/surrogate_A_rayleigh_fit.json")))["density_offset_kg_m3"]
    tD86 = [F(BASE2["D86_F"][PTS_ALL.index(p)]) for p in FIT_PTS]
    tRVP, tRho = BASE2["RVP_psi"] * PSI, density_api(BASE2["API"])
    mix = Mixture(c["names"] + ["mtbe"], c["classes"] + ["Ether"], load_props())
    m = D86Model(mix, load_props(), model="ideal", R=0.0)
    sat = ("n-Paraffin", "iso-Paraffin", "Naphthene")

    def resid(z):
        v, m.R = z[:-1], z[-1]
        E, T = m.run(np.append(v, 0.0), step=1.0)
        cl = class_totals(v, c["classes"])
        r = [sf.W["D86"] * (np.array([T_at(E, T, p) for p in FIT_PTS]) - tD86),
             [sf.W["RVP"] * (rvp_d323(v, c["rho"], c["MW"], c["psat"]) - tRVP),
              sf.W["density"] * (density_15C(v, c["rho"]) + off - tRho),
              sf.W["class"] * (cl.get("Aromatic", 0) - BASE2["aromatics"]),
              sf.W["class"] * (cl.get("Olefin", 0) - BASE2["olefins"]),
              sf.W["class"] * (sum(cl.get(k, 0) for k in sat) - BASE2["saturates"])],
             sf.W["reg"] * (v - v0)]
        return np.concatenate([np.ravel(x) for x in r])

    bounds = [(x, x) if i == ib else (0.0, 25.0) for i, x in enumerate(v0)] + [(0.0, 100.0)]
    res = minimize(lambda z: float(np.sum(resid(z) ** 2)), np.append(v0, 1.0), method="SLSQP", bounds=bounds,
                   constraints=[{"type": "eq", "fun": lambda z: np.sum(z[:-1]) - 100.0}],
                   options={"maxiter": 150, "ftol": 1e-8, "eps": 0.01})
    v = np.clip(res.x[:-1], 0, None); v = 100 * v / v.sum()
    m.R = float(res.x[-1]); E, T = m.run(np.append(v, 0.0), step=1.0)
    fit = {"R": m.R, "D86_resid_C": [round(T_at(E, T, p) - t, 2) for p, t in zip(FIT_PTS, tD86)],
           "RVP_resid_kPa": round(rvp_d323(v, c["rho"], c["MW"], c["psat"]) - tRVP, 2),
           "density_resid": round(density_15C(v, c["rho"]) + off - tRho, 2),
           "classes": {k: round(x, 2) for k, x in class_totals(v, c["classes"]).items()},
           "optimizer": res.message}
    return c, v, fit


def calibrate_without_kopp():
    mix, P, v = cr.setup("A")
    def score(R):
        r = cr.response(R, mix, P, v)
        return sum(((r[k] - TARGET[k]) / SPREAD[k]) ** 2 for k in TARGET)
    grid = [0.0, 0.2, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    s = [score(R) for R in grid]; i = int(np.argmin(s))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    return float(minimize_scalar(score, bounds=(lo, hi), method="bounded", options={"xatol": 0.01}).x)


def measured_curve(d86_F):
    return np.array(PTS_ALL, float), np.array([F(t) for t in d86_F])


def corrected(Es, Ts):
    Em, Tm = measured_curve(BASE2["D86_F"])
    grid = np.arange(Tm[0], Tm[-1] + DT / 2, DT)
    inc = lambda T: T + np.arange(len(T)) * 1e-9
    Emeas = np.interp(grid, Tm, Em); Eb0 = np.interp(grid, inc(Ts[0]), Es[0])
    return grid, [np.maximum.accumulate(Emeas + np.interp(grid, inc(T), E) - Eb0) for E, T in zip(Es, Ts)]


def props(grid, E):
    inc = lambda x: x + np.arange(len(x)) * 1e-9
    r = {f"T{p}": float(np.interp(p, inc(E), grid)) for p in (10, 30, 50, 70, 90)}
    r["E70"] = float(np.interp(70.0, grid, E)); r["E100"] = float(np.interp(100.0, grid, E))
    return r


def main():
    os.makedirs(OUT, exist_ok=True)
    c, v, fit = fit_base2()
    print(f"Base 2 surrogate: R = {fit['R']:.2f}, D86 residuals {fit['D86_resid_C']}, RVP {fit['RVP_resid_kPa']} kPa, "
          f"density {fit['density_resid']} kg/m3, classes {fit['classes']}")
    R_blend = calibrate_without_kopp()
    print(f"Reflux calibrated on reformate without Kopp data: R = {R_blend:.3f}")

    P = load_props(); mix = Mixture(c["names"] + ["mtbe"], c["classes"] + ["Ether"], P)
    m = D86Model(mix, P, model="nrtl", R=R_blend)
    Em, Tm = measured_curve(BASE2["D86_F"]); E5, T5 = measured_curve(FUEL5B["D86_F"])
    meas_b = {f"T{p}": float(np.interp(p, Em, Tm)) for p in (10, 30, 50, 70, 90)}; meas_b["E70"] = float(np.interp(70.0, Tm, Em)); meas_b["E100"] = float(np.interp(100.0, Tm, Em))
    meas_5 = {f"T{p}": float(np.interp(p, E5, T5)) for p in (10, 30, 50, 70, 90)}; meas_5["E70"] = float(np.interp(70.0, T5, E5)); meas_5["E100"] = float(np.interp(100.0, T5, E5))
    rvp_b = mix.rvp(np.append(v, 0.0), "nrtl")
    rows = [["Quantity", "Base 2 measured", "5B measured", "Measured change",
             "Predicted change (11.0 vol%)", "Predicted change (10.3 vol%)", "Error at 11.0 vol%"]]
    pred = {}
    for x in (MTBE_RECIPE, MTBE_ANALYSED):
        Eb, Tb = m.run(np.append(v, 0.0), step=STEP)
        Ex, Tx = m.run(np.append(v * (100 - x) / 100, x), step=STEP)
        grid, (cb, cx) = corrected([Eb, Ex], [Tb, Tx])
        pb, px = props(grid, cb), props(grid, cx)
        d = {k: px[k] - pb[k] for k in pb}
        d["RVP"] = mix.rvp(np.append(v * (100 - x) / 100, x), "nrtl") - rvp_b
        pred[x] = d
    meas_d = {k: meas_5[k] - meas_b[k] for k in meas_b}
    meas_d["RVP"] = (FUEL5B["RVP_psi"] - BASE2["RVP_psi"]) * PSI
    labels = {"RVP": "RVP (kPa)", "T10": "T10 (degC)", "T30": "T30 (degC)", "T50": "T50 (degC)",
              "T70": "T70 (degC)", "T90": "T90 (degC)", "E70": "E70 (vol%)", "E100": "E100 (vol%)"}
    base_vals = dict(meas_b, RVP=BASE2["RVP_psi"] * PSI); b5_vals = dict(meas_5, RVP=FUEL5B["RVP_psi"] * PSI)
    for k, lab in labels.items():
        rows.append([lab, round(base_vals[k], 1), round(b5_vals[k], 1), round(meas_d[k], 2),
                     round(pred[MTBE_RECIPE][k], 2), round(pred[MTBE_ANALYSED][k], 2),
                     round(pred[MTBE_RECIPE][k] - meas_d[k], 2)])
    with open(os.path.join(OUT, "kopp_delta_test.csv"), "w", newline="", encoding="utf-8-sig") as f:
        csv.writer(f).writerows(rows)
    json.dump({"base2_fit": fit, "R_blend_without_kopp": R_blend, "rows": rows,
               "note": "RVP measured by ASTM D5191 (Kopp Table 8); RVP resolution 0.1 psi = 0.7 kPa"},
              open(os.path.join(OUT, "kopp_delta_test.json"), "w"), indent=2, default=str)
    for r in rows: print(" | ".join(str(x) for x in r))


if __name__ == "__main__":
    main()

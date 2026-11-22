"""Part 3, step 2: D86 distillation of reformate-MTBE blends, 0-20 vol% MTBE.

Uses the surrogate and reflux ratio R fitted by surrogate_fit_rayleigh.py. For each liquid
model (nrtl = fitted MTBE parameters, ideal) and each reflux setting (R_blend from
calibrate_reflux.py when available = main case; R fitted on the base curve; R = 0) the D86 curve of every blend is computed and delta-corrected on the
evaporated-volume axis (article Eq. 2):

    E_blend(T) = E_base,measured(T) + [E_blend,model(T) - E_base,model(T)]

evaluated every 0.1 degC between the measured IBP and end point. Local reversals are removed by
a running maximum; the largest removal is reported. T10-T90 are read from the corrected curve,
E70, E100, E150 at fixed temperature, DI = 1.5 T10 + 3 T50 + T90 (degF).

The changes at 10 vol% MTBE are compared with published measurements on MTBE blends
(article Table 5).

Run from the project folder (a few seconds each):
    python src/blend_distillation.py A
    python src/blend_distillation.py B
Outputs in results/distillation/.
"""
import csv, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_models import Mixture, load_props, ROOT
from distillation import D86Model, DI_F

OUT = os.path.join(ROOT, "results", "distillation")
MTBE_PCT = np.arange(0, 21, 1.0)
STEP = 0.5          # vol% per distillation step
DT = 0.1            # degC grid for the delta correction

# Published measured changes (article Table 5): dT10 and dT50 at 10-11 vol% MTBE, E70 rate
MEASURED = {"dT10_10vol_C": [-1.0, -1.5, -0.6], "dT50_10vol_C": [-12.0, -4.5, -7.2],
            "dE70_per_vol_pct": [0.33, 0.55, 0.56, 0.40],
            "sources": "Hamid and Ali [8]; Mustafa [23]; Kopp et al. [22] (11 vol%)"}


def measured_base():
    b = json.load(open(os.path.join(ROOT, "data/processed/base_fuel.json")))["measured"]["D86"]
    E = list(b["evaporated_vol_pct"]) + [b["recovery_vol_pct"]]
    T = list(b["temperature_C"]) + [b["end_point_C"]]
    return np.array(E, float), np.array(T, float)


def strictly_increasing(T):
    return T + np.arange(len(T)) * 1e-9


def corrected_curves(Es, Ts):
    """Es, Ts: lists of model curves, index 0 = base. Returns corrected E(T) on the grid."""
    Em, Tm = measured_base()
    grid = np.arange(Tm[0], Tm[-1] + DT / 2, DT)
    E_meas = np.interp(grid, Tm, Em)
    E_base = np.interp(grid, strictly_increasing(Ts[0]), Es[0])
    out, clip = [], []
    for E, T in zip(Es, Ts):
        Eb = np.interp(grid, strictly_increasing(T), E)
        raw = E_meas + (Eb - E_base)
        mono = np.maximum.accumulate(raw)
        out.append(mono); clip.append(float(np.max(mono - raw)))
    return grid, out, clip


def props_from(grid, E):
    Ti = lambda p: float(np.interp(p, strictly_increasing(E), grid))
    Ei = lambda t: float(np.interp(t, grid, E))
    r = {f"T{p}": Ti(p) for p in (10, 30, 50, 70, 90)}
    r.update({f"E{t}": Ei(t) for t in (70, 100, 150)})
    r["DI_F"] = DI_F(r["T10"], r["T50"], r["T90"])
    return r


def run(which):
    fit = json.load(open(os.path.join(ROOT, "results/surrogate_fit", f"surrogate_{which}_rayleigh_fit.json")))
    sur = json.load(open(os.path.join(ROOT, "data/processed/surrogates.json")))["surrogate_" + which]
    names = [c["component"] for c in sur]; classes = [c["class"] for c in sur]
    v = np.array([fit["fitted_vol_pct"][n] for n in names])
    P = load_props()
    mix = Mixture(names + ["mtbe"], classes + ["Ether"], P)
    R_cal = fit["reflux_ratio_R"]
    settings = [("R_calibrated", R_cal), ("R_0", 0.0)]
    cal = os.path.join(OUT, f"reflux_calibration_{which}.json")
    if os.path.exists(cal):                      # Part 4: R fitted to measured MTBE blends
        settings.insert(0, ("R_blend", json.load(open(cal))["R_blend"]))
    os.makedirs(OUT, exist_ok=True)
    rows, summary = [], {"surrogate": which, "R_settings": dict(settings), "measured_reference": MEASURED, "cases": {}}
    for model in ("nrtl", "ideal"):
        for rname, R in settings:
            m = D86Model(mix, P, model=model, R=R)
            Es, Ts = [], []
            for x in MTBE_PCT:
                E, T = m.run(np.append(v * (100 - x) / v.sum(), x), step=STEP)
                Es.append(E); Ts.append(T)
            grid, Ec, clip = corrected_curves(Es, Ts)
            res = [props_from(grid, E) for E in Ec]
            case = f"{model}_{rname}"
            for x, r, c in zip(MTBE_PCT, res, clip):
                rows.append([case, x] + [round(r[k], 2) for k in ("T10", "T30", "T50", "T70", "T90", "E70", "E100", "E150", "DI_F")] + [round(c, 2)])
            b, t10 = res[0], res[10]
            first = lambda key, ok: next((float(x) for x, r in zip(MTBE_PCT, res) if ok(r[key])), None)
            summary["cases"][case] = {
                "dT10_10vol_C": round(t10["T10"] - b["T10"], 2), "dT50_10vol_C": round(t10["T50"] - b["T50"], 2),
                "dE70_per_vol_pct_0to10": round((t10["E70"] - b["E70"]) / 10, 3),
                "dE70_per_vol_pct_0to15": round((res[15]["E70"] - b["E70"]) / 15, 3),
                "E70_at_20vol": round(res[20]["E70"], 2),
                "first_vol_E70_ge_20": first("E70", lambda e: e >= 20),
                "first_vol_E100_ge_46": first("E100", lambda e: e >= 46),
                "first_vol_T10_le_70": first("T10", lambda t: t <= 70),
                "first_vol_DI_le_1250": first("DI_F", lambda d: d <= 1250),
                "max_running_max_removal_vol_pct": round(max(clip), 2)}
            s = summary["cases"][case]
            print(f"{which} {case:16s} dT10 {s['dT10_10vol_C']:+5.1f}  dT50 {s['dT50_10vol_C']:+5.1f}  "
                  f"dE70/vol {s['dE70_per_vol_pct_0to10']:.2f}  E70@20 {s['E70_at_20vol']:5.1f}  "
                  f"E70>=20 {s['first_vol_E70_ge_20']}  E100>=46 {s['first_vol_E100_ge_46']}  "
                  f"DI<=1250 {s['first_vol_DI_le_1250']}  clip {s['max_running_max_removal_vol_pct']}")
    with open(os.path.join(OUT, f"d86_{which}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "MTBE vol%", "T10", "T30", "T50", "T70", "T90", "E70", "E100", "E150", "DI (F)",
                    "running-max removal (vol%)"])
        w.writerows(rows)
    json.dump(summary, open(os.path.join(OUT, f"d86_{which}_summary.json"), "w"), indent=2)
    print("Measured (Table 5): dT10", MEASURED["dT10_10vol_C"], " dT50", MEASURED["dT50_10vol_C"],
          " dE70/vol", MEASURED["dE70_per_vol_pct"])


if __name__ == "__main__":
    run(sys.argv[1].upper() if len(sys.argv) > 1 else "A")

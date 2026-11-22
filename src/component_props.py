"""Build data/processed/component_props.json: CAS, MW, Tb, liquid density at 15 degC,
vapour pressures at 25, 35, 37.8 and 45 degC, carbon number and UNIFAC groups,
for every surrogate component and MTBE.

Run from the project folder:
    python src/component_props.py
"""
import json, os
from chemicals.identifiers import CAS_from_any, search_chemical
from thermo import Chemical

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T15, T378 = 288.15, 310.928  # 15 degC; 37.8 degC (100 degF)
PSAT_T_C = (25.0, 35.0, 37.8, 45.0)  # temperatures for TVP and RVP

def rho15(ch):
    """Liquid density at 15 degC (kg/m3): DIPPR (Perry 8th ed.) > VDI PPDS > mean of Rackett and
    Hankinson-Thomson COSTALD when no measured correlation is available."""
    vl = ch.VolumeLiquid
    for m in ("DIPPR_PERRY_8E", "VDI_PPDS"):
        if m in vl.all_methods:
            return ch.MW / vl.calculate(T15, m) / 1000, m
    est = [ch.MW / vl.calculate(T15, m) / 1000 for m in ("RACKETT", "HTCOSTALD") if m in vl.all_methods]
    return sum(est) / len(est), "mean(RACKETT,HTCOSTALD)"

def props(name):
    cas = CAS_from_any(name); md = search_chemical(cas)
    ch = Chemical(cas, T=T15)
    rho, src = rho15(ch)
    return {"CAS": cas, "MW": md.MW, "Tb_C": round(ch.Tb - 273.15, 2),
            "rho15_kg_m3": round(rho, 2), "rho15_source": src,
            "psat378_Pa": round(ch.VaporPressure(T378), 1),
            "psat_Pa": {str(t): round(ch.VaporPressure(t + 273.15), 1) for t in PSAT_T_C},
            "n_carbon": int(ch.atoms.get("C", 0)), "UNIFAC_groups": ch.UNIFAC_groups,
            "antoine": antoine(ch)}


def antoine(ch):
    """ln(P/Pa) = A - B/(T/K + C), fitted to the library vapour-pressure curve from 0.6 Tb to
    min(Tb + 150 K, 0.95 Tc); returns A, B, C and the maximum relative deviation."""
    import numpy as np
    from scipy.optimize import curve_fit
    T = np.linspace(0.6 * ch.Tb, min(ch.Tb + 150, 0.95 * ch.Tc), 60)
    lnP = np.log([ch.VaporPressure(t) for t in T])
    f = lambda T, A, B, C: A - B / (T + C)
    p, _ = curve_fit(f, T, lnP, p0=[21, 2800, -50], maxfev=20000)
    dev = float(np.max(abs(np.exp(f(T, *p) - lnP) - 1)))
    return {"A": float(p[0]), "B": float(p[1]), "C": float(p[2]), "T_range_K": [float(T[0]), float(T[-1])],
            "max_rel_dev": round(dev, 5)}

if __name__ == "__main__":
    s = json.load(open(os.path.join(ROOT, "data/processed/surrogates.json")))
    names = sorted({c["component"].lower() for k in s for c in s[k]} | {"mtbe"})
    out = {n: props(n) for n in names}
    json.dump(out, open(os.path.join(ROOT, "data/processed/component_props.json"), "w"), indent=2)
    print(len(out), "components written to data/processed/component_props.json")

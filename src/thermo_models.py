"""Part 2: liquid-phase non-ideality for MTBE blends.

Three liquid models (vapour ideal):
  ideal   Raoult's law (gamma = 1)
  unifac  original UNIFAC (as used for the HYSYS NRTL parameter estimates in the article)
  nrtl    NRTL with MTBE-hydrocarbon parameters regressed on measured binary VLE
          (data/raw/MTBE_binary_parameters_fitted.csv, by hydrocarbon class); hydrocarbon-hydrocarbon
          pairs are taken as ideal (tau = 0)

Functions return true vapour pressure (bubble pressure) and the ASTM D323-basis RVP
(dry vapour pressure at 37.8 degC, vapour/liquid volume = 4).
"""
import csv, json, os
import numpy as np
from thermo.unifac import UNIFAC, UFSG, UFIP

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = 8.314462618
ALPHA = 0.3


def load_props():
    P = json.load(open(os.path.join(ROOT, "data/processed/component_props.json")))
    for p in P.values():
        p["UNIFAC_groups"] = {int(k): v for k, v in p["UNIFAC_groups"].items()}
    return P


def load_nrtl_classes():
    """{class: (b12, b21)} in K, b12 = MTBE->HC, b21 = HC->MTBE."""
    out = {}
    with open(os.path.join(ROOT, "data/raw/MTBE_binary_parameters_fitted.csv")) as f:
        for r in csv.DictReader(f):
            out[r["class"]] = (float(r["NRTL b12 MTBE->HC (K)"]), float(r["NRTL b21 HC->MTBE (K)"]))
    return out


def nrtl_class(cls, n_carbon):
    if cls == "Olefin": return "Olefins"
    if cls == "Aromatic": return "Aromatics"
    if n_carbon == 4: return "C4"
    if n_carbon == 5: return "C5"
    return "C6+ paraffins and naphthenes"


class Mixture:
    """Components in a fixed order; MTBE must be one of them (index i_mtbe)."""

    def __init__(self, names, classes, props):
        self.names = [n.lower() for n in names]
        p = [props[n] for n in self.names]
        self.rho = np.array([x["rho15_kg_m3"] for x in p]); self.MW = np.array([x["MW"] for x in p])
        self.psat = {float(t): np.array([x["psat_Pa"][t] for x in p]) for t in p[0]["psat_Pa"]}
        self.groups = [x["UNIFAC_groups"] for x in p]
        self.i_mtbe = self.names.index("mtbe")
        par = load_nrtl_classes(); n = len(names)
        self.b = np.zeros((n, n))
        for j, (c, x) in enumerate(zip(classes, p)):
            if j == self.i_mtbe: continue
            b12, b21 = par[nrtl_class(c, x["n_carbon"])]
            self.b[self.i_mtbe, j], self.b[j, self.i_mtbe] = b12, b21

    def moles(self, vol_pct):
        v = np.asarray(vol_pct, float) / np.sum(vol_pct)
        return v * self.rho / (self.MW / 1000.0)        # mol per m3 of liquid

    def gammas(self, x, T, model):
        if model == "ideal":
            return np.ones_like(x)
        if model == "unifac":
            return np.array(UNIFAC.from_subgroups(T=T, xs=list(x), chemgroups=self.groups, version=0,
                                                  interaction_data=UFIP, subgroups=UFSG).gammas())
        if model == "nrtl":
            tau = self.b / T; G = np.exp(-ALPHA * tau)
            S = x @ G                                   # S_j = sum_k x_k G_kj
            C = x @ (tau * G)                           # C_j = sum_k x_k tau_kj G_kj
            lng = C / S + (G * (tau - C / S)) @ (x / S)  # standard multicomponent NRTL
            return np.exp(lng)
        raise ValueError(model)

    def tvp(self, vol_pct, T_C, model):
        """Bubble-point (true) vapour pressure, kPa."""
        n = self.moles(vol_pct); x = n / n.sum(); T = T_C + 273.15
        return float(np.sum(x * self.gammas(x, T, model) * self.psat[T_C]) / 1000)

    def rvp(self, vol_pct, model, vl_ratio=4.0, T_C=37.8):
        """Dry vapour pressure at 37.8 degC with V/L = 4 (successive substitution on gamma)."""
        n = self.moles(vol_pct); T = T_C + 273.15; ps = self.psat[T_C]
        x = n / n.sum()
        for _ in range(50):
            k = self.gammas(x, T, model) * ps * vl_ratio / (R * T)
            f = lambda L: np.sum(n / (L + k)) - 1.0
            lo, hi = 1e-9, n.sum()
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if f(mid) > 0: lo = mid
                else: hi = mid
            x_new = n / (0.5 * (lo + hi) + k); x_new /= x_new.sum()
            if np.max(abs(x_new - x)) < 1e-10: break
            x = x_new
        return float(np.sum(x * self.gammas(x, T, model) * ps) / 1000)

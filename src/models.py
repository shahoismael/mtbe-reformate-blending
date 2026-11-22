"""Property models used to fit and check surrogates (open-source, no HYSYS needed).

Units: temperature in degC unless stated, volume fractions in vol%.
Models
  density_15C   ideal volume mixing of COSTALD liquid densities at 15 degC
  rvp_d323      dry vapour pressure at 37.8 degC and vapour/liquid = 4 (ASTM D323 basis), Raoult's law
  tbp_curve     true boiling point curve from component normal boiling points (cumulative liquid volume)
  tbp_to_d86    API Technical Data Book (1994) procedure 3A1.1 (Daubert), inverted
"""
import numpy as np

R = 8.314462618  # J/mol/K
T15 = 288.15
T378 = 310.928

# ---------------------------------------------------------------- density
def density_15C(vol_pct, rho15):
    """Blend density (kg/m3) by ideal volume mixing."""
    v = np.asarray(vol_pct, float) / np.sum(vol_pct)
    return float(np.dot(v, rho15))


# ---------------------------------------------------------------- RVP
def rvp_d323(vol_pct, rho15, MW, psat378, vl_ratio=4.0):
    """Dry vapour pressure (kPa) at 37.8 degC with vapour volume = vl_ratio x liquid volume.

    Liquid moles n_i from the liquid volume; each component splits between liquid (Raoult)
    and an ideal-gas vapour space. Solves sum(x_i) = 1 for the remaining liquid moles L.
    """
    v = np.asarray(vol_pct, float) / np.sum(vol_pct)          # m3 per m3 liquid
    n = v * np.asarray(rho15) / (np.asarray(MW) / 1000.0)     # mol per m3 liquid
    P = np.asarray(psat378, float)                            # Pa
    k = P * vl_ratio / (R * T378)                             # mol/m3 per unit x_i
    f = lambda L: np.sum(n / (L + k)) - 1.0
    lo, hi = 1e-9, n.sum()
    for _ in range(200):                                      # bisection, monotone in L
        mid = 0.5 * (lo + hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    x = n / (0.5 * (lo + hi) + k)
    return float(np.dot(x / x.sum(), P) / 1000.0)


# ---------------------------------------------------------------- TBP and D86
def tbp_curve(vol_pct, nbp_C, points=(0, 10, 30, 50, 70, 90, 100), min_frac=0.05):
    """TBP temperatures (degC) at the given vol% evaporated.

    Components are ordered by boiling point; each is placed at the midpoint of its cumulative
    liquid-volume interval and the curve is interpolated linearly. 0 and 100 vol% are the boiling
    points of the lightest and heaviest components present above min_frac vol%.
    """
    v = np.asarray(vol_pct, float); t = np.asarray(nbp_C, float)
    keep = v > min_frac
    v, t = v[keep], t[keep]
    o = np.argsort(t); v, t = v[o] / v.sum() * 100, t[o]
    cum = np.cumsum(v); mid = cum - v / 2
    xs = np.concatenate([[0], mid, [100]]); ys = np.concatenate([[t[0]], t, [t[-1]]])
    return np.interp(points, xs, ys)


# API TDB 3A1.1 (1994): TBP50 = 0.87180 * D86_50 ** 1.0258 (degF); dTBP = A * dD86 ** B per segment
_SEG = {  # segment: (A, B)
    (90, 100): (0.11798, 1.6606), (70, 90): (3.0419, 0.75497), (50, 70): (2.5282, 0.82002),
    (30, 50): (3.0305, 0.80076), (10, 30): (4.9004, 0.71644), (0, 10): (7.4012, 0.60244)}

def _F(c): return c * 9 / 5 + 32
def _C(f): return (f - 32) * 5 / 9

def d86_to_tbp(d86_C):
    """D86 at 0,10,30,50,70,90,100 vol% -> TBP (degC)."""
    d = [_F(x) for x in d86_C]
    t50 = 0.87180 * d[3] ** 1.0258
    i = {0: 0, 10: 1, 30: 2, 50: 3, 70: 4, 90: 5, 100: 6}
    tbp = {50: t50}
    for lo, hi in [(30, 50), (10, 30), (0, 10)]:
        A, B = _SEG[(lo, hi)]; tbp[lo] = tbp[hi] - A * max(d[i[hi]] - d[i[lo]], 0) ** B
    for lo, hi in [(50, 70), (70, 90), (90, 100)]:
        A, B = _SEG[(lo, hi)]; tbp[hi] = tbp[lo] + A * max(d[i[hi]] - d[i[lo]], 0) ** B
    return [_C(tbp[k]) for k in (0, 10, 30, 50, 70, 90, 100)]

def tbp_to_d86(tbp_C):
    """TBP at 0,10,30,50,70,90,100 vol% -> D86 (degC), exact inverse of d86_to_tbp."""
    t = [_F(x) for x in tbp_C]
    i = {0: 0, 10: 1, 30: 2, 50: 3, 70: 4, 90: 5, 100: 6}
    d = {50: (t[3] / 0.87180) ** (1 / 1.0258)}
    for lo, hi in [(30, 50), (10, 30), (0, 10)]:
        A, B = _SEG[(lo, hi)]; d[lo] = d[hi] - (max(t[i[hi]] - t[i[lo]], 0) / A) ** (1 / B)
    for lo, hi in [(50, 70), (70, 90), (90, 100)]:
        A, B = _SEG[(lo, hi)]; d[hi] = d[lo] + (max(t[i[hi]] - t[i[lo]], 0) / A) ** (1 / B)
    return [_C(d[k]) for k in (0, 10, 30, 50, 70, 90, 100)]

def d86_from_composition(vol_pct, nbp_C):
    return tbp_to_d86(tbp_curve(vol_pct, nbp_C))


def class_totals(vol_pct, classes):
    v = np.asarray(vol_pct, float) / np.sum(vol_pct) * 100
    out = {}
    for c, x in zip(classes, v): out[c] = out.get(c, 0.0) + x
    return out

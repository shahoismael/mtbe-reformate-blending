"""E70 rise per vol% MTBE in measured blends (Table 7, note a).
Hamid and Ali (1995), Table VI; Mustafa (1989), Table 6 (T20 = T10 + reported T20-T10).
Mustafa B (1989) Sci Contrib Oil Gas 12(1):30-35, doi 10.29017/SCOG.12.1.1141. Values checked
against the PDF on 22 Nov 2026. Note: Table 6 blends are MGBT (MTBE + TEL, held at 87 RON;
TEL 0.85 -> 0 cc/USG), base = 70 vol% reformate + 30 vol% light naphtha (Indonesia)."""
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline
def e70(vols, temps, kind):
    v, t = np.asarray(vols, float), np.asarray(temps, float)
    if kind == "linear":
        return float(np.interp(70.0, t, v))
    f = (PchipInterpolator if kind == "pchip" else CubicSpline)(v, t)
    vv = np.linspace(v[0], v[-1], 20001)
    return float(np.interp(70.0, f(vv), vv))
HA_vol = [0, 10, 20, 50, 70, 90, 95]
HA = {0: [29, 45, 58, 108, 147, 191, 206], 5: [27, 45, 58, 104, 144, 189, 205], 10: [27, 44, 56, 96, 147, 188, 204]}
MU = {0: (64.5, 6.5, 90.5, 135.5, 161, 0.8), 10: (63.0, 4.5, 86.0, 135.0, 160, 1.0), 15: (62.0, 4.5, 82.0, 134.0, 160, 1.0)}
for kind in ("linear", "pchip", "cubic"):
    h = {x: e70(HA_vol, HA[x], kind) for x in HA}
    m = {}
    for x, (t10, d, t50, t90, ep, res) in MU.items():
        m[x] = e70([10, 20, 50, 90, 100 - res], [t10, t10 + d, t50, t90, ep], kind)
    print(f"{kind:7s} Hamid&Ali E70 {[round(h[x],1) for x in HA]} rate(10) {(h[10]-h[0])/10:.2f} | "
          f"Mustafa E70 {[round(m[x],1) for x in MU]} rates {(m[10]-m[0])/10:.2f}, {(m[15]-m[0])/15:.2f}")

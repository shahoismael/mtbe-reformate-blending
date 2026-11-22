"""Part 5: uncertainty of the octane-based MTBE bounds (regular, premium, super).

Linear-by-volume blending: ON(x) = (1 - x) ON_base + x bON_MTBE, x = MTBE volume fraction.
The lower bound of a grade is the smallest x that meets both its RON and MON minima.

Uncertain inputs (Monte Carlo, N draws, seed 1):
  base RON    measured 91.0; normal, sigma = reproducibility / 2.77 (ASTM D2699, R = 0.7 assumed
              here; change RON_REPRODUCIBILITY if the lab's value differs)
  base MON    not measured; uniform between MON_LO and MON_HI (article: 79-84; estimates 80.9
              from the Ghosh model and 81.9 from a RON-MON correlation)
  blending RON and MON of MTBE
              empirical: one published measurement in a reformate-rich base is drawn at random
              (article Table 5), then a value uniformly within its reported range

Outputs per grade: median, 5th-95th percentile, P(no MTBE needed), P(bound > oxygen limit),
and Spearman rank correlations of the super-grade bound with each input.

Run from the project folder (a few seconds):
    python src/octane_uncertainty.py
Outputs in results/octane/.
"""
import csv, json, os, sys
import numpy as np
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "octane")

N = 200_000
SEED = 1
RON_BASE = 91.0
RON_REPRODUCIBILITY = 0.7
MON_LO, MON_HI = 79.0, 84.0
O2_LIMIT_VOL = 20.8          # vol% MTBE at 3.7 wt% oxygen (article Section 2.5)

# (source, blending RON range, blending MON range) measured in reformate-rich bases (article Table 5)
BLENDING = [
    ("Ershov et al. [30], 91 wt% reformate", (117, 117), (104, 104)),
    ("Kopp et al. [40], AQIRP Base 2", (118, 118), (103, 103)),
    ("Hamid and Ali [10], reformate", (125, 130), None),
    ("Hamid and Ali [10], Gasoline A", (126, 133), (111, 113)),
    ("Hamid and Ali [10], Gasoline B", (122, 125), (103, 106)),
    ("Mustafa [41], 70% reformate", (122, 123), (121, 121)),
]
GRADES = {"regular": (89, 80), "premium": (92, 82), "super": (95, 85)}


def draw(rng):
    ron0 = rng.normal(RON_BASE, RON_REPRODUCIBILITY / 2.77, N)
    mon0 = rng.uniform(MON_LO, MON_HI, N)
    i_r = rng.integers(0, len(BLENDING), N)
    lo_r = np.array([b[1][0] for b in BLENDING])[i_r]; hi_r = np.array([b[1][1] for b in BLENDING])[i_r]
    bron = rng.uniform(lo_r, hi_r + 1e-12)
    mon_src = [b for b in BLENDING if b[2] is not None]
    i_m = rng.integers(0, len(mon_src), N)
    lo_m = np.array([b[2][0] for b in mon_src])[i_m]; hi_m = np.array([b[2][1] for b in mon_src])[i_m]
    bmon = rng.uniform(lo_m, hi_m + 1e-12)
    return ron0, mon0, bron, bmon


def need(target, base, blend):
    """vol% MTBE to raise base to target (0 if already met)."""
    return np.maximum(0.0, 100.0 * (target - base) / (blend - base))


def main():
    rng = np.random.default_rng(SEED)
    ron0, mon0, bron, bmon = draw(rng)
    os.makedirs(OUT, exist_ok=True)
    summary, rows = {}, [["grade", "median vol%", "p5", "p95", "P(no MTBE needed)", f"P(> {O2_LIMIT_VOL} vol%)",
                          "reference case vol% (MON 82, bRON 119, bMON 103)"]]
    bounds = {}
    for g, (rmin, mmin) in GRADES.items():
        b = np.maximum(need(rmin, ron0, bron), need(mmin, mon0, bmon)); bounds[g] = b
        ref = max(need(rmin, RON_BASE, 119.0), need(mmin, 82.0, 103.0))
        p5, p50, p95 = np.percentile(b, [5, 50, 95])
        summary[g] = {"median": round(p50, 2), "p5": round(p5, 2), "p95": round(p95, 2),
                      "P_no_MTBE": round(float(np.mean(b == 0)), 3),
                      "P_above_O2_limit": round(float(np.mean(b > O2_LIMIT_VOL)), 3),
                      "reference_case": round(float(ref), 2)}
        rows.append([g, round(p50, 1), round(p5, 1), round(p95, 1), summary[g]["P_no_MTBE"],
                     summary[g]["P_above_O2_limit"], round(float(ref), 1)])
        print(f"{g:8s} median {p50:5.1f}  5-95% {p5:5.1f}-{p95:5.1f} vol%  P(none) {summary[g]['P_no_MTBE']:.2f}  "
              f"P(>{O2_LIMIT_VOL}) {summary[g]['P_above_O2_limit']:.3f}  reference {ref:4.1f}")
    sens = {}
    for g in GRADES:
        sens[g] = {k: round(float(spearmanr(x, bounds[g])[0]), 3) if np.std(bounds[g]) > 0 else None
                   for k, x in (("base RON", ron0), ("base MON", mon0), ("blending RON", bron), ("blending MON", bmon))}
    print("Spearman correlation with the super bound:", sens["super"])
    with open(os.path.join(OUT, "octane_bounds.csv"), "w", newline="") as f:
        csv.writer(f).writerows(rows)
    json.dump({"N": N, "seed": SEED, "inputs": {"RON_base": RON_BASE, "RON_reproducibility": RON_REPRODUCIBILITY,
               "MON_base_range": [MON_LO, MON_HI], "blending_sources": BLENDING}, "bounds": summary,
               "spearman": sens}, open(os.path.join(OUT, "octane_summary.json"), "w"), indent=2)


if __name__ == "__main__":
    main()

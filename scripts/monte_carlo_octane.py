"""Monte Carlo propagation of octane inputs (Section 3.3).
Linear-by-volume blending; uniform draws of MTBE blending RON 117-125,
blending MON 103-111 and base MON 79-84; base RON 91 (measured).
Bound of a grade = higher of its RON and MON bounds (vol% MTBE)."""
import numpy as np
rng = np.random.default_rng(1)
N = 200_000
bRON = rng.uniform(117, 125, N)
bMON = rng.uniform(103, 111, N)
MON0 = rng.uniform(79, 84, N)
RON0 = 91.0

def need(target, base, blend):
    return np.maximum(0.0, 100.0 * (target - base) / (blend - base))

grades = {"regular": (89, 80), "premium": (92, 82), "super": (95, 85)}
for g, (ron_min, mon_min) in grades.items():
    b = np.maximum(need(ron_min, RON0, bRON), need(mon_min, MON0, bMON))
    p5, p50, p95 = np.percentile(b, [5, 50, 95])
    print(f"{g:8s} median {p50:5.1f}  5th-95th {p5:5.1f}-{p95:5.1f} vol%  "
          f"P(no MTBE needed) {np.mean(b == 0):.2f}  P(>20.8 vol%) {np.mean(b > 20.8):.3f}")

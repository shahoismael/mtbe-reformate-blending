"""Route 2 (Oil Manager assay + D86 utility): delta correction on the evaporated-volume scale
(Eq. 2) at 0.01 degC intervals with running maximum; T10/T50/T90 by inverse interpolation.
Reproduces Supplementary Table S6. Input: Route2_Assay_OilManager_D86.xlsx, sheet D86_curves."""
import sys, numpy as np, openpyxl
src = sys.argv[1] if len(sys.argv) > 1 else "Route2_Assay_OilManager_D86.xlsx"
fm = np.array([0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90])                     # measured D86, vol%
Tm = np.array([43.2, 65.5, 72.3, 83.4, 93.6, 103.3, 114.4, 126, 138.2, 151, 167.7])  # degC
rows = list(openpyxl.load_workbook(src, data_only=True)["D86_curves"].iter_rows(min_row=6, max_row=36, values_only=True))
pe = np.array([r[0] for r in rows], float)
cur = {x: np.array([r[1 + x] for r in rows], float) for x in range(21)}
Tg = np.arange(43.2, 192.0, 0.01)
F = lambda t: t * 1.8 + 32
def route2(x):
    E = np.interp(Tg, Tm, fm, right=90) + np.interp(Tg, cur[x], pe) - np.interp(Tg, cur[0], pe)
    E = np.maximum.accumulate(E)
    inv = lambda f: float(np.interp(f, E, Tg)) if f < E[-1] else float(Tg[-1])
    T10, T50, T90 = inv(10), inv(50), inv(90 - 1e-9)
    return dict(E70=np.interp(70, Tg, E), E100=np.interp(100, Tg, E), T10=T10, T50=T50, T90=T90,
                DI=1.5 * F(T10) + 3 * F(T50) + F(T90))
res = {x: route2(x) for x in range(21)}
if __name__ == "__main__":
    print("MTBE  E70   E100  T10   T50   T90   DI")
    for x, d in res.items():
        print(f"{x:4d} {d['E70']:5.1f} {d['E100']:5.1f} {d['T10']:5.1f} {d['T50']:5.1f} {d['T90']:5.1f} {d['DI']:5.0f}")

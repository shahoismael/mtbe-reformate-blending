"""Batch (Rayleigh) route: delta correction on the evaporated-volume axis (Eq. 2), Tables S9-S11 and S15.
Run from the package folder on a HYSYS case-study export, as supplied:
    python scripts/route5_rayleigh_delta.py Route5_SurrogateA_Batch_TVP_NRTL_UNIFAC.csv
Stage columns '<T>V - ...' and 'BLENDED - Std Ideal Liq Vol Flow' are read by name; TVP/RVP columns are ignored.
-32767 marks a non-converged stage. By default such states are skipped and responses interpolated linearly.
Surrogate A with Peng-Robinson kij = 0 (Table S11 note): missing stage values are interpolated between neighbouring
states (linear extrapolation from the last two states at 20 vol%) and the 0 vol% state is completed from the fitted
Peng-Robinson case, identical at 0 vol% MTBE:
    python scripts/route5_rayleigh_delta.py Route6_SurrogateA_Batch_TVP_PR_kij0.csv --fill --base Route7b_SurrogateA_Batch_TVP_PR_fitted.csv
Writes <input>_delta.json next to the input file."""
import csv, json, sys, numpy as np

def load(fn):
    rows = list(csv.reader(open(fn))); h = rows[0]
    st = [i for i, x in enumerate(h) if x[:1].isdigit() and x.split(' ')[0].endswith('V')]
    bi = h.index('BLENDED - Std Ideal Liq Vol Flow')
    X, V = [], []
    for r in rows[3:]:
        if not r or not r[0]: continue
        v = np.array([float(r[i]) for i in st])
        X.append(float(r[1])); V.append(np.where(v < -1, np.nan, v / float(r[bi]) * 100))
    return np.array([float(h[i].split('V')[0]) for i in st]), np.array(X), np.array(V)

args = sys.argv[1:]; fn = args[0]; fill = '--fill' in args
Ts, X, V = load(fn)
if '--base' in args:
    _, Xb, Vb = load(args[args.index('--base') + 1]); assert Xb[0] == 0 and not np.isnan(Vb[0]).any()
    V[0] = np.where(np.isnan(V[0]), Vb[0], V[0])
bad = np.isnan(V).any(axis=1)
print('states', len(X), 'with non-converged stages', [float(x) for x in X[bad]], '-> filled' if fill else '-> skipped')
if fill:
    for j in range(V.shape[1]):
        k = ~np.isnan(V[:, j]); xs, ys = X[k], V[k, j]
        for i in np.where(~k)[0]:
            V[i, j] = np.polyval(np.polyfit(xs[-2:], ys[-2:], 1), X[i]) if X[i] > xs.max() else np.interp(X[i], xs, ys)
else:
    X, V = X[~bad], V[~bad]
assert X[0] == 0, 'the 0 vol% MTBE state must be complete (use --base)'
fm = np.array([0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90])
Tm = np.array([43.2, 65.5, 72.3, 83.4, 93.6, 103.3, 114.4, 126, 138.2, 151, 167.7])
F = lambda t: t * 1.8 + 32
Cum = np.cumsum(V, axis=1); Tg = np.arange(43.2, 192, 0.01); Em = np.interp(Tg, Tm, fm, right=np.nan)
res = []
for x, c in zip(X, Cum):
    E = np.maximum.accumulate(np.nan_to_num(Em + np.interp(Tg, Ts, c) - np.interp(Tg, Ts, Cum[0]), nan=0))
    g = lambda t: float(np.interp(t, Tg, E)); inv = lambda f: float(Tg[np.argmax(E >= f)])
    T10, T50, T90 = inv(10), inv(50), inv(90)
    if x == 0: T90 = 167.7
    res.append(dict(x=float(x), E70raw=float(np.interp(70, Ts, c)), E70=g(70), E100=g(100), E150=g(150),
                    T10=T10, T50=T50, T90=T90, DI=1.5 * F(T10) + 3 * F(T50) + F(T90)))
xs = np.array([r['x'] for r in res])
I = {k: np.interp(np.arange(21), xs, [r[k] for r in res]) for k in res[0] if k != 'x'}
def first(k, cond):
    for xx in np.arange(0, 20.01, 0.01):
        if cond(np.interp(xx, xs, [r[k] for r in res])): return round(xx, 2)
print('first met (vol% MTBE): E70>=20', first('E70', lambda v: v >= 20), 'E100>=46', first('E100', lambda v: v >= 46),
      'DI<=1250', first('DI', lambda v: v <= 1250), 'T10<=70', first('T10', lambda v: v <= 70))
for k in range(21): print(k, {q: round(float(I[q][k]), 2) for q in I})
s = np.polyfit(xs[-2:], [r['E70'] for r in res[-2:]], 1)
print('E70 at 20.8 vol% (linear extrapolation from the last two states)', round(float(np.polyval(s, 20.8)), 2))
json.dump({'states': res, 'integer': {q: I[q].tolist() for q in I}}, open(fn.rsplit('.', 1)[0] + '_delta.json', 'w'))

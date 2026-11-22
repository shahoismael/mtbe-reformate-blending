"""Octane and volatility cross-checks reported in Sections 3.2, 3.4 and 3.5 of the article.

1. Ghosh, Hickey and Jaffe (2006), Ind. Eng. Chem. Res. 45, 337-345: Eq. 11 with the
   class-average a-hat parameters (Table 3) and the lump pure-component RON/MON (Table 2),
   applied to surrogate A (Table S1) and its MTBE blends.
2. Stratiev et al. (2026), ChemEngineering 10, 128: RON-MON correlation (Fig. 7) and
   RON-aromatics correlations (Fig. 6) for reformates.
3. Kopp et al. (1993), SAE 930143: AQIRP Base 2 (Table 11) against fuel 5B = 89 vol% Base 2
   + 11 vol% MTBE (recipe, Table 12; Table 13 analysis gives 10.3 lv%); blending values
   back-calculated by volume at 11 vol%. Values checked against the PDF on 22 Nov 2026.
   Kopp V, Bones C, Doerr D, Ho S, et al. (1993) SAE Technical Paper 930143, doi 10.4271/930143.

Run: python scripts/octane_crosschecks.py   (Python 3, no extra packages)
"""

# ---------------------------------------------------------------- 1. Ghosh model
A_HAT = {  # class: (a-hat RON, a-hat MON), Ghosh et al. Table 3
    'n': (2.0559, 0.3092), 'i': (2.0204, 0.4278), 'N': (1.6870, 0.2821),
    'A': (3.3984, 0.4773), 'O': (8.9390, 10.0), 'X': (3.9743, 2.0727)}
K_PN_A, K_PN_B, K_PO_A, K_PO_B = 0.2, 2.4, 0.4, 3.6

# surrogate A (vol%), class, and RON/MON of the matching Ghosh lump
SURROGATE_A = [
    ('n-butane', 'n', 5.0, 94, 89.6), ('i-pentane', 'i', 4.6, 92, 90.3),
    ('n-pentane', 'n', 0.9, 62, 62.6), ('2-methyl-2-butene', 'O', 8.0, 103, 82),
    ('2-methylpentane', 'i', 7.6, 76, 73.9), ('1-hexene', 'O', 0.3, 90, 80),
    ('n-hexane', 'n', 1.1, 24.8, 26), ('benzene', 'A', 1.5, 102.7, 105),
    ('2-methylhexane', 'i', 5.3, 52, 52), ('n-heptane', 'n', 0.3, 0, 0),
    ('2,2,4-trimethylpentane', 'i', 4.8, 105, 98.8), ('methylcyclohexane', 'N', 1.7, 82, 77),
    ('2,5-dimethylhexane', 'i', 4.5, 69, 74.5), ('toluene', 'A', 6.3, 118, 103.5),
    ('1-octene', 'O', 0.3, 90, 77), ('n-octane', 'n', 0.8, -15, -20),
    ('ethylbenzene', 'A', 2.0, 112, 105), ('m-xylene', 'A', 7.6, 112, 105),
    ('2-methyloctane', 'i', 12.5, 15, 22.3), ('n-nonane', 'n', 4.9, -20, -20),
    ('1,2,4-trimethylbenzene', 'A', 11.7, 110, 101), ('1,2,4,5-tetramethylbenzene', 'A', 8.3, 109, 98)]
MTBE = ('MTBE', 'X', 115.2, 97.2)


def ghosh_octane(comp, k):
    """Eq. 11 of Ghosh et al.; k = 0 for RON, 1 for MON; comp = (name, class, vol, RON, MON)."""
    tot = sum(c[2] for c in comp)
    V = [(c[1], c[2] / tot, c[3 + k]) for c in comp]
    VN = sum(v for cl, v, _ in V if cl == 'N')
    VO = sum(v for cl, v, _ in V if cl == 'O')
    I = (K_PN_A * VN + K_PO_A * VO) / (1 + K_PN_B * VN + K_PO_B * VO)
    a = lambda cl: A_HAT[cl][k]
    par = [(cl, v, on) for cl, v, on in V if cl in 'ni']
    num = sum(v * a(cl) * on for cl, v, on in V) + I * sum(v * a(cl) * on for cl, v, on in par)
    den = sum(v * a(cl) for cl, v, _ in V) + I * (sum(v * a(cl) for cl, v, _ in par) - sum(v for _, v, _ in par))
    return num / den


def blend(x):
    return [(n, c, v * (100 - x) / 100, r, m) for n, c, v, r, m in SURROGATE_A] + \
           [(MTBE[0], MTBE[1], x, MTBE[2], MTBE[3])]


print('1. Ghosh et al. (2006) model on surrogate A')
print('   MTBE vol%   RON    MON')
for x in (0, 5, 10, 15, 20):
    b = blend(x)
    print(f'   {x:8d}  {ghosh_octane(b, 0):6.1f} {ghosh_octane(b, 1):6.1f}')
r0, r20 = ghosh_octane(blend(0), 0), ghosh_octane(blend(20), 0)
m0, m20 = ghosh_octane(blend(0), 1), ghosh_octane(blend(20), 1)
print(f'   slope RON {(r20 - r0) / 20:.2f}, MON {(m20 - m0) / 20:.2f} per vol% MTBE (0-20 vol%)')
print('   measured base RON 91.0')

# ---------------------------------------------------------------- 2. Stratiev correlations
RON_BASE = 91.0
print('\n2. Stratiev et al. (2026) correlations')
print(f'   MON = 0.7606 RON + 12.702 (fitted RON 95-101): MON at RON {RON_BASE} = {0.7606 * RON_BASE + 12.702:.1f}, '
      f'at RON 88.5 = {0.7606 * 88.5 + 12.702:.1f}')
AROM_VOL, RHO_AROM, RHO_REF = 34.65, 0.870, 0.7628  # vol%, assumed aromatics density, measured density
arom_wt = AROM_VOL * RHO_AROM / RHO_REF
print(f'   aromatics {AROM_VOL} vol% -> {arom_wt:.1f} wt% (aromatics density {RHO_AROM} g/mL assumed)')
for label, s, c in (('this work (Stratiev)', 0.458, 65.869), ('Moljord et al.', 0.4829, 64.926),
                    ('Bishara et al.', 0.3672, 74.055)):
    ron = s * arom_wt + c
    print(f'   RON from {label:22s} {ron:5.1f}  (measured {RON_BASE}, difference {RON_BASE - ron:+.1f})')

# ---------------------------------------------------------------- 3. Kopp et al. pair
PSI = 6.89476
F_TO_C = lambda f: (f - 32) / 1.8
base = dict(RVP=8.3, T10=125, T30=147, T40=165, T50=189, T90=270, RON=94.7, MON=84.0)
b5 = dict(RVP=8.6, T10=124, T40=157, T50=176, T90=267, RON=97.3, MON=86.1)
x = 0.11
blend_val = lambda p, q: (q - (1 - x) * p) / x
e70_base = 30 + 10 * (158 - base['T30']) / (base['T40'] - base['T30'])   # 70 C = 158 F
e70_5b = 40 + 10 * (158 - b5['T40']) / (b5['T50'] - b5['T40'])
print('\n3. Kopp et al. (1993): Base 2 against 5B (11 vol% MTBE)')
print(f'   blending RON {blend_val(base["RON"], b5["RON"]):.1f}, MON {blend_val(base["MON"], b5["MON"]):.1f}, '
      f'RVP {blend_val(base["RVP"], b5["RVP"]) * PSI:.0f} kPa (+/- {0.1 * PSI / x * 2 ** 0.5:.0f} kPa from 0.1 psi rounding)')
for t in ('T10', 'T50', 'T90'):
    print(f'   d{t} = {(b5[t] - base[t]) / 1.8:+.1f} C')
print(f'   E70: {e70_base:.1f} -> {e70_5b:.1f} vol%, rate {(e70_5b - e70_base) / 11:.2f} vol% per vol% MTBE')
print(f'   base RVP {base["RVP"] * PSI:.1f} kPa, T10 {F_TO_C(base["T10"]):.1f} C, T50 {F_TO_C(base["T50"]):.1f} C')

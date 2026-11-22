"""Octane of reformate-MTBE blends from the PIONA/OS models of Badra et al. (2022),
Fuel 318, 123647, Eq. 6 and Table 7, applied as a labelled extrapolation (Section 2.5).
ON = ON_LBV + phi_MTBE * [beta*phi_iso + gamma*phi_aro + eps*(6 - OS)],
ON_LBV with MTBE RON 118.3 and MON 101.3 (ASTM values quoted by Badra et al.)."""
iso, aro = 0.4049, 0.3465          # Al-Daura reformate, ASTM D6729 (vol fraction)
RON0, RON_MTBE, MON_MTBE = 91.0, 118.3, 101.3
models = {  # (beta, gamma, eps) for RON, then for MON
    "LBV":                ((0, 0, 0),           (0, 0, 0)),
    "LBV + iso + aro":    ((24.2, -11.5, 0),    (31.3, -11.3, 0)),
    "LBV + iso + OS":     ((8.4, 0, 2.653),     (11.3, 0, 2.996)),
    "LBV + iso + aro + OS": ((13.2, -5.3, 1.980), (11.4, -2.0, 2.674)),
}
def need(target, base, slope):
    return max(0.0, 100.0 * (target - base) / slope)
for MON0 in (82, 80, 79):
    OS = RON0 - MON0
    print(f"\nBase MON {MON0} (assumed), OS {OS}")
    for name, (r, m) in models.items():
        sR = (RON_MTBE - RON0) + r[0]*iso + r[1]*aro + r[2]*(6 - OS)
        sM = (MON_MTBE - MON0) + m[0]*iso + m[1]*aro + m[2]*(6 - OS)
        prem = max(need(92, RON0, sR), need(82, MON0, sM))
        sup = max(need(95, RON0, sR), need(85, MON0, sM))
        reg = max(need(89, RON0, sR), need(80, MON0, sM))
        print(f"  {name:22s} bRON {RON0+sR:6.1f} bMON {MON0+sM:6.1f} | "
              f"regular {reg:4.1f}  premium {prem:4.1f}  super {sup:4.1f} vol%")

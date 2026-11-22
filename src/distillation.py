"""Part 3: ASTM D86 distillation model.

The D86 test is modelled as a batch (Rayleigh) distillation at 101.325 kPa in which the flask
vapour passes one equilibrium stage (the flask neck) with internal reflux ratio R before it is
condensed; the thermometer reads the dew-point temperature of the vapour leaving the stage.
R = 0 gives simple Rayleigh distillation. R is calibrated on the measured D86 curve of the
base fuel (d86_calibrate.py), then used unchanged for the MTBE blends.

Liquid model: ideal or NRTL (fitted MTBE-hydrocarbon parameters, thermo_models.py); vapour ideal.
Vapour pressures: Antoine fits stored in data/processed/component_props.json.
"""
import numpy as np
from thermo_models import Mixture

P_ATM = 101325.0


class D86Model:
    def __init__(self, mix: Mixture, props, model="nrtl", R=0.0):
        self.mix, self.model, self.R = mix, model, R
        a = [props[n]["antoine"] for n in mix.names]
        self.A = np.array([x["A"] for x in a]); self.B = np.array([x["B"] for x in a]); self.C = np.array([x["C"] for x in a])
        self.vmol = mix.MW / 1000.0 / mix.rho            # m3 liquid per mol

    def psat(self, T):
        return np.exp(self.A - self.B / (T + self.C))

    def gam(self, x, T):
        return self.mix.gammas(x, T, self.model) if self.model != "ideal" else np.ones_like(x)

    def bubble(self, x, T0=350.0):
        """Bubble temperature (K) and vapour composition of liquid x."""
        T = T0
        for _ in range(60):
            k = self.gam(x, T) * self.psat(T) / P_ATM
            s = np.dot(x, k)
            # d ln s / dT ~ sum x k B/(T+C)^2 / s
            dlns = np.dot(x * k, self.B / (T + self.C) ** 2) / s
            dT = np.log(s) / dlns
            T -= np.clip(dT, -40, 40)
            if abs(dT) < 1e-5: break
        k = self.gam(x, T) * self.psat(T) / P_ATM
        y = x * k
        return T, y / y.sum()

    def dew(self, y, T0=350.0, x0=None):
        """Dew temperature (K) and liquid composition of vapour y."""
        T = T0; x = y.copy() if x0 is None else x0
        for _ in range(40):
            for _ in range(40):
                k = self.gam(x, T) * self.psat(T) / P_ATM
                s = np.sum(y / k)
                dlns = -np.dot(y / k, self.B / (T + self.C) ** 2) / s
                dT = np.log(s) / dlns
                T -= np.clip(dT, -40, 40)
                if abs(dT) < 1e-5: break
            k = self.gam(x, T) * self.psat(T) / P_ATM
            xn = y / k; xn /= xn.sum()
            if np.max(abs(xn - x)) < 1e-8:
                x = xn; break
            x = xn
        return T, x

    def top(self, xB, TB, state):
        """Distillate composition and top (thermometer) temperature for pot liquid xB.

        Stage 1 (flask neck) receives the pot vapour yB and returns liquid x1 at reflux ratio R:
            (R + 1) yB = R x1 + xD,   xD = K(T1) x1   ->   x1_i = (R + 1) yB_i / (R + K_i)
        T1 follows from sum(x1) = 1 and is the thermometer reading.
        """
        TB, yB = self.bubble(xB, TB)
        if self.R <= 0:
            return yB, TB, TB
        R = self.R; T = state.get("T1", TB); x1 = state.get("x1", yB)
        for _ in range(30):
            g = self.gam(x1, T)
            for _ in range(60):
                K = g * self.psat(T) / P_ATM
                s = np.sum((R + 1) * yB / (R + K))
                ds = -np.sum((R + 1) * yB * K * self.B / (T + self.C) ** 2 / (R + K) ** 2)
                dT = (s - 1) / ds
                T -= np.clip(dT, -40, 40)
                if abs(dT) < 1e-6: break
            K = g * self.psat(T) / P_ATM
            xn = (R + 1) * yB / (R + K); xn /= xn.sum()
            if np.max(abs(xn - x1)) < 1e-9:
                x1 = xn; break
            x1 = xn
        xD = K * x1; xD /= xD.sum()
        state.update(T1=T, x1=x1)
        return xD, T, TB

    def run(self, vol_pct, step=0.5, end=95.0):
        """Returns arrays: evaporated vol% and thermometer temperature (degC)."""
        v = np.asarray(vol_pct, float) / np.sum(vol_pct)
        n = v / self.vmol                                  # mol per m3 of initial liquid
        V0 = 1.0; evap = 0.0; TB = 330.0; state = {}
        E, Tt = [], []
        while evap < end:
            x = n / n.sum()
            xD, T1, TB = self.top(x, TB, state)
            E.append(evap); Tt.append(T1 - 273.15)
            dV = step / 100 * V0
            D = dV / np.dot(xD, self.vmol)
            take = np.minimum(D * xD, n)
            n = n - take
            evap += 100 * np.dot(take, self.vmol) / V0
            if n.sum() <= 1e-9: break
        E, Tt = np.array(E), np.array(Tt)
        return E, np.maximum.accumulate(Tt)


# ---------------------------------------------------------------- curve utilities
def T_at(E, T, pct):
    return float(np.interp(pct, E, T))

def E_at(E, T, temp):
    return float(np.interp(temp, T, E))

def summary(E, T):
    out = {f"T{p}": T_at(E, T, p) for p in (0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)}
    out.update({f"E{t}": E_at(E, T, t) for t in (70, 100, 150)})
    return out

def DI_F(T10, T50, T90):
    f = lambda c: c * 9 / 5 + 32
    return 1.5 * f(T10) + 3.0 * f(T50) + f(T90)

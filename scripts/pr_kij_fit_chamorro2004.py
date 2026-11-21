"""Peng-Robinson kij for MTBE + 2,2,4-trimethylpentane and MTBE + toluene (Table S14, C6+ paraffins/naphthenes and
aromatics), regressed on the bubble pressures of Chamorro et al. (2004) at 313.15 K (every third point).
Requires numpy, scipy, thermo 0.6.1. Olefin kij (0.010) is from Bennett et al. (1993); C4 and C5 kij = 0."""
import numpy as np, warnings; warnings.filterwarnings('ignore')
from thermo import ChemicalConstantsPackage, PRMIX, CEOSGas, CEOSLiquid, FlashVL
from scipy.optimize import minimize_scalar
import numpy as np
ch_iso=np.array([[0.0546,15.929],[0.1013,18.467],[0.1509,21.118],[0.1967,23.478],[0.2463,25.974],[0.2996,28.597],[0.3467,30.874],[0.3980,33.326],[0.3981,33.375],[0.4483,35.695],[0.4484,35.672],[0.4984,37.974],[0.4990,37.979],[0.5477,40.187],[0.5493,40.248],[0.5974,42.385],[0.5997,42.476],[0.6471,44.573],[0.6978,46.774],[0.7461,48.890],[0.7958,51.027],[0.8486,53.300],[0.8961,55.350],[0.9467,57.566]])
ch_tol=np.array([[0.0599,11.430],[0.1017,13.883],[0.1507,16.658],[0.1998,19.385],[0.2473,21.956],[0.2991,24.702],[0.3507,27.376],[0.3604,27.855],[0.4006,29.936],[0.4085,30.320],[0.4485,32.357],[0.4580,32.823],[0.4992,34.903],[0.5080,35.322],[0.5494,37.398],[0.5592,37.874],[0.5999,39.884],[0.6105,40.400],[0.6633,42.999],[0.7177,45.688],[0.7748,48.497],[0.8310,51.307],[0.8926,54.396],[0.9716,58.431]])
def flasher(names, kij):
    c, p = ChemicalConstantsPackage.from_IDs(names)
    ek = dict(Tcs=c.Tcs, Pcs=c.Pcs, omegas=c.omegas, kijs=[[0, kij], [kij, 0]])
    return FlashVL(c, p, liquid=CEOSLiquid(PRMIX, ek, HeatCapacityGases=p.HeatCapacityGases),
                   gas=CEOSGas(PRMIX, ek, HeatCapacityGases=p.HeatCapacityGases))
for hc, d, T in (('2,2,4-trimethylpentane', ch_iso[::3], 313.15), ('toluene', ch_tol[::3], 313.15)):
    names = ['methyl tert-butyl ether', hc]
    def err(k):
        F = flasher(names, k)
        return np.mean(np.array([F.flash(T=T, VF=0, zs=[x, 1 - x]).P / 1000 / P - 1 for x, P in d]) ** 2)
    r = minimize_scalar(err, bounds=(-0.05, 0.08), method='bounded')
    print(f"{hc}: kij = 0 rms {np.sqrt(err(0))*100:.2f}%  |  best kij {r.x:.4f} rms {np.sqrt(r.fun)*100:.2f}%")

import numpy as np,warnings,sys,json;warnings.filterwarnings('ignore')
from thermo import Chemical
from thermo.unifac import UNIFAC,UFSG,UFIP
from scipy.optimize import brentq
A={'n-butane':5.0,'isopentane':4.6,'n-pentane':0.9,'2-methyl-2-butene':8.0,'2-methylpentane':7.6,'1-hexene':0.3,'n-hexane':1.1,'benzene':1.5,'2-methylhexane':5.3,'n-heptane':0.3,'2,2,4-trimethylpentane':4.8,'methylcyclohexane':1.7,'2,5-dimethylhexane':4.5,'toluene':6.3,'1-octene':0.3,'n-octane':0.8,'ethylbenzene':2.0,'m-xylene':7.6,'2-methyloctane':12.5,'n-nonane':4.9,'1,2,4-trimethylbenzene':11.7,'1,2,4,5-tetramethylbenzene':8.3}
M='methyl tert-butyl ether'; names=list(A)+[M]
C={k:Chemical(k,T=288.15,P=1e6) for k in names}
vm=np.array([C[k].MW/C[k].rhol for k in names])  # L/mol-ish std liquid molar volume
grp=[C[k].UNIFAC_groups for k in names]
Pc=101.325
def psat(T):
    out=[]
    for k in names:
        try: v=C[k].VaporPressure(T)
        except Exception: v=None
        if v is None or not np.isfinite(v): v=C[k].VaporPressure.extrapolate(T) if hasattr(C[k].VaporPressure,'extrapolate') else 1e8
        if v is None or not np.isfinite(v): v=1e8
        out.append(v)
    return np.array(out)/1000
def gam(x,T,mode):
    if mode=='ideal': return np.ones(len(x))
    return np.array(UNIFAC.from_subgroups(T=T,xs=list(x),chemgroups=grp,version=0,interaction_data=UFIP,subgroups=UFSG).gammas())
def flash(n,T,mode):
    """isothermal flash of moles n at T, P=Pc; returns V moles vector, L moles vector"""
    if n.sum()<1e-12: return np.zeros_like(n),n.copy()
    z=n/n.sum(); ps=psat(T); x=z.copy()
    for it in range(60):
        K=gam(x,T,mode)*ps/Pc
        f=lambda b:(z*(K-1)/(1+b*(K-1))).sum()
        if f(0)<=0: return np.zeros_like(n),n.copy()
        if f(1)>=0: return n.copy(),np.zeros_like(n)
        b=brentq(f,0,1); xn=z/(1+b*(K-1)); xn/=xn.sum()
        if np.abs(xn-x).max()<1e-9: break
        x=0.5*x+0.5*xn
    V=n.sum()*b*K*x; return V,n-V
def tflash_frac(n,phi,mode,Thi):
    """partial condenser: find T<Thi such that vapor moles fraction = 1-phi; returns T,V,L"""
    tot=n.sum()
    g=lambda T:flash(n,T,mode)[0].sum()/tot-(1-phi)
    T=brentq(g,Thi-80,Thi+1e-6) if g(Thi+1e-6)>0 else Thi
    V,L=flash(n,T,mode); return T,V,L
Ts=[40,45,50,55,60,65,70,75,80,85,90,100,110,120,125,130,140,150]
def run(xv,mode,phi=0.0):
    vol={k:A[k]*(100-xv)/100 for k in A}; vol[M]=xv
    n=np.array([vol[k] for k in names])/vm; F=(n*vm).sum()
    L=n.copy(); cum=0; out=[]
    for T in Ts:
        V,L=flash(L,T+273.15,mode)
        if phi>0 and V.sum()>1e-12:
            Tc,D,R=tflash_frac(V,phi,mode,T+273.15); L=L+R; V=D; Tr=Tc-273.15
        else: Tr=T
        cum+=(V*vm).sum()/F*100; out.append((Tr,cum))
    return out
def Tat(out,f):
    T=[a for a,b in out]; c=[b for a,b in out]
    # first point with vapor defines start
    return float(np.interp(f,c,T))
if __name__=='__main__':
    res={}
    for mode in ('unifac','ideal'):
        for phi in (0.0,0.3,0.5):
            o0=run(0,mode,phi); o10=run(10,mode,phi)
            r=dict(T10b=Tat(o0,10),dT10=Tat(o10,10)-Tat(o0,10),dT5=Tat(o10,5)-Tat(o0,5),dT20=Tat(o10,20)-Tat(o0,20),dT50=Tat(o10,50)-Tat(o0,50))
            print(mode,phi,{k:round(v,2) for k,v in r.items()},flush=True)

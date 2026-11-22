"""M2 test: batch (Rayleigh) cascade of surrogate A with and without one rectification stage (partial condenser).
Reproduces the HYSYS batch route outside HYSYS (thermo 0.6.1 / chemicals 1.5.2, ideal vapor, UNIFAC or ideal liquid).
Outputs: reflux.json (evaporated-volume axis, Eq. 2) and reflux_eq1.json (temperature axis, Eq. 1). See Table S16."""
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
def psat(T): return np.array([C[k].VaporPressure(T) for k in names])/1000
def gam(x,T,mode):
    if mode=='ideal': return np.ones(len(x))
    return np.array(UNIFAC.from_subgroups(T=T,xs=list(x),chemgroups=grp,version=0,interaction_data=UFIP,subgroups=UFSG).gammas())
def flash(n,T,mode):
    """isothermal flash of moles n at T, P=Pc; returns V moles vector, L moles vector"""
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
def _summary():
    res={}
    for mode in ('unifac','ideal'):
        for phi in (0.0,0.3,0.5):
            o0=run(0,mode,phi); o10=run(10,mode,phi)
            r=dict(T10b=Tat(o0,10),dT10=Tat(o10,10)-Tat(o0,10),dT5=Tat(o10,5)-Tat(o0,5),dT20=Tat(o10,20)-Tat(o0,20),dT50=Tat(o10,50)-Tat(o0,50))
            print(mode,phi,{k:round(v,2) for k,v in r.items()},flush=True)

# ---- evaporated-volume axis (Eq. 2)
def eq2():
    fm=np.array([0,5,10,20,30,40,50,60,70,80,90]);Tm=np.array([43.2,65.5,72.3,83.4,93.6,103.3,114.4,126,138.2,151,167.7])
    Tg=np.arange(43.2,150,0.01); Em=np.interp(Tg,Tm,fm)
    def curve(o):
        T=np.array([a for a,b in o]);c=np.array([b for a,b in o])
        # prepend zero point at first temperature with vapor minus small
        return T,c
    def E_at(o,t):
        T,c=curve(o); return float(np.interp(t,T,c,left=0))
    res={}
    for mode in ('unifac','ideal'):
        for phi in (0.0,0.3,0.5):
            O={x:run(x,mode,phi) for x in (0,10,15)}
            mono=all(np.diff([a for a,b in O[x]]).min()>-1e-6 for x in O)
            def corr(x):
                E=Em+np.array([E_at(O[x],t)-E_at(O[0],t) for t in Tg]); return np.maximum.accumulate(E)
            E10=corr(10);E15=corr(15)
            T10=lambda E:float(Tg[np.argmax(E>=10)])
            r=dict(mono=mono,E70rate10=(np.interp(70,Tg,E10)-np.interp(70,Tg,Em))/10,E70rate15=(np.interp(70,Tg,E15)-np.interp(70,Tg,Em))/15,
                   dT10_delta=T10(E10)-T10(Em),dT50_delta=float(Tg[np.argmax(E10>=50)]-Tg[np.argmax(Em>=50)]))
            res[f'{mode}_{phi}']=r; print(mode,phi,{k:(round(v,3) if not isinstance(v,bool) else v) for k,v in r.items()},flush=True)
    json.dump(res,open('reflux.json','w'),indent=1)

# ---- temperature axis (Eq. 1)
def eq1():
    fm=np.array([0,5,10,20,30,40,50,60,70,80,90]);Tm=np.array([43.2,65.5,72.3,83.4,93.6,103.3,114.4,126,138.2,151,167.7])
    fg=np.arange(1,90.01,0.01)
    def Tf(o,f):
        T=np.array([a for a,b in o]);c=np.array([b for a,b in o])
        k=c>0; T=T[k];c=c[k]
        # enforce monotone T vs cum
        T=np.maximum.accumulate(T)
        return np.interp(f,c,T)
    out={}
    for mode in ('unifac','ideal'):
        for phi in (0.0,0.3,0.5):
            O={x:run(x,mode,phi) for x in (0,5,10,15,20)}
            base=Tf(O[0],fg); Tmg=np.interp(fg,fm,Tm)
            r={}
            for x in (5,10,15,20):
                Tc=Tmg+Tf(O[x],fg)-base; Tc=np.maximum.accumulate(Tc)
                E70=float(np.interp(70,Tc,fg)); r[x]=dict(E70=round(E70,2),dT10=round(float(np.interp(10,fg,Tc)-72.3),2),dT50=round(float(np.interp(50,fg,Tc)-114.4),2))
            E0=float(np.interp(70,Tmg,fg))
            print(mode,phi,'E70 base %.2f'%E0,r,'rate10 %.3f rate15 %.3f'%((r[10]['E70']-E0)/10,(r[15]['E70']-E0)/15),flush=True)
            out[f'{mode}_{phi}']=r
    json.dump(out,open('reflux_eq1.json','w'),indent=1)

if __name__=='__main__':
    eq2(); eq1()

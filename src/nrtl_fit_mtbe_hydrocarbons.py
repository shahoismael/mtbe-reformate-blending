"""Fit NRTL (alpha=0.3, Bij=0) and check UNIFAC for MTBE-hydrocarbon pairs.
Data: Clark et al. 1997 (Fluid Phase Equilib. 133, 239-246), Chamorro et al. 2004 (Fluid Phase Equilib. 220, 105-112).
Bennett et al. 1993 MTBE+1-heptene: see nrtl_fit_bennett1993.py.
Part 1: ideal-vapour fits (used for C6+ paraffins and aromatics). Part 2: virial-corrected fits (used for C4, C5).
Requires: numpy, scipy, thermo, chemicals."""
import numpy as np,warnings,json;warnings.filterwarnings('ignore')
from scipy.optimize import least_squares
from thermo import Chemical
from thermo.unifac import UNIFAC,UFSG,UFIP
def nrtl(x1,t12,t21,a=0.3):
    x2=1-x1;G12,G21=np.exp(-a*t12),np.exp(-a*t21)
    l1=x2**2*(t21*(G21/(x1+x2*G21))**2+t12*G12/(x2+x1*G12)**2)
    l2=x1**2*(t12*(G12/(x2+x1*G12))**2+t21*G21/(x1+x2*G21)**2)
    return np.exp(l1),np.exp(l2)
ps=lambda n,T:Chemical(n,T=T).Psat/1000
# --- Chamorro 2004, 313.15 K, measured Psat
ch_iso=np.array([[0.0546,15.929],[0.1013,18.467],[0.1509,21.118],[0.1967,23.478],[0.2463,25.974],[0.2996,28.597],[0.3467,30.874],[0.3980,33.326],[0.3981,33.375],[0.4483,35.695],[0.4484,35.672],[0.4984,37.974],[0.4990,37.979],[0.5477,40.187],[0.5493,40.248],[0.5974,42.385],[0.5997,42.476],[0.6471,44.573],[0.6978,46.774],[0.7461,48.890],[0.7958,51.027],[0.8486,53.300],[0.8961,55.350],[0.9467,57.566]])
ch_tol=np.array([[0.0599,11.430],[0.1017,13.883],[0.1507,16.658],[0.1998,19.385],[0.2473,21.956],[0.2991,24.702],[0.3507,27.376],[0.3604,27.855],[0.4006,29.936],[0.4085,30.320],[0.4485,32.357],[0.4580,32.823],[0.4992,34.903],[0.5080,35.322],[0.5494,37.398],[0.5592,37.874],[0.5999,39.884],[0.6105,40.400],[0.6633,42.999],[0.7177,45.688],[0.7748,48.497],[0.8310,51.307],[0.8926,54.396],[0.9716,58.431]])
# --- Clark 1997: x1, P(Pa), T(K)
cl={'n-butane':[[0.0836,230488,298.52],[0.0836,460283,322.89],[0.2124,204081,298.44],[0.2124,408517,322.88],[0.2992,186942,298.46],[0.2992,374760,322.89],[0.4149,164317,298.48],[0.4149,330113,322.90],[0.4978,147528,298.42],[0.4978,298045,322.89]],
'isopentane':[[0.1087,88479,298.41],[0.1087,194574,322.91],[0.2070,84051,298.52],[0.2070,184984,322.90],[0.3031,79233,298.50],[0.3031,175487,322.89],[0.3990,74066,298.53],[0.3990,165281,322.89],[0.4972,68450,298.51],[0.4972,154060,322.89]],
'2,2,4-trimethylpentane':[[0.1235,12257,298.47],[0.1235,31502,322.88],[0.2066,15420,298.48],[0.2066,38548,322.88],[0.2781,17547,298.45],[0.2781,43707,322.87],[0.4054,20951,298.44],[0.4054,52280,322.87],[0.4994,23250,298.42],[0.4994,57821,322.88]],
'o-xylene':[[0.1107,6313,298.47],[0.1107,15772,322.89],[0.2045,9866,298.46],[0.2045,24570,322.88],[0.2929,13061,298.46],[0.2929,32320,322.88],[0.3992,16610,298.45],[0.3992,41132,322.88],[0.5046,19974,298.45],[0.5046,49225,322.88]]}
groups={'MTBE':{1:3,4:1,24:1},'n-butane':{1:2,2:2},'isopentane':{1:3,2:1,3:1},'2,2,4-trimethylpentane':{1:5,2:1,3:1,4:1},'toluene':{9:5,11:1},'o-xylene':{9:4,11:2},'1-heptene':{1:1,2:4,5:1}}
def unif_inf(hc,T):
    G=UNIFAC.from_subgroups(T=T,xs=[1e-6,1-1e-6],chemgroups=[groups['MTBE'],groups[hc]],version=0,interaction_data=UFIP,subgroups=UFSG);return G.gammas()[0]
def unif_P(hc,x,T,p1,p2):
    p1=np.broadcast_to(p1,x.shape);p2=np.broadcast_to(p2,x.shape);out=[]
    for xi,Ti,a,b in zip(x,T,p1,p2):
        g=UNIFAC.from_subgroups(T=Ti,xs=[xi,1-xi],chemgroups=[groups['MTBE'],groups[hc]],version=0,interaction_data=UFIP,subgroups=UFSG).gammas()
        out.append(xi*g[0]*a+(1-xi)*g[1]*b)
    return np.array(out)
R={}
def fit(name,hc,x,P,T,p1,p2,src):
    r=lambda b:(lambda g:(x*g[0]*p1+(1-x)*g[1]*p2)/P-1)(nrtl(x,b[0]/T,b[1]/T))
    f=least_squares(r,[100,-50]);b=f.x
    Pid=x*p1+(1-x)*p2; Pu=unif_P(hc,x,T,p1,p2); Tm=T.mean()
    gi=nrtl(1e-6,b[0]/Tm,b[1]/Tm)[0]
    R[name]=dict(src=src,b12=round(b[0],1),b21=round(b[1],1),A12_cal=round(b[0]*1.98721,1),A21_cal=round(b[1]*1.98721,1),
      MAD_fit=round(np.mean(abs(r(b)))*100,2),MAD_ideal=round(np.mean(abs(Pid/P-1))*100,2),MAD_unifac=round(np.mean(abs(Pu/P-1))*100,2),
      ginf_fit=round(gi,3),ginf_unifac=round(unif_inf(hc,Tm),3),T=round(Tm,1))
T=np.full(len(ch_iso),313.15); fit('isooctane (Chamorro)','2,2,4-trimethylpentane',ch_iso[:,0],ch_iso[:,1],T,59.92,12.973,'Chamorro 2004')
T=np.full(len(ch_tol),313.15); fit('toluene (Chamorro)','toluene',ch_tol[:,0],ch_tol[:,1],T,59.92,7.897,'Chamorro 2004')
for hc,d in cl.items():
    d=np.array(d);x,P,T=d[:,0],d[:,1]/1000,d[:,2]
    p1=np.array([ps('methyl tert-butyl ether',t) for t in T]);p2=np.array([ps(hc,t) for t in T])
    fit(hc+' (Clark)',hc,x,P,T,p1,p2,'Clark 1997')
# Bennett 1-heptene earlier fit at 311 K
R['1-heptene (Bennett)']=dict(src='Bennett 1993',b12=148.0,b21=-89.1,A12_cal=294.1,A21_cal=-177.1,MAD_fit=0.34,MAD_ideal=2.96,MAD_unifac=4.03,ginf_fit=1.13,ginf_unifac=1.39,T=311.2)
for k,v in R.items(): print(k,v)
json.dump(R,open('fits.json','w'),indent=1)
print('Psat check MTBE 313.15 thermo',ps('methyl tert-butyl ether',313.15),'vs 59.92; isooctane',ps('2,2,4-trimethylpentane',313.15),'vs 12.973')

# ---- Part 2: virial vapour correction

from chemicals.virial import BVirial_Tsonopoulos
from scipy.optimize import brentq
def B(n,T):
    c=Chemical(n);return BVirial_Tsonopoulos(T,c.Tc,c.Pc,c.omega)  # m3/mol
def bubbleP(x1,T,b,n1,n2,p1,p2):
    g=nrtl(x1,b[0]/T,b[1]/T);B1,B2=B(n1,T),B(n2,T);R=8.314
    def f(P):
        y1=x1*g[0]*p1/(P*np.exp(B1*(P-p1)*1e3/R/T)); y2=(1-x1)*g[1]*p2/(P*np.exp(B2*(P-p2)*1e3/R/T)); return y1+y2-1
    return brentq(f,1,2000)
for hc,d in cl.items():
    d=np.array(d);x,P,T=d[:,0],d[:,1]/1000,d[:,2]
    p1=np.array([ps('methyl tert-butyl ether',t) for t in T]);p2=np.array([ps(hc,t) for t in T])
    r=lambda b:np.array([bubbleP(xi,Ti,b,'methyl tert-butyl ether',hc,a,c) for xi,Ti,a,c in zip(x,T,p1,p2)])/P-1
    f=least_squares(r,[100,0]);b=f.x;Tm=T.mean()
    print(hc,'virial fit b12 %.1f b21 %.1f MAD %.2f%% ginf %.3f'%(b[0],b[1],np.mean(abs(r(b)))*100,nrtl(1e-6,b[0]/Tm,b[1]/Tm)[0]))
# Chamorro with virial
for nm,hc,dat,p2 in [('iso','2,2,4-trimethylpentane',ch_iso,12.973),('tol','toluene',ch_tol,7.897)]:
    x,P=dat[:,0],dat[:,1];T=313.15
    r=lambda b:np.array([bubbleP(xi,T,b,'methyl tert-butyl ether',hc,59.92,p2) for xi in x])/P-1
    f=least_squares(r,[100,0]);b=f.x
    print(nm,'virial fit b12 %.1f b21 %.1f MAD %.2f%% ginf %.3f'%(b[0],b[1],np.mean(abs(r(b)))*100,nrtl(1e-6,b[0]/T,b[1]/T)[0]))
print('B MTBE 313 cm3/mol',B('methyl tert-butyl ether',313.15)*1e6,'(Chamorro -1331)')

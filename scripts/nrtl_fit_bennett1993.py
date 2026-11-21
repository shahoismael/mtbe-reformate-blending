import numpy as np
from thermo.unifac import UNIFAC, UFSG, UFIP
from scipy.optimize import least_squares
# 311.15 K data: P kPa, y2 (MTBE, vapour), x2 (MTBE, liquid)
D=np.array([[55.45,1.000,1.000],[53.94,.990,.963],[51.77,.970,.900],[49.29,.947,.837],[46.76,.925,.772],[44.31,.898,.711],[42.27,.875,.656],[39.92,.846,.599],[39.89,.849,.599],[37.76,.822,.550],[34.93,.779,.478],[32.00,.728,.406],[29.39,.680,.348],[26.69,.611,.283],[24.53,.549,.226],[21.80,.457,.173],[19.19,.346,.116],[17.79,.273,.084],[15.12,.108,.032],[13.64,0,0]])
T=311.15; P2s,P1s=55.45,13.64  # MTBE, 1-heptene
m=(D[:,1]>0)&(D[:,1]<1); P,y2,x2=D[m,0],D[m,1],D[m,2]
g2=y2*P/(x2*P2s); g1=(1-y2)*P/((1-x2)*P1s)
print('exp gamma MTBE',np.round(g2,3)); print('exp gamma heptene',np.round(g1,3))
# UNIFAC: MTBE = 4 CH3? groups: MTBE: (CH3)3C-O-CH3 -> 3 CH3, 1 C, 1 CH3O ; 1-heptene: CH2=CH (1) + 4 CH2 + 1 CH3
mtbe={1:3,4:1,24:1}   # CH3=1, C=4, CH3O=24
hept={1:1,2:4,5:1}    # CH3, CH2, CH2=CH (5)
def unifac_g(x):
    G=UNIFAC.from_subgroups(T=T,xs=[x,1-x],chemgroups=[mtbe,hept],version=0,interaction_data=UFIP,subgroups=UFSG)
    return G.gammas()
U=np.array([unifac_g(x) for x in x2]); print('UNIFAC gamma MTBE',np.round(U[:,0],3))
Pcalc_U=x2*U[:,0]*P2s+(1-x2)*U[:,1]*P1s; Pid=x2*P2s+(1-x2)*P1s
print('mean abs dev P: ideal %.2f%%, UNIFAC %.2f%%'%(np.mean(abs(Pid/P-1))*100,np.mean(abs(Pcalc_U/P-1))*100))
print('UNIFAC gamma_inf MTBE in heptene',unifac_g(1e-6)[0],'exp approx',g2[-1])
# NRTL fit alpha=0.3
def nrtl(x,t12,t21,a=0.3):
    G12,G21=np.exp(-a*t12),np.exp(-a*t21); x1=x; x2_=1-x
    l1=x2_**2*(t21*(G21/(x1+x2_*G21))**2+t12*G12/(x2_+x1*G12)**2)
    l2=x1**2*(t12*(G12/(x2_+x1*G12))**2+t21*G21/(x1+x2_*G21)**2)
    return np.exp(l1),np.exp(l2)
def res(p):
    gm,gh=nrtl(x2,*p); return (x2*gm*P2s+(1-x2)*gh*P1s)/P-1
f=least_squares(res,[0.1,0.1]); gm,gh=nrtl(x2,*f.x)
print('NRTL tau12,tau21',np.round(f.x,4),'b12,b21 (K)',np.round(f.x*T,1),'MAD P %.2f%%'%(np.mean(abs(res(f.x)))*100),'gamma_inf MTBE',round(nrtl(1e-6,*f.x)[0],3))

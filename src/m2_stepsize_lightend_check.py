import sys,csv,numpy as np
import m2_cascade_core as casc
from m2_cascade_core import *
G5=[40,45,50,55,60,65,70,75,80,85,90,100,110,120,125,130,140,150,160,170,180,190,200]
G25=sorted(set(list(np.arange(40,90.01,2.5))+[100,110,120,125,130,140,150,160,170,180,190,200]))
def runfull(xv,mode,Ts,loss=False):
    vol={k:A[k]*(100-xv)/100 for k in A}; vol[M]=xv
    n=np.array([vol[k] for k in names])/vm; F=(n*vm).sum(); L=n.copy(); out=[]
    ib=[names.index('n-butane')]
    for T in Ts:
        V,L=flash(L,T+273.15,mode)
        if loss: V=V.copy(); V[ib]=0
        out.append((V*vm).sum()/F*100)
    return out
def write(fn,Ts,mode,loss=False):
    rows=[['State','x']+[f'{t:g}V - Phase Std Ideal Liq Vol Flow (Vapour Phase)' for t in Ts]+['BLENDED - Std Ideal Liq Vol Flow'],[''],['']]
    for x in range(21):
        v=runfull(x,mode,Ts,loss); rows.append([f'Case {x}',str(x)]+[f'{a:.5f}' for a in v]+['100'])
        print(fn,x,flush=True)
    csv.writer(open(fn,'w')).writerows(rows)
job=sys.argv[1]
if job=='g5': write('py_uni_g5.csv',G5,'unifac')
if job=='g25': write('py_uni_g25.csv',G25,'unifac')
if job=='loss': write('py_uni_loss.csv',G5,'unifac',True)

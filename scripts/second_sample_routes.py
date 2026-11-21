"""Table S12: route changes of surrogate A applied to the second Al-Daura sample (S2), against the first sample (S1).
Run from the package folder: python scripts/second_sample_routes.py   (writes s2.json)"""
import numpy as np,openpyxl,csv,json
Z=''
F=lambda t:t*1.8+32
# simulated curve sets: dict x -> (T, E) increasing T
def route1():
    ws=openpyxl.load_workbook(Z+'Route1_Surrogate_D86_NRTL.xlsx',data_only=True)['D86_curves']
    rows=list(ws.iter_rows(min_row=5,values_only=True)); hdr=rows[0]; data=[r for r in rows[1:] if r[0] is not None]
    pe=np.array([r[0] for r in data],float); out={}
    for j,x in enumerate(hdr[1:],1):
        if x is None: continue
        T=np.array([r[j] for r in data],float); out[float(x)]=(T,pe)
    return out
def route3():
    ws=openpyxl.load_workbook(Z+'Route3_Surrogate_EFV_NRTL.xlsx',data_only=True)['Raw HYSYS']
    out={}
    for r in ws.iter_rows(min_row=4,values_only=True):
        if not r[0] or r[2] is None or (isinstance(r[4],(int,float)) and r[4]<-1): continue
        x=float(r[1]); out.setdefault(x,[]).append((r[3],100*r[4]/r[6]))
    res={}
    for x,v in out.items():
        v=sorted(v); res[x]=(np.array([a for a,b in v]),np.array([b for a,b in v]))
    return res
def batch(fn,base=None):
    """HYSYS batch export; non-converged stage values (-32767) interpolated between states and,
    if base is given, the 0 vol% state completed from that file (Table S11 note)."""
    def load(f):
        rows=list(csv.reader(open(f))); h=rows[0]
        st=[i for i,x in enumerate(h) if x[:1].isdigit() and x.split(' ')[0].endswith('V')]
        bi=h.index('BLENDED - Std Ideal Liq Vol Flow'); X=[];V=[]
        for r in rows[3:]:
            if not r or not r[0]: continue
            v=np.array([float(r[i]) for i in st]); X.append(float(r[1])); V.append(np.where(v<-1,np.nan,v/float(r[bi])*100))
        return np.array([float(h[i].split('V')[0]) for i in st]),np.array(X),np.array(V)
    Ts,X,V=load(fn)
    if base: V[0]=np.where(np.isnan(V[0]),load(base)[2][0],V[0])
    for j in range(V.shape[1]):
        k=~np.isnan(V[:,j]); xs,ys=X[k],V[k,j]
        for i in np.where(~k)[0]:
            V[i,j]=np.polyval(np.polyfit(xs[-2:],ys[-2:],1),X[i]) if X[i]>xs.max() else np.interp(X[i],xs,ys)
    return {float(x):(Ts,np.cumsum(v)) for x,v in zip(X,V)}
def analyse(sim,Tm,fm,Tend):
    xs=sorted(sim); T0,E0=sim[xs[0]]
    Tg=np.arange(Tm[0],Tm[-1],0.01); Em=np.interp(Tg,Tm,fm); out={}
    for x in xs:
        T,E=sim[x]
        Ec=np.maximum.accumulate(Em+np.interp(Tg,T,E)-np.interp(Tg,T0,E0))
        inv=lambda f: float(Tg[np.argmax(Ec>=f)]) if Ec.max()>=f else np.nan
        T10,T50,T90=inv(10),inv(50),inv(90)
        if x==0: T10,T50,T90=[float(np.interp(f,fm,Tm)) for f in (10,50,90)]
        out[x]=dict(E70=float(np.interp(70,Tg,Ec)),E100=float(np.interp(100,Tg,Ec)),T10=T10,T50=T50,T90=T90,DI=1.5*F(T10)+3*F(T50)+F(T90))
    return out
def bounds(o):
    xs=np.array(sorted(o)); g=lambda k,v: float(np.interp(v,xs,[o[x][k] for x in xs]))
    def first(k,c):
        for v in np.arange(0,20.01,0.01):
            if c(g(k,v)): return int(np.ceil(round(v,2)-1e-9))
        return None
    e0=o[xs[0]]['E70']
    return dict(E70b=first('E70',lambda v:v>=20),E100b=first('E100',lambda v:v>=46),DIb=first('DI',lambda v:v<=1250),T10b=first('T10',lambda v:v<=70),
                E70_0=round(e0,1),E70_10=round(g('E70',10),1),E70_20=round(g('E70',20),1),rate10=round((g('E70',10)-e0)/10,2),E100_0=round(g('E100',0),1),DI_0=round(g('DI',0)),T10_0=round(g('T10',0),1),dT10_10=round(g('T10',10)-g('T10',0),1))
S1=(np.array([43.2,65.5,72.3,83.4,93.6,103.3,114.4,126,138.2,151,167.7]),np.array([0,5,10,20,30,40,50,60,70,80,90.]))
S2=(np.array([43,71,80,86,92,97,103,110,121,144,154.]),np.array([0,10,20,30,40,50,60,70,80,90,95.]))
routes={'Surrogate D86':route1(),'EFV NRTL':route3(),'Batch NRTL':batch(Z+'Route5_SurrogateA_Batch_TVP_NRTL_UNIFAC.csv'),'Batch PR':batch(Z+'Route6_SurrogateA_Batch_TVP_PR_kij0.csv',Z+'Route7b_SurrogateA_Batch_TVP_PR_fitted.csv')}
R={}
for name,sim in routes.items():
    for sname,(Tm,fm) in (('S1',S1),('S2',S2)):
        b=bounds(analyse(sim,Tm,fm,None)); R[(name,sname)]=b; print(name,sname,b)
json.dump({f'{a}|{b}':v for (a,b),v in R.items()},open('s2.json','w'),indent=0)

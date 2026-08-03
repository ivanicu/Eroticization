import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
The killer confound for TEMPO: the onset bins are coarse (0-4, 5-6 ... 19-25, 26+). Two
categories whose population means are both ~15.5 will land in the SAME bin for many people
purely mechanically, while categories 3 years apart rarely share one. That alone makes residual
onset correlation decay with mean-onset distance -- with no developmental tempo at all.

Parametric bootstrap null: generate onsets from person effect + category effect + iid noise
(NO tempo structure, by construction), push them through the SAME binning and the SAME
missingness pattern, and run the identical TEMPO regression. If the null reproduces t=4.6 the
finding is binning. If it gives t~0, tempo is real.
"""
import pandas as pd, numpy as np, re, warnings
from numpy.linalg import lstsq
warnings.filterwarnings('ignore'); rng=np.random.default_rng(31337)
exec(open(round_path('24_attack_rsa.py')).read().split('print("=== is onset a proxy')[0])
catcols=[v for _,v in pairs]
EDGES=[0,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,25.5,99]
MIDS =[2,5.5,7.5,9.5,11.5,13.5,15.5,17.5,22,28]
def binify(x):
    return np.array(MIDS)[np.clip(np.digitize(x,EDGES)-1,0,len(MIDS)-1)]

Oraw=O.copy()                       # binned midpoints, as observed
mask=Oraw.notna().values
mu=np.nanmean(Oraw.values)
pers=np.nanmean(np.where(mask,Oraw.values,np.nan),axis=1)-mu
cat =np.nanmean(np.where(mask,Oraw.values,np.nan),axis=0)-mu
E=Oraw.values-(mu+pers[:,None]+cat[None,:]); sd=np.nanstd(E)
meanons=np.array([Oraw[c].mean() for c in catcols])
print(f"person-effect sd {np.nanstd(pers):.2f}  category-effect sd {np.nanstd(cat):.2f}  residual sd {sd:.2f}")

def tempo_t(values):
    Om=pd.DataFrame(np.where(mask,values,np.nan),columns=catcols)
    Z=Om.sub(Om.mean(axis=1),axis=0); Z=Z.sub(Z.mean(axis=0),axis=1)
    for c in Z.columns:                                    # partial current age, as in iter 3
        m=Z[c].notna()&age.notna()
        if m.sum()>200:
            X=np.c_[np.ones(m.sum()),age[m].values]; b,*_=lstsq(X,Z.loc[m,c].values,rcond=None)
            Z.loc[m,c]=Z.loc[m,c].values-X@b
    CO=np.full((k,k),np.nan); CP=np.full((k,k),np.nan)
    for i in range(k):
        for j in range(i+1,k):
            m=Z.iloc[:,i].notna()&Z.iloc[:,j].notna()&Pres.iloc[:,i].notna()&Pres.iloc[:,j].notna()
            if m.sum()<150: continue
            CO[i,j]=CO[j,i]=np.corrcoef(Z.iloc[:,i][m],Z.iloc[:,j][m])[0,1]
            CP[i,j]=CP[j,i]=np.corrcoef(Pres.iloc[:,i][m],Pres.iloc[:,j][m])[0,1]
    iu=np.triu_indices(k,1); ok=(~np.isnan(CO[iu]))&(~np.isnan(CP[iu]))
    y=CO[iu][ok]; pref=CP[iu][ok]
    T=-np.abs(meanons[iu[0][ok]]-meanons[iu[1][ok]])
    X=np.c_[np.ones(len(y)),pref,T]; b,*_=lstsq(X,y,rcond=None); r=y-X@b
    se=np.sqrt((r**2).sum()/(len(y)-3)*np.diag(np.linalg.pinv(X.T@X)))
    R2=1-(r**2).sum()/((y-y.mean())**2).sum()
    Xb=np.c_[np.ones(len(y)),pref]; bb,*_=lstsq(Xb,y,rcond=None)
    R2b=1-((y-Xb@bb)**2).sum()/((y-y.mean())**2).sum()
    return b[2]/se[2], R2-R2b
t_obs,d_obs=tempo_t(Oraw.values)
print(f"\nOBSERVED   TEMPO t = {t_obs:+.2f}   dR2 = {d_obs:+.4f}")
ts=[];ds=[]
for rep in range(40):
    sim=mu+pers[:,None]+cat[None,:]+rng.normal(0,sd,size=Oraw.shape)
    t,d=tempo_t(binify(sim)); ts.append(t); ds.append(d)
ts=np.array(ts); ds=np.array(ds)
print(f"BINNING NULL (40 draws, no tempo by construction)")
print(f"   TEMPO t  : mean {ts.mean():+.2f}  sd {ts.std():.2f}  max {ts.max():+.2f}")
print(f"   dR2      : mean {ds.mean():+.4f}  sd {ds.std():.4f}")
print(f"\n   z of observed vs binning null : {(t_obs-ts.mean())/ts.std():+.1f}")
print(f"   p (null t >= observed t)      : {np.mean(ts>=t_obs):.3f}")
print(f"\n   VERDICT: {'TEMPO IS A BINNING ARTIFACT' if np.mean(ts>=t_obs)>0.05 else 'tempo survives the binning null'}")

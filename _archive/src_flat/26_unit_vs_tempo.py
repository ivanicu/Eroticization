"""
ITER 4. Two live readings of the onset structure, and they predict DIFFERENT things.

  UNIT ACQUISITION : categories sharing a GCCA coordinate were acquired together.
      predictor = coordinate-loading similarity between categories.
  DEVELOPMENTAL TEMPO : a person's maturational schedule shifts categories that the POPULATION
      places at similar ages together, regardless of any coordinate.
      predictor = -|mean_onset_i - mean_onset_j|  (population arrival-time distance)

Both regressed on within-person onset similarity, with preference similarity controlled.
Prediction matrix: unit -> coordinate term survives, tempo term dies. tempo -> the reverse.
Both -> decompose. Neither -> the iter-3 residual structure is something else again.
"""
import pandas as pd, numpy as np, re, warnings
from numpy.linalg import lstsq, qr, svd
warnings.filterwarnings('ignore'); rng=np.random.default_rng(1123)
exec(open('src/24_attack_rsa.py').read().split('print("=== is onset a proxy')[0])

# strip intensity leakage (iter 3 control) and rebuild the two similarity matrices
Ostr=Ores.copy()
for c in Ores.columns:
    m=Ores[c].notna()&Pres[c].notna()
    X=np.c_[np.ones(m.sum()),Pres.loc[m,c].values]; b,*_=lstsq(X,Ores.loc[m,c].values,rcond=None)
    Ostr.loc[m,c]=Ores.loc[m,c].values-X@b
CO=np.full((k,k),np.nan); CP=np.full((k,k),np.nan)
for i in range(k):
    for j in range(i+1,k):
        m=Ostr.iloc[:,i].notna()&Ostr.iloc[:,j].notna()&Pres.iloc[:,i].notna()&Pres.iloc[:,j].notna()
        if m.sum()<150: continue
        CO[i,j]=CO[j,i]=np.corrcoef(Ostr.iloc[:,i][m],Ostr.iloc[:,j][m])[0,1]
        CP[i,j]=CP[j,i]=np.corrcoef(Pres.iloc[:,i][m],Pres.iloc[:,j][m])[0,1]

# --- category -> GCCA coordinate loading, via the gate column that routes into each block ---
exec(open('src/16_dimensionality.py').read().split("allq=list(B)")[0])
allq=list(B)
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts()
pool=np.array(sorted(nblk[nblk>=8].index)); pm={p:i for i,p in enumerate(pool)}
G=np.load('data/derived/gcca_G.npy')
br=pd.read_csv('data/derived/branching.csv')
gate={int(r.qi):r.gate for _,r in br.iterrows()}
D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
catcols=[v for _,v in pairs]
load=np.full((k,5),np.nan)
for ci,cat in enumerate(catcols):
    qs=[q for q in allq if gate.get(q,'').strip('"')[:len(cat)]==cat or cat in str(gate.get(q,''))]
    if not qs: continue
    vs=[]
    for q in qs:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan); Z[idx]=B[q]['R'][src]
        Z=np.where(np.isnan(Z),np.nanmean(Z,axis=0),Z); Z=Z-Z.mean(0)
        for kk in range(1,6):
            vs.append([kk,np.mean([abs(np.corrcoef(G[:,kk],Z[:,j])[0,1]) for j in range(Z.shape[1])])])
    v=pd.DataFrame(vs,columns=['k','a']).groupby('k').a.mean()
    for kk in range(1,6): load[ci,kk-1]=v.get(kk,np.nan)
have=~np.isnan(load).any(1)
print(f"categories mapped to a GCCA block: {have.sum()}/{k}")

meanons=np.array([O[c].mean() for c in catcols])
iu=np.triu_indices(k,1)
ok=(~np.isnan(CO[iu]))&(~np.isnan(CP[iu]))&have[iu[0]]&have[iu[1]]
def cos(a,b): return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)+1e-12))
COORD=np.array([cos(load[i],load[j]) for i,j in zip(iu[0][ok],iu[1][ok])])
TEMPO=-np.abs(meanons[iu[0][ok]]-meanons[iu[1][ok]])
y=CO[iu][ok]; pref=CP[iu][ok]
print(f"pairs usable: {ok.sum()}\n")

def fit(cols,labels):
    X=np.c_[np.ones(len(y)),pref,*cols]
    b,*_=lstsq(X,y,rcond=None); r=y-X@b
    R2=1-(r**2).sum()/((y-y.mean())**2).sum()
    se=np.sqrt((r**2).sum()/(len(y)-X.shape[1])*np.diag(np.linalg.pinv(X.T@X)))
    return R2,[(l,b[2+i],b[2+i]/se[2+i]) for i,l in enumerate(labels)]
R2b,_=fit([],[])
print(f"baseline (preference only)          R2 = {R2b:.4f}")
for cols,labs in [([COORD],['COORD']),([TEMPO],['TEMPO']),([COORD,TEMPO],['COORD','TEMPO'])]:
    R2,terms=fit(cols,labs)
    s=" ".join(f"{l} b={bb:+.4f} t={t:+.2f}" for l,bb,t in terms)
    print(f"+{'+'.join(labs):14s} R2 = {R2:.4f}  (dR2 {R2-R2b:+.4f})   {s}")
print("\n  correlation between the two predictors: r =", round(float(np.corrcoef(COORD,TEMPO)[0,1]),3))

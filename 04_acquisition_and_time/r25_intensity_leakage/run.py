import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
The confound that manufactures RSA=0.62 for free: ONSET MAY BE A PROXY FOR INTENSITY.
If a person says "13-14" for what they like most and "17-18" for what they like least, then
within-person onset is a monotone transform of within-person preference, and the two
correlation matrices are the same matrix twice.
Test it directly (the diagonal), then partial each category's own preference out of its own
onset and redo the RSA. Also partial out survey column-position distance (recall anchoring).
"""
import pandas as pd, numpy as np, re, warnings
from numpy.linalg import lstsq
warnings.filterwarnings('ignore'); rng=np.random.default_rng(2024)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
KEY={'bondage':'bondage','humiliation':'humiliation','nonconsent':'nonconsent','sadomasochism':'sadomasochism',
 'sensory':'sensory','transformation':'transform','specific roles':'roles','mental alteration':'mentalalteration',
 'pregnancy':'pregnancy','genderplay':'genderplay','exhibitionism':'exhibitionself','multiple partners':'multiplepartners',
 'incest':'incest','bestiality':'bestiality','abnormal body':'abnormalbody','bodily-secretions':'secretions',
 'mythical':'mythical','creepy':'creepy','brutality':'brutality','vore':'vore','clothing':'clothing',
 'body parts':'appearance','gentleness':'gentleness','power dynamics':'powerdynamic','dirtiness':'dirty',
 'eagerness':'eagerness','objects':'objects','toys':'toys'}
pairs=[]
for c in ons:
    lc=c.lower()
    for k_,v in KEY.items():
        if k_ in lc and v in df.columns: pairs.append((c,v)); break
pairs=list(dict.fromkeys(pairs)); k=len(pairs)
names=[re.sub(r'.*interest in |.*sexual interest in |\?.*','',c)[:24] for c,_ in pairs]
colpos={c:i for i,c in enumerate(df.columns)}
O=pd.DataFrame({v:df[c].map(BIN) for c,v in pairs})
P=pd.DataFrame({v:pd.to_numeric(df[v],errors='coerce') for _,v in pairs})
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; age=df['age'].map(AGEMAP)
def dc(M):
    Z=M.copy(); Z=Z.sub(Z.mean(axis=1),axis=0); return Z.sub(Z.mean(axis=0),axis=1)
Ores, Pres = dc(O), dc(P)
for c in Ores.columns:                                   # partial current age out of onset
    m=Ores[c].notna()&age.notna()
    X=np.c_[np.ones(m.sum()),age[m].values]; b,*_=lstsq(X,Ores.loc[m,c].values,rcond=None)
    Ores.loc[m,c]=Ores.loc[m,c].values-X@b

print("=== is onset a proxy for intensity? (within-person, same category) ===")
dg=[]
for i,c in enumerate(Ores.columns):
    m=Ores[c].notna()&Pres[c].notna()
    if m.sum()>200: dg.append(np.corrcoef(Ores.loc[m,c],Pres.loc[m,c])[0,1])
dg=np.array(dg)
print(f"  corr(residual onset, residual preference) same category: mean {dg.mean():+.3f}  median {np.median(dg):+.3f}")
print(f"  categories with |r|>0.15 : {int((np.abs(dg)>0.15).sum())}/{len(dg)}")
print(f"  -> {'CONFOUND IS LIVE' if abs(dg.mean())>0.10 else 'leakage is small but partial it out anyway'}")

# strip each category's own preference from its own onset
Ostr=Ores.copy()
for c in Ores.columns:
    m=Ores[c].notna()&Pres[c].notna()
    X=np.c_[np.ones(m.sum()),Pres.loc[m,c].values]; b,*_=lstsq(X,Ores.loc[m,c].values,rcond=None)
    Ostr.loc[m,c]=Ores.loc[m,c].values-X@b

def rsa(Om):
    CO=np.full((k,k),np.nan); CP=np.full((k,k),np.nan)
    for i in range(k):
        for j in range(i+1,k):
            m=Om.iloc[:,i].notna()&Om.iloc[:,j].notna()&Pres.iloc[:,i].notna()&Pres.iloc[:,j].notna()
            if m.sum()<150: continue
            CO[i,j]=CO[j,i]=np.corrcoef(Om.iloc[:,i][m],Om.iloc[:,j][m])[0,1]
            CP[i,j]=CP[j,i]=np.corrcoef(Pres.iloc[:,i][m],Pres.iloc[:,j][m])[0,1]
    iu=np.triu_indices(k,1); ok=~np.isnan(CO[iu])&~np.isnan(CP[iu])
    return CO,CP,iu,ok
CO,CP,iu,ok=rsa(Ostr)
co,cp=CO[iu][ok],CP[iu][ok]
r_str=np.corrcoef(co,cp)[0,1]
# survey adjacency
pos=np.array([colpos[c] for c,_ in pairs]); Dpos=np.abs(pos[:,None]-pos[None,:]).astype(float)
d=Dpos[iu][ok]
def partial(x,y,z):
    rx=x-np.c_[np.ones(len(z)),z]@lstsq(np.c_[np.ones(len(z)),z],x,rcond=None)[0]
    ry=y-np.c_[np.ones(len(z)),z]@lstsq(np.c_[np.ones(len(z)),z],y,rcond=None)[0]
    return np.corrcoef(rx,ry)[0,1]
r_adj=partial(co,cp,np.c_[d,np.log1p(d)])
nl=[]
for _ in range(2000):
    p=rng.permutation(k); M=CP[np.ix_(p,p)]; v=M[iu][ok]; g=~np.isnan(v)
    nl.append(np.corrcoef(co[g],v[g])[0,1])
nl=np.array(nl)
print(f"\n=== RSA after the attacks ({ok.sum()} pairs) ===")
print(f"  raw (previous iteration)                        : +0.621")
print(f"  own-preference stripped from own onset          : {r_str:+.3f}")
print(f"  + survey column-distance partialled out         : {r_adj:+.3f}")
print(f"  label-permutation null                          : {nl.mean():+.3f} +/- {nl.std():.3f}")
print(f"  p(|null| >= |observed|)                         : {np.mean(np.abs(nl)>=abs(r_str)):.4f}")
print(f"  survives : {'YES' if abs(r_str)>4*nl.std() else 'NO'}")

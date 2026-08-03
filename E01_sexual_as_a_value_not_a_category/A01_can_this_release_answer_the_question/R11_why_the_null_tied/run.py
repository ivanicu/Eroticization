import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
The joint factorization tied its null (1.0x). Two very different readings:
  OVERTURNED  - there is no person-level shared structure (would retract the CCA headline)
  UNVERIFIED  - the instrument had no power to see person structure on top of item structure
P6 says these must not be folded together. Positive control decides it:
does the person side of the factorization carry ANY known person variable?
Also: the loadings printed French text -- quantify multilingual duplicate options.
"""
import pandas as pd, numpy as np, re
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
F=np.load('data/derived/joint_F.npy'); L=np.load('data/derived/joint_L.npy')
items=pd.read_csv('data/derived/joint_items.csv')['item'].values
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
qm=pd.read_csv('data/derived/multiselect_questions.csv')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1000)&(qm.mean_picks>1.5)]
persons=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; p=np.array(sorted(s.person.unique()))
    if len(p)>=1000: persons.append(p)
persons=np.unique(np.concatenate(persons))

print("=== POSITIVE CONTROL: does the joint factorization's PERSON side see anything known? ===")
y=pd.to_numeric(df.loc[persons,'biomale'],errors='coerce').values
m=~np.isnan(y)
for k in range(F.shape[1]):
    r=np.corrcoef(F[m,k],y[m])[0,1]
    print(f"  factor {k+1}:  |r(person score, biomale)| = {abs(r):.3f}")
best=max(abs(np.corrcoef(F[m,k],y[m])[0,1]) for k in range(F.shape[1]))
# multivariate
from numpy.linalg import lstsq
X=np.c_[np.ones(m.sum()),F[m]]
b,*_=lstsq(X,y[m],rcond=None); R2=1-((y[m]-X@b)**2).sum()/((y[m]-y[m].mean())**2).sum()
print(f"  all 8 factors jointly -> R2 on biomale = {R2:.3f}")
print(f"\n  VERDICT: {'UNVERIFIED - instrument blind on the person side, null tie proves nothing' if best<0.15 else 'person side is live; the null tie is informative'}")

print("\n=== multilingual duplicate options ===")
opt=pd.Series([i.split('::',1)[1] for i in items])
FR=r'(é|è|à|ç|ô|û|î)|^(Fessée|Rimmissement|Destination)'
fr=opt.str.contains(FR,regex=True,na=False)
print(f"  options with French orthography: {int(fr.sum())} / {len(opt)}")
print("  examples:", "; ".join(opt[fr].head(8).tolist())[:260])
# do they co-occur inside the SAME block as their English twin?
blk=pd.Series([i.split('::')[0] for i in items])
bad=sorted(set(blk[fr]))
print(f"  blocks containing them: {bad}")
for b in bad[:3]:
    sub=opt[(blk==b)]
    print(f"    block {b}: {len(sub)} options -> {'; '.join(sub.tolist())[:300]}")

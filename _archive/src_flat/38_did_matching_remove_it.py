"""
Self-attack on iter 7's null. Both structural measures live in WITHIN-block profile shape, and
I matched groups on block count. But "porn gave me a wholly new fetish" would most naturally
appear as entering a NEW CATEGORY -- more blocks, or rarer blocks -- which matching deletes.
If induction predicts breadth, the null is scoped to shape and must say so.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
exec(open('src/16_dimensionality.py').read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
IND=[c for c in df.columns if 'induced' in c.lower() and 'fetish' in c.lower()][0]
ORD={'No':0,'Yes: Variations on my current preexisting fetishes':1,
     'Yes: New but still similar to my preexisting fetishes':2,
     'Yes: New and totally different to my preexisting fetishes':3}
g=df[IND].map(ORD).reindex(pool)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').reindex(pool)
breadth=(R>0).sum(axis=1)          # how many CATEGORIES the person has at all
inten  =R[R>0].mean(axis=1)        # how strong, among those they have
rar=(R>0).mul(1/ (R>0).mean(axis=0).clip(lower=.01), axis=1).sum(axis=1)/ (R>0).sum(axis=1).clip(lower=1)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
cov=pd.DataFrame({'age':df['age'].map(AGEMAP),'male':df['biomale'],
                  'porn':pd.to_numeric(df['pornhabit'],errors='coerce')}).reindex(pool)
cov=cov.fillna(cov.median())
def adj(y):
    m=y.notna()
    X=np.c_[np.ones(m.sum()),cov[m].values]; b,*_=lstsq(X,y[m].values,rcond=None)
    o=pd.Series(np.nan,index=y.index); o[m]=y[m].values-X@b; return o
lab=['No','Variations','New but similar','New & totally different']
print("=== does induction predict BREADTH rather than shape? (age/sex/pornhabit adjusted) ===")
rows=[]
for k_ in range(4):
    s=g==k_
    rows.append(dict(answer=lab[k_], n=int(s.sum()),
      n_categories=round(float(breadth[s].mean()),2),
      blocks_entered=round(float(nblk[s].mean()),2),
      intensity=round(float(inten[s].mean()),3),
      rarity_of_own_set=round(float(rar[s].mean()),3)))
print(pd.DataFrame(rows).to_string(index=False))
gg=g.dropna()
for nm,v in [('n_categories',breadth),('blocks_entered',nblk),('intensity',inten),('rarity_of_own_set',rar)]:
    a=adj(v.astype(float)).reindex(gg.index)
    m=a.notna()
    r=stats.spearmanr(gg[m],a[m])
    print(f"  dose-response  {nm:18s} rho={r.statistic:+.4f}  p={r.pvalue:.2e}")
print("\n  raw (unadjusted) breadth spread: "
      f"{breadth[g==0].mean():.2f} -> {breadth[g==3].mean():.2f} categories, "
      f"a {breadth[g==3].mean()-breadth[g==0].mean():+.2f} difference of {len(rate)} possible")

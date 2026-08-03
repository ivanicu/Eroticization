import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Induction predicts breadth at rho .24 while predicting profile SHAPE at ~0. The obvious reverse
reading is nearly forced: someone with 28 interests simply has more candidates to attribute to
porn than someone with 17. Attribution-opportunity predicts breadth explains everything.
Separator: does the RARITY of a person's own set survive controlling breadth? Attribution
opportunity scales with COUNT, not with how unusual the set is.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
IND=[c for c in df.columns if 'induced' in c.lower() and 'fetish' in c.lower()][0]
ORD={'No':0,'Yes: Variations on my current preexisting fetishes':1,
     'Yes: New but still similar to my preexisting fetishes':2,
     'Yes: New and totally different to my preexisting fetishes':3}
g=df[IND].map(ORD).reindex(pool)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').reindex(pool)
H=(R>0); breadth=H.sum(axis=1)
prev=H.mean(axis=0).clip(lower=.01)
rar=H.mul(np.log(1/prev),axis=1).sum(axis=1)/breadth.clip(lower=1)     # mean log-rarity of own set
inten=R[R>0].mean(axis=1)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
cov=pd.DataFrame({'age':df['age'].map(AGEMAP),'male':df['biomale'],
                  'porn':pd.to_numeric(df['pornhabit'],errors='coerce')}).reindex(pool)
cov=cov.fillna(cov.median())
gg=g.dropna()
def partial(y,extra):
    m=y.notna()&gg.notna().reindex(y.index).fillna(False)
    X=np.c_[np.ones(m.sum()),cov[m].values]+0.0
    if extra is not None: X=np.c_[X,extra[m].values.reshape(m.sum(),-1)]
    b,*_=lstsq(X,y[m].values,rcond=None); res=y[m].values-X@b
    return stats.spearmanr(gg.reindex(y.index)[m],res)
for nm,y,ex in [('rarity | age,sex,porn',rar,None),
                ('rarity | + breadth',rar,breadth.astype(float)),
                ('rarity | + breadth + breadth^2',rar,pd.DataFrame({'b':breadth,'b2':breadth**2})),
                ('intensity | + breadth',inten,breadth.astype(float)),
                ('breadth | age,sex,porn',breadth.astype(float),None)]:
    r=partial(y,ex)
    print(f"  {nm:32s} rho={r.statistic:+.4f}  p={r.pvalue:.2e}")
print("\n  attribution-opportunity predicts: rarity | breadth -> 0")
print("  domain-expansion predicts:        rarity survives, because the set gets unusual, not just long")

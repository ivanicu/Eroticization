"""
ITER 6. "Do you feel as though you've INDUCED new fetishes into yourself through porn?" -- an
ordinal scale of ACQUIRED erotic value, n=15,502, graded exactly along Ivan's phase-3 axis
(variations -> similar -> totally different).

Temporal separator, because self-report alone cannot carry this:
  if porn induced an interest, that interest's onset should POSTDATE porn onset.
  if the report is retrospective narrative, the ordering should not move with the answer.
Porn onset is asked separately ("At what age did you begin watching porn or reading erotic
content at least semiregularly").
Second discriminator: narrative predicts a UNIFORM shift of all a person's onsets; real
induction predicts a CONCENTRATED one -- a few categories far after porn onset, the rest before.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(1597)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
PBIN={'11-12yo':11.5,'13-14yo':13.5,'15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
IND=[c for c in df.columns if 'induced' in c.lower() and 'fetish' in c.lower()][0]
PORN=[c for c in df.columns if 'begin watching porn' in c.lower()][0]
print(df[IND].value_counts().to_string())
ORD={'No':0,'Yes: Variations on my current preexisting fetishes':1,
     'Yes: New but still similar to my preexisting fetishes':2,
     'Yes: New and totally different to my preexisting fetishes':3}
g=df[IND].map(ORD)
porn=df[PORN].map(PBIN)
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300 and c!=PORN]
O=pd.DataFrame({c:df[c].map(BIN) for c in ons})
gap=O.sub(porn,axis=0)                      # interest onset minus porn onset
n=O.notna().sum(1)
ok=g.notna()&porn.notna()&(n>=5)
print(f"\nusable: {ok.sum():,}   (>=5 onsets, porn onset present, induction answered)")

AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; age=df['age'].map(AGEMAP)
def adj(y):
    """control porn-onset age and current age -- later porn onset mechanically shrinks the gap"""
    m=y.notna()&ok&porn.notna()&age.notna()
    X=np.c_[np.ones(m.sum()),porn[m].values,age[m].values,n[m].values]
    b,*_=lstsq(X,y[m].values,rcond=None); r=pd.Series(np.nan,index=df.index); r[m]=y[m].values-X@b
    return r
mean_gap=adj(gap.mean(axis=1)); max_gap=adj(gap.max(axis=1))
frac_after=adj((gap>0).sum(axis=1)/n)
conc=adj(gap.max(axis=1)-gap.mean(axis=1))     # concentration: how far the latest sits above own mean

print("\n=== mean interest-onset MINUS porn-onset, by induction answer (age/porn-onset/n adjusted) ===")
lab={0:'No',1:'Variations',2:'New but similar',3:'New & totally different'}
rows=[]
for k_ in [0,1,2,3]:
    s=ok&(g==k_)
    rows.append(dict(answer=lab[k_], n=int(s.sum()),
        mean_gap=round(float(mean_gap[s].mean()),3),
        frac_onsets_after_porn=round(float(frac_after[s].mean()),4),
        max_gap=round(float(max_gap[s].mean()),3),
        concentration=round(float(conc[s].mean()),3)))
T=pd.DataFrame(rows); print(T.to_string(index=False))
gg=g[ok].values
for nm,v in [('mean_gap',mean_gap),('frac_after',frac_after),('max_gap',max_gap),('concentration',conc)]:
    x=v[ok].values; m=~np.isnan(x)
    r=stats.spearmanr(gg[m],x[m])
    print(f"  dose-response over the 4 ordered answers: {nm:14s} rho={r.statistic:+.4f}  p={r.pvalue:.2e}")
print("\n  narrative predicts a UNIFORM shift -> mean_gap moves, concentration does not.")
print("  induction predicts a CONCENTRATED shift -> concentration moves at least as much.")

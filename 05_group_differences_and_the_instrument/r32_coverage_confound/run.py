import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
PORNHABIT reorganises the grammar more than sex does (deficit .2285 vs .0965). Before believing
that, the mechanical explanation: DIFFERENTIAL MISSINGNESS. The survey is a gated tree, so a
group that entered FEWER blocks has its loadings estimated on sparser data with more
mean-imputation, and congruence drops for reasons that have nothing to do with erotic grammar.
Test: does each split's deficit track the between-group gap in blocks entered? Then re-run the
big splits with block-count MATCHED between groups.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore'); rng=np.random.default_rng(271828)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
pool_all=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool_all).fillna(0)
def loadings(people,K=5):
    ppl=np.array(sorted(set(people)&set(pool_all)))
    if len(ppl)<600: return None
    pm={p:i for i,p in enumerate(ppl)}; cols=[]
    for q in allq:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(ppl),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.0)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(ppl)),COV.loc[ppl].values]
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    return svd(Z,full_matrices=False)[2][:K]
def cong(a,b):
    if a is None or b is None: return np.nan
    Qa,_=qr(a.T,mode='reduced'); Qb,_=qr(b.T,mode='reduced')
    return float(np.mean(svd(Qa.T@Qb,compute_uv=False)))
T=pd.read_csv('data/derived/deficit_reference.csv')
def median_split(col):
    s=pd.to_numeric(df[col],errors='coerce'); m=s.median()
    return set(df.index[s>m]), set(df.index[s<m])
MOD=[c for c in df.columns if 'type of erotic content you prefer tends to be:' in c][0]
ANI=[c for c in df.columns if 'type of erotic content you prefer tends to be more' in c][0]
REL=[c for c in df.columns if 'preferred relationship style' in c][0]
SPL={'MODALITY written/visual':(set(df.index[df[MOD].isin(['Mostly written','Entirely written'])]),
                                set(df.index[df[MOD].isin(['Mostly visual','Entirely visual'])])),
 'ANIMATION drawn/live':(set(df.index[df[ANI].isin(['Mostly animated/drawn','Entirely animated/drawn'])]),
                         set(df.index[df[ANI].isin(['Mostly live action vid/photos','Entirely live action vid/photos'])])),
 'SEX male/female':(set(df.index[df.biomale==1]),set(df.index[df.biomale==0])),
 'OPENNESS hi/lo':median_split('opennessvariable'),'NEUROTICISM hi/lo':median_split('neuroticismvariable'),
 'EXTROVERSION hi/lo':median_split('extroversionvariable'),'POWERLESSNESS hi/lo':median_split('powerlessnessvariable'),
 'PORNHABIT hi/lo':median_split('pornhabit'),
 'MONOGAMY yes/no':(set(df.index[df[REL]=='Monogamous']),set(df.index[df[REL]=='Not monogamous']))}
gap={}
for name,(g1,g2) in SPL.items():
    a=nblk.reindex(sorted(set(g1)&set(pool_all))).mean(); b=nblk.reindex(sorted(set(g2)&set(pool_all))).mean()
    gap[name]=abs(a-b)
T['blockgap']=T.split.map(gap)
print("does the deficit just track differential coverage of the gated tree?")
print(T[['split','deficit','blockgap']].sort_values('deficit').to_string(index=False))
print(f"\n  corr(deficit, |blocks-entered gap|) = {T.deficit.corr(T.blockgap):+.3f}   (n={len(T)} splits)")

print("\n=== block-count MATCHED re-run of the two biggest splits ===")
for name in ['PORNHABIT hi/lo','SEX male/female']:
    g1,g2=SPL[name]
    g1=np.array(sorted(set(g1)&set(pool_all))); g2=np.array(sorted(set(g2)&set(pool_all)))
    # match on blocks entered by stratified subsampling
    s1=nblk.reindex(g1); s2=nblk.reindex(g2); keep1=[];keep2=[]
    for v in sorted(set(s1.unique())|set(s2.unique())):
        a=g1[s1.values==v]; b=g2[s2.values==v]; m=min(len(a),len(b))
        if m: keep1+=list(rng.choice(a,m,replace=False)); keep2+=list(rng.choice(b,m,replace=False))
    k1,k2=np.array(keep1),np.array(keep2)
    c=cong(loadings(k1),loadings(k2))
    ceil=[cong(loadings(p[:len(k1)]),loadings(p[len(k1):len(k1)+len(k2)])) for p in [rng.permutation(pool_all) for _ in range(6)]]
    raw=float(T[T.split==name].deficit.iloc[0])
    print(f"  {name}: n {len(k1):,} vs {len(k2):,} (mean blocks {nblk.reindex(k1).mean():.1f} vs {nblk.reindex(k2).mean():.1f})")
    print(f"     matched deficit = {np.nanmean(ceil)-c:.4f}   unmatched was {raw:.4f}"
          f"   -> {100*(np.nanmean(ceil)-c)/raw:.0f}% of it survives")

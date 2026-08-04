"""
E01 A07 R01 -- do the subspace-congruence comparisons survive a K sweep?

R14-R15 established that this project's CCA machinery is monotone in K everywhere: the coordinate
count rose 4->16 as K went 4->24, and the top canonical correlation rose 0.304->0.461 with it.
Three surviving headlines are computed with the SAME machinery at a fixed, never-swept K=5:
    consumption -> coordinates   0.0439  (triple-matched)
    modality written/visual      0.0334  (block-matched)
    sex, the reference           0.0927  (block-matched)
    neuroticism, the null-ish    0.0121  (triple-matched)
If K contaminates deficits the way it contaminated counts, all four move together and the
COMPARISONS are what matter, not the magnitudes.

ESTIMAND        the subspace-congruence deficit for each split, and -- the quantity the claims
                actually rest on -- the ORDERING of the four splits, as a function of K.
IDENTIFICATION  magnitudes are identified only at fixed K (R15). The ORDERING may be identified
                even when magnitudes are not; testing that is the round.
SCOPE           population BKS public n=15,503 · instrument none · baseline matched random-split
                ceiling recomputed inside every cell · regime correlations attenuated ~25%.
WORLDS          A  deficits are real properties: magnitudes may scale with K, ORDERING is stable
                   (sex > pornhabit > modality > neuroticism at every K)
                B  deficits are K artifacts: the ordering compresses or flips
                PREDICTION MATRIX
                            K small        K large
                  A         same order     same order, magnitudes scaled
                  B         some order     a different order
KILL            PRE-REGISTERED: if the rank order of the four splits changes at any K in the swept
                range, the congruence COMPARISONS are withdrawn from README.md, not just rescaled.
                If the order holds at every K, the comparisons stand and only the magnitudes are
                declared K-relative.
POSITIVE CTRL   sex must remain the largest deficit at every K -- it is the largest documented
                effect in this dataset (published d=0.62 on the pain gap). If sex is NOT largest
                somewhere, the instrument is broken at that K rather than informative.
                Fails at g=0: a split of the pool at random must give deficit ~0 at every K.
NEGATIVE CTRL   random split at matched group sizes, recomputed per cell -- this IS the ceiling.
SHAM            a coin-flip variable with the same marginal as pornhabit.
PLACEBO         split on row index parity: must be ~0 at every K.
NOISE FLOOR     ceiling sd across 5 random splits per cell.
MULTIPLICITY    6 K values x 6 splits x 3 seeds = 108 congruence computations, all reported.
SPECIFICATION   K in {2,3,5,8,12,16} x split in {sex, pornhabit, modality, neuroticism, sham,
                placebo} x seed in 3.
SEEDS           3, verified to change the random-split draw.
IMPOSSIBLE      independently replicated · causally identified · cross-dataset -- one release.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import pandas as pd, numpy as np, warnings, hashlib, itertools
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
Rt=df[rate].apply(pd.to_numeric,errors='coerce')
theta=((Rt>0).sum(1)/Rt.notna().sum(1).clip(lower=1)).reindex(pool)
def loadings(people,K):
    ppl=np.array(sorted(set(people)&set(pool)))
    if len(ppl)<500: return None
    pm={p:i for i,p in enumerate(ppl)}; cols=[]
    for q in allq:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(ppl),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(ppl)),COV.loc[ppl].values]
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    return svd(Z,full_matrices=False)[2][:K]
def cong(a,b):
    if a is None or b is None: return np.nan
    Qa,_=qr(a.T,mode='reduced'); Qb,_=qr(b.T,mode='reduced')
    return float(np.mean(svd(Qa.T@Qb,compute_uv=False)))
male=df['biomale'].reindex(pool); ph=pd.to_numeric(df['pornhabit'],errors='coerce').reindex(pool)
neu=pd.to_numeric(df['neuroticismvariable'],errors='coerce').reindex(pool)
MOD=[c for c in df.columns if 'type of erotic content you prefer tends to be:' in c][0]
mod=df[MOD].reindex(pool)
tq=pd.qcut(theta,6,labels=False,duplicates='drop')
def matched(g1,g2,rng,use_theta_sex):
    i1,i2=pd.Index(g1),pd.Index(g2)
    def key(i):
        parts=[nblk.reindex(i).astype(int)]
        if use_theta_sex: parts+=[tq.reindex(i).fillna(-1).astype(int),male.reindex(i).fillna(-1).astype(int)]
        return pd.Series(list(zip(*[p.values for p in parts])),index=i)
    k1,k2=key(i1),key(i2); a=[];b=[]
    for v in set(k1)|set(k2):
        x=k1.index[k1==v].values; y=k2.index[k2==v].values; m=min(len(x),len(y))
        if m: a+=list(rng.choice(x,m,replace=False)); b+=list(rng.choice(y,m,replace=False))
    return np.array(a),np.array(b)
rows=[]
for K,seed in itertools.product([2,3,5,8,12,16],[7,17,27]):
    rng=np.random.default_rng(seed)
    sham=pd.Series(rng.random(len(pool))<(ph>ph.median()).mean(),index=pool)
    SPL={'sex':(pool[(male==1).values],pool[(male==0).values],False),
         'pornhabit':(pool[(ph>ph.median()).values],pool[(ph<ph.median()).values],True),
         'modality':(pool[mod.isin(['Mostly written','Entirely written']).values],
                     pool[mod.isin(['Mostly visual','Entirely visual']).values],False),
         'neuroticism':(pool[(neu>neu.median()).values],pool[(neu<neu.median()).values],True),
         'sham_coinflip':(pool[sham.values],pool[~sham.values],True),
         'placebo_rowparity':(pool[np.arange(len(pool))%2==0],pool[np.arange(len(pool))%2==1],False)}
    for name,(g1,g2,mt) in SPL.items():
        a,b=matched(g1,g2,rng,mt)
        if len(a)<500: continue
        c=cong(loadings(a,K),loadings(b,K))
        ceil=[cong(loadings(p[:len(a)],K),loadings(p[len(a):len(a)+len(b)],K)) for p in [rng.permutation(pool) for _ in range(5)]]
        rows.append(dict(K=K,seed=seed,split=name,n=len(a),deficit=float(np.nanmean(ceil))-c,
                         ceiling=float(np.nanmean(ceil)),ceil_sd=float(np.nanstd(ceil))))
G=pd.DataFrame(rows); G.to_csv(OUT/'k_sweep.csv',index=False)
P=G.pivot_table(index='split',columns='K',values='deficit',aggfunc='median')
print("=== congruence deficit by K (median over 3 seeds) ===")
print(P.round(4).to_string())
print("\n=== the quantity the claims rest on: RANK ORDER at each K (1 = largest deficit) ===")
Rk=P.rank(ascending=False,axis=0).astype(int)
print(Rk.to_string())
print("\n=== POSITIVE CONTROL: is sex largest at every K? ===")
print("   ", {int(k):("YES" if Rk.loc['sex',k]==1 else f"NO (rank {Rk.loc['sex',k]})") for k in P.columns})
print("\n=== fails at g=0 ===")
print(f"    placebo (row parity) deficit across K: {P.loc['placebo_rowparity'].abs().max():.4f} max")
print(f"    sham (coin flip)     deficit across K: {P.loc['sham_coinflip'].abs().max():.4f} max")
print(f"\n  ceiling by K: {P.columns.tolist()} -> {[round(float(G[G.K==k].ceiling.median()),3) for k in P.columns]}")
print(f"  ceiling sd (noise floor): {G.ceil_sd.median():.4f}")
real=['sex','pornhabit','modality','neuroticism']
orders=[tuple(P[k].reindex(real).sort_values(ascending=False).index) for k in P.columns]
stable=len(set(orders))==1
print(f"\nMULTIPLICITY  {len(G)} congruence computations (6 K x 6 splits x 3 seeds), all reported")
print(f"  distinct orderings of the four real splits across K: {len(set(orders))}")
for k,o in zip(P.columns,orders): print(f"    K={k:2d}: {' > '.join(o)}")
print("\nPRE-REGISTERED KILL, evaluated:")
print("  -> ordering STABLE at every K : comparisons stand, magnitudes are K-relative" if stable
      else "  -> ordering CHANGES across K : congruence COMPARISONS WITHDRAWN")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")

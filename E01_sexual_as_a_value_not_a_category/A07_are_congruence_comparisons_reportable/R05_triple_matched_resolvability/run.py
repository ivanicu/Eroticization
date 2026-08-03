"""
E01 A07 R05 -- is the consumption headline reportable at the matching level it was published at?

A07 was opened mid-project and has rounds but no decision statement, which the E/A/R structure
forbids: an arc that cannot name the decision it made safe should not have been opened. The
measurement needed to write one is missing.

#11 published consumption -> coordinates at 0.0439, TRIPLE-matched (block count + theta + sex).
#35 measured resolvability for consumption at 0.0951 -- but that is the BLOCK-ONLY matched value.
The published number has never had its seed spread measured, and it is half the size, which is
exactly where #34 found modality failing.

ESTIMAND        the consumption congruence deficit at each matching level, with its across-seed
                spread, and whether |effect| > 2 x spread at the published level.
IDENTIFICATION  identified; this re-runs the published computation with the seed varied.
WORLDS          A  reportable: ratio > 2 at triple matching
                B  unreportable like modality: matching shrinks the effect below its own noise
KILL (CONDITIONAL) gate: sex must be resolvable at the SAME matching level (positive control) AND
                   the placebo must be unresolvable. Otherwise UNVERIFIED.
                   then: consumption ratio > 2 at triple matching -> A07 closes with a reportable
                         comparison; ratio <= 2 -> the published 0.0439 joins modality as
                         unresolvable and A07 closes with almost nothing reportable
POSITIVE CTRL   sex at the same matching level -- and note that sex-matching makes the sex split
                degenerate, so the positive control is sex at BLOCK+THETA matching, stated rather
                than silently substituted.
NEGATIVE CTRL   row-parity placebo at every matching level.
NOISE FLOOR     across-seed spread, 6 seeds.
MULTIPLICITY    3 matching levels x 3 splits x 6 seeds, all reported.
SEEDS           6.
IMPOSSIBLE      sex at triple matching -- matching on sex removes the split. Reported as N/A.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
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
tq=pd.qcut(theta,6,labels=False,duplicates='drop')
male=df['biomale'].reindex(pool); ph=pd.to_numeric(df['pornhabit'],errors='coerce').reindex(pool)
def loadings(people,K=5):
    ppl=np.array(sorted(set(people)&set(pool)))
    if len(ppl)<400: return None
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
def key(i,level):
    parts=[nblk.reindex(i).astype(int)]
    if level in ('block+theta','triple'): parts.append(tq.reindex(i).fillna(-1).astype(int))
    if level=='triple': parts.append(male.reindex(i).fillna(-1).astype(int))
    return pd.Series(list(zip(*[p.values for p in parts])),index=i)
def deficit(split,level,seed):
    rng=np.random.default_rng(seed)
    if split=='consumption': g1,g2=pool[(ph>ph.median()).values],pool[(ph<ph.median()).values]
    elif split=='sex': g1,g2=pool[(male==1).values],pool[(male==0).values]
    else: g1,g2=pool[np.arange(len(pool))%2==0],pool[np.arange(len(pool))%2==1]
    i1,i2=pd.Index(sorted(set(g1))),pd.Index(sorted(set(g2)))
    k1,k2=key(i1,level),key(i2,level); a=[];b=[]
    for v in set(k1)|set(k2):
        x=k1.index[k1==v].values; y=k2.index[k2==v].values; m=min(len(x),len(y))
        if m: a+=list(rng.choice(x,m,replace=False)); b+=list(rng.choice(y,m,replace=False))
    a,b=np.array(a),np.array(b)
    if min(len(a),len(b))<400: return np.nan
    c=cong(loadings(a),loadings(b))
    ceil=[cong(loadings(p[:len(a)]),loadings(p[len(a):len(a)+len(b)])) for p in [rng.permutation(pool) for _ in range(3)]]
    return float(np.nanmean(ceil))-c
rows=[]
for split,level,seed in itertools.product(['consumption','sex','placebo'],
                                          ['block','block+theta','triple'],[1,2,3,4,5,6]):
    if split=='sex' and level=='triple': continue     # degenerate: matching on sex removes the split
    d=deficit(split,level,seed)
    if np.isfinite(d): rows.append(dict(split=split,level=level,seed=seed,deficit=d))
G=pd.DataFrame(rows); G.to_csv(OUT/'triple_matched.csv',index=False)
S=G.groupby(['split','level']).deficit.agg(['median','min','max'])
S['spread']=S['max']-S['min']; S['ratio']=S['median'].abs()/S['spread'].clip(lower=1e-9)
S['resolvable']=S['median'].abs()>2*S['spread']
print("=== deficit, seed spread and resolvability by matching level (6 seeds) ===")
print(S.round(4).to_string())
print("\n  N/A: sex at triple matching -- matching on sex removes the split")
try:
    cons=S.loc[('consumption','triple')]; sexc=S.loc[('sex','block+theta')]; plc=S.loc[('placebo','triple')]
except KeyError:
    print("  cells missing"); sys.exit()
gate_pos=bool(sexc['resolvable']); gate_neg=bool(not plc['resolvable'])
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  sex resolvable at block+theta (pos ctrl) : {'PASS' if gate_pos else 'FAIL'} (ratio {sexc['ratio']:.2f})")
print(f"  placebo NOT resolvable at triple         : {'PASS' if gate_neg else 'FAIL'} (ratio {plc['ratio']:.2f})")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED")
elif cons['ratio']>2: print(f"  -> REPORTABLE : consumption at triple matching {cons['median']:.4f} +/- {cons['spread']:.4f}, ratio {cons['ratio']:.2f}")
else: print(f"  -> UNRESOLVABLE : consumption {cons['median']:.4f} against spread {cons['spread']:.4f}, ratio {cons['ratio']:.2f} -- joins modality")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")

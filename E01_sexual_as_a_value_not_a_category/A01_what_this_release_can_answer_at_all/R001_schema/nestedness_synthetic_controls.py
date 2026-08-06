"""
E01 A05 R13 -- validate the nestedness measure with the control that killed the concentration one.

#46 killed "breadth has no shape" by generating synthetic populations with KNOWN answers: a
base-rate population showed -1.6% concentration where it should have shown 0, and a genuinely
concentrated population showed the same -1.6%. The measure could not tell them apart, and its
published value sat inside its own Jensen bias.

Nestedness is now one of A05's three remaining supports (published: containment 0.7938 against a
size-matched base-rate null of 0.7278, excess +0.0660 = 24.2% of the distance from chance to
perfect). It has never had the same treatment. It differs structurally from the concentration
measure -- one null draw per pair rather than an average of six -- so the Jensen mechanism should
NOT apply, and this round is a test of that reasoning rather than a fishing expedition.

ESTIMAND        the containment excess over a size-matched base-rate null, on three populations
                whose true answer is known or observed.
IDENTIFICATION  identified: the synthetic populations have their nestedness set by construction.
WORLDS          A  the measure is sound: base-rate synthetic ~0, nested synthetic large, real +0.066
                B  biased like the concentration measure: base-rate synthetic departs from 0
KILL (CONDITIONAL) gate: the NESTED synthetic must show a clearly larger excess than the base-rate
                   synthetic; otherwise the measure cannot grade nestedness and no verdict follows.
                   then: |base-rate synthetic excess| < 0.01 -> measure SOUND, the 24.2% stands
                         |base-rate synthetic excess| > 0.03 -> BIASED like #46, support withdrawn
                         otherwise -> UNVERIFIED
POSITIVE CTRL   a synthetically NESTED population: every person's set is a prefix of a fixed
                popularity ordering, so containment is perfect by construction.
NEGATIVE CTRL   a base-rate population: sets drawn from base rates at the observed sizes.
NOISE FLOOR     across-seed spread, 5 seeds.
MULTIPLICITY    3 populations x 5 seeds, all reported.
SEEDS           5.
IMPOSSIBLE      nothing -- the controls have analytic answers, which is why this check is cheap.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
H=(R>0).astype(float).where(R.notna()); H=H[H.notna().sum(1)>=40]
V=(H.fillna(0).values>0.5); K=V.shape[1]; base=V.mean(0); p=base/base.sum(); sizes=V.sum(1)
print(f"people {len(V):,}  categories {K}  median set size {int(np.median(sizes))}")
order=np.argsort(-base)                      # popularity ordering, for the nested control
def make(kind,rng):
    M=np.zeros_like(V)
    for i,k_ in enumerate(sizes):
        k_=int(min(k_,K))
        if kind=='nested': M[i,order[:k_]]=True                     # perfect Guttman by construction
        elif kind=='baserate': M[i,rng.choice(K,k_,replace=False,p=p)]=True
    return M
def nestedness(M,rng,npairs=15000):
    sz=M.sum(1); n=len(M); obs=[];nul=[]
    def rnd(k_):
        m=np.zeros(K,bool); m[rng.choice(K,min(int(k_),K),replace=False,p=p)]=True; return m
    tries=0
    while len(obs)<npairs and tries<npairs*4:
        tries+=1
        i,j=rng.integers(0,n,2)
        if sz[i]<=sz[j] or sz[j]<5: continue
        obs.append((M[i]&M[j]).sum()/sz[j])
        a,b=rnd(sz[i]),rnd(sz[j])
        nul.append((a&b).sum()/max(b.sum(),1))
    return float(np.mean(obs)),float(np.mean(nul))
rows=[]
for seed in (1,2,3,4,5):
    rng=np.random.default_rng(seed)
    for name,M in [('real',V),('NESTED [pos ctrl]',make('nested',rng)),('base-rate [neg ctrl]',make('baserate',rng))]:
        o,n=nestedness(M,rng)
        rows.append(dict(pop=name,seed=seed,observed=o,null=n,excess=o-n,
                         pct_of_gap=100*(o-n)/max(1-n,1e-9)))
G=pd.DataFrame(rows); G.to_csv(OUT/'nestedness_controls.csv',index=False)
S=G.groupby('pop')[['observed','null','excess','pct_of_gap']].agg(['median','min','max'])
print("\n=== containment, null, and excess by population (5 seeds) ===")
print(G.groupby('pop')[['observed','null','excess','pct_of_gap']].median().round(4).to_string())
print("\n  seed spread of the excess:")
print(G.groupby('pop').excess.agg(lambda s: s.max()-s.min()).round(4).to_string())
med=G.groupby('pop').excess.median()
nested=float(med['NESTED [pos ctrl]']); baser=float(med['base-rate [neg ctrl]']); real=float(med['real'])
gate=nested>baser+0.05
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  NESTED synthetic clearly exceeds base-rate synthetic : {'PASS' if gate else 'FAIL'} ({nested:+.4f} vs {baser:+.4f})")
if not gate: print("  -> gate FAILED : the measure cannot grade nestedness, UNVERIFIED")
elif abs(baser)<0.01: print(f"  -> MEASURE SOUND : base-rate synthetic excess {baser:+.4f}; the real excess {real:+.4f} (published +0.0660) stands")
elif abs(baser)>0.03: print(f"  -> BIASED like #46 : base-rate synthetic excess {baser:+.4f}; support withdrawn")
else: print(f"  -> UNVERIFIED : base-rate synthetic excess {baser:+.4f}")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")

"""
E01 A08 R05 -- sort the surviving headlines into the two classes #49 and #50 created.

#49: the cross-domain grammar exists as a correlation (CCA 0.198) and carries NO predictable
variance (pairwise R2 -0.002). #50: the onset RSA (+0.599) DOES carry it (31% of its same-domain
ceiling). Every remaining headline is a correlation or a congruence and none has been sorted.

Common pipeline for each claim: the published statistic, the held-out predictive R2, the SAME-DOMAIN
CEILING measured in the same round (per #50's procedural fix -- compute the ceiling first, then
express everything as a fraction of it), and the permuted null.

  nestedness       published: containment excess +0.0655, 24.0% of chance->perfect
                   prediction: can population base rates alone predict WHICH categories a person
                   endorses, held out?  ceiling: another random person's set with the same count.
  sex deficit      published: congruence deficit 0.093
                   prediction: does sex predict a person's within-block profile?
                   ceiling: profile-half -> profile-half.

ESTIMAND        per claim, held-out predictive R2 as a FRACTION of its own same-domain ceiling.
IDENTIFICATION  identified; every quantity is out-of-sample.
WORLDS          A  the claim's structure travels: R2 a substantial fraction of the ceiling
                B  thin direction like #49: R2 at the null despite a large published statistic
KILL (CONDITIONAL) gate -- CEILING FIRST, per #50: each claim's ceiling must exceed its own null by
                   3x. A claim whose ceiling is at its null cannot be sorted and is reported as
                   UNSORTABLE rather than assigned a class.
                   then: R2 > 25% of ceiling -> CARRIES PREDICTABLE VARIANCE
                         R2 <= null           -> THIN DIRECTION
                         between              -> partial, reported as the fraction
POSITIVE CTRL   each claim's own same-domain ceiling.
NEGATIVE CTRL   permuted persons, per claim.
SEEDS           5.
MULTIPLICITY    2 claims x 3 quantities x 5 seeds, all reported, plus the two already sorted.
IMPOSSIBLE      the maturational schedule -- its statistic is a within-person RANK agreement whose
                predictive analogue is a ranking task, not an R2, so it is not commensurable with
                this pipeline. Listed UNSORTED rather than forced.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq, svd
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
H=(R>0).astype(float).where(R.notna()); H=H[H.notna().sum(1)>=40]
V=H.fillna(0).values>0.5
def r2(X,Y,rng,permute=False):
    n=len(X); idx=rng.permutation(n); tr,te=idx[:int(.7*n)],idx[int(.7*n):]
    Yt=Y[rng.permutation(n)] if permute else Y
    A=np.c_[np.ones(len(tr)),X[tr]]; b,*_=lstsq(A,Yt[tr],rcond=None)
    P=np.c_[np.ones(len(te)),X[te]]@b
    return 1-((Yt[te]-P)**2).sum()/max((Yt[te]**2).sum(),1e-12)
rows=[]
for seed in (1,2,3,4,5):
    rng=np.random.default_rng(seed)
    # --- NESTEDNESS: does the population ordering predict which categories a person holds? ---
    Y=V.astype(float); Yc=Y-Y.mean(0)
    base=np.tile(Y.mean(0),(len(Y),1))                      # population base rates, same for all
    cnt=Y.sum(1,keepdims=True)
    Xpop=np.c_[base[:,:1]*0+cnt, cnt**2]                    # only person-level info: their count
    cols=rng.permutation(Y.shape[1]); h=Y.shape[1]//2
    rows.append(dict(seed=seed,claim='nestedness',
        pred=r2(Xpop,Yc,rng), ceiling=r2(Yc[:,cols[:h]],Yc[:,cols[h:2*h]],rng),
        null=r2(Xpop,Yc,rng,permute=True)))
    # --- SEX DEFICIT: does sex predict a person's profile? ---
    idx=H.index
    sex=pd.to_numeric(df['biomale'],errors='coerce').reindex(idx).fillna(0.5).values.reshape(-1,1)
    rows.append(dict(seed=seed,claim='sex deficit',
        pred=r2(sex,Yc,rng), ceiling=r2(Yc[:,cols[:h]],Yc[:,cols[h:2*h]],rng),
        null=r2(sex,Yc,rng,permute=True)))
G=pd.DataFrame(rows); G.to_csv(OUT/'sort.csv',index=False)
S=G.groupby('claim')[['pred','ceiling','null']].median().round(4)
S['pct_of_ceiling']=(100*S['pred']/S['ceiling']).round(1)
S['ceiling_over_null']=(S['ceiling']/S['null'].abs().clip(lower=1e-9)).round(1)
print("=== CEILING FIRST (per #50), then everything as a fraction of it ===")
print(S.to_string())
print("\n=== the two already sorted, for the table ===")
print("  cross-domain CCA (#49)   published 0.198   prediction -0.0022 vs null -0.0110   THIN DIRECTION")
print("  onset RSA        (#50)   published 0.599   prediction +0.0136, 31% of ceiling   CARRIES VARIANCE")
print("\n  UNSORTED: the maturational schedule -- a within-person RANK agreement, not commensurable")
print("            with an R2 pipeline. Listed rather than forced.")
for c in S.index:
    r=S.loc[c]
    gate=r['ceiling_over_null']>3
    tag='UNSORTABLE (ceiling at its null)' if not gate else (
        'CARRIES PREDICTABLE VARIANCE' if r['pct_of_ceiling']>25 else
        ('THIN DIRECTION' if r['pred']<=abs(r['null']) else f"partial, {r['pct_of_ceiling']:.0f}% of ceiling"))
    print(f"\n  {c:14s} ceiling {r['ceiling']:+.4f} ({r['ceiling_over_null']:.1f}x its null) · prediction {r['pred']:+.4f} -> {tag}")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")

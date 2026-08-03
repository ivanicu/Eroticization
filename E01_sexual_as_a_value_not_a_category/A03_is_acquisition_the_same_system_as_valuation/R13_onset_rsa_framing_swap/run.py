"""
E01 A03 R13 -- the onset RSA under the framing that broke the central claim.

#49 showed the cross-domain grammar exists as a CORRELATION (CCA 0.198, tier-ordered) and carries
almost no PREDICTABLE VARIANCE (pairwise R2 -0.002). Every number was right; the framing answered
"is there a shared direction" while I read it as "does A tell you about B".

The onset RSA has the same structure: a correlation between two similarity matrices, +0.599 after
stripping intensity leakage and recall anchoring, z about 11 against a category-permutation null. It
has never been asked the predictable-variance question, and it is the load-bearing evidence for
A03's "acquired together tracks liked together".

  published framing : corr(onset-similarity matrix, preference-similarity matrix) over 344 pairs
  this framing      : PERSON-LEVEL prediction. Can a held-out person's PREFERENCE profile be
                      predicted from their ONSET profile? Different unit (person, not pair),
                      different estimand family, different failure modes.

ESTIMAND        held-out R2 predicting a person's 27-category preference residuals from their
                27-category onset residuals.
IDENTIFICATION  identified; out-of-sample by construction.
WORLDS          A  the RSA reflects structure that travels: person-level prediction works
                B  correlation in a thin direction, like #49: prediction ~ 0 despite RSA 0.599
KILL (CONDITIONAL) gate: the positive control (predict preference from a DISJOINT HALF of the
                   preference columns) must work; otherwise the pipeline cannot predict profiles at
                   all and no comparison follows.
                   then: onset->preference R2 > 3x |null| -> the RSA survives the framing swap
                         R2 at the null                    -> same shape as #49, and A03's evidence
                                                              is a correlation in a thin direction
POSITIVE CTRL   preference half -> preference half, same pipeline.
NEGATIVE CTRL   permuted persons.
SEEDS           5.
MULTIPLICITY    3 predictor sets x 5 seeds, all reported.
IMPOSSIBLE      a person-level ground truth for "acquisition structure" -- onset is self-report.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('24_attack_rsa.py')).read().split('print("=== is onset a proxy')[0])
# Ores / Pres are the double-centred, age-adjusted onset and preference residual frames
On=Ores.copy(); Pr=Pres.copy()
keep=(On.notna().sum(1)>=8)&(Pr.notna().sum(1)>=8)
On=On[keep]; Pr=Pr[keep]
On=On.fillna(0.0); Pr=Pr.fillna(0.0)
print(f"people with >=8 onsets and >=8 preferences: {len(On):,}   categories: {On.shape[1]}")
def r2(X,Y,rng,permute=False):
    n=len(X); idx=rng.permutation(n); tr,te=idx[:int(.7*n)],idx[int(.7*n):]
    Yt=Y[rng.permutation(n)] if permute else Y
    A=np.c_[np.ones(len(tr)),X[tr]]
    b,*_=lstsq(A,Yt[tr],rcond=None)
    P=np.c_[np.ones(len(te)),X[te]]@b
    return 1-((Yt[te]-P)**2).sum()/max((Yt[te]**2).sum(),1e-12)
Xo=On.values; Yp=Pr.values
rows=[]
for seed in (1,2,3,4,5):
    rng=np.random.default_rng(seed)
    cols=rng.permutation(Yp.shape[1]); h=Yp.shape[1]//2
    rows.append(dict(seed=seed,
        onset_to_pref=r2(Xo,Yp,rng),
        onset_to_pref_null=r2(Xo,Yp,rng,permute=True),
        pref_to_pref_POSCTRL=r2(Yp[:,cols[:h]],Yp[:,cols[h:2*h]],rng),
        onset_to_onset=r2(Xo[:,cols[:h]] if h<=Xo.shape[1] else Xo, Xo[:,cols[h:2*h]] if 2*h<=Xo.shape[1] else Xo, rng)))
G=pd.DataFrame(rows); G.to_csv(OUT/'framing_swap.csv',index=False)
M=G.drop(columns='seed').median()
S=G.drop(columns='seed').agg(lambda s: s.max()-s.min())
print("\n=== held-out person-level prediction R2 (5 seeds) ===")
for k in M.index: print(f"  {k:24s} {M[k]:+.4f}   seed spread {S[k]:.4f}")
print(f"\n  published RSA framing: corr(onset-similarity, preference-similarity) = +0.599, z about 11")
gate=M['pref_to_pref_POSCTRL']>0.05
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  positive control: preference half -> preference half works : {'PASS' if gate else 'FAIL'} ({M['pref_to_pref_POSCTRL']:+.4f})")
if not gate: print("  -> gate FAILED : the pipeline cannot predict profiles at all, UNVERIFIED")
else:
    eff=M['onset_to_pref']; nul=abs(M['onset_to_pref_null'])
    print(f"  onset -> preference {eff:+.4f} against |null| {nul:.4f}")
    if eff>3*nul and eff>0.01: print("  -> SURVIVES THE FRAMING SWAP : onset structure predicts preference at the person level")
    else: print("  -> SAME SHAPE AS #49 : the RSA is a correlation in a direction that carries little predictable variance")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")

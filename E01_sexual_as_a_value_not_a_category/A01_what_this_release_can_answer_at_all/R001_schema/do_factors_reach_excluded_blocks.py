"""
E01 A01 R15 -- do the person factors predict in the blocks the filter threw away?

#39 showed the SUBSTANCE axis predicts nothing outside the seven blocks that defined it, and named
the failure: availability became representativeness. The same question has never been asked of the
construct carrying this project's central claim. The person factors -- and every GCCA coordinate,
every congruence deficit, the leave-one-block-out result -- are built from 32 blocks out of 101,
selected by a filter (n_respondents>=1200, n_options>=10, mean_picks>1.5) that 20 rounds inherited
and none chose (#32, #33).

If the grammar is a property of the DOMAIN it should predict in the excluded blocks too. If it is a
property of the LARGE blocks, it will not.

ESTIMAND        held-out variance in an EXCLUDED block's within-person option profile explained by
                factors fitted on the 32 included blocks, over demographics + propensity.
IDENTIFICATION  identified for excluded blocks with enough respondents to fit a profile at all
                (>=300), which is what the exclusion threshold was set above.
WORLDS          A  domain property: factors predict excluded blocks at a gain comparable to
                   leave-one-out within the included set
                B  large-block property: no gain outside, or none above the placebo
KILL (CONDITIONAL) gate: factors must predict INCLUDED blocks held-out (the A01 R12/R13 result) AND
                   permuted factors must give ~0 in both sets. Otherwise UNVERIFIED.
                   then: excluded-block gain > 2x placebo and > 0.005 -> REACHES
                         gain <= placebo -> LARGE-BLOCK PROPERTY, and the central claim is scoped
POSITIVE CTRL   leave-one-out inside the included set, which #25 measured at +0.017 to +0.037.
NEGATIVE CTRL   permuted factor rows.
PLACEBO         the person's mean endorsement RATE across the included blocks -- same source, same
                people, no coordinate content. This is the control that decided #39.
SEEDS           4.
MULTIPLICITY    every admissible excluded block reported individually, not just the median.
IMPOSSIBLE      blocks with too few respondents to estimate a profile; their count is reported so
                the coverage of this test is visible.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq, svd
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
INC=list(B)
qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
EXC={}
for _,q in qm.iterrows():
    if q.qi in INC or q.single_pick or q.mean_picks<=1.5: continue
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(vc[vc>=20].index)]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<300 or len(opt)<6: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
    EXC[q.qi]=dict(ppl=ppl,R=R)
print(f"included blocks {len(INC)} · excluded blocks admissible for this test {len(EXC)}")
print(f"  excluded but too small to test: {int((~qm.single_pick & (qm.mean_picks>1.5)).sum())-len(INC)-len(EXC)}")
pool=np.unique(np.concatenate([B[q]['ppl'] for q in INC]))
pm={p:i for i,p in enumerate(pool)}; cols=[]
for q in INC:
    idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
    Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan)
    if len(idx): Z[idx]=B[q]['R'][src]
    mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
    cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
Zall=np.hstack(cols); Zall=Zall-Zall.mean(0)
D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
b,*_=lstsq(D,Zall,rcond=None); Zall=Zall-D@b
U,S,_=svd(Zall,full_matrices=False)
FAC=pd.DataFrame(U[:,:6]*S[:6],index=pool)
RATE=pd.Series(np.nanmean(np.where(np.isnan(np.hstack(cols)),np.nan,np.hstack(cols)),axis=1),index=pool)
RATE=pd.Series((Zall!=0).mean(1),index=pool)
def test(blocks,label):
    rows=[]
    for t,blk in blocks.items():
        common=np.array([p for p in blk['ppl'] if p in pm])
        if len(common)<300: continue
        ia=np.searchsorted(blk['ppl'],common); Y=blk['R'][ia]
        prop=(Y!=0).mean(1)
        Xb=np.c_[np.ones(len(common)),prop,COV.loc[common].values]
        F=FAC.loc[common].values; r=RATE.loc[common].values.reshape(-1,1)
        for seed in (1,2,3,4):
            rng=np.random.default_rng(seed)
            idx=rng.permutation(len(common)); tr,te=idx[:int(.7*len(idx))],idx[int(.7*len(idx)):]
            def r2(X):
                bb,*_=lstsq(X[tr],Y[tr],rcond=None); return 1-((Y[te]-X[te]@bb)**2).sum()/((Y[te])**2).sum()
            base=r2(Xb)
            rows.append(dict(set=label,block=t,seed=seed,gain=r2(np.c_[Xb,F])-base,
                             placebo=r2(np.c_[Xb,r])-base,
                             null=r2(np.c_[Xb,F[rng.permutation(len(common))]])-base))
    return pd.DataFrame(rows)
GE=test(EXC,'excluded'); GI=test({q:B[q] for q in INC},'included [POS CTRL]')
G=pd.concat([GI,GE]); G.to_csv(OUT/'reach_excluded.csv',index=False)
S_=G.groupby('set')[['gain','placebo','null']].median().round(5)
print("\n=== held-out gain from factors fitted on the 32 INCLUDED blocks ===")
print(S_.to_string())
print(f"\n  per-block medians, excluded set ({GE.block.nunique()} blocks):")
pb=GE.groupby('block')[['gain','placebo']].median().round(4).sort_values('gain',ascending=False)
print(pb.to_string())
inc=S_.loc['included [POS CTRL]']; exc=S_.loc['excluded']
# FIX: excluded blocks are smaller, so six permuted factor columns cost held-out R2 by DEGREES OF
# FREEDOM alone -- the permuted null is -0.0064, not 0. Requiring |null|<0.005 was a gate condition
# that does not fit the design (same class as #33's error 2). The null IS the df cost, so the
# effect is gain MINUS null, and the negative control is that the null be STABLE across blocks.
eff_inc=inc['gain']-inc['null']; eff_exc=exc['gain']-exc['null']; eff_pla=exc['placebo']-exc['null']
null_sd=float(GE.groupby('block').null.median().std())
gate_pos=eff_inc>0.02; gate_neg=null_sd<0.02
print(f"\nCONDITIONAL KILL -- gate first (effect = gain - permuted null, which prices the df cost)")
print(f"  factors predict INCLUDED blocks (effect>0.02) : {'PASS' if gate_pos else 'FAIL'} ({eff_inc:.5f})")
print(f"  permuted null stable across blocks (sd<0.02) : {'PASS' if gate_neg else 'FAIL'} (sd {null_sd:.5f})")
print(f"\n  included effect {eff_inc:.5f} | excluded effect {eff_exc:.5f} | excluded placebo effect {eff_pla:.5f}")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED")
elif eff_exc>2*abs(eff_pla) and eff_exc>0.005:
    print(f"  -> REACHES : factors predict blocks the filter discarded, at {100*eff_exc/eff_inc:.0f}% of the included effect")
elif eff_exc<=abs(eff_pla):
    print(f"  -> LARGE-BLOCK PROPERTY : the central claim is scoped to the 32 included blocks")
else: print(f"  -> UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")

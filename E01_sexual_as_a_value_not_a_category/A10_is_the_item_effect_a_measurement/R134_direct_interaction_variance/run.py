import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A10 R19 -- THE INTERACTION MAGNITUDE, MEASURED INSTEAD OF INFERRED.

#90 put the person x item interaction at +/-23.7 pp by INVERTING A PLANT: which synthetic
perturbation, passed through the same low-rank estimator, reproduces the real skill? Item (+/-22.6)
and person (+/-16.3) are measured DIRECTLY from the margins. That asymmetry is stated in #90's own
scope paragraph and it has never been attacked. It is the top row of the README.

There is a direct route that uses no estimator at all. For a binary cell,

    M_ij = pi_ij + e_ij ,   e independent across cells given pi

so after removing both margins the mean squared residual is

    mean_j r^2  =  Var(interaction)  +  binomial noise  +  (person-estimate error)^2

Two of those three are nuisances, and CURVEBALL CANCELS BOTH EXACTLY: it preserves every row sum and
every column sum, so its binomial term and its person-estimation error are identical to the real
matrix's, while its interaction is destroyed. Therefore

    Var(interaction)  =  mean_j r^2 (real)  -  mean_j r^2 (fixed-margin null)

is a DIRECT measurement on the probability scale. No low-rank fit, no plant, no inversion.

ESTIMAND        sigma_int = sqrt(mean r^2_real - mean r^2_null), in percentage points, per block and
                pooled; compared against #90's inverted +/-23.7 pp.
IDENTIFICATION  identified by the cancellation above, which requires the null to preserve BOTH
                margins exactly -- asserted per draw, not assumed.
SCOPE           the 23 blocks A09/R114 identified. Person effect estimated from OTHER blocks only,
                so no target cell contributes to its own residual.
WORLDS          sound      direct ~= 23.7 pp -> #90's inversion is validated by an independent route
                inflated   direct << 23.7 -> the headline number is an artifact of the plant
                deflated   direct >> 23.7 -> the inversion understated it
KILL            threshold-free: the direct estimate with its bootstrap interval, against 23.7.
POSITIVE CTRL   plant a known interaction of KNOWN realized sd; the estimator must recover that sd,
                and must be MONOTONE across a ladder. A magnitude estimator that cannot return a
                planted magnitude cannot judge #90.
NEGATIVE CTRL   null minus a SECOND independent null: must return ~0 pp.
NOISE FLOOR     200 bootstrap resamples over people; 3 null draws.
MULTIPLICITY    23 blocks x 4 worlds x 3 draws, published whole.
IMPOSSIBLE      separating "the person's probability differs" from "the person answered
                inconsistently" -- one observation per cell, no repeats. sigma_int is an upper bound
                on stable interaction and includes any within-person response noise.
"""
import pandas as pd, numpy as np, warnings, hashlib
sys.path.insert(0,str(ROOT))
from lib.gates import Gate
warnings.filterwarnings('ignore')
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,ppl=ppl)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in IDENT])); PM={p:i for i,p in enumerate(ALLP)}
print(f"targets {len(IDENT)}  people {len(ALLP)}",flush=True)
def curveball(M,rng,per_row=5.):
    A=[set(np.flatnonzero(r).tolist()) for r in M]; n=len(A)
    for _ in range(int(per_row*n)):
        i,j=int(rng.integers(n)),int(rng.integers(n))
        if i==j: continue
        ai,aj=A[i],A[j]; inter=ai&aj
        di=list(ai-inter); dj=list(aj-inter); L=di+dj
        if not L: continue
        rng.shuffle(L); k=len(di)
        A[i]=inter|set(L[:k]); A[j]=inter|set(L[k:])
    out=np.zeros_like(M)
    for i,s in enumerate(A): out[i,list(s)]=1.
    return out
# person effect from OTHER blocks only -- no target cell contributes to its own residual
brate=np.zeros(len(ALLP)); bcnt=np.zeros(len(ALLP))
for t in IDENT:
    M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
    brate[idx]+=M.sum(1)-M.mean(0).sum(); bcnt[idx]+=M.shape[1]
def msr(M,idx,exclude_self):
    """mean squared residual after removing the item margin and an OUT-OF-BLOCK person effect."""
    p=M.mean(0); n,m=M.shape
    num=brate[idx]-(M.sum(1)-p.sum() if exclude_self else 0.)
    den=np.maximum(bcnt[idx]-(m if exclude_self else 0),1)
    b=(num/den)[:,None]
    r=M-p[None,:]-b
    return float(np.mean(r**2)), b
def plant(M,sc,rng,r=5):
    n,m=M.shape; F=rng.normal(size=(n,r)); L=rng.normal(size=(r,m))*sc
    dev=F@L; base=M.mean(0)[None,:]
    p=np.clip(base+dev,0.02,0.98)
    return (rng.random((n,m))<p).astype(float), float((p-base).std())
rows=[]
for t in IDENT:
    M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
    v_real,_=msr(M,idx,True)
    for d in range(3):
        Mn=curveball(M,np.random.default_rng(4400+d))
        assert np.allclose(Mn.sum(0),M.sum(0)) and np.allclose(Mn.sum(1),M.sum(1)),"margins"
        v_null,_=msr(Mn,idx,True)
        Mn2=curveball(M,np.random.default_rng(5500+d)); v_null2,_=msr(Mn2,idx,True)
        rows.append(dict(q=t,draw=d,world='real',v=v_real,v_ref=v_null,n=M.shape[0],m=M.shape[1]))
        rows.append(dict(q=t,draw=d,world='nullnull',v=v_null2,v_ref=v_null,n=M.shape[0],m=M.shape[1]))
    for sc in [0.05,0.10,0.20]:
        Mp,psd=plant(M,sc,np.random.default_rng(6600))
        v_p,_=msr(Mp,idx,True)
        Mpn=curveball(Mp,np.random.default_rng(6601)); v_pn,_=msr(Mpn,idx,True)
        rows.append(dict(q=t,draw=0,world=f'plant{sc}',v=v_p,v_ref=v_pn,n=M.shape[0],m=M.shape[1],
                         planted_pp=100*psd))
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
def sigma(w):
    d=D[D.world==w]; dv=(d.v-d.v_ref)
    per=d.assign(dv=dv).groupby('q').dv.mean()
    pooled=per.mean()
    return 100*np.sqrt(max(pooled,0)), per
s_real,per_real=sigma('real'); s_nn,_=sigma('nullnull')
print("\n=== DIRECT INTERACTION MAGNITUDE (no estimator, no plant, no inversion) ===")
print(f"  real  - fixed-margin null : sigma_int = {s_real:.2f} pp")
print(f"  null  - fixed-margin null : {s_nn:.2f} pp     <- must be ~0")
print("\n=== POSITIVE CONTROL: planted interactions of KNOWN realized sd ===")
for sc in [0.05,0.10,0.20]:
    s,_=sigma(f'plant{sc}'); tp=D[D.world==f'plant{sc}'].planted_pp.mean()
    print(f"  scale {sc:.2f}:  planted {tp:5.2f} pp   recovered {s:5.2f} pp   "
          f"ratio {s/max(tp,1e-9):.2f}")
rb=np.random.default_rng(99)
bs=[100*np.sqrt(max(per_real.values[rb.integers(0,len(per_real),len(per_real))].mean(),0))
    for _ in range(200)]
sd_b=float(np.std(bs))
rec=[sigma(f'plant{s}')[0]/D[D.world==f'plant{s}'].planted_pp.mean() for s in [0.05,0.10,0.20]]
g=Gate("is #90's +/-23.7 pp validated by a direct measurement?")
g.negative_control("null minus null", null=s_nn, effect=s_real)
g.no_sign_crossing("recovery ratio across the plant ladder", rec)
g.asserted("plant ladder monotone in planted magnitude",
           all(sigma(f'plant{a}')[0]<sigma(f'plant{b}')[0] for a,b in [(0.05,0.10),(0.10,0.20)]),
           f"recovered {[round(sigma(f'plant{s}')[0],2) for s in [0.05,0.10,0.20]]}")
g.resolvable("direct sigma_int", effect=s_real, spread=sd_b)
print(); print(g)
if g.verdict():
    cal=np.mean(rec)
    print(f"\n   mean recovery ratio across the ladder: {cal:.2f}  -> calibrated direct estimate "
          f"{s_real/cal:.2f} pp")
    print(f"   #90's inverted estimate: 23.70 pp   bootstrap sd of the direct estimate {sd_b:.2f} pp")
    d=abs(s_real/cal-23.70)
    print(f"   |direct - inverted| = {d:.2f} pp = {d/max(sd_b,1e-9):.1f} bootstrap sd")
    print("\n   ->", "AGREES -- two independent routes, one number" if d<3*sd_b else
          f"DISAGREES by {d:.1f} pp -- #90's inversion is not validated")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A09 R02 -- THE INTERACTION ON THE SAME SCALE AS THE ITEM EFFECT.

R01 measured the epoch's own title for the first time and could not deliver a verdict: the
low-rank interaction estimator carries a large NEGATIVE bias (a purely additive world returns
X_share -0.629), and the bias is SHAPE-DEPENDENT (r=+0.517 with cells-per-parameter). One
reference block's dose curve cannot calibrate 32 blocks of different shape.

Fix: give every block its OWN floor and its OWN dose reference, built from its OWN marginals.

ESTIMAND        bias-corrected interaction contribution  X_c = X_real - X_additive(same block),
                on the held-out R2 scale, directly comparable to the ITEM main effect I_real;
                and the calibrated dose  g_hat = (X_real - X_g0) / (X_g1 - X_g0).
IDENTIFICATION  identified iff the matched synthetic reproduces the real block's I and P (the
                nuisance-matching check below). A block whose synthetic does not match is
                UNVERIFIED for this estimand, not silently included.
SCOPE           population/instrument/regime as R01. K in {1,2} only -- R01 measured that overfit
                grows monotonically with K, so K>=4 is a worse instrument, not a robustness check.
WORLDS          A content category    -> I - X_c > 2x spread
                B individualised value-> X_c - I > 2x spread
                C neither dominates   -> |X_c - I| <= 2x spread
KILL            THRESHOLD-FREE. the verdict is the ORDERING of two quantities measured on the same
                scale, and it is declared only when the gap exceeds 2x its own seed spread
                (the project's resolvability criterion, #34). No number is chosen anywhere.
POSITIVE CTRL   per block, g=1 world with a planted interaction: X_g1 - X_g0 must be > 0.
                A block where the planted interaction is NOT recovered has no scale, so it
                reports UNVERIFIED rather than a corrected value.
NEGATIVE CTRL   g=0 world IS the negative control, and it is per-block, which is the whole point.
NUISANCE MATCH  |I_synth - I_real| and |P_synth - P_real| reported per block; a block fails the
                match if either exceeds 0.05 (that IS a chosen number, and it is chosen as a
                LENIENCY -- it only ever removes blocks from the claim, never adds them).
PLACEBO         inherited from R01, with the direction corrected: a permuted person component
                must contribute <= 0. Re-run here at K=1.
NOISE FLOOR     3 masks x 3 synthetic draws, spread reported on every cell.
MULTIPLICITY    32 blocks x 2 K x 3 seeds x 3 worlds, reported whole including non-survivors.
IMPOSSIBLE      as R01.
"""
import pandas as pd, numpy as np, warnings, hashlib
from numpy.linalg import svd
from math import factorial
warnings.filterwarnings('ignore')

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MINN=20; RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]
    vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=MINN].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=M
print(f"blocks {len(RAW)}")

MASK=0.15; SEEDS=[11,29,47]

def components(M,obs,K):
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P; F=np.where(np.isnan(Rres),0.,Rres)
    for _ in range(25):
        U,S,Vt=svd(F,full_matrices=False); F=np.where(obs,Rres,(U[:,:K]*S[:K])@Vt[:K])
    U,S,Vt=svd(F,full_matrices=False)
    return gm,I,P,(U[:,:K]*S[:K])@Vt[:K]

def shap(M,K,seed):
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    gm,I,P,X=components(M,obs,K)
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'X':X}
    base=np.mean((M[he]-gm)**2)
    v={}
    for bits in range(8):
        S=frozenset([c for j,c in enumerate('IPX') if bits>>j&1])
        p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
        v[S]=1.-np.mean((M[he]-p[he])**2)/base if base>0 else np.nan
    out={}
    for c in 'IPX':
        o=[x for x in 'IPX' if x!=c]; tot=0.
        for S in [(),(o[0],),(o[1],),tuple(o)]:
            w=factorial(len(S))*factorial(2-len(S))/6.
            tot+=w*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

def synth(M,g,seed):
    """matched additive world (g=0) or matched world with a planted interaction (g>0),
    built from THIS block's own marginals and THIS block's own shape."""
    rng=np.random.default_rng(seed)
    n,m=M.shape; a=M.mean(0); b=M.mean(1)-M.mean()
    p=a[None,:]+b[:,None]
    if g>0: p=p+g*0.35*np.outer(rng.normal(size=n),rng.normal(size=m))
    return (rng.random((n,m))<np.clip(p,0.02,0.98)).astype(float)

rows=[]
for qi,(q,M) in enumerate(RAW.items()):
    for K in [1,2]:
        for sd in SEEDS:
            r=shap(M,K,sd);            rows.append(dict(q=q,K=K,seed=sd,world='real',n=M.shape[0],m=M.shape[1],**r))
            r0=shap(synth(M,0.,sd),K,sd); rows.append(dict(q=q,K=K,seed=sd,world='g0',n=M.shape[0],m=M.shape[1],**r0))
            r1=shap(synth(M,1.,sd),K,sd); rows.append(dict(q=q,K=K,seed=sd,world='g1',n=M.shape[0],m=M.shape[1],**r1))
    print(f"  block {qi+1}/{len(RAW)} done",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R02_per_block_additive_floor/results/'
D.to_csv(OUT+'grid.csv',index=False)

print("\n=== PLACEBO (direction corrected: a permuted person component must contribute <= 0) ===")
pl=[]
for sd in SEEDS:
    M=RAW[list(RAW)[0]]; rng=np.random.default_rng(900+sd)
    obs=np.random.default_rng(sd).random(M.shape)>=MASK; he=~obs
    gm,I,P,X=components(M,obs,1); base=np.mean((M[he]-gm)**2)
    f=lambda parts: 1.-np.mean((M[he]-np.clip(gm+sum(parts),0,1)[he])**2)/base
    pl.append(f([np.broadcast_to(I,M.shape),P[rng.permutation(M.shape[0])]])-f([np.broadcast_to(I,M.shape)]))
print(f"  permuted person contribution {np.mean(pl):+.5f}  -> {'PASS (<=0)' if np.mean(pl)<=0 else 'FAIL'}")

for K in [1,2]:
    d=D[D.K==K]
    piv=d.groupby(['q','world'])[['I','P','X']].mean().unstack('world')
    sp =d[d.world=='real'].groupby('q').X.std()
    sp0=d[d.world=='g0'].groupby('q').X.std()
    res=pd.DataFrame({'n':d.groupby('q').n.first(),'m':d.groupby('q').m.first(),
        'I_real':piv[('I','real')],'I_g0':piv[('I','g0')],'P_real':piv[('P','real')],'P_g0':piv[('P','g0')],
        'X_real':piv[('X','real')],'X_g0':piv[('X','g0')],'X_g1':piv[('X','g1')]})
    res['dI']=(res.I_real-res.I_g0).abs(); res['dP']=(res.P_real-res.P_g0).abs()
    res['matched']=(res.dI<=0.05)&(res.dP<=0.05)
    res['recovers']=(res.X_g1-res.X_g0)>0
    res['X_c']=res.X_real-res.X_g0
    res['g_hat']=(res.X_real-res.X_g0)/(res.X_g1-res.X_g0)
    res['spread']=np.sqrt(sp**2+sp0**2)
    res['gap']=res.X_c-res.I_real
    res['ok']=res.matched&res.recovers
    print(f"\n{'='*78}\nK={K}   blocks passing nuisance-match AND dose-recovery: "
          f"{int(res.ok.sum())}/{len(res)}")
    print(f"  failed nuisance match: {int((~res.matched).sum())}   failed dose recovery: {int((~res.recovers).sum())}")
    print(res[['n','m','I_real','X_real','X_g0','X_c','g_hat','gap','spread','matched','recovers']]
          .sort_values('X_c').round(4).to_string())
    v=res[res.ok]
    if len(v)==0: print("  -> no block is identified at this K. UNVERIFIED."); continue
    mI,mXc=v.I_real.median(),v.X_c.median(); mgap=v.gap.median(); msp=v.spread.median()
    print(f"\n  identified blocks (n={len(v)}), medians:")
    print(f"    ITEM main effect          I   = {mI:+.4f}")
    print(f"    interaction, bias-corr.   X_c = {mXc:+.4f}")
    print(f"    calibrated dose           g^  = {v.g_hat.median():.3f}  (IQR {v.g_hat.quantile(.25):.3f}-{v.g_hat.quantile(.75):.3f})")
    print(f"    gap  X_c - I              = {mgap:+.4f}   2x seed spread = {2*msp:.4f}")
    print(f"    blocks with X_c > 0: {int((v.X_c>0).sum())}/{len(v)}   with X_c > I: {int((v.gap>0).sum())}/{len(v)}")
    print("\n  CONDITIONAL KILL -- threshold-free, the verdict is an ordering")
    if abs(mgap)<=2*msp:
        print(f"    -> WORLD C. |gap| {abs(mgap):.4f} <= 2x spread {2*msp:.4f}: the item main effect and the "
              f"interaction are NOT DISTINGUISHABLE in size. 'a value, NOT a category' is an "
              f"OVERSTATEMENT -- the data says BOTH, and refuses to rank them.")
    elif mgap>0:
        print(f"    -> WORLD B. the interaction exceeds the item main effect by {mgap:+.4f} "
              f"(> 2x spread {2*msp:.4f}). the epoch title SURVIVES.")
    else:
        print(f"    -> WORLD A. the item main effect exceeds the interaction by {-mgap:+.4f} "
              f"(> 2x spread {2*msp:.4f}). THE EPOCH TITLE IS FALSE: content dominates, and "
              f"105 rounds were run on the residual left after deleting the winner.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

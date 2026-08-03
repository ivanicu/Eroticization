import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
E01 A09 R01 -- THE EPOCH'S OWN TITLE, MEASURED FOR THE FIRST TIME.

The epoch is called "sexual as a value, not a category". That is a claim about the RELATIVE SIZE
of two things: the ITEM main effect (content is content -- an option is erotic, and people differ
mainly in threshold) and the PERSON x ITEM INTERACTION (the same content carries different value
for different people -- an individualised readout weight). 105 rounds have been run inside this
epoch and not one of them measured it, because the shared loader computes

    R = M - M.mean(0) ; R = R - R.mean(1)

on line 1 -- THE ITEM MAIN EFFECT IS DELETED BEFORE ANYTHING IS EVER LOOKED AT. Every claim in
this project is a claim about the interaction, made after the rival was removed from the data.

ESTIMAND        share of held-out predictable variance in the RAW binary endorsement matrix
                attributable to each of {ITEM main effect, PERSON main effect, INTERACTION},
                as a Shapley decomposition over all 6 orderings (the decomposition is
                path-dependent; the point estimate alone is not reportable).
IDENTIFICATION  within a block the matrix is COMPLETE for people who entered it, so masking
                cells at random is a genuine held-out design, not an imputation of the gate.
SCOPE           population: people who entered a block and endorsed >=1 option in it (people
                who entered and endorsed nothing are absent from the long table -- a floor
                selection, stated not corrected). instrument: binary endorsement, not a rating.
                baseline: grand mean of training cells. regime: 32 blocks, >=8 options, >=1200 people.
WORLDS          A  content category   -> ITEM share dominates, interaction near its null floor
                B  individualised value -> INTERACTION share >= ITEM share
                C  recursive mix      -> both large, neither dominant
KILL            pre-registered, and the threshold is DERIVED not chosen: the interaction is real
                only if its share exceeds the ADDITIVE synthetic control's interaction share by
                more than 2x the seed spread (the project's own resolvability criterion, #34).
                The A-vs-B verdict is then the ORDERING of item vs interaction share, reported
                with the ordering's own seed spread.
POSITIVE CTRL   graded synthetic worlds at interaction strength g in {0, 0.25, 0.5, 1.0} on the
                real matrix's own marginals. Must be monotone in g AND must return the additive
                floor at g=0 (a control that passes at g=0 is a control that cannot fail).
NEGATIVE CTRL   shuffle option labels independently within each person: destroys the item main
                effect and the interaction, preserves person breadth exactly. ITEM share -> 0.
                World it excludes: "the item share is an artefact of row sums".
PLACEBO         a component built from a permuted person index must contribute exactly 0.
NOISE FLOOR     measured: 3 independent masks per cell, spread reported everywhere.
MULTIPLICITY    32 blocks x K in {1,2,4,8} x 3 seeds x 8 subset models, reported whole.
SPECIFICATION   K sweep is the rank axis; the 6 orderings are the attribution axis; both published
                including the cells that disagree.
IMPOSSIBLE      causally identified (no intervention on a survey), temporally resolved (no
                timestamps), independently replicated (one release), construct validated (no
                external gold standard for "erotic value").
"""
import pandas as pd, numpy as np, warnings
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
    RAW[q.qi]=M                                        # NOT centred. that is the whole point.
print(f"blocks {len(RAW)}   raw (uncentred) matrices, shapes "
      f"{min(m.shape[0] for m in RAW.values())}-{max(m.shape[0] for m in RAW.values())} x "
      f"{min(m.shape[1] for m in RAW.values())}-{max(m.shape[1] for m in RAW.values())}")

MASK_FRAC=0.15; SEEDS=[11,29,47]; KS=[1,2,4,8]

def components(M,obs,K):
    """three components computed from TRAINING cells only, in a fixed construction order.
    the order of CONSTRUCTION is fixed; the order of ATTRIBUTION is swept by Shapley below."""
    T=np.where(obs,M,np.nan)
    gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm)
    I=(cm-gm)[None,:]                                   # ITEM main effect
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm)
    P=rm[:,None]                                        # PERSON main effect
    Rres=T1-P
    F=np.where(np.isnan(Rres),0.,Rres)                  # soft-impute on the double-centred residual
    for _ in range(25):
        U,S,Vt=svd(F,full_matrices=False)
        L=(U[:,:K]*S[:K])@Vt[:K]
        F=np.where(obs,Rres,L)
    U,S,Vt=svd(F,full_matrices=False)
    X=(U[:,:K]*S[:K])@Vt[:K]                            # INTERACTION
    return gm,I,P,X

def r2(M,he,gm,parts):
    p=np.clip(gm+sum(parts) if parts else np.full(M.shape,gm),0.,1.)
    if np.isscalar(p): p=np.full(M.shape,p)
    base=np.mean((M[he]-gm)**2)
    return 1.-np.mean((M[he]-p[he])**2)/base if base>0 else np.nan

SH=[]
for i,c in enumerate('IPX'):
    w=[]
    for S in range(4):
        for mask in range(8):
            pass
    SH.append(c)

def shapley(v):
    """v: dict from frozenset(subset of 'IPX') -> R2.  exact Shapley over the 3 components."""
    C='IPX'; out={}
    for c in C:
        tot=0.
        others=[x for x in C if x!=c]
        for r in range(3):
            for S in ([()] if r==0 else ([(others[0],),(others[1],)] if r==1 else [tuple(others)])):
                w=factorial(len(S))*factorial(3-len(S)-1)/factorial(3)
                tot+=w*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

def orderings(v):
    """all 6 orderings -> marginal contribution of each component. the path-dependence itself."""
    C='IPX'; rows=[]
    import itertools
    for perm in itertools.permutations(C):
        cur=(); d={}
        for c in perm:
            d[c]=v[frozenset(cur+(c,))]-v[frozenset(cur)]; cur=cur+(c,)
        rows.append(dict(order=''.join(perm),**d))
    return pd.DataFrame(rows)

def run_matrix(M,K,seed):
    rng=np.random.default_rng(seed)
    obs=rng.random(M.shape)>=MASK_FRAC; he=~obs
    gm,I,P,X=components(M,obs,K)
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'X':X}
    v={}
    for bits in range(8):
        S=frozenset([c for j,c in enumerate('IPX') if bits>>j&1])
        v[S]=r2(M,he,gm,[comp[c] for c in S])
    return v

# ---------------- controls, built on a real block's own marginals ----------------
ref=sorted(RAW,key=lambda q:RAW[q].size)[len(RAW)//2]
Mref=RAW[ref]; n,m=Mref.shape
def synth(g,seed):
    rng=np.random.default_rng(1000+seed)
    a=Mref.mean(0); b=Mref.mean(1); b=b-b.mean()
    u=rng.normal(size=n); v=rng.normal(size=m)
    p=np.clip(a[None,:]+b[:,None]+g*0.35*np.outer(u,v),0.02,0.98)
    return (rng.random(p.shape)<p).astype(float)

print("\n=== POSITIVE CONTROL (graded): synthetic worlds on the reference block's marginals ===")
print(f"reference block {ref}  shape {n}x{m}   K=4")
ctrl=[]
for g in [0.,0.25,0.5,1.0]:
    for sd in SEEDS:
        sh=shapley(run_matrix(synth(g,sd),4,sd))
        t=sum(max(x,0) for x in sh.values())
        ctrl.append(dict(g=g,seed=sd,**{c:sh[c] for c in 'IPX'},
                         X_share=sh['X']/t if t>0 else np.nan))
C=pd.DataFrame(ctrl); cg=C.groupby('g').agg(['mean','std'])
print(C.groupby('g')[['I','P','X','X_share']].mean().round(4).to_string())
FLOOR=float(C[C.g==0].X_share.mean()); FLOOR_SD=float(C[C.g==0].X_share.std())
print(f"\nadditive-world interaction share (the FLOOR)  {FLOOR:+.4f}  seed sd {FLOOR_SD:.4f}")
mono=list(C.groupby('g').X_share.mean())
print(f"monotone in g: {all(mono[i]<mono[i+1] for i in range(3))}   {[round(x,4) for x in mono]}")
print(f"fails at g=0 (share below 2x its own sd): {abs(FLOOR)<2*max(FLOOR_SD,1e-9) or FLOOR<0.02}")

print("\n=== NEGATIVE CONTROL: option labels shuffled within each person ===")
neg=[]
for sd in SEEDS:
    rng=np.random.default_rng(500+sd)
    Ms=np.array([rng.permutation(row) for row in Mref])
    sh=shapley(run_matrix(Ms,4,sd)); t=sum(max(x,0) for x in sh.values())
    neg.append(dict(seed=sd,**{c:sh[c] for c in 'IPX'},I_share=sh['I']/t if t>0 else np.nan))
N=pd.DataFrame(neg); print(N.round(4).to_string(index=False))
print(f"item share after within-person shuffle: {N.I_share.mean():+.4f} (was the item effect real? "
      f"{'YES' if N.I_share.mean()<0.1 else 'NO -- row sums alone reproduce it'})")

print("\n=== PLACEBO: person component built from a permuted person index ===")
pl=[]
for sd in SEEDS:
    rng=np.random.default_rng(900+sd)
    obs=np.random.default_rng(sd).random(Mref.shape)>=MASK_FRAC; he=~obs
    gm,I,P,X=components(Mref,obs,4)
    Pp=P[rng.permutation(n)]
    pl.append(r2(Mref,he,gm,[np.broadcast_to(I,Mref.shape),Pp])-
              r2(Mref,he,gm,[np.broadcast_to(I,Mref.shape)]))
print(f"contribution of a permuted person effect: {np.mean(pl):+.5f} (must be ~0, is it? "
      f"{'YES' if abs(np.mean(pl))<0.002 else 'NO'})")

# ---------------- the estimand, over the whole grid ----------------
print("\n=== THE GRID: 32 blocks x K in {1,2,4,8} x 3 seeds ===")
rows=[]; ordrows=[]
for q,M in RAW.items():
    for K in KS:
        for sd in SEEDS:
            v=run_matrix(M,K,sd); sh=shapley(v)
            t=sum(max(x,0) for x in sh.values())
            rows.append(dict(q=q,K=K,seed=sd,n=M.shape[0],m=M.shape[1],full=v[frozenset('IPX')],
                             **{c:sh[c] for c in 'IPX'},
                             **{c+'_share':(sh[c]/t if t>0 else np.nan) for c in 'IPX'}))
            if K==4: ordrows.append(orderings(v).assign(q=q,seed=sd))
G=pd.DataFrame(rows); O=pd.concat(ordrows)
G.to_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/'
         'R01_item_person_interaction/results/grid.csv',index=False)

print("\nheld-out R2 of the full model, and the three shares, by K (mean over 32 blocks x 3 seeds)")
print(G.groupby('K')[['full','I_share','P_share','X_share']].mean().round(4).to_string())

print("\nper-block at K=4, sorted by interaction share (ALL blocks, none dropped)")
pb=G[G.K==4].groupby('q').agg(n=('n','first'),m=('m','first'),full=('full','mean'),
        I=('I_share','mean'),P=('P_share','mean'),X=('X_share','mean'),
        Xsd=('X_share','std')).sort_values('X')
print(pb.round(4).to_string())

print("\nPATH DEPENDENCE -- marginal contribution of each component by ordering (K=4, all blocks)")
print(O.groupby('order')[['I','P','X']].mean().round(4).to_string())
print(f"item marginal ranges [{O.groupby('order').I.mean().min():.4f}, "
      f"{O.groupby('order').I.mean().max():.4f}] across orderings -- "
      f"this is why a single decomposition number is not reportable")

med=G[G.K==4].groupby('q')[['I_share','P_share','X_share']].mean()
mi,mp,mx=med.I_share.median(),med.P_share.median(),med.X_share.median()
sdx=G[G.K==4].groupby('q').X_share.std().median()
print("\n" + "="*72)
print("CONDITIONAL KILL -- gates first, threshold only if they pass")
g1=all(mono[i]<mono[i+1] for i in range(3)); print(f"  (a) positive control monotone in g      : {'PASS' if g1 else 'FAIL'}")
g2=FLOOR<0.02; print(f"  (b) positive control fails at g=0       : {'PASS' if g2 else 'FAIL'} (floor {FLOOR:+.4f})")
g3=N.I_share.mean()<0.1; print(f"  (c) negative control kills item share   : {'PASS' if g3 else 'FAIL'} ({N.I_share.mean():+.4f})")
g4=abs(np.mean(pl))<0.002; print(f"  (d) placebo returns zero                : {'PASS' if g4 else 'FAIL'}")
if g1 and g2 and g3 and g4:
    real=(mx-FLOOR)>2*max(sdx,FLOOR_SD)
    print(f"\n  median block, K=4:  ITEM {mi:.1%}   PERSON {mp:.1%}   INTERACTION {mx:.1%}")
    print(f"  interaction above the additive floor by {mx-FLOOR:+.4f}, "
          f"2x seed spread = {2*max(sdx,FLOOR_SD):.4f}  -> {'RESOLVABLE' if real else 'UNRESOLVABLE'}")
    if not real: print("  -> INTERACTION UNVERIFIED. no A-vs-B verdict is licensed.")
    elif mx>=mi:  print("  -> WORLD B: the interaction carries at least as much as the item. "
                        "'a value, not a category' SURVIVES.")
    elif mi>=0.60 and mx<=0.15: print("  -> WORLD A: content dominates. THE EPOCH TITLE IS FALSE.")
    else: print(f"  -> WORLD C: neither dominates (item {mi:.1%} vs interaction {mx:.1%}). "
                f"the epoch title is an OVERSTATEMENT, not a finding.")
else:
    print("\n  -> gates failed. UNVERIFIED, and that is not an acquittal.")
import hashlib
print(f"\nartifact sha1 {hashlib.sha1(G.to_csv(index=False).encode()).hexdigest()[:12]}")

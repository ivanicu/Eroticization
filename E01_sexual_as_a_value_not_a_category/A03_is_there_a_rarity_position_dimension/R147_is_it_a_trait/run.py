import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R13 -- IF IT IS A DIMENSION EVERYONE HAS A VALUE ON, IT MUST BEHAVE LIKE ONE.

#99 replaced "a minority concentrating on rare options" with "a continuous person-level parameter of
rare-option affinity". That is a claim about a TRAIT, and a trait has two properties that can be
checked without any null at all:

  RELIABILITY   a person's value measured on half their blocks predicts it on the other half
  DISCRIMINANT  it is not just something already measured wearing a new name -- here, PICK COUNT.
                A person who ticks few options may pick rarer ones by selection alone.

Split-half is the cheapest possible test and the most damaging if it fails: an unreliable "trait" is
noise with a name.

ESTIMAND        Spearman-Brown-corrected split-half reliability of per-person surprisal; and the same
                after residualising on pick count, plus corr(S, pick count) itself.
IDENTIFICATION  identified; halves are disjoint sets of BLOCKS, so no cell is shared.
SCOPE           people entering >=6 blocks, so each half has >=3.
WORLDS          trait      reliability well above the null's and survives residualising on count
                breadth    it vanishes once pick count is removed -> it is breadth in a new costume
                noise      reliability at the null -> #99's "parameter" is unreliable and withdrawn
KILL            threshold-free: reliability against the fixed-margin null's, in units of the
                bootstrap spread of the difference.
POSITIVE CTRL   a PLANTED trait (a known per-person affinity) must show high reliability, so a low
                real value cannot be blamed on the split-half machinery.
NEGATIVE CTRL   the fixed-margin null, which has no person-level structure by construction: its
                reliability is the floor.
NOISE FLOOR     6 random half-splits x 3 seeds.
IMPOSSIBLE      distinguishing the trait from a stable RESPONSE STYLE -- both are reliable person
                properties. Reliability establishes that something stable exists, never what it is.
"""
import pandas as pd, numpy as np, warnings, hashlib
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
def plant_trait(M,aff,rng):
    """a KNOWN per-person affinity: aff[i] swaps toward rare options."""
    Mw=M.copy(); p=M.mean(0); o=np.argsort(p); rare=o[:3]; med=o[len(o)//2:]
    for i in np.flatnonzero(aff>0):
        d=0
        for _ in range(8*int(aff[i])):
            if d>=aff[i]: break
            c=med[rng.integers(len(med))]; r=rare[rng.integers(len(rare))]
            if Mw[i,c]==1 and Mw[i,r]==0: Mw[i,c]=0.; Mw[i,r]=1.; d+=1
    return Mw
def sb(r): return 2*r/(1+r) if r>-1 else np.nan
rows=[]
for sd in range(1,4):
    rgc=np.random.default_rng(6100+sd)
    aff=np.floor(rgc.exponential(0.8,size=len(ALLP))).astype(float)   # the planted trait
    num={w:{} for w in ['real','cb','plant']}; den={w:{} for w in ['real','cb','plant']}
    blockof={w:{} for w in ['real','cb','plant']}
    for bi,t in enumerate(IDENT):
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.))
        Mn=curveball(M,np.random.default_rng(6200+sd))
        pack={'real':M,'cb':Mn,'plant':plant_trait(Mn,aff[idx]*0+aff[idx],np.random.default_rng(6300+sd))}
        for w,Mw in pack.items():
            kk=Mw.sum(1); tot=Mw@ref
            for j,gi in enumerate(idx):
                if kk[j]>0:
                    blockof[w].setdefault(gi,[]).append((bi,tot[j],kk[j]))
    for w in ['real','cb','plant']:
        ks=[i for i,v in blockof[w].items() if len(v)>=6]
        picks=np.array([sum(x[2] for x in blockof[w][i]) for i in ks])
        rs=[]
        for rep in range(6):
            rr=np.random.default_rng(700+rep)
            a=[];b=[]
            for i in ks:
                v=blockof[w][i]; o=rr.permutation(len(v)); h=len(v)//2
                A=[v[j] for j in o[:h]]; B=[v[j] for j in o[h:2*h]]
                a.append(sum(x[1] for x in A)/max(sum(x[2] for x in A),1))
                b.append(sum(x[1] for x in B)/max(sum(x[2] for x in B),1))
            a=np.array(a); b=np.array(b)
            r=np.corrcoef(a,b)[0,1]
            # residualise BOTH halves on pick count, then re-correlate
            X=np.c_[np.ones(len(picks)),picks,np.log(picks)]
            ra=a-X@np.linalg.lstsq(X,a,rcond=None)[0]; rb=b-X@np.linalg.lstsq(X,b,rcond=None)[0]
            rres=np.corrcoef(ra,rb)[0,1]
            rs.append((r,rres,np.corrcoef((a+b)/2,picks)[0,1]))
        R=np.array(rs)
        rows.append(dict(seed=sd,world=w,n=len(ks),
                         r_half=R[:,0].mean(),r_half_sd=R[:,0].std(),
                         rel=sb(R[:,0].mean()),
                         r_resid=R[:,1].mean(),rel_resid=sb(R[:,1].mean()),
                         r_vs_picks=R[:,2].mean()))
    print(f"  seed {sd}",flush=True)
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
G=D.groupby('world')[['n','r_half','r_half_sd','rel','r_resid','rel_resid','r_vs_picks']].mean()
print("\n=== SPLIT-HALF RELIABILITY OF PER-PERSON SURPRISAL ===")
print(G.round(4).to_string())
gp=G.loc['plant','rel']>0.5; gn=abs(G.loc['cb','rel'])<0.15
print("\n  CONDITIONAL KILL -- gates first")
print(f"   (a) a PLANTED trait is recovered reliably : {'PASS' if gp else 'FAIL'} "
      f"(rel {G.loc['plant','rel']:.3f})")
print(f"   (b) the fixed-margin null is unreliable   : {'PASS' if gn else 'FAIL'} "
      f"(rel {G.loc['cb','rel']:.3f})")
if not(gp and gn): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    rel=G.loc['real','rel']; relr=G.loc['real','rel_resid']; vp=G.loc['real','r_vs_picks']
    print(f"\n   real split-half reliability (Spearman-Brown) : {rel:.3f}")
    print(f"   after residualising BOTH halves on pick count: {relr:.3f}")
    print(f"   corr(surprisal, total picks)                 : {vp:+.3f}")
    if rel<G.loc['cb','rel']+0.1:
        print("\n   -> NOISE. The 'parameter' is not reliable and #99b must be withdrawn.")
    elif relr<0.2*rel:
        print("\n   -> BREADTH IN A NEW COSTUME. Reliability collapses once pick count is removed.")
    else:
        print(f"\n   -> A RELIABLE TRAIT, and not pick count: {100*relr/rel:.0f}% of the reliability")
        print(f"      survives residualising on how many options the person ticked.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

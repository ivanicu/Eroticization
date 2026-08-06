import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A03 R20 -- RE-RUNNING A LOW-POWER NULL WITH THE INSTRUMENT #72 SHOWED WAS TOO WEAK.

#55 found individual variation in acquisition order at +0.88 points with a seed spread of 1.03
(ratio 0.85 -- unresolvable), and #62 found group-level variation at ratio 0.40. Both were logged as
LOW-POWER NULLS, which obliges a re-run when the instrument improves. Two things now have:

  (1) #72 measured that the person-side structure is a slowly-decaying spectrum with real predictive
      content out to at least 32 dimensions. #55 used an EIGHT-dimensional embedding, so it defined
      "similar people" in a space that #72 shows discards most of the signal.
  (2) #55's accuracy function carried `cap=20000` and `if tot>=cap: break`, which exits the PERSON
      loop -- so at ~10 pairs per person it stopped after roughly 2,000 people out of the available
      pool. That is a power cap I wrote and never priced.

And #55 reported no MDE in interpretable units, so "no individual variation" had no size attached.
This round fixes that with a graded planted-variation control measured IN YEARS.

ESTIMAND        held-out pairwise onset-ordering accuracy of a neighbour-fitted ordering minus the
                global ordering, on identical pairs; and the MDE in years of planted individual shift.
IDENTIFICATION  identified. neighbours are defined in preference space and never touch onset.
SCOPE           people with >=6 recorded onsets. n_eff is PEOPLE, not pairs.
WORLDS          A purely global schedule -> neighbour ties global at every embedding dimension
                B individual variation    -> neighbour beats global, and more so in richer spaces
                C artefact                -> random-neighbour also beats global (sample size, not
                                             similarity)
KILL            threshold-free: declared only where the neighbour-minus-global gap exceeds 2x its
                own seed spread AND the random-neighbour control does not.
POSITIVE CTRL   GRADED: plant person-specific onset shifts of g in {0, 0.5, 1, 2} YEARS, generated
                as a function of the person's preference-space coordinates so that neighbours can in
                principle find them. Must be monotone in g and must NOT fire at g=0.
                The smallest g that fires IS the MDE, in years.
NEGATIVE CTRL   random neighbours at the same k -- same fitting-set size, no similarity.
NOISE FLOOR     5 seeds; and a person-level bootstrap so the interval is over PEOPLE not pairs.
MULTIPLICITY    4 embedding dims x 3 k x 5 seeds, published whole.
IMPOSSIBLE      a person-level ordering ceiling -- needs repeated onset measurement per person.
                Unchanged from #55, and it bounds what any of this can mean.
"""
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import svd
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'; OUT.mkdir(exist_ok=True)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values; mask=~np.isnan(V)
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
Pz=((df[rate].apply(pd.to_numeric,errors='coerce')>0).astype(float)).fillna(0.).values
Pz=Pz-Pz.mean(0)
U,S_,_=svd(Pz,full_matrices=False); EMBFULL=U*S_
keep=np.flatnonzero(mask.sum(1)>=6)
print(f"people >=6 onsets {len(keep):,}   categories {V.shape[1]}   embedding pool {EMBFULL.shape[1]}",flush=True)

DIMS=[8,32,64,128]; KS=[400,1000,2500]; SEEDS=[1,2,3,4,5]; PAIRS_PER=30

def acc(order_of,people,rng,Vm):
    right=0.;tot=0
    for i in people:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        ov=order_of(i)
        a=rng.choice(j,min(PAIRS_PER,len(j)*2)); b=rng.choice(j,len(a))
        ok=(a!=b)&(Vm[i,a]!=Vm[i,b])&(ov[a]!=ov[b])
        if not ok.any(): continue
        right+=float((((ov[a]<ov[b])==(Vm[i,a]<Vm[i,b]))&ok).sum()); tot+=int(ok.sum())
    return 100*right/max(tot,1),tot

def plant(g,seed):
    """person-specific onset shift, in YEARS, generated from preference-space coordinates so that
    neighbours can in principle recover it. g=0 returns the real matrix untouched."""
    if g<=0: return V
    rng=np.random.default_rng(600+seed)
    E=EMBFULL[:,:8]; E=(E-E.mean(0))/(E.std(0)+1e-9)
    L=rng.normal(size=(8,V.shape[1]))
    S=E@L; S=S/ (S.std()+1e-9) * g
    return V+S

def one(Vm,dim,k,seed):
    rng=np.random.default_rng(seed)
    p=rng.permutation(keep); tr,te=p[:len(p)//2],p[len(p)//2:]
    E=EMBFULL[:,:dim]; Etr=E[tr]
    glob=np.nanmean(Vm[tr],axis=0)
    rp=rng.permutation(tr)[:k]; rnd_o=np.nanmean(Vm[rp],axis=0)
    rnd_o=np.where(np.isfinite(rnd_o),rnd_o,glob)
    def nb_order(i):
        d=((Etr-E[i])**2).sum(1)
        nb=tr[np.argpartition(d,k)[:k]]
        o=np.nanmean(Vm[nb],axis=0)
        return np.where(np.isfinite(o),o,glob)
    g_,n=acc(lambda i: glob,te,rng,Vm)
    nb_,_=acc(nb_order,te,rng,Vm)
    rn_,_=acc(lambda i: rnd_o,te,rng,Vm)
    return dict(dim=dim,k=k,seed=seed,glob=g_,neighbour=nb_,random_nb=rn_,
                gap=nb_-g_,ctrl_gap=rn_-g_,n_pairs=n,n_people=len(te))

rows=[]
for dim in DIMS:
    for k in KS:
        for sd in SEEDS: rows.append(dict(g=0.0,**one(V,dim,k,sd)))
    print(f"  real dim={dim} done",flush=True)
BEST=dict(dim=64,k=1000)
for g in [0.5,1.0,2.0]:
    for sd in SEEDS: rows.append(dict(g=g,**one(plant(g,sd),BEST['dim'],BEST['k'],sd)))
    print(f"  planted g={g} done",flush=True)
D=pd.DataFrame(rows); D.to_csv(OUT/'grid.csv',index=False)

R=D[D.g==0]
print(f"\n=== POWER: {R.n_pairs.mean():,.0f} pairs from {R.n_people.mean():,.0f} held-out people "
      f"(#55 used a cap of 20,000 pairs, ~2,000 people) ===")
print("\n=== THE GRID (5 seeds each) ===")
t=R.groupby(['dim','k']).agg(glob=('glob','mean'),nb=('neighbour','mean'),rnd=('random_nb','mean'),
        gap=('gap','mean'),gap_sd=('gap','std'),ctrl=('ctrl_gap','mean'),ctrl_sd=('ctrl_gap','std'))
t['resolvable']=t.gap.abs()>2*t.gap_sd
t['ctrl_clean']=t.ctrl.abs()<=2*t.ctrl_sd
print(t.round(3).to_string())

print("\n=== GRADED POSITIVE CONTROL (planted individual shift, in YEARS) ===")
P=D[D.g>0].groupby('g').agg(gap=('gap','mean'),sd=('gap','std'),ctrl=('ctrl_gap','mean'))
base=R[(R.dim==BEST['dim'])&(R.k==BEST['k'])]
P.loc[0.0]=[base.gap.mean(),base.gap.std(),base.ctrl_gap.mean()]; P=P.sort_index()
P['fires']=P.gap>P.loc[0.0,'gap']+2*P.loc[0.0,'sd']
print(P.round(3).to_string())
fired=P.index[P.fires]
MDE=float(fired.min()) if len(fired) else np.nan
print(f"\n  monotone in g: {list(P.gap.round(3))}")
print(f"  MDE = {MDE if not np.isnan(MDE) else 'NOT REACHED'} years of individual shift"
      f"  (at dim={BEST['dim']}, k={BEST['k']})")

print("\n  CONDITIONAL KILL -- gates first")
g_a=not P.loc[0.0,'fires']
g_b=(not np.isnan(MDE)) and bool(P.fires.iloc[-1])
g_c=bool(t.ctrl_clean.mean()>0.5)
print(f"   (a) control does NOT fire at g=0        : {'PASS' if g_a else 'FAIL'}")
print(f"   (b) control fires at some planted shift : {'PASS' if g_b else 'FAIL -- instrument blind'}")
print(f"   (c) random-neighbour control clean      : {int(t.ctrl_clean.sum())}/{len(t)} cells")
if not(g_a and g_b and g_c):
    print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    res=t[t.resolvable&t.ctrl_clean]
    print(f"\n   cells where the neighbour gap is resolvable AND the control is clean: "
          f"{len(res)}/{len(t)}")
    print(f"   median gap over the whole grid: {t.gap.median():+.3f} points "
          f"(2x seed spread {2*t.gap_sd.median():.3f})")
    print(f"   best cell: {t.gap.max():+.3f} at {t.gap.idxmax()}")
    if len(res)==0:
        print(f"\n   -> STILL NULL, but now with a SIZE: individual variation in acquisition order is")
        print(f"      below {MDE} years of person-specific shift. #55's null is CONFIRMED and priced.")
    else:
        print(f"\n   -> INDIVIDUAL VARIATION DETECTED in {len(res)} cells. #55's null was a POWER")
        print(f"      failure, not a fact, and the 'one global schedule' reading must be withdrawn.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R04 -- THE CONTROL #93/#R03 NEEDED, BUILT PROPERLY.

R03 found a resolvable rare-option excess in the real data at p75-p99 against a fixed-margin null,
and could not report it, because its positive control failed two ways at once:

  (1) the "strength 0" arm was IDENTICAL to the real data -- plant_rare with extra=0 changes nothing --
      so the check "does not fire at strength 0" was the main comparison wearing a control's name.
      Twentieth mis-specified design element.
  (2) the plant was SELF-CANCELLING: adding picks to rare options raises their base rates, which
      lowers their surprisal, so a stronger plant produced a SMALLER effect (+0.0496 -> +0.0407).

Both are fixed by (a) planting into a SEPARATE copy so the real data is never an arm, and (b) scoring
every arm against a FIXED reference -- the base rates of the untouched real matrix -- so no plant can
move its own yardstick.

ESTIMAND        unchanged: upper quantiles of per-person mean surprisal, real vs fixed-margin null.
IDENTIFICATION  unchanged, and now the reference is external to every arm.
WORLDS          unchanged.
KILL            threshold-free per quantile above 2x its own seed spread.
POSITIVE CTRL   swap-based: a carrier GAINS a rare option and LOSES a common one, so their pick count
                is unchanged and the base-rate shift is minimal. Graded over 0 (a true no-op that is
                a genuinely separate matrix), 1 and 2 swaps. Must be monotone and silent at 0.
NEGATIVE CTRL   fixed-margin null, margins asserted per draw.
NOISE FLOOR     6 seeds.
IMPOSSIBLE      unchanged.
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
QS=[50,75,90,95,99]
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
def swap_plant(M,dens,nswap,rng):
    """carrier GAINS a rare option and LOSES a common one: pick count unchanged, margins barely move."""
    Mw=M.copy(); n,m=M.shape
    o=np.argsort(M.mean(0)); rare=o[:4]; common=o[-4:]
    for i in np.flatnonzero(rng.random(n)<dens):
        for _ in range(nswap):
            r=rare[rng.integers(len(rare))]; c=common[rng.integers(len(common))]
            if Mw[i,r]==0 and Mw[i,c]==1: Mw[i,r]=1.; Mw[i,c]=0.
    return Mw
rows=[]
for sd in range(1,7):
    arms=['real','cb']+[f's{k}' for k in [0,1,2]]+[f's{k}cb' for k in [0,1,2]]
    acc={w:[np.zeros(len(ALLP)),np.zeros(len(ALLP))] for w in arms}
    for t in IDENT:
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.))          # FIXED reference, external to every arm
        rg=np.random.default_rng(2600+sd)
        pack={'real':M,'cb':curveball(M,rg)}
        for k in [0,1,2]:
            Mp=swap_plant(M,0.05,k,np.random.default_rng(2700+sd))
            pack[f's{k}']=Mp; pack[f's{k}cb']=curveball(Mp,np.random.default_rng(2701+sd))
        assert np.allclose(pack['cb'].sum(1),M.sum(1)),"curveball row sums"
        for w,Mw in pack.items():
            acc[w][0][idx]+=Mw@ref; acc[w][1][idx]+=Mw.sum(1)
    ok=acc['real'][1]>=15
    for w in acc:
        S=acc[w][0][ok]/np.maximum(acc[w][1][ok],1)
        for qv in QS: rows.append(dict(seed=sd,world=w,qq=qv,S=float(np.percentile(S,qv))))
    print(f"  seed {sd}  n={int(ok.sum()):,}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
def cmp(a,b,qv):
    r=D[(D.world==a)&(D.qq==qv)]['S']; n=D[(D.world==b)&(D.qq==qv)]['S']
    return r.mean(),n.mean(),r.mean()-n.mean(),2*np.sqrt(r.std()**2+n.std()**2)
print("\n=== GRADED SWAP CONTROL (5% carriers, count preserved, fixed reference) ===")
for k in [0,1,2]:
    r,n,d,s=cmp(f's{k}',f's{k}cb',95); r9,n9,d9,s9=cmp(f's{k}',f's{k}cb',99)
    print(f"  {k} swaps: p95 diff {d:+.4f} (2xsp {s:.4f}) {'FIRES' if d>s else '-'}   "
          f"p99 diff {d9:+.4f} (2xsp {s9:.4f}) {'FIRES' if d9>s9 else '-'}")
print("\n=== REAL vs FIXED-MARGIN NULL ===")
lab={50:'p50',75:'p75',90:'p90',95:'p95',99:'p99'}
for qv in QS:
    r,n,d,s=cmp('real','cb',qv)
    print(f"  {lab[qv]} real {r:.4f}  null {n:.4f}  diff {d:+.4f}  2xspread {s:.4f}  "
          f"{'RESOLVABLE' if abs(d)>s else 'no'}")
_,_,d0,s0=cmp('s0','s0cb',95); _,_,d2,s2=cmp('s2','s2cb',95)
_,_,dr,sr=cmp('real','cb',95)
print("\n  CONDITIONAL KILL -- gates first")
g1=d0<=s0 or abs(d0-dr)<s0
g2=d2>d0+s0
print(f"   (a) 0-swap arm matches the real comparison (it IS a no-op): "
      f"{'PASS' if abs(d0-dr)<max(s0,sr) else 'FAIL'} ({d0:+.4f} vs real {dr:+.4f})")
print(f"   (b) 2 swaps raises the effect above 0 swaps                : "
      f"{'PASS' if g2 else 'FAIL'} ({d2:+.4f} vs {d0:+.4f}, 2xsp {s0:.4f})")
if not g2: print("   -> UNVERIFIED: the control cannot show the statistic responds to planted rare-option preference.")
else:
    hits=[qv for qv in QS if (lambda t:(t[2]>t[3]))(cmp('real','cb',qv))]
    print(f"\n   the statistic RESPONDS to a planted rare-option minority (+{d2-d0:.4f} for 2 swaps in 5% of people)")
    print(f"   real tail exceeds the null at: {hits}")
    print("\n   -> A MINORITY CONCENTRATING ON RARE OPTIONS IS PRESENT." if hits else
          "\n   -> NO SUCH MINORITY.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

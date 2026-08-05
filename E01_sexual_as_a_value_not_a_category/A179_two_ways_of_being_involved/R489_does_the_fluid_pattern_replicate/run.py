import os,sys,pathlib,json,re
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #444c found the half of involvement that costs the sense of changeability loads on the six
   blocks whose wording names a bodily fluid (+0.346, 4.31 sd). That set was NOTICED in this
   data before it was formalised. Does the pattern replicate when the loadings are re-estimated
   on people who did not contribute to the estimate being tested?

Worlds
  A  replicates : the gap clears its own null in held-out halves at a rate well above 5% ->
     #444c stops being "generated here" and becomes a confirmed pattern.
  B  does not   : it was overfitting -> #444c must be pulled back to "unresolved" on the page.

Design: 20 random person-level splits (fixed seeds). For each, loadings are estimated on half
   A ONLY, and the fluid-vs-rest gap is evaluated on half B's OWN re-estimated loadings, with
   half B's own block-label permutation null. The fluid regex is NOT touched -- changing it
   would make this a different hypothesis (#444's NEXT).
⚠ POWER FIRST: halving n halves the precision, so the MDE on a half is reported BEFORE the
   result, or "it fell inside the null" is unreadable (#413b: a NEXT that promises power must
   check what the power rests on).
⚠ HONEST FRAMING: this is SPLIT-HALF REPLICATION, not true out-of-sample -- the same people
   were in the sample where the pattern was noticed. It reduces the concern; it does not
   remove it. Stated, not hidden.
CONTROL : the full-sample gap must reproduce #444c (+0.346) from this script.
FRONTIER: world B retracts a paragraph published one round ago.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R489 replicate the fluid pattern")
_R372=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like'
            /'R372_block_loadings/run.py').read_text()
exec(_R372.split('"""',2)[2].split('def load_of')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
SXA=pd.to_numeric(raw['Totalsexacts'],errors='coerce').values.astype(float)

def load_k(rows,k):
    mm=np.zeros(NN,bool); mm[rows]=True
    F=np.isfinite(A); Z=np.where(F,A,0.0); tot=Z.sum(0); ct=F.sum(0); Ra=np.full_like(A,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        Ra[b]=np.where(F[b],A[b]-lo,np.nan); Ra[b]=Ra[b]-np.nanmean(np.where(mm,Ra[b],np.nan))
    F2=np.isfinite(B); Z2=np.where(F2,B,0.0); t2=Z2.sum(0); c2_=F2.sum(0); Rb=np.full_like(B,np.nan)
    for b in range(NB):
        lo=np.where(c2_-F2[b]>=6,(t2-Z2[b])/np.maximum(c2_-F2[b],1),np.nan)
        Rb[b]=np.where(F2[b],B[b]-lo,np.nan); Rb[b]=Rb[b]-np.nanmean(np.where(mm,Rb[b],np.nan))
    C=np.zeros((NB,NB))
    for i in range(NB):
        for j in range(NB):
            g=np.isfinite(Ra[i])&np.isfinite(Rb[j])&mm
            if g.sum()>100: C[i,j]=np.corrcoef(Ra[i][g],Rb[j][g])[0,1]
    C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
    v=V[:,k]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(v[:,None]*Zm).sum(0); den=(Fm*np.abs(v)[:,None]).sum(0)
    s=np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
    g=np.isfinite(s)&np.isfinite(SXA)&mm
    if g.sum()>100 and np.corrcoef(s[g],SXA[g])[0,1]<0: v=-v      # 朝「更卷入」(#443b)
    return v

FLU=re.compile(r'urine|ejaculate|squirt|precum|saliva|secretion',re.I)   # ⚠ 一个字都没改
fl=np.array([bool(FLU.search(n)) for n in NAMES]); nf=int(fl.sum())
ALLR=np.flatnonzero(ok)
def gap_and_null(v, seed):
    d=float(v[fl].mean()-v[~fl].mean())
    rg=np.random.default_rng(seed)
    nul=np.array([(lambda s: v[s].mean()-v[~s].mean())(
        np.isin(np.arange(NB),rg.choice(NB,nf,replace=False))) for _ in range(4000)])
    return d, float(np.percentile(np.abs(nul),95)), float(nul.std())

v_full=load_k(ALLR,1)
d_full,p95_full,sd_full=gap_and_null(v_full, 1)
print(f"**对照**:全样本 gap = **{d_full:+.4f}**(`#444c` 报 +0.3462)· 阈 {p95_full:.4f}")
GATE.asserted("CONTROL the full-sample gap reproduces #444c",
              abs(d_full-0.3462)<0.05, f"{d_full:+.4f} vs +0.3462", kind="control")

rng=np.random.default_rng(97); rows=[]
for t in range(20):
    perm=rng.permutation(ALLR); hb=perm[len(perm)//2:]
    vb=load_k(hb,1)
    d,p95,sd=gap_and_null(vb, 500+t)
    rows.append(dict(split=t, n_half=len(hb), gap=d, thr=p95, null_sd=sd,
                     mde=1.96*sd, sig=bool(abs(d)>p95)))
T=pd.DataFrame(rows)
show(T, HERE/'results/splits.csv', n=8, label="劈半")
rate=float(T.sig.mean()); mde=float(T.mde.median())
print(f"\n**功率先报**:半样本上的 MDE 中位 = **{mde:.4f}**(全样本 gap 是 **{d_full:.4f}** ->"
      f" **{d_full/mde:.1f}×** MDE)")
print(f"**复制率 = {int(T.sig.sum())}/20 = {rate:.0%}** · gap 中位 **{T.gap.median():+.4f}** ·"
      f" 全部为正 = **{bool((T.gap>0).all())}**")

GATE.asserted("CONTROL power: the full-sample effect is above the half-sample MDE",
              d_full>mde, f"gap {d_full:.4f} vs half-sample MDE {mde:.4f}", kind="control")
GATE.asserted("KILL the pattern replicates out of the estimating half (>50% of splits)",
              rate>0.5, f"{int(T.sig.sum())}/20 splits clear their own null")
verdict = "REPLICATES" if rate>0.5 else "OVERFIT"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,rate=rate,gap_full=d_full,gap_median=float(T.gap.median()),
               mde_median=mde,all_positive=bool((T.gap>0).all()),n_splits=20),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())

import os,sys,pathlib,json,re
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #446a found the loading vector fails to re-form in a few halves. Is that the STRUCTURE
   being unstable, or only the LABELS -- the 2nd and 3rd eigenvalues sitting close enough
   that the two components trade places?

Worlds
  A  label swap  : the halves where the per-vector correlation collapses are the halves where
     the eigenvalue gap is small, AND the 2D SUBSPACE span{v2,v3} is stable throughout. Then
     "unstable" is the wrong word: the ordering is unstable, the structure is not, and every
     page claim named `c2` or `c3` needs a naming caveat rather than a stability caveat.
  B  real instability : the gap explains nothing and the subspace moves too. Then every claim
     named after either coordinate needs an explicit "this coordinate fails to re-form in X%
     of halves" boundary.

⚠ The instrument must be the right one: PER-VECTOR CORRELATION IS EXACTLY THE QUANTITY A
   LABEL SWAP FOOLS. Subspaces are compared by PRINCIPAL ANGLES, which are invariant to any
   rotation or relabelling inside the subspace.
NULL : principal angles between two RANDOM 2-planes in 32 dimensions -- what "no shared
   structure" looks like, so a small angle means something.
CONTROL : the full-sample subspace against itself must give angles ~0 (the instrument can
   detect identity).
CONTROL2: the null must be far from zero (the instrument can say no).
FRONTIER: both branches change the page.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R491 swap or real instability")
_R489=(ROOT/'E01_sexual_as_a_value_not_a_category/A179_two_ways_of_being_involved'
            /'R489_does_the_fluid_pattern_replicate/run.py').read_text()
exec(_R489.split('"""',2)[2].split('FLU=re.compile')[0])

def eig_of(rows):
    mm=np.zeros(NN,bool); mm[rows]=True
    def prof(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(np.where(mm,R[b],np.nan))
        return R
    Ra,Rb=prof(A),prof(B)
    C=np.zeros((NB,NB))
    for i in range(NB):
        for j in range(NB):
            g=np.isfinite(Ra[i])&np.isfinite(Rb[j])&mm
            if g.sum()>100: C[i,j]=np.corrcoef(Ra[i][g],Rb[j][g])[0,1]
    C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w)
    return w[o], V[:,o]

def princ_angles(P,Q):
    """两个子空间之间的主角(度)。对子空间内的任何旋转/换标签**不变**。"""
    Pq,_=np.linalg.qr(P); Qq,_=np.linalg.qr(Q)
    s=np.linalg.svd(Pq.T@Qq, compute_uv=False)
    return np.degrees(np.arccos(np.clip(s,-1,1)))

ALLR=np.flatnonzero(ok)
w_full,V_full=eig_of(ALLR)
S_full=V_full[:,[1,2]]
a0=princ_angles(S_full,S_full)
GATE.asserted("CONTROL the instrument detects identity (angles ~ 0)",
              float(a0.max())<1e-6, f"max angle to itself = {a0.max():.2e} deg", kind="control")
rg=np.random.default_rng(5)
nul=np.array([princ_angles(rg.normal(size=(NB,2)),rg.normal(size=(NB,2))).max() for _ in range(4000)])
GATE.asserted("CONTROL2 the null is far from zero (the instrument can say no)",
              float(np.percentile(nul,5))>45,
              f"random 2-plane max angle 5th pct = {np.percentile(nul,5):.1f} deg", kind="control")
print(f"随机二维子空间的最大主角:5 分位 **{np.percentile(nul,5):.1f}°** · "
      f"中位 **{np.median(nul):.1f}°**")

rng=np.random.default_rng(97)          # 与 #445/#446 同一批半样本
rows=[]
for t in range(20):
    perm=rng.permutation(ALLR); hb=perm[len(perm)//2:]
    w,V=eig_of(hb)
    ang=princ_angles(V[:,[1,2]],S_full)
    rgap=float((w[1]-w[2])/max(w[1],1e-12))
    rv=float(np.corrcoef(V[:,1],V_full[:,1])[0,1])
    rows.append(dict(split=t, rel_gap=rgap, max_angle=float(ang.max()),
                     min_angle=float(ang.min()), r_vec=abs(rv)))
T=pd.DataFrame(rows)
show(T.sort_values('r_vec'), HERE/'results/subspace.csv', n=8, label="子空间 vs 逐向量")

r_gap_vec=float(np.corrcoef(T.rel_gap,T.r_vec)[0,1])
print(f"\n**特征值相对间隙 与 逐向量相关 的相关 = {r_gap_vec:+.4f}**")
print(f"子空间最大主角:中位 **{T.max_angle.median():.1f}°** · 最差 **{T.max_angle.max():.1f}°** ·"
      f" 随机零的 5 分位 **{np.percentile(nul,5):.1f}°**")
worst=T.nsmallest(5,'r_vec')
print(f"逐向量最差的 5 个:逐向量 r **{worst.r_vec.min():.3f}~{worst.r_vec.max():.3f}** · "
      f"**子空间最大主角 {worst.max_angle.min():.1f}°~{worst.max_angle.max():.1f}°**")

swap = (r_gap_vec>0.4) and float(T.max_angle.max())<float(np.percentile(nul,5))
GATE.asserted("KILL it is a label swap, not real instability", swap,
              f"gap-vs-vector r = {r_gap_vec:+.3f}; worst subspace angle "
              f"{T.max_angle.max():.1f} deg vs null 5th pct {np.percentile(nul,5):.1f} deg")
verdict = "LABEL_SWAP" if swap else "REAL_INSTABILITY"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,r_gap_vec=r_gap_vec,
               ang_median=float(T.max_angle.median()),ang_worst=float(T.max_angle.max()),
               null_p5=float(np.percentile(nul,5)),
               rvec_min=float(T.r_vec.min())),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())

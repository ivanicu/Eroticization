import os,sys,pathlib,json,re
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #447c found one direction of the (c2, c3) pair re-forms in every half and the other fails
   in about one half in twenty. Principal angles carry no labels, so which direction is the
   stable one was left open. What is it -- and is it the fluid direction?

Worlds
  A  the stable direction IS the fluid direction -> #446b's "content is sturdier" gets its
     mechanism: the content is sturdy BECAUSE it happens to sit on the well-conditioned
     direction of the pair, and the two facts merge into one.
  B  it is NOT -> then the fluid effect rides the FRAGILE direction and replicates anyway,
     which is STRONGER than A: the content would be independent of the coordinate rather than
     explained by it, and the page has to be upgraded, not merged.

PRE-REGISTERED PREDICTION: **A**. ⚠ Four of my six predictions this session were wrong
(#433b #434 #437 #442 vs #439 #445 right). Written to be killed.

Method: for each half, the SVD of (half subspace)^T (full subspace) gives paired directions
ordered by preservation. The FIRST pair spans the smallest principal angle -- that is the
stable direction; the SECOND is the fragile one. Both are projected back into block space.
⚠ Singular vectors carry a sign/rotation freedom inside the subspace, so every recovered
direction is anchored the same way before any comparison (#443b; five burns).
CONTROL : the stable directions must agree with each other far more than random block vectors
          do -- otherwise "stable" was not recovered and nothing below means anything.
CONTROL2: the fluid regex is untouched, character for character (#444's NEXT).
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R492 what is the stable direction")
_R489=(ROOT/'E01_sexual_as_a_value_not_a_category/A179_two_ways_of_being_involved'
            /'R489_does_the_fluid_pattern_replicate/run.py').read_text()
exec(_R489.split('"""',2)[2].split('FLU=re.compile')[0])
FLU=re.compile(r'urine|ejaculate|squirt|precum|saliva|secretion',re.I)   # ⚠ 一字未改
fl=np.array([bool(FLU.search(n)) for n in NAMES]); nf=int(fl.sum())

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

ALLR=np.flatnonzero(ok)
w_full,V_full=eig_of(ALLR); S_full=V_full[:,[1,2]]
Qf,_=np.linalg.qr(S_full)

def paired_dirs(rows):
    _,V=eig_of(rows); Ph,_=np.linalg.qr(V[:,[1,2]])
    U,s,Wt=np.linalg.svd(Ph.T@Qf)
    full_side=Qf@Wt.T                 # 全样本侧的配对方向,按保存度排序
    return full_side[:,0], full_side[:,1], s   # 稳 · 脆 · 奇异值

def anchor(v, ref):                    # 子空间内的符号自由度(#443b)
    return -v if float(v@ref)<0 else v

rng=np.random.default_rng(97)          # 与 #445–#447 同一批半样本
ST=[]; FR=[]; SV=[]
for t in range(20):
    perm=rng.permutation(ALLR); hb=perm[len(perm)//2:]
    a,b,s=paired_dirs(hb); ST.append(a); FR.append(b); SV.append(s)
ST=np.array(ST); FR=np.array(FR); SV=np.array(SV)
ref=ST[0].copy()
ST=np.array([anchor(v,ref) for v in ST])
ref2=FR[0].copy(); FR=np.array([anchor(v,ref2) for v in FR])

pw_st=np.array([abs(float(np.corrcoef(ST[i],ST[j])[0,1])) for i in range(20) for j in range(i+1,20)])
pw_fr=np.array([abs(float(np.corrcoef(FR[i],FR[j])[0,1])) for i in range(20) for j in range(i+1,20)])
rg=np.random.default_rng(11)
nul=np.array([abs(float(np.corrcoef(rg.normal(size=NB),rg.normal(size=NB))[0,1])) for _ in range(4000)])
print(f"稳方向 两两 |r| 中位 **{np.median(pw_st):.4f}** · 脆方向 **{np.median(pw_fr):.4f}** · "
      f"随机零 95 分位 **{np.percentile(nul,95):.4f}**")
GATE.asserted("CONTROL the stable directions really do agree with each other",
              float(np.median(pw_st))>float(np.percentile(nul,95)),
              f"median |r| {np.median(pw_st):.4f} vs null p95 {np.percentile(nul,95):.4f}", kind="control")
GATE.asserted("CONTROL2 the fluid regex is untouched", nf==6,
              f"{nf} blocks matched by wording", kind="control")

def gapz(v, seed):
    d=float(v[fl].mean()-v[~fl].mean())
    r2=np.random.default_rng(seed)
    nl=np.array([(lambda s: v[s].mean()-v[~s].mean())(
        np.isin(np.arange(NB),r2.choice(NB,nf,replace=False))) for _ in range(4000)])
    return abs(d)/max(float(nl.std()),1e-12)

zs=np.array([gapz(v,800+i) for i,v in enumerate(ST)])
zf=np.array([gapz(v,900+i) for i,v in enumerate(FR)])
T=pd.DataFrame([dict(q='稳方向 体液 |z| 中位',v=float(np.median(zs))),
                dict(q='稳方向 |z| 区间 2.5%',v=float(np.percentile(zs,2.5))),
                dict(q='稳方向 |z| 区间 97.5%',v=float(np.percentile(zs,97.5))),
                dict(q='脆方向 体液 |z| 中位',v=float(np.median(zf))),
                dict(q='脆方向 |z| 区间 2.5%',v=float(np.percentile(zf,2.5))),
                dict(q='脆方向 |z| 区间 97.5%',v=float(np.percentile(zf,97.5))),
                dict(q='稳方向两两|r| 中位',v=float(np.median(pw_st))),
                dict(q='脆方向两两|r| 中位',v=float(np.median(pw_fr)))])
show(T, HERE/'results/directions.csv', n=8, label="稳 vs 脆")

isA = float(np.median(zs))>float(np.median(zf))
GATE.asserted("KILL the stable direction is the fluid direction (预注册预测 A)", isA,
              f"stable fluid |z| median {np.median(zs):.2f} vs fragile {np.median(zf):.2f}")
verdict = "STABLE_IS_FLUID" if isA else "FLUID_RIDES_THE_FRAGILE_ONE"
print(f"\n判决 = {verdict}   (预注册预测 = A / STABLE_IS_FLUID)")
json.dump(dict(verdict=verdict,z_stable=float(np.median(zs)),z_fragile=float(np.median(zf)),
               pw_stable=float(np.median(pw_st)),pw_fragile=float(np.median(pw_fr)),
               null_p95=float(np.percentile(nul,95)),prediction="A"),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())

import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #443c found two halves of "being more involved" with different consequences. What are they
   made of? R372 did the block loadings for the third component only; the second component's
   loadings have never been computed, and #443c says that is the half that reaches agency.

⚠⚠ R372's own opening is the constraint on this round, quoted from it:
   "a loading ordering does NOT constitute a name for a construct -- #201 and #202 died exactly
   that way, telling a story off an ordering and then being killed by the data."
   **So the product here is a DESCRIPTIVE TABLE, not a naming.** Nothing below names anything.

Worlds
  A  same blocks, different weights : the two halves load on the same blocks -> "two kinds of
     involvement" narrows to "two weightings of one content set".
  B  different blocks               : then the two halves are two content sets, and what is in
     them is the first thing on this page that could be pointed at -- BY A LATER ROUND, not
     this one.

Method: R372's `load_of` spliced and generalised to any component; person-level bootstrap
  (300), each draw sign-aligned to the full-sample reference (#368a).
Orientation: both components turned to face MORE INVOLVEMENT via the count anchor, exactly as
  #443b did, so the two loading vectors are comparable.
CONTROL : the third component's loadings must reproduce R372's (r ~ 1 against its saved csv).
CONTROL2: the bootstrap must be able to say "no" -- report how many of 32 have intervals that
  DO contain zero.
CLOSURE (it describes; it decides nothing).
"""
import numpy as np, pandas as pd, warnings, re
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R488 what are the halves made of")
_R372=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like'
            /'R372_block_loadings/run.py').read_text()
exec(_R372.split('"""',2)[2].split('def load_of')[0])

def load_k(rows,k,ref=None):
    mm=np.zeros(NN,bool); mm[rows]=True
    def prof_(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(np.where(mm,R[b],np.nan))
        return R
    Ra,Rb=prof_(A),prof_(B)
    C=np.zeros((NB,NB))
    for i in range(NB):
        for j in range(NB):
            g=np.isfinite(Ra[i])&np.isfinite(Rb[j])&mm
            if g.sum()>200: C[i,j]=np.corrcoef(Ra[i][g],Rb[j][g])[0,1]
    C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
    v=V[:,k]
    if ref is not None and float(v@ref)<0: v=-v
    return v

ALLR=np.flatnonzero(ok)
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
SXA=pd.to_numeric(raw['Totalsexacts'],errors='coerce').values.astype(float)
def score_k(v):
    mm=np.zeros(NN,bool); mm[ALLR]=True
    def prof_(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo_=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo_,np.nan); R[b]=R[b]-np.nanmean(np.where(mm,R[b],np.nan))
        return R
    Ra,Rb=prof_(A),prof_(B)
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(v[:,None]*Zm).sum(0); den=(Fm*np.abs(v)[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)

REF={}
for k,nm in ((1,'c2'),(2,'c3')):
    v=load_k(ALLR,k); s=score_k(v)
    g=np.isfinite(s)&np.isfinite(SXA)&ok
    if np.corrcoef(s[g],SXA[g])[0,1]<0: v=-v      # 朝「更卷入」(#443b)
    s=score_k(v); g=np.isfinite(s)&np.isfinite(SXA)&ok
    REF[nm]=v
    print(f"{nm}(定向后)vs 性行为计数 = **{np.corrcoef(s[g],SXA[g])[0,1]:+.4f}**")

NBOOT=300; rg=np.random.default_rng(2929)
BS={nm:np.array([load_k(ALLR[rg.integers(0,len(ALLR),len(ALLR))],k,ref=REF[nm])
                 for _ in range(NBOOT)]) for k,nm in ((1,'c2'),(2,'c3'))}
rows=[]
for i in range(NB):
    r=dict(block=i, name=NAMES[i])
    for nm in ('c2','c3'):
        lo,hi=np.percentile(BS[nm][:,i],[2.5,97.5])
        r[f'{nm}_load']=float(REF[nm][i]); r[f'{nm}_lo']=float(lo); r[f'{nm}_hi']=float(hi)
        r[f'{nm}_sig']=bool((lo>0)or(hi<0))
    rows.append(r)
T=pd.DataFrame(rows)
T['gap']=(T.c2_load-T.c3_load).abs()
T.to_csv(HERE/'results/loadings_both.csv',index=False)

r_shape=float(np.corrcoef(REF['c2'],REF['c3'])[0,1])
print(f"\n**两条载荷向量的相关 = {r_shape:+.4f}**(32 个块)")
print(f"区间不含零:c2 **{int(T.c2_sig.sum())}/32** · c3 **{int(T.c3_sig.sum())}/32**")
print(f"⚠ 含零的:c2 **{32-int(T.c2_sig.sum())}** · c3 **{32-int(T.c3_sig.sum())}** -> **自助说得出「不」**")

print("\n**差别最大的 8 块(描述性,不是命名 —— `R372` 的警告)**")
show(T.sort_values('gap',ascending=False)[['block','c2_load','c3_load','gap','c2_sig','c3_sig','name']],
     HERE/'results/biggest_gap.csv', n=8, label="载荷差")

try:
    old=pd.read_csv('E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/'
                    'R372_block_loadings/results/loadings.csv').sort_values('v_block')
    r_old=float(np.corrcoef(np.abs(REF['c3']),np.abs(old.v_load.values))[0,1])
except Exception as e:
    r_old=float('nan')
GATE.asserted("CONTROL the third component's loadings reproduce R372's (up to sign)",
              (r_old==r_old) and abs(r_old)>0.95, f"|loading| corr vs R372 = {r_old:+.4f}", kind="control")
GATE.asserted("CONTROL2 the bootstrap can say no",
              int(T.c2_sig.sum())<NB and int(T.c3_sig.sum())<NB,
              f"intervals containing zero: c2 {NB-int(T.c2_sig.sum())}, c3 {NB-int(T.c3_sig.sum())}",
              kind="control")
diff = abs(r_shape)<0.5
GATE.asserted("KILL the two halves load on different blocks", diff,
              f"loading-vector correlation = {r_shape:+.4f}")
verdict = "DIFFERENT_BLOCKS(⚠ 见 #444b:该判据不可能失败)" if diff else "SAME_BLOCKS_DIFFERENT_WEIGHTS"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,r_shape=r_shape,c2_sig=int(T.c2_sig.sum()),
               c3_sig=int(T.c3_sig.sum()),NB=int(NB),r_vs_R372=r_old),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())

# ---------------------------------------------------------------- #444b/#444c(发布前追加)
# ⚠ #444b: the KILL above compared the two loading vectors' correlation. They are eigenvectors
# of the SAME symmetric matrix, so v2 . v3 = 0 EXACTLY -- the test could not fail (§9's
# arithmetic trap). Kept on the record, replaced by two tests that CAN fail.
print(f"\n⚠ **#444b 正交是构造的**:v2·v3 = {float(REF['c2']@REF['c3']):+.2e} -> 上面那个 KILL 不可能失败。")
rng2=np.random.default_rng(7); ov=[]
for k in (4,6,8,10,12):
    a=set(np.argsort(-np.abs(REF['c2']))[:k]); b=set(np.argsort(-np.abs(REF['c3']))[:k])
    nul=[]
    for _ in range(4000):
        X=rng2.normal(size=(NB,2)); Q,_=np.linalg.qr(X)
        nul.append(len(set(np.argsort(-np.abs(Q[:,0]))[:k])&set(np.argsort(-np.abs(Q[:,1]))[:k])))
    nul=np.array(nul)
    ov.append(dict(k=k,overlap=len(a&b),null_mean=float(nul.mean()),
                   null_p95=float(np.percentile(nul,95))))
OV=pd.DataFrame(ov); show(OV, HERE/'results/topk_overlap.csv', n=6, label="高载荷块重叠")
GATE.asserted("KILL(重做,可失败):高载荷块集合的重叠越出随机正交零",
              bool((OV.overlap>OV.null_p95).any()),
              f"overlaps {list(OV.overlap)} vs null p95 {list(OV.null_p95)}")

# #444c: a content test defined by WORDING ONLY (never by the loadings)
FLU=re.compile(r'urine|ejaculate|squirt|precum|saliva|secretion',re.I)
fl=np.array([bool(FLU.search(n)) for n in NAMES])
rows2=[]
for cn in ('c2','c3'):
    v=REF[cn]; d=float(v[fl].mean()-v[~fl].mean())
    nul=np.array([(lambda s: v[s].mean()-v[~s].mean())(
        np.isin(np.arange(NB),rng2.choice(NB,int(fl.sum()),replace=False))) for _ in range(20000)])
    rows2.append(dict(comp=cn, fluid_mean=float(v[fl].mean()), rest_mean=float(v[~fl].mean()),
                      diff=d, z=d/float(nul.std()), p95=float(np.percentile(np.abs(nul),95)),
                      sig=bool(abs(d)>np.percentile(np.abs(nul),95))))
F=pd.DataFrame(rows2); show(F, HERE/'results/fluid_test.csv', n=4, label="体液块")
GATE.asserted("CONTROL3 the fluid set is defined by wording and includes a non-fitting member",
              bool(fl.sum()==6 and REF['c2'][[i for i in range(NB) if 'secretions' in NAMES[i]][0]]<0.1),
              "the general 'bodily secretions' block is in the set and does NOT fit", kind="control")
print(f"\n⚠ **命名边界**(`R372` 的开头):一个载荷排序**不构成**一个构念的名字 —— 本轮不命名。")

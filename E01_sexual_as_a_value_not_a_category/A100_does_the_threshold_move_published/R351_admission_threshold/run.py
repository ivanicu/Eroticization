import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A100 R351 -- 那个块对准入门槛,动不动**已发表的**数

`#305a`:`Cmat` 里 `mm.sum() > 200` 这个块对准入门槛,能把 `c1` 的重现度从 0.68 撬到 0.96。
**那它动不动 `c3` 本身、动不动公开页面上的 `c3 ↔ 羞耻 = −0.1286`?**
若动,`#231` 起的每一个 c 系数都有一个**未列出的旋钮**。

ESTIMAND        全样本下门槛 ∈ {0, 100, 200(发布口径), 400, 800},每次重估 `c1/c2/c3`,报
                ① 与门槛 200 版本的 |cos| ② 各自 ↔ 羞耻的相关 ③ 被置零的格子比例。
KILL            **若 `c3 ↔ 羞耻` 跨门槛可分辨地变 -> 进 `CALIBER.md`,且公开页面那个数要带上它;
                若四个门槛给同一个数 -> 这条旋钮出局,而 `#305a` 的撬动要归因于 n 小,不是门槛本身
                —— 那也是一个结果。**
POSITIVE CTRL   门槛 0(全部格子)必须与门槛 100 几乎一样(否则是别的东西在动)。
NEGATIVE CTRL   **`S` 不经过 `C`** —— 它必须跨门槛**逐位不变**。这是一个内建的、免费的负对照:
                若 `S` 也动了,那动的不是门槛。
IMPOSSIBLE      门槛只影响 `C` 的**哪些格子被填**;它不影响块内剖面本身。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def fit_apply')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); m=np.zeros(NN,bool); m[ALLR]=True
def prof_(X):
    F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0)
    R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(np.where(m,R[b],np.nan))
    return R
Ra,Rb=prof_(A),prof_(B)
CNT=np.zeros((NB,NB)); RAW=np.full((NB,NB),np.nan)
for i in range(NB):
    for j in range(NB):
        mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
        CNT[i,j]=mm.sum()
        if mm.sum()>2: RAW[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
def cs_for(thr,ref=None):
    C=np.where((CNT>thr)&np.isfinite(RAW),RAW,0.0); C=(C+C.T)/2
    w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
    # ⚠ 特征向量的符号是**任意的**(R210:73 的老坑)。不对齐的话,门槛 400 会印出一个
    #    「符号翻转」的相关,而同时 |cos| 还有 0.9466 —— 那两件事不可能同时为真,
    #    因为 |cos| 取了绝对值。对齐到发布口径(门槛 200)那一版。
    if ref is not None:
        for k in range(V.shape[1]):
            if float(V[:,k]@ref[:,k])<0: V[:,k]=-V[:,k]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0); cs=[]
    for k in range(3):
        num=(V[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(V[:,k])[:,None]).sum(0)
        cs.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    return V,cs,float(((CNT<=thr).sum())/(NB*NB))
def cor(u,v):
    k=np.isfinite(u)&np.isfinite(v)&m
    return float(np.corrcoef(u[k],v[k])[0,1]) if k.sum()>200 else np.nan
# 负对照:S 不经过 C
cvS=np.zeros(NN); psS=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    v=np.where(n>0,(M@rr)/np.maximum(n,1),np.nan); g=np.isfinite(v)
    cvS[ppl[g]]+=1; psS[ppl[g]]+=v[g]
S=np.where(cvS>=8,psS/np.maximum(cvS,1),np.nan)
S_sh=cor(S,sh)
print(f"块对计数:中位 **{int(np.median(CNT)):,}** · 最小 **{int(CNT.min()):,}** · "
      f"发布门槛 200 排掉 **{100*((CNT<=200).sum())/(NB*NB):.1f}%** 的格子\n")
V200,_,_=cs_for(200)
rows=[]
for thr in (0,100,200,400,800):
    V,cs,z=cs_for(thr,ref=V200)
    co=[float(V[:,k]@V200[:,k]) for k in range(3)]   # 已对齐,所以直接读带号的内积
    sc=[cor(cs[k],sh) for k in range(3)]
    rows.append(dict(thr=thr,zero=z,cos1=co[0],cos2=co[1],cos3=co[2],
                     sh1=sc[0],sh2=sc[1],sh3=sc[2],s_shame=cor(S,sh)))
    print(f"  门槛 {thr:>4}:置零 {100*z:>4.1f}% · 与 200 版的 |cos| "
          f"c1 **{co[0]:.4f}** c2 **{co[1]:.4f}** c3 **{co[2]:.4f}** · "
          f"↔羞耻 c1 {sc[0]:+.4f} c2 {sc[1]:+.4f} **c3 {sc[2]:+.4f}**")
T=pd.DataFrame(rows); check_columns(T,'R351')
T.to_csv(pathlib.Path(__file__).parent/'results'/'threshold.csv',index=False)
rng3=float(T.sh3.max()-T.sh3.min()); rngC=float(T.cos3.min())
print(f"\n★ `c3 ↔ 羞耻` 跨五个门槛:**{T.sh3.min():+.4f} … {T.sh3.max():+.4f}**(极差 **{rng3:.4f}**)")
print(f"★ `c3` 方向与发布版的 |cos| 最小:**{rngC:.4f}**")
print(f"⚠ 负对照 `S ↔ 羞耻` 跨门槛:{' · '.join(f'{v:+.6f}' for v in T.s_shame)} "
      f"(极差 {float(T.s_shame.max()-T.s_shame.min()):.2e})")
p0=T[T.thr==0].iloc[0]; p1=T[T.thr==100].iloc[0]
gg=Gate('那个准入门槛动不动已发表的数')
r200=float(T[T.thr==200].iloc[0].sh3)
for _,r in T.iterrows():
    gg.sign_flip_needs_direction_change(f'⚠ guard 20:门槛 {int(r.thr)} 的方向一致性 vs 派生量符号',
                                        r.cos3,r200,r.sh3)
USE=T[T.cos3.abs()>=0.95]
gg.asserted('★ 在方向还认得出来的区间(|cos| >= 0.95)内,`c3 ↔ 羞耻` 的极差',
            float(USE.sh3.max()-USE.sh3.min())<0.036,
            f"门槛 {list(USE.thr.astype(int))} -> 极差 **{float(USE.sh3.max()-USE.sh3.min()):.4f}**"
            f"(`#348` 的嵌套噪声 0.0359)")
gg.asserted('★ 正对照:门槛 0 与门槛 100 必须几乎一样',
            abs(p0.sh3-p1.sh3)<0.01 and abs(p0.cos3-p1.cos3)<0.02,
            f"c3↔羞耻 {p0.sh3:+.4f} vs {p1.sh3:+.4f} · |cos| {p0.cos3:.4f} vs {p1.cos3:.4f}")
gg.asserted('★ 负对照:`S` 不经过 `C`,必须跨门槛逐位不变',
            float(T.s_shame.max()-T.s_shame.min())<1e-12,
            f"S↔羞耻 极差 {float(T.s_shame.max()-T.s_shame.min()):.2e} —— 若不为零,动的不是门槛")
gg.asserted('★ 注册的 kill:`c3 ↔ 羞耻` 跨门槛是否可分辨地变(以 `#348` 的嵌套噪声 0.036 为尺)',
            rng3>0.036,
            f"极差 **{rng3:.4f}**;`#348` 里 c3 的嵌套-样本内差是 0.0359 —— "
            f"比它大就说明这条旋钮不比「换一批人」小")
gg.asserted('⚠ 方向本身动了多少',True,
            ' · '.join(f"门槛 {int(r.thr)} |cos|={r.cos3:.4f}" for _,r in T.iterrows()))
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

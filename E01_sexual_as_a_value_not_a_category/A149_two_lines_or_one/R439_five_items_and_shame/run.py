import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A149 R439 -- 「常规的东西对他不色」这条线,和羞耻那条线,是同一条吗

页面上现在有两条线:
**羞耻线**:羞耻 -> 少做(`#384b`);做过 -> 觉得改不掉(`#386a`);**羞耻完全不预测「能不能改」**(`#384a`)。
**可塑性线**:常规的东西对他越不色 -> 越觉得能停下来(`#392c` 翻正方向后的 `#390a`/`#393a`)。

**两个活着的世界:**
**A 同源** —— 常规不色的人也更羞耻(两者都来自「我的性是偏的」这个感觉)-> `corr(五题, 羞耻)` **正**且可观;
**B 正交** —— 两条线**互不相干**,而那本身是页面上一个**结构性**的发现 -> **≈ 0** 且 MDE 够小。

ESTIMAND        `羞耻 ~ 五题 + S + c3⁻ + 类别数`,主量 = **五题的系数**。
判据(**先标支**,`#379c`)
                【两支】guard 26 用 **MDE 扫描**;负对照;offset 零非退化;
                        **阳性参照**:`c3⁻ ↔ 羞耻`(页面已知 +0.129)必须在同一模型里开火。
                【非零支】越阈 -> 世界 A(同源),并报符号。
                【零支】未越阈**且 MDE < 0.05** -> 世界 B(正交),**这个零可发布**。
⚠ 零的种类     `offset_control`:五题与 `S`、类别数都相关(`S` 本身就是稀有度)->
                **这个零不该是零**;零 = **控制 `S`/`c3⁻`/类别数之后**、用 `lib.nulls.perm_in`
                在掩码内打乱五题分数的分布。
⚠ 新基础设施   **零的构造全部从 `lib.nulls` 取,不手写**(`#394e` 的可证伪检验,本轮是第一次)。
⚠ `#392e`      五题分数进模型前,先打印**取值集合、众数、与方向已知锚的相关**。
IMPOSSIBLE      ① 五题与 `S` 相关 -> 控制 `S` 可能过度控制(两者都含「常规 vs 冷门」)->
                   **原始与控制后都报**;② 全自报,同源方差。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, sham_control, controls as null_controls
nc=null_controls(); nc.pop('_detail')
print(f"⚠ `lib.nulls.controls()` 用前自检:**{'全部通过' if all(nc.values()) else nc}**\n")
assert all(nc.values()), nc
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
INT=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i]) for i in range(5)]),1)

# ---- ⚠ `#392e` 的规矩:进模型前先看它自己 ----
u=np.unique(V[np.isfinite(V)]); anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
ga=np.isfinite(INT)&np.isfinite(anc)
print(f"⚠ **`#392e`:五题分数进模型前先看它自己**")
print(f"   原始题的取值集合 = {u.tolist()} · 众数 = **{float(pd.Series(V[np.isfinite(V)]).mode().iloc[0]):g}**")
print(f"   与方向已知的锚 `Totalsexacts` 相关 = **{np.corrcoef(INT[ga],anc[ga])[0,1]:+.4f}** "
      f"-> **分数越高 = 常规的东西越不色**(`#392c` 已定,这里复核)\n")

M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(INT)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def fit(y,xs,g=None):
    g=M if g is None else g; k=int(g.sum()); gg=g&np.isfinite(y)
    for v in xs: gg=gg&np.isfinite(v)
    k=int(gg.sum())
    X=np.column_stack([np.ones(k)]+[z(v,gg) for v in xs]); yy=z(y,gg)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
CTRL=[S,C3,ncat]
bRaw,_=fit(sh,[INT]); bAdj,sAdj=fit(sh,[INT]+CTRL)
print(f"n=**{n:,}** · corr(五题, `S`) = **{np.corrcoef(INT[M],S[M])[0,1]:+.4f}** · "
      f"corr(五题, 类别数) = **{np.corrcoef(INT[M],ncat[M])[0,1]:+.4f}**")
print(f"五题 -> 羞耻:原始 **{bRaw[0]:+.4f}** · 控制 `S`/`c3⁻`/类别数后 **{bAdj[0]:+.4f}** (se {sAdj[0]:.4f})")
NP_=400
NUL=np.array([fit(perm_in(sh,M,7700+s),[INT]+CTRL)[0][0] for s in range(NP_)])   # ★ lib.nulls
THR=float(np.percentile(np.abs(NUL),95))
print(f"\n⚠ offset 零(**`lib.nulls.perm_in`,掩码内打乱**;这个零不该是零 —— "
      f"五题与 `S`/类别数都相关):**{NUL.mean():+.5f} ± {NUL.std():.5f}** · |值| 95 分位 **{THR:.5f}**")
print(f"   实测 **{bAdj[0]:+.4f}** -> **{(bAdj[0]-NUL.mean())/max(NUL.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(bAdj[0])>THR else '**未越阈**'}")
bPOS,sPOS=fit(sh,[C3,S,ncat])
print(f"\n★ 阳性参照(**真实数据**,页面已知 `c3⁻ ↔ 羞耻` ≈ +0.129):**{bPOS[0]:+.4f}** vs 阈 {THR:.4f} -> "
      f"{'**开火**' if abs(bPOS[0])>THR else '**不开火**'}")
# ⚠ **第一版的负对照是错的设计,而它当场就 FAIL 了 —— 按构造 5% 的时候会。**
# 它是**一次抽样**,而阈**就是同一个分布的 95 分位** ——
# **一个从建阈的同一个分布里抽出来的负对照,不可能提供信息:它 5% 的时候「失败」,而那正是阈的定义。**
# 正确的负对照是**越阈率**:多次独立抽样越过阈的比例应当 ≈ 5%。
NEGN=200
negs=np.array([fit(perm_in(sh,M,50000+s),[INT]+CTRL)[0][0] for s in range(NEGN)])
rate=float((np.abs(negs)>THR).mean())
print(f"  负对照(**越阈率**,`perm_in` 打乱羞耻 {NEGN} 次):**{100*rate:.1f}%** "
      f"(阈按定义应给 ≈5%;合格区间 1–12%)")
bNEG=[float(np.median(np.abs(negs)))]
print(f"\nguard 26 的正对照 = **MDE 扫描**(`#384d`),每级 30 次:")
MDE=None
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(3500+int(gg*100)*53+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(INT,M)+rg.standard_normal(n)
        if abs(fit(y,[INT]+CTRL)[0][0])>THR: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
MEANINGFUL=0.05
print(f"   **MDE = {MDE_:.2f}** · 有意义 **{MEANINGFUL:.2f}**")
pd.DataFrame([dict(v_raw=bRaw[0],v_adj=bAdj[0],v_se=sAdj[0],v_thr=THR,v_mde=MDE_,
                   v_pos=bPOS[0],v_n=n)]).to_csv(pathlib.Path(__file__).parent/'results'/'five_shame.csv',index=False)
g=Gate('「常规不色」这条线和羞耻那条线,是同一条吗')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 有意义幅度',MDE_,MEANINGFUL,True,what='MDE 扫描 80% 检出')
g.asserted('★【两支】阳性参照:`c3⁻ ↔ 羞耻` 必须在同一模型里开火',abs(bPOS[0])>THR,
           f"{bPOS[0]:+.4f} vs {THR:.4f}",kind='control')
g.asserted('★【两支】负对照:**越阈率**必须 ≈5%(单次抽样的负对照是空的 —— 阈就是它的 95 分位)',
           0.01<=rate<=0.12,f"越阈率 {100*rate:.1f}% · 中位 |b| {bNEG[0]:.4f} vs 阈 {THR:.4f}",kind='control')
g.asserted('★【两支】offset 零非退化',NUL.std()>0,f"{NUL.mean():+.5f} ± {NUL.std():.5f}",kind='control')
if MDE_<=MEANINGFUL and abs(bPOS[0])>THR and 0.01<=rate<=0.12:
    if abs(bAdj[0])>THR:
        g.asserted('★【非零支】越阈 -> 世界 A(同源)',True,
                   f"{bAdj[0]:+.4f} · 符号 {'正' if bAdj[0]>0 else '负'}")
    else:
        g.asserted('★【零支】未越阈且 MDE < 有意义 -> 世界 B(两条线正交)',MDE_<MEANINGFUL,
                   f"{bAdj[0]:+.4f} vs {THR:.4f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

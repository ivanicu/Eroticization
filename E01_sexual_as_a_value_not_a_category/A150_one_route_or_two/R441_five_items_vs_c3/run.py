import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A150 R441 -- 五题(常规不色)与 `c3⁻`(页面上最大的羞耻路径),是同一件事吗

`#395a`:五题分数预测羞耻 **+0.0657**。`c3⁻` 预测羞耻 **+0.129**(页面)。
**而 `corr(五题, c3⁻)` 从没被估过** —— 所以这**不是** `#396b` 那一类空分离器:**这条腿是新的。**

两个活着的世界:
**A 同一件事的两种测法** -> 相关高,且放进同一模型时**互相吃掉** -> **页面上那两条路要合并**;
**B 两个不同的东西** -> 相关低,且**各自存活** -> **页面上多一条独立的羞耻路径**。

⚠ **`c3⁻` 是特征向量,符号是任意的**(guard 24,`#368a` 是本项目第四次栽在这上面)——
**进模型前先定向**,并按 `#392e` 打印它的取值范围与它与方向已知锚的相关。

ESTIMAND        ① `corr(五题, c3⁻)`;② `羞耻 ~ 五题 + c3⁻ + S + 类别数`,两者的**偏系数**。
判据(**先标支**,`#379c`)
                【两支】guard 24(`c3⁻` 已定向)· guard 26(**MDE 扫描**)·
                        负对照用**越阈率**(`#395b`)· offset 零非退化。
                【非零支】**两者都越阈** -> 世界 B(两条独立的路);
                          **只剩一个** -> 世界 A(合并),并报是哪一个活下来。
                【零支】两者都塌 -> 共线到无法分辨,那时报「这份数据分不开」。
⚠ 零的种类     `offset_control`:两者都与 `S`/类别数相关 -> **这个零不该是零**;
                零 = **控制之后**、用 `lib.nulls.perm_in` 在掩码内打乱的分布(`#394e` 第二轮)。
IMPOSSIBLE      ① `c3⁻` 是在**同一份数据**上估的坐标(`#347`:坐标本身会过拟合)-> 绝对值是上界;
                ② 五题是**评分**,`c3⁻` 是**块剖面的特征向量投影** —— 量纲不同,**只比是否存活,不比大小**;
                ③ 「不是同一件事」不等于「因果上独立」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
print(f"⚠ `lib.nulls.controls()` 用前自检:**全部通过**\n")
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
INT=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i]) for i in range(5)]),1)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(INT)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)

# ---- ⚠ `#392e` + guard 24:两个量进模型前,各自先看清楚 ----
gA=Gate('两个量进模型前先各自看清楚')
anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
gI=M&np.isfinite(anc)
print("⚠ **`#392e`:进模型前先看它自己**")
print(f"   五题分数:范围 [{np.nanmin(INT[M]):+.3f}, {np.nanmax(INT[M]):+.3f}] · "
      f"与方向已知锚 `Totalsexacts` 相关 **{np.corrcoef(INT[gI],anc[gI])[0,1]:+.4f}** "
      f"-> **分数越高 = 常规越不色**(`#392c`)")
print(f"   `c3⁻`:范围 [{np.nanmin(C3[M]):+.3f}, {np.nanmax(C3[M]):+.3f}] · "
      f"⚠ **它是特征向量,符号任意** -> guard 24 定向:")
_an=gA.eigenvector_is_anchored('★ `c3⁻` 已对着羞耻定向(约定:指向更多羞耻)',C3,sh,'羞耻')
print(gA)
assert _an, "c3⁻ 未定向 —— 不许在未定向的分量上写任何标签(#368a)"

print(f"\nn=**{n:,}**")
R=float(np.corrcoef(INT[M],C3[M])[0,1])
print(f"① **`corr(五题, c3⁻)` = {R:+.4f}**  ← 这条腿**从没被估过**,所以它是新信息(不是 `#396b` 那一类)")
def fit(y,xs):
    X=np.column_stack([np.ones(n)]+[z(v,M) for v in xs]); yy=z(y,M)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(n-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
bI,sI=fit(sh,[INT,S,ncat]); bC,sC=fit(sh,[C3,S,ncat]); bJ,sJ=fit(sh,[INT,C3,S,ncat])
print(f"\n② 单独放:五题 **{bI[0]:+.4f}** (se {sI[0]:.4f}) · `c3⁻` **{bC[0]:+.4f}** (se {sC[0]:.4f})")
print(f"   一起放:五题 **{bJ[0]:+.4f}** (se {sJ[0]:.4f}) · `c3⁻` **{bJ[1]:+.4f}** (se {sJ[1]:.4f})")
NP_=400
nI=[];nC=[]
for s_ in range(NP_):
    b2,_=fit(perm_in(sh,M,8800+s_),[INT,C3,S,ncat]); nI.append(b2[0]); nC.append(b2[1])
nI=np.array(nI); nC=np.array(nC)
TI=float(np.percentile(np.abs(nI),95)); TC=float(np.percentile(np.abs(nC),95))
print(f"\n⚠ offset 零(**`lib.nulls.perm_in`** 掩码内打乱羞耻;这个零不该是零 —— 两者都与 `S`/类别数相关):")
print(f"   五题阈 **{TI:.4f}** · `c3⁻` 阈 **{TC:.4f}**")
OI=abs(bJ[0])>TI; OC=abs(bJ[1])>TC
print(f"   -> 五题 {'**越阈**' if OI else '**塌掉**'} · `c3⁻` {'**越阈**' if OC else '**塌掉**'}")
NEGN=200
negs=np.array([fit(perm_in(sh,M,60000+s),[INT,C3,S,ncat])[0][0] for s in range(NEGN)])
rate=float((np.abs(negs)>TI).mean())
print(f"\n负对照(**越阈率**,`#395b` —— 单次抽样的负对照是空的):**{100*rate:.1f}%**(合格 1–12%)")
print(f"\nguard 26 的正对照 = **MDE 扫描**,每级 30 次(种在五题上):")
MDE=None
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(2600+int(gg*100)*59+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(INT,M)+rg.standard_normal(n)
        if abs(fit(y,[INT,C3,S,ncat])[0][0])>TI: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 = 单独放时的 **{max(abs(bI[0]),abs(bC[0])):.4f}**")
pd.DataFrame([dict(v_what='五题',v_alone=bI[0],v_joint=bJ[0],v_thr=TI,v_over=bool(OI)),
              dict(v_what='c3⁻',v_alone=bC[0],v_joint=bJ[1],v_thr=TC,v_over=bool(OC))]).to_csv(
    pathlib.Path(__file__).parent/'results'/'two_routes.csv',index=False)
g=Gate('五题与 c3⁻ 是同一件事吗')
g.eigenvector_is_anchored('★【两支】guard 24:`c3⁻` 已定向',C3,sh,'羞耻')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 单独放时的效应',MDE_,max(abs(float(bI[0])),abs(float(bC[0]))),True,
    what='MDE 扫描 80% 检出')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.asserted('★【两支】offset 零非退化',nI.std()>0 and nC.std()>0,
           f"五题 {nI.std():.4f} · c3⁻ {nC.std():.4f}",kind='control')
if 0.01<=rate<=0.12 and MDE_<=max(abs(bI[0]),abs(bC[0])):
    g.asserted('★【非零支】两者**都**越阈 -> 世界 B(两条独立的路)',OI and OC,
               f"五题 {bJ[0]:+.4f}/{TI:.4f} · c3⁻ {bJ[1]:+.4f}/{TC:.4f} · "
               f"corr = {R:+.4f} -> **世界 {'B 两条独立的路' if OI and OC else ('A 合并' if OI!=OC else '共线到分不开')}**")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

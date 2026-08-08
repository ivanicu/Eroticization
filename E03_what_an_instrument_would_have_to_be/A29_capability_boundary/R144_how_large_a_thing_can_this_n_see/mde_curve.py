"""E03·A29·R144 —— 这个 n 能看见多大的东西

**类型:FRONTIER(它改的是能力边界 —— §0 的第 ④ 个更新目标)。A28 关弧,A29 开弧。**
**心理学的那一句(本轮要判的):社会侧这条线上,「测不到」到底是「没有」,还是「这个 n 看不见」?**

## 仪器(硬规则②/④)
唯一的仪器是 SCCS `barry1977agents` 的五个管教块。**没有第二具仪器**
—— `#700` 已枚举 D-PLACE 全部 12 个档案,同时有管教与性规范的只此一具。

## G1 ESTIMAND
统计量沿用 `#699`:**六对上「预测符号 × 观测带号ρ」的均值**(语义载荷,不重新拟合)。
**MDE(n) = 该 n 下 halo 零的 95% 分位** —— 即「观测要多大才刚好越过」。
**逐 n 报:n ∈ {40, 50, 61, 75, 97}。**
## ⑤ 最强混淆(`#701` 预注册)
**零的宽度不只随 n 变,也随「五块本来有多像」变** ⇒
**同时报每个 n 上「打乱前 vs 打乱后」的块间平均 |ρ|**,**不许把宽度全归给 n**。
## ⑧ 判据(`#701` 写死)
**若 n=61 的 MDE 高于本页面社会侧所有已报效应 ⇒ 记「社会侧的分辨率天花板已经到了」,写进「做不到的」。**
## ④ 正对照
**`#699` 的 +0.3243(n=97)必须落在曲线的「可检测」一侧**,否则曲线算错了,当场停。
## IMPOSSIBLE(不写 planned)
SCCS 只有 186 个社会,五块交集 97 ⇒ **n 不可能再大**;
**本轮量的是分辨率,不是任何声明的真伪**;因果:横断面民族志编码。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
Dd=pd.read_csv(B+"data.csv")
W=Dd.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
BL={"以身作则":[429,430,431,432],"讲道理":[437,438,439,440],"体罚":[453,454,455,456],
    "放任":[465,466,467,468],"疼爱":[469,470,471,472]}
LOAD={"体罚":+1,"讲道理":+1,"以身作则":0,"放任":-1,"疼爱":-1}
S=pd.DataFrame({k:W[[f"SCCS{i}" for i in v]].mean(axis=1) for k,v in BL.items()}).dropna()
print(f"五块皆有分的社会 **{len(S)}** 个 —— 这是 n 的上限,SCCS 给不出更多")
def sp(a,b): return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def nrm(x,y):
    r=sp(x,y)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(np.asarray(x,float)); ys=np.sort(np.asarray(y,float)); ys=ys if r>0 else ys[::-1]
    c=sp(xs,ys)
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
cols=list(S.columns)
def stat(M):
    vals=[]
    for a,b in combinations(cols,2):
        p=LOAD[a]*LOAD[b]
        if p==0: continue
        v=nrm(M[a],M[b])
        if np.isfinite(v): vals.append(p*v)
    return float(np.mean(vals)) if vals else np.nan
def absmean(M):
    vals=[abs(nrm(M[a],M[b])) for a,b in combinations(cols,2)]
    vals=[v for v in vals if np.isfinite(v)]
    return float(np.mean(vals)) if vals else np.nan
obs_full=stat(S)
print(f"观测(n={len(S)})= **{obs_full:+.4f}**(`#699` 报 +0.3243)\n")
rng=np.random.default_rng(20260806)
rows=[]
for n in (40,50,61,75,97):
    if n>len(S): continue
    nul=[]; sim_pre=[]; sim_post=[]
    reps = 60 if n<97 else 1
    per  = 400 if n<97 else 2000
    for _ in range(reps):
        idx=rng.choice(len(S),size=n,replace=False) if n<len(S) else np.arange(len(S))
        sub=S.iloc[idx]
        sim_pre.append(absmean(sub))
        A=sub.to_numpy(float)
        for _ in range(per//reps if reps>1 else per):
            Ap=A.copy()
            for i in range(Ap.shape[0]): Ap[i,:]=rng.permutation(Ap[i,:])
            P=pd.DataFrame(Ap,index=sub.index,columns=cols)
            v=stat(P)
            if np.isfinite(v): nul.append(v)
            if len(sim_post)<200: sim_post.append(absmean(P))
    nul=np.array(nul); q=float(np.quantile(nul,.95))
    rows.append((n,q,float(np.median(nul)),float(np.mean(sim_pre)),float(np.mean(sim_post))))
    print(f"  n={n:>3d}  **MDE = halo 零 95% 分位 = {q:+.4f}**  零中位 {np.median(nul):+.4f}  "
          f"块间平均|ρ| 打乱前 {np.mean(sim_pre):.3f} → 打乱后 {np.mean(sim_post):.3f}")
mde61=[r[1] for r in rows if r[0]==61][0]
PAGE={"#699 轴统计量(n=97)":0.3243,"#700 轴×谴责(n=75)":0.3631,"#701 轴×谴责(n=61)":0.2549,
      "#701 轴×罕见(n=61)":0.2769,"#698 跨手段|ρ|中位(n=97)":0.3241}
print(f"\n③ 判据:n=61 的 MDE = **{mde61:.4f}**,对照页面上社会侧已报效应:")
below=[k for k,v in PAGE.items() if v<mde61]
for k,v in PAGE.items(): print(f"   {k:26s} {v:+.4f}  {'⛔ 低于 MDE' if v<mde61 else '✅ 可检测'}")
G=Gate("这个 n 能看见多大的东西")
p1=G.positive_control("#699 的 +0.3243(n=97)必须落在可检测一侧",
                      planted=0.3243,floor=[r[1] for r in rows if r[0]==97][0],spread=0.01)
if p1:
    v=(f"**社会侧的分辨率天花板已经到了:n=61 的 MDE 是 {mde61:.4f},"
       f"页面上 {len(below)}/{len(PAGE)} 条社会侧效应低于它 ⇒ 在 61 个社会上,"
       f"这条线报的多数数字都在自己的分辨率以下**" if below else
       f"**还有余量:n=61 的 MDE 是 {mde61:.4f},页面上社会侧效应全部高于它**")
else: v="UNVERIFIED —— 正对照失败,曲线算错了"
print(f"\n{v}"); print(G)
json.dump(dict(n_max=int(len(S)),obs_full=obs_full,
               curve=[{"n":n,"mde":q,"null_median":m,"sim_pre":a,"sim_post":b} for n,q,m,a,b in rows],
               page_effects=PAGE,below_mde=below,verdict=v,unchallenged=True),
          open(OUT/"mde_curve.json","w"),indent=1,ensure_ascii=False)

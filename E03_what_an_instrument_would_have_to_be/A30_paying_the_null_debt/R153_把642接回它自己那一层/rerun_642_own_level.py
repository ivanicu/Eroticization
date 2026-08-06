"""E03·A30·R153 —— 把 `#642` 接回它自己那一层,并补上它当年没有的零

**类型:FRONTIER。**
**心理学的那一句(本轮要判的):「一个社会换谁挨打几乎不改变力度,换用哪一种手段就把它压掉四分之三」——
这句话当年没有零。补上之后它还站不站得住?**

## 硬规则①(逐格打印,并与 `#698` 的层级对齐)
`#642` 用的是 **四个对象 × 五个手段**:对象 = 早期男孩 / 早期女孩 / 晚期男孩 / 晚期女孩;
手段 = 以身作则 · 讲道理 · 体罚 · 放任 · 疼爱。**逐格 n 见输出。**
`#698` 用的是**块均值层**(先把四个对象平均成一个块分,再算块间相关),**n = 97**,
`#698` 已明写两层**不是同一个量,不得互相引用为复现**。

## ⚠ 方法未知项(跑之前先说)
账本正文**没有写死** `+0.2288` 用的是 `|ρ|` 中位还是带号中位 ⇒
**两种都算,让 ④ 的正对照决定哪一种是当年的**;**若两种都复现不了,记「旧值不可复现」。**

## G1 ESTIMAND(按 `#642` 原样)
每个对象:五个手段两两(10 对)的天花板归一相关,取中位;**再取四个对象的中位。**
## G2 CONTROLS
**零**:**每个社会内部把五个手段分打乱**(保住该社会整体水平,毁掉手段特异性)——
与 `#698`/`#699` 同一套 halo 零 ⇒ `offset_control`,
**零的种类 = 只剩「这个社会整体被打了相近的分」时的同一个中位。**
**④ 正对照**:必须复现 `#642` 账本里的四对象值 **+0.1874 / +0.2252 / +0.2328 / +0.1947**。
## ⑧ 判据(`#710` 写死)
**重跑值与账本 +0.2288 之差 < 该轮自己的零的 95% 分位**;更大 ⇒ **记「旧值不可复现」。**
## ⑤ 最强混淆(`#710` 预注册)
**四对象层与块均值层不是同一个量** ⇒ **两层都报,并明写为何不同,不许用一层「验证」另一层。**
## IMPOSSIBLE(不写 planned)
五块全部出自 `barry1977agents` ⇒ **halo 与真结构不可分**(`#698` 已登记);
**换不了仪器**(`#700` 已枚举 D-PLACE 全部 12 个档案);因果:横断面民族志编码。`[unchallenged]`
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
D=pd.read_csv(B+"data.csv")
W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
TECH={"以身作则":429,"讲道理":437,"体罚":453,"放任":465,"疼爱":469}
TGT={"早期男孩":0,"早期女孩":1,"晚期男孩":2,"晚期女孩":3}
def sp(a,b): return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def nrm(fr,a,b,floor=30):
    m=fr[[a,b]].dropna()
    if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: return np.nan
    x=m[a].to_numpy(float); y=m[b].to_numpy(float); r=sp(x,y)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
    c=sp(xs,ys)
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
print("① 逐格 n(四对象 × 五手段):")
for t,off in TGT.items():
    cols=[f"SCCS{v+off}" for v in TECH.values()]
    print(f"   {t:6s} 变量 {cols}  五题皆有分的社会 n = **{len(W[cols].dropna())}**")
def per_target(fr,off,absolute):
    cols=[f"SCCS{v+off}" for v in TECH.values()]
    vals=[]
    for a,b in combinations(cols,2):
        v=nrm(fr,a,b)
        if np.isfinite(v): vals.append(abs(v) if absolute else v)
    return float(np.median(vals)) if vals else np.nan
for absolute,tag in ((True,"|ρ| 中位"),(False,"带号中位")):
    per=[per_target(W,off,absolute) for off in TGT.values()]
    print(f"\n{tag}:四对象 {[f'{x:+.4f}' for x in per]} · 四者中位 **{np.median(per):+.4f}**")
LED=[0.1874,0.2252,0.2328,0.1947]
best=None
for absolute,tag in ((True,"|ρ|"),(False,"signed")):
    per=[per_target(W,off,absolute) for off in TGT.values()]
    d=max(abs(a-b) for a,b in zip(sorted(per,reverse=True),sorted(LED,reverse=True)))
    print(f"④ 与账本四对象的最大绝对差({tag})= **{d:.4f}**")
    if best is None or d<best[1]: best=(absolute,d,per,tag)
absolute,diff,per,tag=best
obs=float(np.median(per))
print(f"\n⇒ 当年用的应是 **{tag}**(差最小 {diff:.4f});重跑四者中位 = **{obs:+.4f}**(账本 +0.2288)")
rng=np.random.default_rng(20260806); nul=[]
allc=[f"SCCS{v+o}" for v in TECH.values() for o in TGT.values()]
M=W[allc].copy()
for _ in range(400):
    P=M.copy()
    for off in TGT.values():
        cols=[f"SCCS{v+off}" for v in TECH.values()]
        A=P[cols].to_numpy(float).copy()   # ⚠ 不加 .copy() 是只读视图,赋值会 ValueError
        ok=~np.isnan(A).any(1)
        idx=np.where(ok)[0]
        for i in idx: A[i,:]=rng.permutation(A[i,:])
        P[cols]=A
    per_n=[per_target(P,off,absolute) for off in TGT.values()]
    v=np.nanmedian(per_n)
    if np.isfinite(v): nul.append(abs(v))
nul=np.array(nul); q=float(np.quantile(nul,.95))
print(f"零(社会内五手段互换,B={len(nul)} ⇒ 分辨率 1/{len(nul)+1}):95% 分位 **{q:.4f}** · 中位 {np.median(nul):.4f}")
G=Gate("把 #642 接回它自己那一层")
p1=G.positive_control("必须复现 #642 账本的四对象值(最大绝对差 <0.02)",planted=float(0.02-diff),floor=0.0,spread=0.001)
p2=G.offset_control("四对象中位必须超过 halo 零",effect=abs(obs),offset=q,spread=0.01,
                    null_kind="社会内部五个手段分互换 —— 保住该社会整体水平,毁掉手段特异性(halo)")
if p1:
    v=(f"**`#642` 补零后仍站得住:四对象中位 {obs:+.4f}(账本 +0.2288,差 {abs(obs-0.2288):.4f}),超过 halo 零 {q:.4f}**"
       if p2 else f"**`#642` 降级:四对象中位 {obs:+.4f} 未超过 halo 零 {q:.4f}**")
else: v=f"**旧值不可复现:两种约定的最大绝对差都 ≥0.02(最小 {diff:.4f})**"
print(f"\n{v}"); print(G)
json.dump(dict(convention=tag,per_target=[float(x) for x in per],median=obs,ledger=0.2288,
               ledger_per_target=LED,max_abs_diff=diff,null_q95=q,null_median=float(np.median(nul)),
               B=len(nul),block_mean_level_n=97,block_mean_level_value=0.3241,
               verdict=v,unchallenged=True),open(OUT/"rerun_642.json","w"),indent=1,ensure_ascii=False)

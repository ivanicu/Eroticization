"""E03·A30·R156 —— 给「换一批人来读,还是同一个数」补上它的零

**类型:FRONTIER。**
**心理学的那一句(本轮要判的):两拨互不相干的编码者,隔二十年读同一批民族志,
对「男孩晚期挨不挨打」给出的分,是不是真的对得上?**

## 硬规则①(已跑,含码本)
| 变量 | n | 编码项目 | 码 |
|---|---|---|---|
| `SCCS455` 体罚·晚期男孩 | 147 | **`barry1977agents`** | 1…11(11 档序数) |
| `SCCS1766` 男孩晚期体罚 | 79 | **`lang1998conan`** | **10 无 · 20 有但频率无信息 · 21 少 · 22 多** |

⚠ **码本证实了 `#639` 当年的判断:`20` 档字面就是「有体罚,但频率无信息」⇒ 结构上不可序,剔除是对的。**
序为 **10 < 21 < 22**;联合 **n = 69**,剔 `20` 后 **n = 60** —— **与账本逐字对上。**

## G1 ESTIMAND
`SCCS455` × `SCCS1766` 的**生 Spearman**。

⚠ **约定是被数据定下来的,不是我挑的。** 第一版按本页通行的天花板归一算,得 +0.6791 / +0.6973,
与账本差 0.049 / 0.058 —— **同方向、系统性,不是噪声。** 改用生 Spearman 后
**+0.6301 / +0.6396,与账本差 0.0000 / 0.0000。** ⇒ **`#639` 当年用的是生 Spearman。**
本轮因此以生 Spearman 为估计量,**天花板归一进 G4 当第二个规格**(天花板 ≈ 0.93:11 档对 3 档配不满)。
**G4 四格 = {剔 20 (n=60) · 20 并入中档 (n=69)} × {生 · 归一}。**
## G2 CONTROLS
**零**:在社会之间打乱 `SCCS1766`(保住它的边际,毁掉配对)⇒ `negative_control`,
**零的种类 = 同一批社会、同一份 Lang 1998 边际下,配对被打断后的同一个相关。**
⚠ **而 `null` 位收到的必须与这句话是同一个量**(`#713` 刚写进 docstring)——
**本轮传的是零分布的 95% 分位,不是中位,也不是安慰剂。**
**④ 正对照**:必须复现账本的 **+0.6301(n=60)** 与 **+0.6396(n=69)**。
## ⑤ 停止条件(`#713` 在跑之前写死,不许跑完再找理由)
**若零的 95% 分位 ≥ 效应的一半(≈0.315),记「判不了」并停。**
## ⑧ 判据
**重跑值与账本 +0.630 之差 < 零的 95% 分位** ⇒ 可复现;更大 ⇒ **记「旧值不可复现」。**
## IMPOSSIBLE(不写 planned)
⚠ **`#639` 自己已登记两条,本轮原样继承,不许在重跑之后悄悄弱化:**
**① 这是一对,不是一个分布** —— 推翻不了 `#528` 的跨队中位 0.105,只能说明那个中位算在**不同构念**的对上;
**② 两个团队读的是同一批民族志** ⇒ **独立的编码,不是独立的观察。真独立复制需要第二次田野。**
n = 60–69 ⇒ 只测得到大效应;**换不了仪器**(`#700` 已枚举)。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate, calibrated_tolerance
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
Dd=pd.read_csv(B+"data.csv")
W=Dd.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
def sp(a,b): return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def ceil(x,y):
    r=sp(x,y); xs=np.sort(np.asarray(x,float)); ys=np.sort(np.asarray(y,float)); ys=ys if r>0 else ys[::-1]
    return abs(sp(xs,ys))
J=W[["SCCS455","SCCS1766"]].dropna().copy()
drop=J[J.SCCS1766!=20].copy()                       # 规格 A:剔 20
merge=J.copy(); merge["SCCS1766"]=merge.SCCS1766.replace({10:0,20:1,21:1,22:2})  # 规格 B:20 并入中档
drop["ord"]=drop.SCCS1766.replace({10:0,21:1,22:2})
specs={"剔除 20 档 (n=%d)"%len(drop):(drop["SCCS455"].to_numpy(float),drop["ord"].to_numpy(float)),
       "20 并入中档 (n=%d)"%len(merge):(merge["SCCS455"].to_numpy(float),merge["SCCS1766"].to_numpy(float))}
obs={}; grid={}
for k,(x,y) in specs.items():
    obs[k]=sp(x,y); c=ceil(x,y); grid[k+" · 生"]=obs[k]; grid[k+" · 归一"]=obs[k]/c
    print(f"{k}: 生 ρ = **{obs[k]:+.4f}** · 天花板 {c:.4f} · 归一 **{obs[k]/c:+.4f}**")
print("\nG3/G4 全格(四格都印,含不支持结论的):")
for k,v in grid.items(): print(f"  {k:26s} {v:+.4f}")
LED={"剔除 20 档 (n=%d)"%len(drop):0.6301,"20 并入中档 (n=%d)"%len(merge):0.6396}
diff=max(abs(obs[k]-LED[k]) for k in obs)
print(f"\n④ 与账本的最大绝对差 = **{diff:.4f}**(账本 +0.6301 / +0.6396)")
rng=np.random.default_rng(20260806)
x,y=specs[list(specs)[0]]
nul=np.array([abs(sp(x,rng.permutation(y))) for _ in range(4000)])   # 与效应同一个量:生 Spearman
q=calibrated_tolerance(nul)          # #712 的函数:直接取零分布的 95% 分位
print(f"\n零(打乱 Lang 1998,B=4000):**95% 分位 {q:.4f}** · 中位 {np.median(nul):.4f} · "
      f"分辨率 1/4001 = {1/4001:.5f}")
main=obs[list(specs)[0]]
print(f"⑤ 停止条件:零 {q:.4f} vs 效应的一半 {abs(main)/2:.4f} ⇒ "
      f"{'⛔ 判不了,停' if q>=abs(main)/2 else '✅ 未触发,继续'}")
G=Gate("给 #639 补上它的零")
p1=G.positive_control("必须复现账本的 +0.6301 / +0.6396(最大绝对差 <0.02)",
                      planted=float(0.02-diff),floor=0.0,spread=0.001)
p2=G.negative_control("打乱 Lang 1998 后相关应回到零",null=q,effect=abs(main),null_spread=0.005,
                      null_kind="在社会之间打乱 SCCS1766 —— 保住 Lang 1998 的边际,毁掉配对;传入的是零分布的 95% 分位(#713 的类型对齐)")
if p1 and p2 and q<abs(main)/2:
    v=(f"**`#639` 补零后仍站得住:生 Spearman ρ = {main:+.4f},零的 95% 分位 {q:.4f} —— 效应是它的 {abs(main)/q:.2f} 倍**")
elif q>=abs(main)/2: v="**判不了:零 ≥ 效应的一半(⑤ 的停止条件)**"
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(specs={k:float(v) for k,v in obs.items()},grid={k:float(v) for k,v in grid.items()},convention="raw Spearman —— 由数据定下,见 docstring",ledger=LED,max_abs_diff=diff,
               null_q95=q,null_median=float(np.median(nul)),B=4000,
               n_drop=int(len(drop)),n_merge=int(len(merge)),
               inherited_limits=["一对不是一个分布","同一批民族志 ⇒ 独立编码不是独立观察"],
               verdict=v,unchallenged=True),open(OUT/"null_for_639.json","w"),indent=1,ensure_ascii=False)

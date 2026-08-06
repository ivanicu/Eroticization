"""E03·A34·R179 —— 天花板归一与衰减校正,修的是不是同一件事

**类型:FRONTIER。这是 `#735`② —— 而它的判据必须在复合之前写下。**

**心理学的那一句(它决定这一页上每一个归一数怎么读):
这一页所有的「紧密度」都除过一次天花板。若天花板本身已经吃掉了一部分测量误差,
那么这一页的每一个数都被抬高过,而我从没检验过这件事。**

## 缺口
`#735` 用**衰减校正**得到 0.30;页上的 **0.4154** 用的是**天花板归一**。
`#735`② 写下:**两者看起来正交(一个修边际不对称,一个修测量误差),而我没有验证过。**
⚠ **实数据回答不了这个问题:只有 4 个题、6 对,n 太小。** ⇒ **合成世界**(`realstat` 攻击阶梯第 4 级)。

## G1 ESTIMAND
**估计量的偏差**:`E[估计] − ρ_true`,在 (ρ_true × 信度 × 边际偏斜) 的网格上,
对四个估计量各算一次:**生 · 只归一 · 只衰减校正 · 两者复合**。
## W1 / W2 / W3
| 世界 | 复合后的偏差 | 读法 |
|---|---|---|
| **W1 正交** | ≈ 0(而单用任一个都偏低) | **两者可以复合**,页上的数可以再乘一次校正 |
| **W2 重复校正** | **显著 > 0** | **不许复合** —— 而且**只归一那一列若也 > 0,页上每一个归一数都被抬高过** |
| **W3 天花板本身就够** | 只归一 ≈ 0,复合 > 0 | 天花板已经吃掉了测量误差,**衰减校正在这具数据上是多余的** |

⚠ **W2 与 W3 的正结果都不利于这一页** —— 它们说页上的归一数偏高。**这是本轮设计成能出的结果。**

## G2 CONTROLS
**④ 正对照**:**边际不同 + 信度 = 1** 时,`只归一` 的偏差必须**小于** `生` —— 这是天花板归一**自称**在做的事。
⚠ **第一版的正对照写错了**(预设 Spearman 能从 4 档数据取回 ρ_true),而 **4 档离散化本身衰减 12–13%**;
⚠ **第一版的设计更错**:两个题用了同一组切点 ⇒ **天花板恒等于 1 ⇒ 归一从没被启动过。**
⚠ **且必须在 g=0 时失败**:**ρ_true = 0** 时四个估计量都必须 ≈ 0 ——
**一个在真值为零时仍返回正数的校正,是在无中生有。**
**零** = `negative_control`,**零的种类 = ρ_true 设为 0 的同一套合成流程 ——
保住样本量、档数、边际偏斜与噪声,只把真实关联抽掉。**
## G3:ρ_true ∈ {0, .2, .4, .6} × 信度 ∈ {0.7, 0.85, 1.0} × 偏斜 ∈ {对称, 中度, 强} = **36 格全报**。
## G4:档数 4(与 GSS 同)· 每格 200 次重复 · 报中位偏差与 5–95% 区间。
## ⑤ 停止条件(跑之前写死)
- **正对照不过(信度=1、对称时任一估计量偏差 ≥0.02)⇒ UNVERIFIED 并停。**
- **ρ_true=0 时任一估计量的中位 ≥0.05 ⇒ 该估计量记「无中生有」,单列并停用。**
- **复合的中位偏差 > +0.10(相对 ρ_true)⇒ 判 W2,页上写明不许复合。**
- **只归一那一列的中位偏差 > +0.10 ⇒ 页上每一个归一数都要带「在有测量误差时偏高」这个限定。**
## IMPOSSIBLE(不写 planned)
⚠ **本轮结构上一具仪器也没有用,而这不是遗漏:换不了仪器,因为它问的是估计量本身的性质,
不是任何一份数据的性质。** 把它跑在 GSS 上反而会把「估计量偏多少」与「GSS 的真值是多少」混在一起,
而后者没有已知答案 —— **这正是必须用合成世界的理由。**
**这是合成世界,不是 GSS** ⇒ 它能说的是**估计量的性质**,**不能**说 GSS 的真实值是多少;
真实数据的噪声结构未必是这里的加性高斯 ⇒ **结论的 scope 是「在这一类噪声下」。**
`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def ceil_of(a,b):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    r=sp(a,b); return abs(sp(x,y if r>0 else y[::-1]))
# ⚠ 第一版给**两个题用同一组切点** ⇒ 两列边际相同 ⇒ 排序后逐位单调 ⇒ **天花板恒等于 1**
#   ⇒ **天花板归一从没被启动过,而本轮的全部目的就是测它。**
#   修法:**每一格给两个题不同的切点** —— 这才是 GSS 的真实处境
#   (`homosex` 有 5 档且边际与 `premarsx` 差很远)。
CUT={"边际相同":([-0.67,0,0.67],[-0.67,0,0.67]),
     "边际中度不同":([-0.67,0,0.67],[-0.2,0.5,1.2]),
     "边际强烈不同":([-1.0,-0.3,0.3],[0.4,1.0,1.7])}
def gen(rho,rel,cutpair,n,rng):
    """潜变量 -> 加噪(信度 rel)-> 按 cuts 切成 4 档。"""
    L=rng.multivariate_normal([0,0],[[1,rho],[rho,1]],size=n)
    if rel<1:
        s=np.sqrt((1-rel)/rel)
        L=L+rng.normal(0,s,size=L.shape)
    ca,cb=cutpair
    return np.digitize(L[:,0],ca),np.digitize(L[:,1],cb)
def estimators(a,b,rel):
    r=sp(a,b); c=ceil_of(a,b)
    return dict(生=r, 只归一=(r/c if c>1e-9 else np.nan),
                只衰减=r/rel, 复合=((r/c)/rel if c>1e-9 else np.nan))
NAMES=["生","只归一","只衰减","复合"]
rng=np.random.default_rng(20260806); N=3000; REP=200
grid={}
for rho in (0.0,0.2,0.4,0.6):
    for rel in (0.7,0.85,1.0):
        for sk,cuts in CUT.items():
            acc={k:[] for k in NAMES}
            for _ in range(REP):
                a,b=gen(rho,rel,cuts,N,rng)
                for k,v in estimators(a,b,rel).items(): acc[k].append(v)
            grid[(rho,rel,sk)]={k:float(np.nanmedian(v)) for k,v in acc.items()}
print(f"合成世界:n={N} · 每格 {REP} 次 · 4 档 · {len(grid)} 格全报\n")
# ⚠ 第一版的正对照预设「信度=1、边际相同时 Spearman 应等于 ρ_true」—— **那是假的**:
#   **4 档离散化本身就衰减 12–13%**,而四个估计量没有一个是设计来修它的。
#   修法:正对照改成检验**天花板归一自己宣称的性质** ——
#   **当两个题的边际不同、且信度=1 时,`只归一` 的偏差必须小于 `生`。**
#   若不小于,**天花板归一就没有做它自称在做的事**,后面的复合问题不必谈。
print("=== ④ 正对照(改):边际不同 + 信度=1 时,`只归一` 的偏差必须小于 `生` ===")
pc=True
for rho in (0.2,0.4,0.6):
    for sk in ("边际中度不同","边际强烈不同"):
        d=grid[(rho,1.0,sk)]
        br,bn=abs(d["生"]-rho),abs(d["只归一"]-rho)
        ok=bn<br; pc&=ok
        print(f"  ρ_true={rho} {sk}: 生 {d['生']:+.3f}(|偏|{br:.3f}) · 只归一 {d['只归一']:+.3f}(|偏|{bn:.3f})  {'✅' if ok else '⛔'}")
print("  ⚠ 离散化本身的衰减(信度=1、边际相同格)= "
      + " · ".join(f"ρ={r}: {grid[(r,1.0,'边际相同')]['生']-r:+.3f}" for r in (0.2,0.4,0.6))
      + "  —— **四个估计量都不修它,如实登记。**")
print("\n=== 零(g=0):ρ_true=0 时四个估计量都必须 ≈ 0 ===")
zero_bad=[]
for rel in (0.7,0.85,1.0):
    for sk in CUT:
        d=grid[(0.0,rel,sk)]
        for k in NAMES:
            if abs(d[k])>=0.05: zero_bad.append((rel,sk,k,d[k]))
        if (rel,sk)in((0.7,"边际强烈不同"),(1.0,"边际相同")):
            print(f"  信度{rel} 偏斜{sk}: "+" · ".join(f"{k} {d[k]:+.4f}" for k in NAMES))
print(f"  ⇒ 在 ρ_true=0 的 9 格 × 4 个估计量里,中位 ≥0.05 的:**{len(zero_bad)}** "
      f"{'✅ 没有一个无中生有' if not zero_bad else zero_bad[:4]}")
print("\n=== G3 全格:相对偏差 (估计 − ρ_true)/ρ_true,ρ_true>0 的 27 格 ===")
print(f"{'ρ':>5s}{'信度':>6s}{'偏斜':>6s}" + "".join(f"{k:>10s}" for k in NAMES))
rel_bias={k:[] for k in NAMES}
for rho in (0.2,0.4,0.6):
    for rel in (0.7,0.85,1.0):
        for sk in CUT:
            d=grid[(rho,rel,sk)]
            row=[(d[k]-rho)/rho for k in NAMES]
            for k,v in zip(NAMES,row): rel_bias[k].append(v)
            print(f"{rho:>5.1f}{rel:>6.2f}{sk:>6s}" + "".join(f"{v:>+10.3f}" for v in row))
print("\n=== 中位相对偏差(27 格)===")
for k in NAMES: print(f"  {k:6s} **{np.median(rel_bias[k]):+.3f}**  (5–95% [{np.quantile(rel_bias[k],.05):+.3f}, {np.quantile(rel_bias[k],.95):+.3f}])")
mn=float(np.median(rel_bias["只归一"])); mc=float(np.median(rel_bias["复合"]))
G=Gate("天花板归一与衰减校正修的是不是同一件事")
p1=G.positive_control("信度=1 且对称时四个估计量都必须 ≈ ρ_true",planted=1.0 if pc else 0.0,floor=0.0,spread=0.1)
p2=G.negative_control("ρ_true=0 时估计量必须回到零",
    null=float(max(abs(grid[(0.0,r,s)][k]) for r in (0.7,0.85,1.0) for s in CUT for k in NAMES)),
    effect=0.4,null_spread=0.005,
    null_kind="ρ_true 设为 0 的同一套合成流程 —— 保住样本量、档数、边际偏斜与噪声,只把真实关联抽掉")
if not p1: v="**UNVERIFIED:正对照没过**"
elif mc>0.10 and mn>0.10: v=f"**W2+:复合偏 {mc:+.1%},而且**只归一本身就偏 {mn:+.1%}** ⇒ 不许复合,且页上每个归一数都要带限定**"
elif mc>0.10: v=f"**W2:复合偏 {mc:+.1%} ⇒ 不许复合**"
elif abs(mn)<=0.10 and abs(mc)<=0.10: v=f"**W1:两者正交(只归一 {mn:+.1%} · 复合 {mc:+.1%})⇒ 可以复合**"
else: v=f"**判不了:只归一 {mn:+.1%} · 复合 {mc:+.1%}**"
print(f"\n{v}"); print(G)
json.dump(dict(grid={f"{a}|{b}|{c}":grid[(a,b,c)] for a,b,c in grid},
  median_rel_bias={k:float(np.median(rel_bias[k])) for k in NAMES},
  zero_violations=len(zero_bad),verdict=v,unchallenged=True),open(OUT/"two.json","w"),indent=1,ensure_ascii=False)

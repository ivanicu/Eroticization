"""E03·A30·R161 —— 那些比值的分母,自己稳不稳

**类型:FRONTIER。**

**心理学的那一句(元层,但它决定这一页上二十多句话能不能说):
我报的每一个「效应是零的 N 倍」,分母都是估的。那个估计自己抖多少?**

## 为什么现在做(`#718`①,已被推后三轮)
量出来了:**页上 28 行带 ⟨比值⟩ 的,25 行从未声明分母是全枚举还是抽样估的。**
而 `#716` 实测过一次:同一个 MFQ k=6 的零,抽 1,500 给 **0.0162**、抽 3,000 给 **0.0313**、
**全枚举 0.0299** —— **分母抖了 1.9 倍**,而它是每一个比值的分母。

## W1 / W2(预测矩阵)—— 这两个世界给的下一步完全不同
| | B=4000 置换零的 q95 跨种子极差 | 读法 | 下一步 |
|---|---|---|---|
| **W1 普遍不稳** | ≳10% | 这一页的比值**广泛不可靠** | 全部重估,页上比值改成区间 |
| **W2 只是组合抽样特有** | ≲3% | 不稳来自**结构化块空间的组合抽样**,不是置换 | 债已还清,只需标注哪几格是组合抽的 |

**W1 是我不高兴的那个** —— 它会把这一页二十多行的比值一起降级。

## G1 ESTIMAND
同一个置换零的 **95% 分位,跨 5 个种子的极差 ÷ 中位**(相对抖动)。
## G2 CONTROLS
**④ 正对照**:同一具装置在 **B 很小(200)** 时必须**明显更抖** —— 若 B=200 也稳,
说明这个量根本测不到抖动,**记「仪器不灵」并停**(不是「零很稳」)。
**对照臂**:把 `#716` 的组合抽样按同样口径量一次(抽 1,500 vs 全枚举),**两个臂并排**。
## ⑤ 停止条件(跑之前写死)
- B=200 的相对抖动 **不大于** B=4000 的 ⇒ **记「仪器不灵」并停。**
- B=4000 的相对抖动 **≥10% ⇒ 判 W1**,页上比值全部要改区间。
- **≤3% ⇒ 判 W2。** 落在 3–10% ⇒ **记「判不了」,报区间。**
## IMPOSSIBLE(不写 planned)
只在 `#714` 这一具装置上量(SCCS n=60 的秩相关置换)—— **换不了仪器**:
别的行的原始数据与统计量各不相同,一轮里重估不完;本轮判的是**这一类零**,不是每一行。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
B0=pathlib.Path("data/external/dplace/repo/datasets/SCCS/")
W=pd.read_csv(B0/"data.csv").pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
J=W[["SCCS455","SCCS1766"]].dropna(); J=J[J.SCCS1766!=20].copy()
J["ord"]=J.SCCS1766.replace({10:0,21:1,22:2})
x=J["SCCS455"].to_numpy(float); y=J["ord"].to_numpy(float)
print(f"硬规则①:`#714` 的装置 —— SCCS455 × SCCS1766(剔 20 档),**n = {len(J)}**,生 Spearman")
print(f"  实测效应 {sp(x,y):+.4f}(账本 +0.6301)\n")
print("置换零的 95% 分位,跨 5 个种子:")
res={}
for Bn in (200,4000,20000):
    q=[]
    for s in (1,2,3,4,5):
        r=np.random.default_rng(s)
        q.append(float(np.quantile(np.abs([sp(x,r.permutation(y)) for _ in range(Bn)]),0.95)))
    rng_=max(q)-min(q); rel=rng_/np.median(q)
    res[Bn]=dict(q=q,spread=rng_,rel=rel)
    print(f"  B={Bn:>5,}  "+" ".join(f"{v:.4f}" for v in q)+
          f"   极差 {rng_:.4f} · **相对抖动 {rel*100:5.1f}%**")
rel4=res[4000]["rel"]; rel200=res[200]["rel"]
print(f"\n④ 正对照:B=200 的抖动 {rel200*100:.1f}% {'>' if rel200>rel4 else '≤'} B=4000 的 {rel4*100:.1f}%  "
      f"{'✅ 仪器测得到抖动' if rel200>rel4 else '⛔ ⑤ 触发:仪器不灵,停'}")
# 对照臂:#716 的组合抽样,同样口径
print("\n对照臂 —— `#716` 的组合抽样(MFQ k=6 从 593,770 块里抽),同样跨 5 个种子:")
import pyreadstat
ITEM={"emotionally":0,"weak":0,"cruel":0,"compassion":0,"animal":0,"kill":0,"treated":1,"unfairly":1,
 "rights":1,"fairly":1,"justice":1,"rich":1,"lovecountry":2,"betray":2,"loyalty":2,"history":2,"family":2,
 "team":2,"respect":3,"traditions":3,"chaos":3,"kidrespect":3,"sexroles":3,"soldier":3,
 "decency":4,"disgusting":4,"god":4,"harmlessdg":4,"unnatural":4,"chastity":4}
mfq,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
X=mfq[list(ITEM)].dropna(); it=list(ITEM); n=len(it); M=np.eye(n)
for a in range(n):
    for b in range(a+1,n):
        xa=X[it[a]].to_numpy(float); yb=X[it[b]].to_numpy(float)
        r=sp(xa,yb); xs=np.sort(xa); ys=np.sort(yb); ys=ys if r>0 else ys[::-1]
        M[a,b]=M[b,a]=r/abs(sp(xs,ys))
truth={frozenset(i for i,d in ITEM.items() if d==D) for D in set(ITEM.values())}
allb=[c for c in itertools.combinations(range(n),6) if frozenset(it[j] for j in c) not in truth]
Vfull=np.array([min(M[a,b] for a,b in itertools.combinations(c,2)) for c in allb])
exact=float(np.quantile(Vfull,0.95))
comb={}
for Bn in (1500,3000):
    q=[]
    for s in (1,2,3,4,5):
        r=np.random.default_rng(s); q.append(float(np.quantile(Vfull[r.choice(len(allb),Bn,replace=False)],0.95)))
    comb[Bn]=dict(q=q,rel=(max(q)-min(q))/abs(np.median(q)))
    print(f"  抽 {Bn:,}  "+" ".join(f"{v:.4f}" for v in q)+
          f"   极差 {max(q)-min(q):.4f} · **相对抖动 {comb[Bn]['rel']*100:5.1f}%**(全枚举精确值 {exact:.4f})")
G=Gate("那些比值的分母自己稳不稳")
p1=G.positive_control("B=200 必须比 B=4000 抖得明显(否则仪器测不到抖动)",
    planted=float(rel200-rel4),floor=0.0,spread=0.002)
p2=G.negative_control("B=4000 的置换零应当稳",null=float(rel4),effect=0.10,null_spread=0.005,
    null_kind="同一具置换装置、同一个 n=60 的配对,只换随机种子 —— 测的是分母自己的抽样误差")
if not p1: v="**判不了:仪器测不到抖动**"
elif rel4>=0.10: v=f"**W1:B=4000 的置换零相对抖动 {rel4*100:.1f}% ≥ 10% ⇒ 页上比值要改区间**"
elif rel4<=0.03: v=(f"**W2:B=4000 的置换零只抖 {rel4*100:.1f}%,而组合抽样抖 {comb[1500]['rel']*100:.1f}% "
                    f"⇒ 不稳来自结构化块空间的组合抽样,不是置换零**")
else: v=f"**判不了({rel4*100:.1f}% 落在 3–10%)—— 报区间**"
print(f"\n{v}"); print(G)
json.dump(dict(permutation={str(k):v for k,v in res.items()},combinatorial={str(k):v for k,v in comb.items()},
   exact=exact,verdict=v,unchallenged=True),open(OUT/"stability.json","w"),indent=1,ensure_ascii=False)

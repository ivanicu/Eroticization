"""E03·A32·R166 —— 十个品格里,性克制排第几

**类型:FRONTIER。这是 `#723`③ 要求先找、再跑的那个对象层缺口。**

**心理学的那一句:一个社会教孩子的十件事里,有九件是对四种孩子一视同仁地教的。
性克制是不是那个例外 —— 唯一一件「看你是男孩还是女孩、几岁」才决定的?**

## 缺口从哪来(不是我挑的,是页面自己写着的)
`#617` 确立:**一个社会不决定暴力对不对,它决定对谁做才不对。**
`#618` 记的是:**同一个问题问到「性」时测不到** —— 三个搜索版本没有一个通过自己的正对照,
诚实状态是「没测到,不是不存在」,而失败的原因是**「对谁」写在变量之间,不在任何一个之内**。
⇒ **而 `#717` 恰好证明了怎么用这种结构**(5 手段 × 4 对象)。**同一把工具,指向 `#618` 留下的洞。**

## 池(硬规则①先印)
`barry1976traits` 44 个变量里,**「品格 × 四对象」的完整族恰好 10 个**,而
**`Sexual Restraint`(SCCS330–333)是其中之一**,`Self-restraint`(SCCS326–329)是它的**近邻对照**
—— 同一件事减去性。**零池 = 这 40 个变量,k=4,C(40,4)=91,390,可全枚举。**

## W1 / W2
| | 性克制四对象的最弱一环 | 读法 |
|---|---|---|
| **W1 它也一视同仁** | 与其余九个品格同档 | 性只是第十件品格,**页上「性特殊」的说法要收窄** |
| **W2 它是那个例外** | **明显最低** | **一个社会对性的要求,是按「对谁」分档的** —— `#617` 的形状延伸到性 |

⚠ **两个结果我都不高兴:** W1 收窄「性特殊」;W2 则给 `#717`「接缝在手段上,不在孩子身上」开一个例外,
**而那句话就在页面上。**

## G1 ESTIMAND
每个品格四对象块的**最弱一环(天花板归一 · 最优符号)**,以及 **它在十个品格里的名次**。
## G2 CONTROLS
**零** = 同池同 k 全枚举(排除 10 个真块),**`negative_control`,零的种类 = 同一编码团队
barry1976traits 的同一个 40 变量池、同样 k=4、同样取最优符号,只打散「哪四个算一族」。**
**④ 正对照**:十个品格块必须**全部高于零的中位**(否则「品格 × 对象」这个结构本身不成立,后面不用看)。
**跨队臂(硬规则④)**:`SCCS596/597`(**Whyte 1978,另一个团队**)直接编「有没有双重标准」。
**若 W2 成立,双重标准 = Yes 的社会,其男孩/女孩性克制之差应更大。**
## ⑤ 停止条件(跑之前写死)
- 有任一品格块低于零的中位 ⇒ **结构不成立,记「判不了」并停。**
- **性克制若不是十个里最低的三个之一 ⇒ W2 不成立,如实报它的名次,不许改口径去找。**
- 跨队臂 **n < 30 ⇒ 记「判不了」,不进结论**(`#641` 的地板写在每一对上,不是联合 n)。
## IMPOSSIBLE(不写 planned)
**换不了仪器**:只有这一个团队编过「品格 × 四对象」的网格(`#700` 已枚举);
跨队臂只有 Whyte 一个二值/三值变量,**不能替代复制**。横断面,无因果。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, re
from lib.blocks import pairmat, opt_batch, weakest_optimal, weakest_greedy
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
V=pd.read_csv(B+"variables.csv",low_memory=False); Dd=pd.read_csv(B+"data.csv")
W=Dd.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
P=V[V.source.astype(str).str.contains("barry1976traits",na=False)]
fam={}
for _,r in P.iterrows():
    m=re.match(r'(.+?):\s*(Early|Late)\s+(Boy|Girl)s?$',str(r.title))
    if m: fam.setdefault(m.group(1),[]).append((r.id,f"{m.group(2)} {m.group(3)}"))
FAM={k:[x[0] for x in sorted(v,key=lambda y:y[1])] for k,v in fam.items() if len(v)==4}
POOL=[v for L in FAM.values() for v in L]
print(f"=== 硬规则①:`barry1976traits` 池 {len(P)} 变量 · 完整「品格 × 四对象」族 **{len(FAM)}** · 池 {len(POOL)} 变量 ===")
for k,L in FAM.items():
    print(f"  {k:18s} {L}  联合 n={len(W[L].dropna()):3d}  逐变量 n="+"/".join(str(int(W[c].notna().sum())) for c in L))
M=pairmat(W,POOL,floor=30); IDX={c:i for i,c in enumerate(POOL)}
obs={k:weakest_optimal(M,[IDX[c] for c in L]) for k,L in FAM.items()}
grd={k:weakest_greedy(M,[IDX[c] for c in L]) for k,L in FAM.items()}
truth={frozenset(L) for L in FAM.values()}
allb=[c for c in itertools.combinations(POOL,4) if frozenset(c) not in truth]
arr=np.array([[IDX[c] for c in b] for b in allb])
V0=opt_batch(M,arr); V0=V0[np.isfinite(V0)]
q50,q95,q99=(float(np.quantile(V0,x)) for x in (0.50,0.95,0.99))
print(f"\n=== 零:同池同 k **全枚举 {len(allb):,} 块**(排除 10 真块)===")
print(f"  中位 {q50:+.4f} · **95% 分位 {q95:+.4f}** · 99% 分位 {q99:+.4f}")
order=sorted(obs,key=lambda k:-obs[k])
print(f"\n=== G3 十个品格全报,按最弱一环排 ===")
print(f"{'名次':>3s} {'品格':20s}{'贪心':>9s}{'最优符号':>10s}{'÷零95%':>9s}")
for i,k in enumerate(order,1):
    star=" ★" if k=="Sexual Restraint" else ("  ←近邻对照" if k=="Self-restraint" else "")
    print(f"{i:>3d} {k:20s}{grd[k]:>+9.4f}{obs[k]:>+10.4f}{obs[k]/q95:>9.2f}{star}")
rank=order.index("Sexual Restraint")+1
below=[k for k in obs if obs[k]<=q50]
G=Gate("十个品格里性克制排第几")
p1=G.positive_control("十个品格块必须全部高于零的中位(否则结构不成立)",
    planted=float(min(obs.values())-q50),floor=0.0,spread=0.005)
p2=G.negative_control("同池随机四变量应低于真块",null=q95,effect=float(np.median(list(obs.values()))),
    null_spread=0.005,null_kind="同一编码团队 barry1976traits 的同一个 40 变量池、同样 k=4、同样最优符号,只打散「哪四个算一族」")
# 跨队臂
DS={}
for v,nm in (("SCCS596","婚前性·无双重标准"),("SCCS597","婚外性·无双重标准")):
    if v not in W.columns: continue
    J=W[[v]+FAM["Sexual Restraint"]].dropna()
    if len(J)<30: DS[nm]=dict(n=int(len(J)),verdict="判不了(n<30)"); continue
    boy=J[[FAM["Sexual Restraint"][0],FAM["Sexual Restraint"][2]]].mean(axis=1)
    girl=J[[FAM["Sexual Restraint"][1],FAM["Sexual Restraint"][3]]].mean(axis=1)
    gap=(girl-boy).abs()
    yes=gap[J[v]==1]; no=gap[J[v]!=1]
    DS[nm]=dict(n=int(len(J)),n_yes=int(len(yes)),n_no=int(len(no)),
                gap_yes=float(yes.median()) if len(yes) else None,gap_no=float(no.median()) if len(no) else None)
print(f"\n=== 跨队臂:Whyte 1978「有没有双重标准」× Barry 1976 的男女差 ===")
for nm,d in DS.items():
    if "verdict" in d: print(f"  {nm}: n={d['n']} ⇒ **{d['verdict']}**"); continue
    print(f"  {nm}: 联合 n={d['n']}(有双标 {d['n_yes']} · 无双标 {d['n_no']})· "
          f"|女−男| 中位 有双标 **{d['gap_yes']:.3f}** vs 无双标 **{d['gap_no']:.3f}**"
          f"{'  ⚠ 某一臂 <30,只作观察' if min(d['n_yes'],d['n_no'])<30 else ''}")
if not p1: v="**判不了:有品格块低于零的中位,结构不成立**"
elif rank>=8: v=f"**W2:性克制在十个品格里排第 {rank}(倒数第 {11-rank})—— 它是那个例外**"
else: v=f"**W1 不成立、W2 也不成立:性克制排第 {rank},不在最低的三个里 —— 如实报名次,不改口径**"
print(f"\n{v}"); print(G)
json.dump(dict(pool=POOL,weakest_optimal=obs,weakest_greedy=grd,null_median=q50,null_q95=q95,null_q99=q99,
   n_blocks=int(V0.size),rank_sexual=rank,order=order,double_standard=DS,verdict=v,unchallenged=True),
   open(OUT/"rank.json","w"),indent=1,ensure_ascii=False)

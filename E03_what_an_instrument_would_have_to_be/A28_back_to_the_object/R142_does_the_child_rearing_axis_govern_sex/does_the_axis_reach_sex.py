"""E03·A28·R142 —— 管孩子的那条轴,管不管性

**类型:FRONTIER。这是把 E02 的对象与社会侧仅存的结构连上的那一步。**
**心理学的那一句(本轮要判的):一个社会管孩子用「施加」还是「温暖」,
和它拿一件性的事怎么办,是不是同一条轴?**

## 硬规则①(已跑,从对象读出)
⚠ **我的关键词检索是坏仪器**:「sexual」命中 248 个变量,而绝大多数是 **sexual division of labor**
(性别分工),不是性规范。**按来源这条可判定的规则收窄**:
`frayser1985varieties`(45)与 `broude1976cross`(13)才是性规范编码,**两者都不是 `barry1977agents`**。

| 变量 | 码 | 联合 n | 处置 |
|---|---|---|---|
| **`SCCS165`** 婚前性态度·女 | **1=Expected … 6=Strongly disapproved**,干净序数 | **75** | ✅ **唯一可用** |
| `SCCS169` 婚外性 | 「单一标准/双重标准·仅丈夫/…」——**类型不是强度** | 67 | ⛔ **不可序,淘汰**(与 `#641` 丢 `SCCS620` 同理) |
| `SCCS961` 婚前性限制 | 码 4「仅男性可」夹在宽严之间 | **32** | ⛔ 不可序 + n 太小,淘汰 |
| `SCCS177` 同性性行为频率 | 二值 absent/present | 45 | ⛔ 是**频率**不是**规范**,本轮不用 |

## G1 ESTIMAND(方向在跑之前写死)
`施加轴 = z(体罚) + z(讲道理) − z(放任) − z(疼爱)`(语义载荷取自 `#699`,不重新拟合)。
**主量 = 天花板归一 `ρ(施加轴, SCCS165)`,n = 75。**
**W1 预测:正**(管教越施加 ⇒ 对婚前性越谴责);**W2:两条独立的轴 ⇒ 零附近。**
## G2 CONTROLS
**松零**:在社会之间打乱 `SCCS165` ⇒ `negative_control`。
**紧零 = halo**:社会内部把五个块分互换后重算轴(`SCCS165` 不动)——
  问的是「**具体的手段图样**是否重要,还是只要这个社会整体水平」⇒ `offset_control`,
  **零的种类 = 只剩社会整体水平时的同一个相关。**
**⑤ 跨项目基线(预注册)**:`SCCS165` 出自 `broude1976cross`,施加轴出自 `barry1977agents`
⇒ **跨项目对本来就更低**(`#528`/`#665` 量到的跨团队中位约 **0.105**)。
**必须把这条基线一起报,不许把跨项目的低值读成「没关系」。**
## G3/G4:四个成分块各自 × `SCCS165` 也全部照登。
## 仪器(硬规则②/④)—— 换不了仪器,而这是枚举出来的,不是断言
逐个查 D-PLACE 全部 12 个档案的 `variables.csv`(管教类 / 性规范类 命中数):
`SCCS` **57 / 30 ← 只有它两者都有** · `EA` 0 / 1 · `WNAI` 0 / 0 · `Binford` 0 / 0 ·
`Jenkins`·`Kreft`·`MODIS`·`TEOW`·`GMTED2010`·`GSHHS`·`ecoClimate` 全 0 / 0。
**⇒ 要同时有「管教手段」与「性规范」,这个星球上公开的跨文化档案里只此一具。没有第二具仪器。**
(而 `#672` 那次闸也这样拦过我一次,**当时第二具仪器是存在的** —— 所以这条每次都要真去查。)

## IMPOSSIBLE(不写 planned)
**n = 75 ⇒ 只能测到中等以上的效应**;`SCCS169`/`SCCS961` 码不可序 ⇒ **本站点问不了婚外性与婚前限制**;
**跨项目 ⇒ halo 在这里不再抬高,但衰减无法量化**;因果:横断面民族志编码。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
Dd=pd.read_csv(B+"data.csv")
W=Dd.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
BL={"以身作则":[429,430,431,432],"讲道理":[437,438,439,440],"体罚":[453,454,455,456],
    "放任":[465,466,467,468],"疼爱":[469,470,471,472]}
LOAD={"体罚":+1,"讲道理":+1,"以身作则":0,"放任":-1,"疼爱":-1}
S=pd.DataFrame({k:W[[f"SCCS{i}" for i in v]].mean(axis=1) for k,v in BL.items()}).dropna()
Z=(S-S.mean())/S.std()
axis=sum(LOAD[k]*Z[k] for k in BL if LOAD[k]!=0)
norm=W["SCCS165"]
J=pd.concat([axis.rename("axis"),norm.rename("norm")],axis=1).dropna()
print(f"联合 n = **{len(J)}**")
def sp(a,b): return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def nrm(x,y):
    r=sp(x,y)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(np.asarray(x,float)); ys=np.sort(np.asarray(y,float)); ys=ys if r>0 else ys[::-1]
    c=sp(xs,ys)
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
obs=nrm(J.axis.to_numpy(),J.norm.to_numpy())
print(f"主量 归一 ρ(施加轴, 婚前性谴责) = **{obs:+.4f}**(W1 预测为正)")
rng=np.random.default_rng(20260806)
loose=np.array([abs(nrm(J.axis.to_numpy(),rng.permutation(J.norm.to_numpy()))) for _ in range(3000)])
qL=float(np.quantile(loose,.95))
V=S.to_numpy(float); cols=list(S.columns); tight=[]
for _ in range(1500):
    A=V.copy()
    for i in range(A.shape[0]): A[i,:]=rng.permutation(A[i,:])
    P=pd.DataFrame(A,index=S.index,columns=cols); Zp=(P-P.mean())/P.std()
    ax=sum(LOAD[k]*Zp[k] for k in BL if LOAD[k]!=0)
    j=pd.concat([ax.rename("a"),norm.rename("n")],axis=1).dropna()
    v=nrm(j.a.to_numpy(),j.n.to_numpy())
    if np.isfinite(v): tight.append(abs(v))
tight=np.array(tight); qT=float(np.quantile(tight,.95))
print(f"松零 95% 分位 **{qL:.4f}**(中位 {np.median(loose):.4f})")
print(f"紧零 halo 95% 分位 **{qT:.4f}**(中位 {np.median(tight):.4f})")
print(f"⑤ 跨项目基线(`#528`/`#665` 跨团队中位)≈ **0.105** —— 本轮是跨项目对,**低值不得读成「没关系」**")
print(f"\nG4 四个成分块各自 × 婚前性谴责(全部照登):")
grid={}
for k in ("体罚","讲道理","放任","疼爱","以身作则"):
    jj=pd.concat([Z[k].rename("a"),norm.rename("n")],axis=1).dropna()
    v=nrm(jj.a.to_numpy(),jj.n.to_numpy()); grid[k]=v
    print(f"   {k:6s} **{v:+.4f}**  n={len(jj)}")
G=Gate("管孩子的轴管不管性")
p1=G.positive_control("施加轴内部必须仍成立(沿用 #699:六对预测 6/6、统计量 +0.3243 远超两零)",
                      planted=0.3243,floor=0.1231,spread=0.01)
p2=G.negative_control("打乱婚前性谴责后,相关应回到零",null=float(np.median(loose)),effect=abs(obs),
                      null_spread=0.01,null_kind="在社会之间打乱 SCCS165 —— 保边际,毁配对")
p3=G.offset_control("必须超过「只剩社会整体水平」时的同一个相关",effect=abs(obs),offset=qT,spread=0.01,
                    null_kind="社会内部五块互换后重算轴 —— 保住该社会整体水平,毁掉手段图样(halo)")
if p1 and p2:
    v=(f"**是同一条轴:归一 ρ = {obs:+.4f},方向与预注册一致,且同时超过松零 {qL:.4f} 与 halo 紧零 {qT:.4f}**"
       if (p3 and obs>0) else
       f"**不是同一条轴(或本设计判不出):ρ = {obs:+.4f},松零 {qL:.4f} / halo 紧零 {qT:.4f} "
       f"—— ⚠ 而这是跨项目对,基线本就约 0.105,**低值不得读成「没关系」**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(J)),obs=obs,loose_q95=qL,loose_median=float(np.median(loose)),
               tight_q95=qT,tight_median=float(np.median(tight)),cross_project_baseline=0.105,
               grid={k:(None if not np.isfinite(x) else float(x)) for k,x in grid.items()},
               verdict=v,unchallenged=True),open(OUT/"axis_reaches_sex.json","w"),indent=1,ensure_ascii=False)

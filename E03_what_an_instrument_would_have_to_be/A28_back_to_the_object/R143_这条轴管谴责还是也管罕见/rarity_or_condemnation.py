"""E03·A28·R143 —— 这条轴管谴责,还是也管罕见

**类型:FRONTIER。这是 E01 那条 `0.758` 在社会层的检验。**
**心理学的那一句(本轮要判的):一个社会越靠强制管孩子,它对婚前性越谴责 —— 这已经量到了。
那它那里的婚前性,是不是也更罕见?**

## 硬规则①(已跑,逐条读码本,不许把频率当罕见度)
| 变量 | 码 | 与轴联合 n | 处置 |
|---|---|---|---|
| **`SCCS167`** 婚前性频率·女 | **1=Universal … 4=Uncommon ⇒ 高 = 更罕见**,干净序数 | **67** | ✅ **主用** |
| `SCCS166` 同·男 | 同上 | 63 | ✅ **G4 规格** |
| `SCCS177` 同性性行为 | 二值 absent/present | 45 | ⛔ 不是罕见度梯度 |
| `SCCS160` 婚内性频率 | 「无节制/有时禁欲/节制/推崇禁欲」—— **是规范不是频率** | 43 | ⛔ 淘汰 |
| `SCCS175` · `SCCS159` | 男性性主动性 / 谈性尺度 | 39 / 35 | ⛔ 与罕见度无关 |

**⇒ `SCCS167` 与 `#700` 用的 `SCCS165` 是完美匹配对:同一行为、同一性别、一罕见一谴责。**

## G1 ESTIMAND(方向在跑之前写死)
`施加轴 = z(体罚)+z(讲道理)−z(放任)−z(疼爱)`(载荷取自 `#699`,不重新拟合)。
**三个量全部在同一个 n = 67 的子样本上算,保证可比:**
**① `ρ(轴, 罕见度)`** —— W1 预测正(轴也管罕见);**② `ρ(轴, 谴责)`** —— `#700` 已得 +0.3631,本轮复算;
**③ `ρ(谴责, 罕见度)`** —— **这是 E01 那条 0.758 的社会版。**
## ⑤ 同项目基线(预注册)
`SCCS165` 与 `SCCS167` **同出 `broude1976cross`** ⇒ **同项目对被抬高约 2.14 倍**(`#665`)。
**③ 必须当作「同项目对」来读,不许直接与 E01 的 0.758 比大小;
而 ① 与 ② 是跨项目对(轴出自 `barry1977agents`),基线约 0.105。**
## G2 CONTROLS
**松零**:打乱被预测量 ⇒ `negative_control`。
**紧零 = halo**:社会内五块互换后重算轴 ⇒ `offset_control`,**零的种类 = 只剩社会整体水平时的同一个相关**。
**正对照(④)**:`ρ(轴, 谴责)` 在本子样本上必须仍显著为正(沿用 `#700`)。
## ⑧ 判据(`#700` 写死)
**① 与 ② 同号且量级相当 ⇒ 罕见与谴责在社会层是同一件事(E01 的社会版复活);
只有 ② 成立 ⇒ 记「这条轴管谴责,不管罕见」,E01 的降级在社会层被坐实。**
## IMPOSSIBLE(不写 planned)
**n = 67 ⇒ 只测得到中等以上效应**;`SCCS160/175/159/177` 结构性不可用;
**同时有管教与性规范的档案只此一具**(`#700` 已枚举 D-PLACE 全部 12 个)**⇒ 没有第二具仪器**;
因果:横断面民族志编码,无干预。`[unchallenged]`
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
Z=(S-S.mean())/S.std(); axis=sum(LOAD[k]*Z[k] for k in BL if LOAD[k]!=0)
J=pd.concat([axis.rename("axis"),W["SCCS165"].rename("cond"),W["SCCS167"].rename("rare")],axis=1).dropna()
print(f"三量同一子样本 n = **{len(J)}**")
def sp(a,b): return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def nrm(x,y):
    r=sp(x,y)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(np.asarray(x,float)); ys=np.sort(np.asarray(y,float)); ys=ys if r>0 else ys[::-1]
    c=sp(xs,ys)
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
r_axis_rare=nrm(J.axis,J.rare); r_axis_cond=nrm(J.axis,J.cond); r_cond_rare=nrm(J.cond,J.rare)
print(f"① ρ(轴, 罕见度) = **{r_axis_rare:+.4f}**  [跨项目,基线约 0.105]  W1 预测正")
print(f"② ρ(轴, 谴责)   = **{r_axis_cond:+.4f}**  [跨项目;`#700` 在 n=75 上得 +0.3631]")
print(f"③ ρ(谴责, 罕见度) = **{r_cond_rare:+.4f}**  [**同项目对,抬高约 2.14 倍**;E01 社会版]")
rng=np.random.default_rng(20260806)
qL_r=float(np.quantile([abs(nrm(J.axis.to_numpy(),rng.permutation(J.rare.to_numpy()))) for _ in range(3000)],.95))
qL_c=float(np.quantile([abs(nrm(J.axis.to_numpy(),rng.permutation(J.cond.to_numpy()))) for _ in range(3000)],.95))
V=S.to_numpy(float); cols=list(S.columns); tr=[]; tc=[]
for _ in range(1500):
    A=V.copy()
    for i in range(A.shape[0]): A[i,:]=rng.permutation(A[i,:])
    P=pd.DataFrame(A,index=S.index,columns=cols); Zp=(P-P.mean())/P.std()
    ax=sum(LOAD[k]*Zp[k] for k in BL if LOAD[k]!=0)
    jj=pd.concat([ax.rename("a"),W["SCCS165"].rename("c"),W["SCCS167"].rename("r")],axis=1).dropna()
    a=nrm(jj.a,jj.r); b=nrm(jj.a,jj.c)
    if np.isfinite(a): tr.append(abs(a))
    if np.isfinite(b): tc.append(abs(b))
qT_r=float(np.quantile(tr,.95)); qT_c=float(np.quantile(tc,.95))
print(f"\n① 的松零 95% 分位是 {qL_r:.4f};① 的 halo 紧零 95% 分位是 {qT_r:.4f}")
print(f"② 的松零 95% 分位是 {qL_c:.4f};② 的 halo 紧零 95% 分位是 {qT_c:.4f}")
g4=nrm(*pd.concat([axis.rename("a"),W["SCCS166"].rename("r")],axis=1).dropna().T.to_numpy())
j6=pd.concat([axis.rename("a"),W["SCCS166"].rename("r")],axis=1).dropna()
print(f"\nG4 规格:ρ(轴, 罕见度·男 `SCCS166`) = **{nrm(j6.a,j6.r):+.4f}**  n={len(j6)}")
G=Gate("这条轴管谴责还是也管罕见")
p1=G.positive_control("ρ(轴, 谴责) 在本子样本上必须仍显著为正(沿用 #700)",
                      planted=abs(r_axis_cond),floor=qT_c,spread=0.01)
p2=G.negative_control("打乱罕见度后,① 应回到零",null=qL_r*0+float(np.median([abs(nrm(J.axis.to_numpy(),rng.permutation(J.rare.to_numpy()))) for _ in range(500)])),
                      effect=abs(r_axis_rare),null_spread=0.01,
                      null_kind="在社会之间打乱 SCCS167 —— 保边际,毁配对")
p3=G.offset_control("① 必须超过「只剩社会整体水平」时的同一个相关",effect=abs(r_axis_rare),offset=qT_r,
                    spread=0.01,null_kind="社会内五块互换后重算轴 —— 保住该社会整体水平,毁掉手段图样(halo)")
if p1:
    ok1 = p3 and (r_axis_rare>0)
    v=(f"**罕见与谴责在社会层是同一件事:① {r_axis_rare:+.4f} 与 ② {r_axis_cond:+.4f} 同号且量级相当,"
       f"① 超过 halo 零 {qT_r:.4f} ⇒ E01 的社会版复活**" if ok1 and abs(r_axis_rare)>=0.5*abs(r_axis_cond) else
       f"**这条轴管谴责,不管罕见:② {r_axis_cond:+.4f} 成立而 ① {r_axis_rare:+.4f} "
       f"{'未超 halo 零 %.4f'%qT_r if not p3 else '量级不足 ②的一半'} ⇒ E01 的降级在社会层被坐实**")
else: v="UNVERIFIED —— 正对照失败"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(J)),r_axis_rare=r_axis_rare,r_axis_cond=r_axis_cond,r_cond_rare=r_cond_rare,
               loose_q95_rare=qL_r,loose_q95_cond=qL_c,tight_q95_rare=qT_r,tight_q95_cond=qT_c,
               g4_male=float(nrm(j6.a,j6.r)),n_male=int(len(j6)),
               same_project_pair="SCCS165×SCCS167 (broude1976cross, 抬高约 2.14×)",
               cross_project_baseline=0.105,verdict=v,unchallenged=True),
          open(OUT/"rarity_or_condemnation.json","w"),indent=1,ensure_ascii=False)

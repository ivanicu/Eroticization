"""E03·A28·R140 —— 跨手段那 +0.2288,撑不撑得住同一个 halo 零

**类型:FRONTIER。** `#697` 换的方向:halo 在手段块**内部**赢了,
**那么块与块之间的 +0.2288,是这份数据里唯一可能不被 halo 解释的社会侧结构。**
**心理学的那一句(本轮要判的):「一个社会用哪一种手段管孩子」,是不是也只是
「这个编码者给这个社会打了相近的分」?**

## 硬规则①(已跑,从对象读出)
`barry1977agents` 共 **144** 个变量,其中管教手段是 **10 个块 × 4 个对象**
(早期男孩 / 早期女孩 / 晚期男孩 / 晚期女孩)。`#642` 用的五个块与真实非缺失 n:

| 块 | 变量 | n |
|---|---|---|
| 以身作则 Use of Example | `SCCS429–432` | 152–155 |
| 讲道理 Lecturing | `SCCS437–440` | 134–141 |
| 体罚 Corporal Punishment | `SCCS453–456` | 140–147 |
| 放任 Permissiveness | `SCCS465–468` | 167–169 |
| 疼爱 Affection | `SCCS469–472` | 154–155 |

**⚠ ⑤ 预注册的混淆是活的:五个块全部出自 `barry1977agents` ⇒ 块间与块内共享同一批编码者,
halo 可能跨块。这条写进「做不到的」,不当作已控制。**

## G1 ESTIMAND
每个社会每个块的分 = 该块四个对象的均值。
**主量 = 十对块间「天花板归一相关」的中位**,**|ρ| 与带号两种都报**(G4;`#642` 的符号是混的)。
## G2 CONTROLS —— 两个零
**松零**:五个块分各自独立在社会之间打乱 ⇒ `negative_control`(该是零)。
**紧零 = halo(本轮的判据所在)**:每个社会内部把它的五个块分在五个块之间打乱,
保住该社会整体水平,毁掉手段特异性 ⇒ `offset_control`,
**零的种类 = 只剩「这个社会整体被打了相近的分」时的同一个中位。**
## ⑧ 判据(`#697` 在跑之前写死)
**跨手段中位必须超过 halo 零的 95% 分位;不超 ⇒ 社会侧这一整条线都是 halo,如实记。**
## IMPOSSIBLE(不写 planned)
**五块同一编码项目 ⇒ halo 与真结构在本站点不可分**;每对 n ≤ 169;
**这些编码只有 D-PLACE 一处 ⇒ 换不了仪器**;因果:横断面民族志编码,无干预。`[unchallenged]`
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
BLOCKS={"以身作则":[429,430,431,432],"讲道理":[437,438,439,440],"体罚":[453,454,455,456],
        "放任":[465,466,467,468],"疼爱":[469,470,471,472]}
S=pd.DataFrame({k:W[[f"SCCS{i}" for i in v]].mean(axis=1) for k,v in BLOCKS.items()})
S=S.dropna()
print(f"五块皆有分的社会 **{len(S)}** 个")
def sp(a,b): return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def med(fr,absolute=True,floor=30):
    v=[]
    for a,b in combinations(fr.columns,2):
        m=fr[[a,b]].dropna()
        if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: continue
        x=m[a].to_numpy(float); y=m[b].to_numpy(float); r=sp(x,y)
        if not np.isfinite(r) or abs(r)<1e-12: continue
        xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        c=sp(xs,ys)
        if np.isfinite(c) and abs(c)>1e-9:
            u=r/abs(c); v.append(abs(u) if absolute else u)
    return (float(np.median(v)) if v else np.nan), len(v)
obs_abs,np_=med(S,True); obs_sgn,_=med(S,False)
print(f"观测:十对块间归一相关 **|ρ| 中位 {obs_abs:+.4f}** · **带号中位 {obs_sgn:+.4f}**(可用对 {np_})")
print(f"   (账本 `#642` 记 +0.2288 —— 那是四对象各自算完再取中位,层级不同,不直接可比)")
rng=np.random.default_rng(20260806); loose=[]; tight=[]
V=S.to_numpy(float)
for _ in range(2000):
    P=pd.DataFrame({c:rng.permutation(S[c].to_numpy()) for c in S.columns},index=S.index)
    m,_=med(P,True)
    if np.isfinite(m): loose.append(m)
for _ in range(2000):
    A=V.copy()
    for i in range(A.shape[0]): A[i,:]=rng.permutation(A[i,:])
    P=pd.DataFrame(A,index=S.index,columns=S.columns)
    m,_=med(P,True)
    if np.isfinite(m): tight.append(m)
loose=np.array(loose); tight=np.array(tight)
qL=float(np.quantile(loose,.95)); qT=float(np.quantile(tight,.95))
print(f"\n松零(五块各自独立打乱)95% 分位 **{qL:.4f}**(中位 {np.median(loose):.4f})")
print(f"紧零 halo(社会内五块互换)95% 分位 **{qT:.4f}**(中位 {np.median(tight):.4f})")
G=Gate("跨手段那一条撑不撑得住 halo")
p1=G.positive_control("跨手段 |ρ| 中位必须远超松零(仪器能开火)",planted=obs_abs,floor=qL,spread=0.01)
p2=G.negative_control("五块各自独立打乱后,中位应回到零",null=float(np.median(loose)),effect=obs_abs,
                      null_spread=0.01,null_kind="五个块分各自独立在社会之间打乱 —— 保每块边际,毁全部配对")
p3=G.offset_control("跨手段中位必须超过「只剩社会整体水平」时的同一个中位",
                    effect=obs_abs,offset=qT,spread=0.01,
                    null_kind="社会内部把五个块分在五个块之间打乱 —— 保住该社会整体水平,毁掉手段特异性(halo)")
if p1 and p2:
    v=(f"**跨手段结构撑住了 halo:|ρ| 中位 {obs_abs:+.4f} 超过 halo 零的 95% 分位 {qT:.4f} "
       f"⇒ 这是这份数据里唯一不被 halo 解释的社会侧结构**" if p3 else
       f"**也是 halo:|ρ| 中位 {obs_abs:+.4f} 不超过 halo 零的 95% 分位 {qT:.4f} "
       f"⇒ 社会侧这一整条线都是 halo,如实记**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n_soc=int(len(S)),obs_abs=obs_abs,obs_signed=obs_sgn,n_pairs=np_,
               loose_q95=qL,loose_median=float(np.median(loose)),
               tight_q95=qT,tight_median=float(np.median(tight)),
               blocks={k:[f"SCCS{i}" for i in v] for k,v in BLOCKS.items()},
               verdict=v,unchallenged=True),open(OUT/"across_technique.json","w"),indent=1,ensure_ascii=False)

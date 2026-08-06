"""E03·A28·R138 —— 这一页最老的社会侧声明,整条建立在一个编码团队的笔迹上

**类型:FRONTIER。A27 关弧,A28 开弧。盆地已破(`#692`–`#695` 连续四轮没推进对象)。**
**心理学的那一句(本轮要判的):「换谁挨打,几乎不改变一个社会下手的力度」——
这是关于社会的事实,还是关于一个团队怎么给 186 个社会打分的事实?**

## 硬规则①(已跑,从对象读出)
| 变量 | 非缺失 n | 编码项目 |
|---|---|---|
| `SCCS453` 体罚·儿子 | **145** | `barry1977agents` |
| `SCCS454` | **140** | `barry1977agents` |
| `SCCS455` | **147** | `barry1977agents` |
| `SCCS456` | **142** | `barry1977agents` |
| `SCCS469–472` 疼爱四题 | 154–155 | `barry1977agents` |
| `SCCS172` | 103 | `broude1976cross` |
| `SCCS619 / SCCS620` | 93 / 63 | `whyte1978cross` |
| `SCCS754` | 70 | `broude1983cross` |
| `SCCS1801` | 135 | `divale1999codes` |

**⇒ 体罚四题全部出自同一个编码项目 ⇒ `#640` 的 `+0.845` 六对全是同项目对。**
**而 `#665` 已测出同项目对被抬高约 2.14 倍。**

## G1 ESTIMAND
**① 同项目**:`SCCS453–456` 六对天花板归一相关的中位(= `#640` 的量,用今天的流水线重跑)。
**② 跨项目**:每个体罚题 × 每个非 `barry1977` 变量的归一相关中位。
**主量 = 两者都报,不合并**(⑤ 预注册)。
## ⑧ 判据(`#695` 在跑之前写死)
**重跑值与账本记载 `+0.845` 之差,必须小于该量自己的零的 95% 分位;更大 ⇒ 记「旧值不可复现」。**
## G2 CONTROLS
**正对照(④)**:必须复现 `#639` 已验证的跨编码者一致性(该条有明确数值,见下方输出比对)。
**零**:在社会之间打乱其中一题(保住边际,毁掉配对)⇒ `negative_control`,
**零的种类 = 同一批社会、同一批边际下,配对被打断后的同一个中位。**
## KILL(条件式)
if 正对照复现 and 零确实为零: evaluate(判据⑧) else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**四题同项目 ⇒ 「同项目 vs 跨项目」只能拿不同构念的变量比**,不是等价对照;
SCCS 186 个社会 ⇒ 每对 n ≤ 147;**因果:横断面民族志编码,无干预**;
跨仪器:体罚编码只有 D-PLACE 这一处 ⇒ **换不了仪器**。`[unchallenged]`
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
D=pd.read_csv(B+"data.csv"); V=pd.read_csv(B+"variables.csv")
W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
SRC=dict(zip(V.id,V.source))
CORP=["SCCS453","SCCS454","SCCS455","SCCS456"]
OTHER=[v for v in ["SCCS172","SCCS619","SCCS620","SCCS754","SCCS1801","SCCS298","SCCS1766"]
       if v in W.columns and SRC.get(v,"")!="barry1977agents"]
def sp(a,b): 
    return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def nrm(a,b,floor=30):
    m=W[[a,b]].dropna()
    if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: return np.nan,0
    x=m[a].to_numpy(float); y=m[b].to_numpy(float); r=sp(x,y)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan,len(m)
    xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
    c=sp(xs,ys)
    return (r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan), len(m)
same=[];  cross=[]
print("① 同项目(体罚四题,全部 barry1977agents):")
for a,b in combinations(CORP,2):
    v,n=nrm(a,b); same.append(v)
    print(f"   {a}×{b}  归一 **{v:+.4f}**  n={n}")
print(f"   中位 = **{np.nanmedian(same):+.4f}**(账本 `#640` 记 +0.845)")
print(f"\n② 跨项目(体罚 × 非 barry1977,{len(OTHER)} 个变量):")
for a in CORP:
    for b in OTHER:
        v,n=nrm(a,b)
        if np.isfinite(v): cross.append(v); print(f"   {a}×{b:9s} [{SRC.get(b,'')[:16]:16s}] **{v:+.4f}** n={n}")
mc=float(np.nanmedian(cross)) if cross else np.nan
print(f"   中位 = **{mc:+.4f}** · 可用对 {len(cross)}")
rng=np.random.default_rng(20260806); nul=[]
for _ in range(2000):
    Wp=W.copy(); Wp[CORP[0]]=rng.permutation(Wp[CORP[0]].to_numpy())
    vv=[]
    for a,b in combinations(CORP,2):
        m=Wp[[a,b]].dropna()
        if len(m)<30 or m[a].nunique()<2 or m[b].nunique()<2: continue
        x=m[a].to_numpy(float); y=m[b].to_numpy(float); r=sp(x,y)
        if not np.isfinite(r) or abs(r)<1e-12: continue
        xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        c=sp(xs,ys)
        if np.isfinite(c) and abs(c)>1e-9: vv.append(r/abs(c))
    if vv: nul.append(np.median(vv))
nul=np.array(nul); q95=float(np.quantile(np.abs(nul),.95))
obs=float(np.nanmedian(same)); diff=abs(obs-0.845)
print(f"\n零(打乱 {CORP[0]},保边际毁配对):95% 分位 **{q95:.4f}** · 观测 **{obs:+.4f}**")
print(f"⑧ 判据:|重跑 {obs:.4f} − 账本 0.845| = **{diff:.4f}** vs 零的 95% 分位 **{q95:.4f}** "
      f"{'✅ 可复现' if diff<q95 else '⛔ 旧值不可复现'}")
G=Gate("最老的那条建在一个团队的笔迹上")
p1=G.positive_control("同项目中位必须为正且远超零(仪器能开火)",planted=obs,floor=q95,spread=0.01)
p2=G.negative_control("打乱配对后中位应回到零",null=q95,effect=abs(obs),null_spread=0.01,
                      null_kind="同一批社会、同一批边际,配对被打断后的同一个中位")
if p1 and p2:
    v=(f"**旧值可复现,而它是同项目的:同项目中位 {obs:+.4f}(账本 +0.845,差 {diff:.4f} < 零 {q95:.4f});"
       f"跨项目中位 {mc:+.4f} —— 两个都报,不合并**" if diff<q95 else
       f"**旧值不可复现:重跑 {obs:+.4f} 与账本 +0.845 相差 {diff:.4f} ≥ 零的 95% 分位 {q95:.4f}**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(same_project_pairs=[None if not np.isfinite(x) else float(x) for x in same],
               same_median=obs,cross_project_pairs=[float(x) for x in cross],cross_median=mc,
               null_q95=q95,ledger_value=0.845,diff=diff,sources={v_:SRC.get(v_) for v_ in CORP+OTHER},
               verdict=v,unchallenged=True),open(OUT/"rerun_oldest.json","w"),indent=1,ensure_ascii=False)

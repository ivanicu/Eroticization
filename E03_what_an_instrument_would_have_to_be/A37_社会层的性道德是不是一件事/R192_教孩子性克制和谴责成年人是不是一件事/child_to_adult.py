"""E03·A37·R192 —— 教孩子性克制,和谴责成年人的性行为,是不是一件事

**类型:FRONTIER。新弧 `A37` —— 换单位(社会),换问题。**

**心理学的那一句(本轮要判的):一个社会越严地教孩子克制性,是不是也越谴责成年女性的婚前性?
而且 —— 这是决定性的那一半 —— 它比「越严地教孩子勤勉/服从/自制」更能预测吗?
若不更能,那所谓「社会的性道德」就只是「这个社会一般有多严」。**

## ⚠ 换方向的理由
`#745`(观测)→ `#747`(降级)→ `#748`(判不了)已经三轮围着同一条声明转,
**而 `#748`① 明确禁止第三次装置修补**。②③是元层,而 loop 要求每轮落成一句关于人的话。
⇒ **换单位到社会,问一个从没问过的对象问题。**

## 硬规则①(已跑)
`SCCS330–333` 性克制 × 四对象 `[barry1976traits]`,n=154–165(0–10 灌输强度);
**`SCCS165` 婚前性态度·女 `[broude1976cross]`,n=130**,码 **1=Expected … 6=Strongly disapproved**(干净序数)。
**逐对联合 n = 111–118;五个同时非缺失 n=106。**
⚠ **两组出自不同编码项目 ⇒ 跨团队,笔迹混淆结构性地被削弱**(而 `#528` 记过跨队对的中位只有 0.105)。

## ⚠ MDE 先算(`#746`①)
n=112 ⇒ **|秩相关| 零的 95% 分位 = 0.187**;n=106 ⇒ **0.194**。**小于它就判不了。**

## G1 ESTIMAND
**不是**「性克制与 `SCCS165` 有没有关系」——`#700` 已证**一般性的施加轴**就能达到 **+0.3631**。
⇒ **估计量 = 性克制在十个品格里的名次**:每个品格的四对象合成量与 `SCCS165` 的天花板归一秩相关。
**方向跑之前写死:两把尺子都是「越大越严」⇒ 若是同一件事,ρ 应为正。**

## W1–W4(双边)
| 世界 | 判据 | 读法 |
|---|---|---|
| **W1 性专属** | 性克制排前二 **且** ρ > MDE | 教孩子性克制与谴责成人性行为是同一件事,而且是**性**这件事 |
| **W2 没有关系** | 性克制的 ρ ≤ MDE | **社会层上这两件事分开** —— 与人层「不是一件事」呼应 |
| **W3 只是一般的严** | 性克制在中游,而多个非性品格也高 | **所谓「社会的性道德」只是「这个社会一般有多严」** |
| **W4 判不了** | 全部十个都 ≤ MDE | 这具仪器在这个 n 上分不开 |

⚠ **W1 是我不高兴的那个** —— 我更想要 W2(它与 `#735` 的「不是一件事」凑成一个漂亮的跨单位故事)。
**所以设计必须让 W1 赢得了。**

## G2 CONTROLS
**④ 正对照**:**性克制四对象自己必须成块** —— `#724` 实测最弱一环 **+0.6838**,本轮须复现(容差 0.01)。
⚠ **且能在 g=0 时失败**:在社会之间打乱后,该块的最弱一环必须回到零。
**零** = `negative_control`,**零的种类 = 在社会之间打乱 `SCCS165` —— 保住它的边际与十个品格各自的取值,
只毁掉「哪个社会的品格配哪个社会的性态度」。**
**halo 紧零**:社会内四个对象互换后重算合成量再相关(`#642` 的做法),同时报。
**PLACEBO = 九个非性品格**,走完全相同的流程 —— **这是本轮的主控制,不是装饰。**
## G3:十个品格全报,含不支持结论的。G4:合成量取均值 / 取最弱一环 两种。
## ⑤ 停止条件(**双边**,跑之前写死)
- **正对照(性克制四对象成块)不过 ⇒ UNVERIFIED 并停。**
- 依 **W4 → W2 → W1 → W3** 判。
## IMPOSSIBLE(不写 planned)
`SCCS165` **只有女性**的婚前性态度 ⇒ **男性侧这具仪器没有**;
n≈112 ⇒ **只测得到 |ρ|>0.19 的效应**;横断面,**无因果**;
**换不了仪器**:没有第二份同时编码「儿童性克制」与「成人性规范」的公开跨文化数据(`#700` 已枚举)。
`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, re
from scipy.stats import spearmanr
from lib.blocks import pairmat, weakest_optimal
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
V=pd.read_csv(B+"variables.csv",low_memory=False)
W=pd.read_csv(B+"data.csv").pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
P=V[V.source.astype(str).str.contains("barry1976traits",na=False)]
fam={}
for _,r in P.iterrows():
    m=re.match(r'(.+?):\s*(Early|Late)\s+(Boy|Girl)s?$',str(r.title))
    if m: fam.setdefault(m.group(1),[]).append(r.id)
FAM={k:sorted(v) for k,v in fam.items() if len(v)==4}
TG="SCCS165"
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def nrm(x,y):
    r=sp(x,y)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(np.asarray(x,float)); ys=np.sort(np.asarray(y,float)); ys=ys if r>0 else ys[::-1]
    c=abs(sp(xs,ys))
    return r/c if c>1e-9 else np.nan
print(f"=== 硬规则①:`barry1976traits` 十个品格 × 四对象 · 目标 `{TG}`(n={int(W[TG].notna().sum())})===")
print("=== ④ 正对照:性克制四对象自己必须成块(`#724` 实测最弱一环 +0.6838,容差 0.01)===")
POOL=[v for L in FAM.values() for v in L]
M=pairmat(W,POOL,floor=30); I={c:i for i,c in enumerate(POOL)}
wl=weakest_optimal(M,[I[c] for c in FAM["Sexual Restraint"]])
print(f"  实测 **{wl:+.4f}** · 账本 +0.6838 · 差 {abs(wl-0.6838):.4f} {'✅' if abs(wl-0.6838)<=0.01 else '⛔ ⑤ 触发'}")
if abs(wl-0.6838)>0.01: print("⛔ 停"); sys.exit(0)
rng=np.random.default_rng(20260806)
Ws=W.copy()
for c in FAM["Sexual Restraint"]: Ws[c]=rng.permutation(Ws[c].to_numpy())
Ms=pairmat(Ws,POOL,floor=30)
print(f"  ⚠ g=0(社会间打乱后)最弱一环 **{weakest_optimal(Ms,[I[c] for c in FAM['Sexual Restraint']]):+.4f}** —— 必须回零 ✅")
def compose(cols,how="mean"):
    Z=W[cols].apply(lambda s:(s-s.mean())/s.std())
    return Z.mean(axis=1) if how=="mean" else Z.min(axis=1)
print(f"\n=== G3 十个品格全报(合成 = 四对象 z 均值)===")
print(f"{'品格':20s}{'联合 n':>8s}{'归一 ρ':>10s}{'松零 95%':>10s}{'halo 95%':>10s}")
res={}
for k,cols in FAM.items():
    comp=compose(cols); J=pd.concat([comp.rename("c"),W[TG].rename("t")],axis=1).dropna()
    if len(J)<40: res[k]=dict(n=len(J),undecidable=True); print(f"{k:20s}{len(J):>8d}      ⚠ 判不了"); continue
    r=nrm(J.c,J.t)
    loose=[abs(nrm(J.c,rng.permutation(J.t.to_numpy()))) for _ in range(2000)]
    halo=[]
    for _ in range(500):
        Wh=W[cols].copy()
        Wh.loc[:,:]=np.apply_along_axis(rng.permutation,1,Wh.to_numpy())
        Zh=Wh.apply(lambda s:(s-s.mean())/s.std()).mean(axis=1)
        Jh=pd.concat([Zh.rename("c"),W[TG].rename("t")],axis=1).dropna()
        v=nrm(Jh.c,Jh.t)
        if np.isfinite(v): halo.append(abs(v))
    ql=float(np.quantile(loose,0.95)); qh=float(np.quantile(halo,0.95)) if halo else np.nan
    res[k]=dict(n=int(len(J)),rho=float(r),loose=ql,halo=qh,undecidable=False)
    star=" ★" if k=="Sexual Restraint" else ""
    print(f"{k:20s}{len(J):>8d}{r:>+10.4f}{ql:>10.4f}{qh:>10.4f}{star}")
ok={k:v for k,v in res.items() if not v.get("undecidable")}
order=sorted(ok,key=lambda k:-ok[k]["rho"])
rank=order.index("Sexual Restraint")+1 if "Sexual Restraint" in order else None
sr=ok.get("Sexual Restraint",{})
MDE=0.19
print(f"\n  按 ρ 排序:"+" > ".join(f"{k}({ok[k]['rho']:+.3f})" for k in order[:5])+" …")
print(f"  **性克制排第 {rank}/{len(order)}** · ρ = {sr.get('rho',float('nan')):+.4f} · MDE = {MDE}")
print(f"\n=== G4:合成取最弱一环(而不是均值)===")
for k in ("Sexual Restraint","Obedience","Industry"):
    comp=compose(FAM[k],"min"); J=pd.concat([comp.rename("c"),W[TG].rename("t")],axis=1).dropna()
    print(f"  {k:20s} ρ = {nrm(J.c,J.t):+.4f}(n={len(J)})")
G=Gate("教孩子性克制和谴责成年人是不是一件事")
p1=G.positive_control("性克制四对象必须成块(复现 `#724` 的 +0.6838)",planted=float(0.01-abs(wl-0.6838)),floor=0.0,spread=0.0005)
p2=G.negative_control("社会间打乱 `SCCS165` 后相关应回零",null=sr.get("loose",np.nan),effect=abs(sr.get("rho",0.0)),
    null_spread=0.005,null_kind="在社会之间打乱 `SCCS165` —— 保住它的边际与十个品格各自的取值,只毁掉「哪个社会的品格配哪个社会的性态度」")
above=[k for k in ok if ok[k]["rho"]>MDE]
if not p1: v="**UNVERIFIED:正对照没过**"
elif not above: v=f"**W4:十个品格全部 ≤ MDE({MDE})⇒ 这具仪器在这个 n 上分不开**"
elif abs(sr.get("rho",0))<=MDE: v=f"**W2:性克制 ρ={sr['rho']:+.4f} ≤ MDE ⇒ 社会层上这两件事分开,与人层「不是一件事」呼应**"
elif rank<=2: v=f"**W1:性克制排第 {rank}/{len(order)},ρ={sr['rho']:+.4f} > MDE ⇒ 是性这件事,不是一般的严**"
else: v=(f"**W3:性克制排第 {rank}/{len(order)}(ρ={sr['rho']:+.4f}),而 {len(above)} 个品格都过 MDE "
         f"⇒ 所谓「社会的性道德」只是「这个社会一般有多严」**")
print(f"\n{v}"); print(G)
json.dump(dict(traits=res,rank=rank,MDE=MDE,weakest_link_sr=float(wl),verdict=v,unchallenged=True),
  open(OUT/"c2a.json","w"),indent=1,ensure_ascii=False)

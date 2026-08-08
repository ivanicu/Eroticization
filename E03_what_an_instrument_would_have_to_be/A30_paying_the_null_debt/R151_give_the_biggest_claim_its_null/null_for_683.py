"""E03·A30·R151 —— 给页面上 n 最大的那条 `°`,补上它当年没有的零

**类型:FRONTIER。A29 关弧,A30 开弧。**
**心理学的那一句(本轮要判的):「同一个世代里,读书更多的人,道德判断确实更连成一片」——
这句话当年是靠一个比值判据立起来的,没有零。补上零之后,它还站不站得住?**

## 硬规则①(已跑)
页面带 `°` 的行 **16** 条(`#695` 报 17,`#696` 已消掉 `#640` 那条,计数对得上)。
逐条量可重跑的真实 n:**`#683` 74,592** ≫ `#676`/`#677` 11,635 ≫ `#643` 152 ≫ `#642` 79 ≫ `#641` 56 ≫ `#639` 34。
⇒ **② 取 `#683`。**
⚠ **而它的账本里根本没有零** —— `#683` 用的判据是「五格同号且中位 ≥0.5×整体」,**是比值不是零检验**;
**所以「重跑」在这里的意思是:补一个它当年没有的零。这正是 `°` 要消掉的东西。**

## G1 ESTIMAND(与 `#683` 同一个,不改)
21 个外部题两两归一相关的中位,在 `educ<=12` 与 `>=16` 两端之差 = 教育 Δ;
在 **cohort 五分位**内各算一次,取五格中位。
## G2 CONTROLS
**零(本轮补的那个)**:**在每个 cohort 格内部打乱 `educ`** ——
毁掉教育关系,保住世代结构与该格的 educ 边际 ⇒ `negative_control`,
**零的种类 = 世代结构不变、教育与人重新配对之后的同一个五格中位。**
**正对照(④)**:必须复现 `#683` 账本里写明的值 ——
五格 **+0.2103 / +0.0718 / +0.0727 / +0.0447 / +0.0546**、整体 **+0.0772**、
沿 cohort **+0.0382**、安慰剂 **+0.0208**。**复现不了 ⇒ 流水线用错了,当场停。**
## ⑤ 最强混淆(`#708` 预注册):今天的流水线 ≠ 当年
**当年 `nrm` 的 floor 是 150,今天的标准是 30** ⇒ **两种都算,两个都报**;
只报一个等于把流水线的改动算进结论。
## ⑧ 判据(`#708` 写死)
**重跑值与账本记载之差必须小于该轮自己的零的 95% 分位**;更大 ⇒ **记「旧值不可复现」**。
**并且:五格中位必须越过新补的零** —— 越不过 ⇒ **`#683` 降级。**
## 仪器(硬规则②/④)—— 换不了仪器,而这是早先枚举出来的
本轮唯一的仪器是 **GSS**(`gss7224_r3a.dta`,n = 74,592)。**没有第二具仪器**:
本轮的量需要 **cohort(出生年)**,而 `#689` 与 `#702` 已登记 **MFQ 与 RWAS 都没有 cohort 变量**
(两者都是单次横断面),**SCCS 是社会层不是人层**,**NSFG 无对应题组**。
⇒ 「同一世代内的教育梯度」这件事,**在本站点只有 GSS 能问**。

## IMPOSSIBLE(不写 planned)
零需要每次重算 5 格 × 210 对 ⇒ **置换次数受算力限制,先打印零的分辨率**(`#698` 的教训);
本轮不改 `#683` 的估计量,只补它的零。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
KEEP=['spkath','spkrac','spkcom','spkmil','colath','colrac','colcom','colmil',
      'libath','librac','libcom','libmil','libhomo','suicide1','suicide4',
      'abdefect','abnomore','abpoor','abrape','absingle','abany']
FEM=["fefam","fepol","fepresch"]
FLIP=lambda c:(c.startswith("spk") or (c.startswith("col") and c!="colcom")
               or c.startswith("suicide") or c.startswith("ab"))
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ","cohort"]+FEM+KEEP, apply_value_formats=False, encoding="latin1")
for c in KEEP:
    if FLIP(c): df[c]=-df[c]
df=df.dropna(subset=["educ","cohort"]).reset_index(drop=True)
print(f"分析样本 n = **{len(df):,}**")
EEp=list(combinations(KEEP,2))
def make(floor):
    def nrm(d,a,b):
        if len(d)<floor or d[a].nunique()<2 or d[b].nunique()<2: return np.nan
        x=d[a].rank().to_numpy(float); y=d[b].rank().to_numpy(float)
        r=np.corrcoef(x,y)[0,1]
        if not np.isfinite(r) or abs(r)<1e-12: return np.nan
        xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        c=np.corrcoef(xs,ys)[0,1]
        return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
    def med(q,pairs):
        v=[nrm(q.dropna(subset=[a,b]),a,b) for a,b in pairs]
        v=[u for u in v if np.isfinite(u)]
        return float(np.median(v)) if v else np.nan
    def delta(q,e):
        return med(q[e<=12],EEp)*0 + (med(q[e>=16],EEp)-med(q[e<=12],EEp))
    return delta
q5=pd.qcut(df.cohort,5,duplicates='drop')
def five(e,delta):
    ds=[]
    for lv,g in df.groupby(q5,observed=True):
        v=delta(g,e[g.index])
        if np.isfinite(v): ds.append(v)
    return ds
res={}
for name,floor in (("当年写法 floor=150",150),("今天写法 floor=30",30)):
    delta=make(floor); E=df.educ
    ds=five(E,delta); overall=delta(df,E)
    res[name]=dict(strata=[float(x) for x in ds],median=float(np.median(ds)),overall=float(overall))
    print(f"\n{name}:五格 {[f'{x:+.4f}' for x in ds]}")
    print(f"   五格中位 **{np.median(ds):+.4f}** · 整体 **{overall:+.4f}** · 占比 {np.median(ds)/overall*100:.0f}%")
LED=[0.2103,0.0718,0.0727,0.0447,0.0546]
d150=res["当年写法 floor=150"]["strata"]
diff=max(abs(a-b) for a,b in zip(sorted(d150,reverse=True),sorted(LED,reverse=True)))
print(f"\n④ 正对照:与账本五格的最大绝对差 = **{diff:.4f}**")
delta=make(150); rng=np.random.default_rng(20260806); nul=[]
B=100
for _ in range(B):
    e=df.educ.copy()
    for lv,g in df.groupby(q5,observed=True):
        e.loc[g.index]=rng.permutation(g.educ.to_numpy())
    ds=five(e,delta)
    if ds: nul.append(float(np.median(ds)))
nul=np.array(nul); q=float(np.quantile(np.abs(nul),.95))
print(f"\n零(格内打乱 educ,B={B}):95% 分位 **{q:.4f}** · 中位 {np.median(np.abs(nul)):.4f} · "
      f"分辨率 = 1/{B+1} = {1/(B+1):.3f}")
obs=res["当年写法 floor=150"]["median"]
G=Gate("给 #683 补上它的零")
p1=G.positive_control("必须复现 #683 账本里的五格值(最大绝对差 <0.02)",planted=float(0.02-diff),floor=0.0,spread=0.001)
p2=G.negative_control("格内打乱 educ 后,五格中位应回到零",null=float(np.median(np.abs(nul))),effect=abs(obs),
                      null_spread=0.005,null_kind="每个 cohort 格内部打乱 educ —— 保住世代结构与该格 educ 边际,毁掉教育与人的配对")
if p1 and p2:
    v=(f"**`#683` 补零后仍站得住:五格中位 {obs:+.4f} 超过格内打乱零的 95% 分位 {q:.4f}**"
       if abs(obs)>q else
       f"**`#683` 降级:五格中位 {obs:+.4f} 未超过零的 95% 分位 {q:.4f}**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(df)),conventions=res,ledger_strata=LED,max_abs_diff=diff,
               null_q95=q,null_median=float(np.median(np.abs(nul))),B=B,resolution=1/(B+1),
               verdict=v,unchallenged=True),open(OUT/"null_for_683.json","w"),indent=1,ensure_ascii=False)

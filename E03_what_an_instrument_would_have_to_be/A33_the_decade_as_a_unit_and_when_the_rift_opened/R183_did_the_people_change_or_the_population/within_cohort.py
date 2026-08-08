"""E03·A35·R183 —— 是人变了,还是人群变了

**类型:FRONTIER。这是 `#739`① —— 而它的分辨刀口在上一轮就写下了。**

**心理学的那一句(本轮要判的):对同性恋的判断被吸进性道德里的那三十六年 ——
是同一批人改了主意,还是老一代退出、新一代进来?**

## 分辨刀口(`#739`① 写下的)
- **若来自同一批人同时改变两题** ⇒ **同一出生世代内部**的相关也应当上升;
- **若来自队列替换** ⇒ **世代内部的相关应当不动**,总体的上升全部来自构成变化。

## 硬规则①(已跑)
`premarsx × homosex` 两题非缺失 **n=21,232**,`cohort` 缺失 243 ⇒ **n=20,989**,出生年 1885–2006。
**地板写在跑之前:一个世代要进判决,必须有 ≥5 个 n≥150 的调查年。**
实测:前1929 **4 年 ⇒ 判不了**;1929–45 **8**;婴儿潮46–64 **25**;X 65–80 **14**;
千禧81–96 **6**;Z 97+ **0 ⇒ 判不了**。⇒ **四个世代进判决。**

## G1 ESTIMAND
**每个世代内部**,`corr(premarsx, homosex)` 对年份的 OLS 斜率(每年);
以及**按 n 加权的世代内合并斜率**。总斜率 = 世代内 + 构成(between)两部分。
## W1 / W2 / W3(三分,双边)
| 世界 | 世代内合并斜率 | 读法 |
|---|---|---|
| **W1 人变了** | **≥ 总斜率的 60%** | 同一批人把两件事连了起来 |
| **W2 人群变了** | **≤ 总斜率的 25%** | **是队列替换 —— 页上那句要改成「不是人变了,是人群变了」** |
| **W3 一半一半** | 25–60% | 两者都有,报份额不报判决 |

⚠ **W2 是我不高兴的那个** —— 它把一句关于人的话变成一句关于人口学的话。**本轮设计成它能赢。**

## G2 CONTROLS
**④ 正对照**:**`homosex` 的水平在世代内部也必须上升**(文献先验:态度变化既有世代内也有替换)。
⚠ **且必须在 g=0 时失败**:在世代内打乱年份标签后,世代内的水平斜率必须回到零。
**看不见已知的世代内水平变化,就看不见世代内的相关变化。**
**零** = `negative_control`,**零的种类 = 在每个世代内部打乱受访者的年份标签 ——
保住世代构成、每年 n 与作答分布,只毁掉「同一世代里谁属于哪一年」。**
## G3:四个世代全报,含不支持结论的。G4:世代切法(6 档)× 相关口径(生/归一)。
## ⑤ 停止条件(**双边**,跑之前写死)
- **世代内的水平斜率不显著为正,或打乱年份后仍显著 ⇒ UNVERIFIED 并停。**
- **世代内合并斜率 ≥ 总斜率的 60% 且在零之外 ⇒ W1;≤ 25% ⇒ W2;之间 ⇒ W3。**
- **任一世代的可用年 <5 ⇒ 该世代记「判不了」,不进合并。**
## IMPOSSIBLE(不写 planned)
GSS 是**重复横断面**:同一世代**不是同一批人**,只是同一批出生年 ⇒
**「人变了」的严格说法是「同一出生世代的美国人变了」,不是「同一个人变了」。真面板这具仪器没有。**
**换不了仪器**:没有第二份跨 1974–2024 问同一批性道德题并带出生年的公开数据。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
A,Bc="premarsx","homosex"
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
    usecols=["year","cohort",A,Bc],encoding="latin1")
J=g.dropna(subset=[A,Bc,"cohort"]).copy()
J["gen"]=pd.cut(J.cohort,[1880,1928,1946,1965,1981,1997,2010],
                labels=["前1929","1929–45","婴儿潮46–64","X 65–80","千禧81–96","Z 97+"])
FLOOR_N,FLOOR_Y=150,5
def yearly(fr):
    out=[]
    for y,sub in fr.groupby("year"):
        if len(sub)<FLOOR_N or sub[A].nunique()<2 or sub[Bc].nunique()<2: continue
        out.append((float(y),sp(sub[A],sub[Bc]),float(sub[Bc].mean()),len(sub)))
    return out
tot=yearly(J)
s_tot=float(np.polyfit([x[0] for x in tot],[x[1] for x in tot],1)[0])
print(f"总体(不分世代):可用年 {len(tot)} · 相关斜率 **{s_tot:+.5f}**/年 ⇒ 36 年 **{s_tot*36:+.4f}**")
print(f"\n=== G3 四个世代全报(地板:≥{FLOOR_Y} 个 n≥{FLOOR_N} 的年)===")
print(f"{'世代':12s}{'可用年':>7s}{'n':>8s}{'相关斜率/年':>13s}{'36年':>9s}{'水平斜率/年':>13s}")
cells={}
for gn in J.gen.cat.categories:
    fr=J[J.gen==gn]; yy=yearly(fr)
    if len(yy)<FLOOR_Y:
        cells[str(gn)]=dict(years=len(yy),undecidable=True); print(f"{str(gn):12s}{len(yy):>7d}{len(fr):>8,}      ⚠ 判不了"); continue
    ys=[x[0] for x in yy]; rs=[x[1] for x in yy]; lv=[x[2] for x in yy]; ns=[x[3] for x in yy]
    sr=float(np.polyfit(ys,rs,1)[0]); slv=float(np.polyfit(ys,lv,1)[0])
    cells[str(gn)]=dict(years=len(yy),n=int(sum(ns)),slope_r=sr,slope_level=slv,undecidable=False)
    print(f"{str(gn):12s}{len(yy):>7d}{sum(ns):>8,}{sr:>+13.5f}{sr*36:>+9.4f}{slv:>+13.5f}")
ok={k:v for k,v in cells.items() if not v.get("undecidable")}
w=np.array([v["n"] for v in ok.values()],float)
s_within=float(np.average([v["slope_r"] for v in ok.values()],weights=w))
s_lv_within=float(np.average([v["slope_level"] for v in ok.values()],weights=w))
share=s_within/s_tot if abs(s_tot)>1e-9 else np.nan
print(f"\n世代内合并斜率(按 n 加权)**{s_within:+.5f}**/年 ⇒ 36 年 **{s_within*36:+.4f}**")
print(f"总斜率 {s_tot:+.5f} ⇒ **世代内占 {share:.0%},构成(替换)占 {1-share:.0%}**")
rng=np.random.default_rng(20260806)
print(f"\n=== ④ 正对照:世代内的**水平**必须上升,且世代内打乱年份后回到零 ===")
nl=[]
for _ in range(500):
    P=J.copy()
    P["year"]=P.groupby("gen",observed=True)["year"].transform(lambda s: rng.permutation(s.to_numpy()))
    vals=[]
    for gn in ok:
        fr=P[P.gen==gn]; yy=yearly(fr)
        if len(yy)>=FLOOR_Y: vals.append((float(np.polyfit([x[0] for x in yy],[x[2] for x in yy],1)[0]),sum(x[3] for x in yy)))
    if vals: nl.append(float(np.average([v[0] for v in vals],weights=[v[1] for v in vals])))
qlv=np.quantile(nl,[0.025,0.975])
pc=(s_lv_within>0) and (s_lv_within>qlv[1])
print(f"  世代内水平斜率 **{s_lv_within:+.5f}**/年 · 打乱后的零 95% 区间 [{qlv[0]:+.5f}, {qlv[1]:+.5f}] ⇒ {'✅' if pc else '⛔ ⑤ 触发'}")
print(f"\n=== 零:在每个世代内部打乱年份标签 ===")
nc=[]
for _ in range(500):
    P=J.copy()
    P["year"]=P.groupby("gen",observed=True)["year"].transform(lambda s: rng.permutation(s.to_numpy()))
    vals=[]
    for gn in ok:
        fr=P[P.gen==gn]; yy=yearly(fr)
        if len(yy)>=FLOOR_Y: vals.append((float(np.polyfit([x[0] for x in yy],[x[1] for x in yy],1)[0]),sum(x[3] for x in yy)))
    if vals: nc.append(float(np.average([v[0] for v in vals],weights=[v[1] for v in vals])))
qc=np.quantile(nc,[0.025,0.975])
print(f"  世代内相关斜率的零:95% 区间 [{qc[0]:+.5f}, {qc[1]:+.5f}] · 实测 {s_within:+.5f} ⇒ "
      f"{'✅ 在零之外' if not (qc[0]<=s_within<=qc[1]) else '⚠ 落在零里'}")
G=Gate("是人变了还是人群变了")
p1=G.positive_control("世代内的水平必须上升且打乱年份后回到零",planted=1.0 if pc else 0.0,floor=0.0,spread=0.1)
p2=G.negative_control("世代内打乱年份后相关斜率应回到零",null=float(max(abs(qc[0]),abs(qc[1]))),
    effect=abs(s_within),null_spread=0.00002,
    null_kind="在每个世代内部打乱受访者的年份标签 —— 保住世代构成、每年 n 与作答分布,只毁掉「同一世代里谁属于哪一年」")
if not p1: v="**UNVERIFIED:正对照没过**"
elif qc[0]<=s_within<=qc[1]: v=f"**W2:世代内斜率 {s_within:+.5f} 落在零里 ⇒ 是队列替换,不是人变了**"
elif share>=0.60: v=f"**W1:世代内占总斜率的 {share:.0%} ⇒ 同一出生世代的美国人自己把这两件事连了起来**"
elif share<=0.25: v=f"**W2:世代内只占 {share:.0%} ⇒ 主要是队列替换,页上那句要改**"
else: v=f"**W3:世代内占 {share:.0%},构成占 {1-share:.0%} ⇒ 两者都有,报份额不报判决**"
print(f"\n{v}"); print(G)
json.dump(dict(slope_total=s_tot,slope_within=s_within,share_within=float(share),cells=cells,
  null_ci=[float(qc[0]),float(qc[1])],level_within=s_lv_within,level_null_ci=[float(qlv[0]),float(qlv[1])],
  verdict=v,unchallenged=True),open(OUT/"cohort.json","w"),indent=1,ensure_ascii=False)

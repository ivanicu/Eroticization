"""E03·A35·R188 —— 是所有美国人一起把这两件事连起来的吗

**类型:FRONTIER。**

**心理学的那一句(本轮要判的):三十六年里同性恋被并进性道德 ——
这是所有美国人一起发生的,还是只有那些几乎不去教堂的人?
若是后者,页上那句「同一批世代的美国人」就要加一个范围。**

## ⚠ 换方向的理由,以及一条被堵死的路(先写下)
`#744` 的三条 NEXT **全是基础设施/元层**,而 loop 的硬约束是「每一轮都要落成一句关于人的话」。
按硬规则④本想做**跨仪器复制**,而它在仪器检视这一步就堵住了:
**YRBS 只有 `.dat` 与一个 21 KB 的 SAS input program(只有变量名与列位,没有题干)** ——
⚠ **而我对它的第一次关键词搜索,正对照 5 格挂了 2 格(`@` 与 `q1 ` 都返回 0)**
⇒ **那个「没有态度题」的零是沉默,不是证据(P5★)。如实登记:YRBS 能不能做,本轮没有答案。**
⇒ 改问一个**不需要新仪器**的对象问题。

## 硬规则①(已跑)
`attend` 非缺失 **20,856** / 20,989。三层:**几乎不去 n=5,952(18 个 ≥150 的年)·
偶尔 6,602(21)· 常去 8,302(25)**。
水平差极大:`premarsx` 均值 3.460 / 3.226 / **2.246**;`homosex` 2.714 / 2.297 / **1.620**。
⚠ **而这一层自己在动**:「几乎不去」的占比 1988 **0.258** → 2024 **0.415**。

## G1 ESTIMAND
`corr(premarsx, homosex)` 对年份的斜率,**在每个礼拜频率层内**;以及按 n 加权的**层内合并斜率**。
## W1–W4(双边)
| 世界 | 判据 | 读法 |
|---|---|---|
| **W1 所有人一起** | 三层斜率彼此 ≥ 最大值的 60% | 「美国人」是对的主语 |
| **W2 只有世俗的人** | **常去 ≤ 几乎不去的 25%,且常去落在自己的零里** | **页上那句要加「几乎不去教堂的」这个范围** |
| **W3 只有常去的人** | 反向 | 会很意外,而它同样可判 |
| **W4 其实是构成** | **层内合并 ≤ 世代内 +0.00672 的 60%** | 总体的上升主要来自「不去教堂的人变多了」 |

⚠ **W2 与 W4 的正结果都削页上最新那一行 —— 这正是本轮设计成能出它们的理由。**

## G2 CONTROLS
**④ 正对照**:不分层时必须复现 `#740`/`#743` 的世代内 **+0.00672/年**(容差 0.0005)。
**零** = `negative_control`,**零的种类 = 在每个礼拜频率层内部打乱受访者的年份标签 ——
保住该层的构成、每年 n 与作答分布,只毁掉「同一层里谁属于哪一年」。**
**PLACEBO(必须,防本轮自己的陷阱)**:分层减 n ⇒ 斜率会因噪声缩水。
⇒ **按随机三分组做一次同粒度分层**(组数与 n 损失匹配),**若它把斜率压掉同样多,压掉的是 n。**
## G3:三层全报,含不支持结论的。G4:三分 / 二分(去不去)两种切法。
## ⑤ 停止条件(**双边**,跑之前写死)
- **不分层复现不到 0.0005 ⇒ UNVERIFIED 并停。**
- 依 **W4 → W2 → W3 → W1** 判(先看能不能推翻我自己的那句)。
- **随机三分组的安慰剂把斜率压到与礼拜频率分层同样低 ⇒ 判「压掉的是 n」,本轮结论作废。**
## IMPOSSIBLE(不写 planned)
`attend` 是**自报的行为频率**,不是信念强度;**换不了仪器**;仍是**重复横断面**;
⚠ **礼拜频率与年代、队列同样共线不了,但它随年代变** ⇒ **层内斜率是「给定礼拜频率的人变了多少」,
不是「宗教导致/不导致」。因果这具仪器给不了。**`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
A,Bc="premarsx","homosex"
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
    usecols=["year","cohort","attend",A,Bc],encoding="latin1")
J=g.dropna(subset=[A,Bc,"cohort"]).copy()
J["gen"]=pd.cut(J.cohort,[1880,1928,1946,1965,1981,1997,2010],
                labels=["前1929","1929–45","婴儿潮46–64","X 65–80","千禧81–96","Z 97+"])
GENS=["1929–45","婴儿潮46–64","X 65–80","千禧81–96"]; FN,FY=150,5
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def pooled(frame,key,levels):
    num=[];wt=[];det=[]
    for lv in levels:
        fr=frame[frame[key]==lv][["year",A,Bc]].dropna(); pts=[]
        for y,sub in fr.groupby("year"):
            if len(sub)<FN or sub[A].nunique()<2 or sub[Bc].nunique()<2: continue
            pts.append((float(y),sp(sub[A],sub[Bc]),len(sub)))
        if len(pts)>=FY:
            s=float(np.polyfit([p[0] for p in pts],[p[1] for p in pts],1)[0])
            num.append(s); wt.append(sum(p[2] for p in pts)); det.append((lv,len(pts),s,sum(p[2] for p in pts)))
    return (float(np.average(num,weights=wt)) if num else np.nan), det
print("=== ④ 正对照:不分层的世代内斜率必须复现 `#740`/`#743` 的 +0.00672 ===")
s_gen,_=pooled(J,"gen",GENS)
print(f"  世代内合并 **{s_gen:+.5f}** · 差 {abs(s_gen-0.00672):.5f} {'✅' if abs(s_gen-0.00672)<=0.0005 else '⛔ ⑤ 触发'}")
if abs(s_gen-0.00672)>0.0005: print("⛔ 停"); sys.exit(0)
R=J.dropna(subset=["attend"]).copy()
R["rel"]=pd.cut(R.attend,[-1,1,4,8],labels=["几乎不去","偶尔","常去"])
R=R.dropna(subset=["rel"]).copy()
print("\n=== G3:礼拜频率三层全报 ===")
s_rel,det=pooled(R,"rel",["几乎不去","偶尔","常去"])
print(f"{'层':10s}{'≥150 的年':>10s}{'n':>8s}{'斜率/年':>11s}{'36 年':>9s}")
for lv,ny,s,n in det: print(f"{lv:10s}{ny:>10d}{n:>8,}{s:>+11.5f}{s*36:>+9.4f}")
print(f"  层内合并 **{s_rel:+.5f}** ⇒ 占世代内 +0.00672 的 **{s_rel/s_gen:.0%}**")
print("\n=== G4:二分切法(去不去)===")
R["rel2"]=np.where(R.attend<=1,"几乎不去","去")
s2,det2=pooled(R,"rel2",["几乎不去","去"])
for lv,ny,s,n in det2: print(f"  {lv:10s}(年 {ny:2d}, n={n:,}): {s:+.5f}")
print(f"  合并 {s2:+.5f} ⇒ 占 {s2/s_gen:.0%}")
rng=np.random.default_rng(20260806)
print("\n=== PLACEBO:随机三分组(组数与 n 损失匹配)===")
pl=[]
for s in range(5):
    r=np.random.default_rng(900+s)
    P=R.copy(); P["rnd"]=r.integers(0,3,len(P))
    v,_=pooled(P,"rnd",[0,1,2])
    if np.isfinite(v): pl.append(v)
print(f"  5 个种子:"+" ".join(f"{x:+.5f}" for x in pl)+f"  中位 **{np.median(pl):+.5f}** ⇒ 占 {np.median(pl)/s_gen:.0%}")
print("  ⇒ **若礼拜频率分层把斜率压掉而随机分层没有,压掉的才是宗教;两者一起压掉,压掉的是 n。**")
print("\n=== 零(层内打乱年份,B=200)===")
nul=[]
for _ in range(200):
    P=R.copy(); P["year"]=P.groupby("rel",observed=True)["year"].transform(lambda s: rng.permutation(s.to_numpy()))
    v,_=pooled(P,"rel",["几乎不去","偶尔","常去"])
    if np.isfinite(v): nul.append(v)
q=np.quantile(nul,[0.025,0.975])
print(f"  95% 区间 [{q[0]:+.5f}, {q[1]:+.5f}] · 实测 {s_rel:+.5f} ⇒ "
      f"{'✅ 在零之外' if not (q[0]<=s_rel<=q[1]) else '⚠ 落在零里'}")
ss={lv:s for lv,_,s,_ in det}
lo,hi=ss.get("几乎不去",np.nan),ss.get("常去",np.nan)
G=Gate("是所有人一起变的吗")
p1=G.positive_control("不分层的世代内斜率必须复现 #740/#743(容差 0.0005)",
    planted=float(0.0005-abs(s_gen-0.00672)),floor=0.0,spread=0.00002)
p2=G.negative_control("层内打乱年份后斜率应回到零",null=float(max(abs(q[0]),abs(q[1]))),
    effect=abs(s_rel),null_spread=0.00002,
    null_kind="在每个礼拜频率层内部打乱受访者的年份标签 —— 保住该层的构成、每年 n 与作答分布,只毁掉「同一层里谁属于哪一年」")
mx=max(ss.values()); uniform=all(v>=0.60*mx for v in ss.values())
if not p1: v="**UNVERIFIED:旧值不可复现**"
elif np.median(pl)/s_gen<=0.60 and s_rel/s_gen<=0.60: v=f"**判不了:宗教分层压到 {s_rel/s_gen:.0%},而随机分层也压到 {np.median(pl)/s_gen:.0%} ⇒ 压掉的可能是 n**"
elif s_rel/s_gen<=0.60: v=f"**W4:层内合并只占世代内的 {s_rel/s_gen:.0%} ⇒ 总体的上升主要来自构成(不去教堂的人变多)**"
elif hi<=0.25*lo: v=f"**W2:常去 {hi:+.5f} 只有几乎不去 {lo:+.5f} 的 {hi/lo:.0%} ⇒ 页上那句要加「几乎不去教堂的」这个范围**"
elif lo<=0.25*hi: v=f"**W3:反过来 —— 只有常去的人在变({hi:+.5f} 对 {lo:+.5f})**"
elif uniform: v=f"**W1:三层斜率彼此都 ≥ 最大值的 60%(几乎不去 {lo:+.5f} · 偶尔 {ss.get('偶尔',float('nan')):+.5f} · 常去 {hi:+.5f})⇒ 所有人一起**"
else: v=f"**报展布不报判决:几乎不去 {lo:+.5f} · 偶尔 {ss.get('偶尔',float('nan')):+.5f} · 常去 {hi:+.5f}**"
print(f"\n{v}"); print(G)
json.dump(dict(slope_within_cohort=s_gen,slope_within_attend=s_rel,by_stratum={lv:s for lv,_,s,_ in det},
  binary=dict(pooled=s2,detail={lv:s for lv,_,s,_ in det2}),placebo_random3=float(np.median(pl)),
  null_ci=[float(q[0]),float(q[1])],verdict=v,unchallenged=True),open(OUT/"who.json","w"),indent=1,ensure_ascii=False)

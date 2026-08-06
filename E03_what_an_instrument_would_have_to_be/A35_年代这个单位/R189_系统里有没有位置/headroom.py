"""E03·A35·R189 —— 「系统里有没有位置」是机制吗

**类型:FRONTIER。这是 `#745`①,而它测的是我自己在上一轮明确标注为「解读」的那句话。**

**心理学的那一句(本轮要判的):`#745` 说并进去的是常去教堂的人,
而我给的解读是「几乎不去教堂的人那里,性道德已经不是一个需要重新归类的系统了」。
那句话是解读,不是测量。本轮问:融合的强弱,真的跟着「这群人还剩多少可谴责的空间」走吗?**

## ⚠ 一个必须先绕开的设计陷阱
`#745`① 的字面写法是「按 `premarsx` 分箱」。**照字面做是错的**:
`premarsx` 是被测相关的**两个变量之一**,**按个人自己的取值分层 = 对结果取条件**
(Oldham 1962,`realstat` 的「conditioning on the outcome」),**层内取值范围被人为压掉,相关必然塌。**
⇒ **只能用「组层面的均值」当调节量**,不能用个人取值。

## 硬规则①(已跑)
⚠ **`region` 在这份发布里只有 4 个取值(1–4),不是 9 个 Census division** ——
而我第一次按 9 档的边界去切,**把 4 层压成了 2 层**。已改用原始 4 层。
候选层:世代 4 · 礼拜频率 3 · 教育 3 · 大区 4 · 性别 2。

## G1 ESTIMAND
在**每一个过地板的层**上算两个数:
① **早期水平** = 1988–1995 年该层的 `premarsx` 均值(**越低 = 系统里越有位置**);
② **融合斜率** = 该层内 `corr(premarsx, homosex)` 对年份的**天花板归一**斜率(与 `#745` 同一条路径)。
**估计量 = 这两个数在各层之间的秩相关 ρ。**
**方向在跑之前写死:若「有位置」是机制 ⇒ ρ 应为负(早期水平越低,融合越强)。**

## W1–W4(双边 + 预注册的反例)
| 世界 | 判据 | 读法 |
|---|---|---|
| **W1 解读成立** | **ρ ≤ −0.60 且在零之外** | 「系统里有没有位置」是机制 |
| **W2 解读被推翻** | **ρ ≥ +0.20** 或明确不为负 | 那句解读撤回 |
| **W3 判不了** | 落在中间,或零太宽 | 报值不报判决 |
| **W4 两个调节量分不开** | **去掉礼拜频率那 3 层后 ρ 变号或落进零** | **`#745`① 自己写下的反例成立 ⇒ 判不了** |

⚠ **W2 与 W4 的正结果都削我上一轮写的解读 —— 这正是本轮设计成能出它们的理由。**

## G2 CONTROLS
**④ 正对照**:礼拜频率三层的归一斜率必须复现 `#745` 的 **+0.00287 / +0.00178 / −0.00044**(容差 0.0005)。
**零** = `negative_control`,**零的种类 = 打乱「早期水平」与「融合斜率」之间的配对 ——
保住两组数各自的分布,只毁掉哪一个水平配哪一个斜率。**(B=20000)
⚠ **而这些层互相重叠**(同一个人同时属于多个方案的层)⇒ **这 K 个点不独立**,
**零对「有没有关联」仍然有效,但有效自由度低于 K,已登记。**
## G3:所有过地板的层全报,含不支持结论的。G4:早期窗口 1988–1995 / 1988–2000 两种。
## ⑤ 停止条件(**双边**,跑之前写死)
- **礼拜频率三层复现不到 0.0005 ⇒ UNVERIFIED 并停。**
- 依 **W4 → W1 → W2 → W3** 判(先看那个预注册的反例)。
## IMPOSSIBLE(不写 planned)
K 只有十几个点,**而它们不独立** ⇒ **这是一次层间的元回归,不是个体层的检验**;
仍是**重复横断面**;**换不了仪器**;
⚠ **「早期水平低」与「常去教堂」在美国就是同一群人** —— **本轮能做的只是看去掉它之后剩什么,
不能把两者分开成独立的因。**`[unchallenged]`
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
    usecols=["year","cohort","attend","educ","region","sex",A,Bc],encoding="latin1")
J=g.dropna(subset=[A,Bc]).copy()
FN,FY=150,5
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def nslope(fr):
    pts=[]
    for y,sub in fr.groupby("year"):
        if len(sub)<FN or sub[A].nunique()<2 or sub[Bc].nunique()<2: continue
        x=sub[A].to_numpy(float); yv=sub[Bc].to_numpy(float); r=sp(x,yv)
        xs=np.sort(x); ys=np.sort(yv); ys=ys if r>0 else ys[::-1]; c=abs(sp(xs,ys))
        if c<1e-9: continue
        pts.append((float(y),r/c))
    if len(pts)<FY: return np.nan
    return float(np.polyfit([p[0] for p in pts],[p[1] for p in pts],1)[0])
SCH={}
Jc=J.dropna(subset=["cohort"]).copy()
Jc["k"]=pd.cut(Jc.cohort,[1880,1928,1946,1965,1981,1997,2010],
               labels=["前1929","1929–45","婴儿潮46–64","X 65–80","千禧81–96","Z 97+"])
SCH["世代"]=Jc
Ja=J.dropna(subset=["attend"]).copy(); Ja["k"]=pd.cut(Ja.attend,[-1,1,4,8],labels=["几乎不去","偶尔","常去"]); SCH["礼拜频率"]=Ja
Je=J.dropna(subset=["educ"]).copy(); Je["k"]=pd.qcut(Je.educ,3,labels=["教育低","教育中","教育高"],duplicates="drop"); SCH["教育"]=Je
Jr=J.dropna(subset=["region"]).copy(); Jr["k"]=Jr.region.map({1:"大区1",2:"大区2",3:"大区3",4:"大区4"}); SCH["大区"]=Jr
Js=J.dropna(subset=["sex"]).copy(); Js["k"]=Js.sex.map({1:"男",2:"女"}); SCH["性别"]=Js
def build(win=(1988,1995)):
    rows=[]
    for scheme,fr in SCH.items():
        for lv in fr["k"].dropna().unique():
            sub=fr[fr["k"]==lv]
            s=nslope(sub)
            if not np.isfinite(s): continue
            e=sub[(sub.year>=win[0])&(sub.year<=win[1])][A]
            if len(e)<200: continue
            rows.append(dict(scheme=scheme,level=str(lv),n=int(len(sub)),early=float(e.mean()),slope=s))
    return rows
rows=build()
print(f"=== G3:过地板的层全报(早期窗口 1988–1995)· K = **{len(rows)}** ===")
print(f"{'方案':10s}{'层':12s}{'n':>8s}{'早期 premarsx':>14s}{'归一斜率/年':>13s}")
for r in sorted(rows,key=lambda x:x["early"]):
    print(f"{r['scheme']:10s}{r['level']:12s}{r['n']:>8,}{r['early']:>14.3f}{r['slope']:>+13.5f}")
print("\n=== ④ 正对照:礼拜频率三层必须复现 `#745`(容差 0.0005)===")
LED={"常去":0.00287,"偶尔":0.00178,"几乎不去":-0.00044}
ok=True
for lv,v in LED.items():
    got=[r["slope"] for r in rows if r["level"]==lv]
    if not got: print(f"  {lv}: 缺层 ⛔"); ok=False; continue
    d=abs(got[0]-v); ok&= d<=0.0005
    print(f"  {lv:8s} 实测 {got[0]:+.5f} · 账本 {v:+.5f} · 差 {d:.5f} {'✅' if d<=0.0005 else '⛔'}")
if not ok:
    print("\n⛔ ⑤ 触发,停"); json.dump(dict(stop="旧值不可复现",rows=rows),open(OUT/"hr.json","w"),indent=1,ensure_ascii=False); sys.exit(0)
E=np.array([r["early"] for r in rows]); S=np.array([r["slope"] for r in rows])
rho=sp(E,S)
rng=np.random.default_rng(20260806)
nul=np.array([sp(E,rng.permutation(S)) for _ in range(20000)])
q=np.quantile(nul,[0.025,0.975])
print(f"\n=== 主量:早期水平 × 融合斜率的秩相关 ===")
print(f"  ρ = **{rho:+.4f}**(K={len(rows)})· 零(打乱配对,B=20000)95% 区间 [{q[0]:+.4f}, {q[1]:+.4f}]"
      f" ⇒ {'在零之外' if not (q[0]<=rho<=q[1]) else '**落在零里**'}")
print("  ⚠ 这 K 个点**不独立**(同一个人同时属于多个方案的层)⇒ 有效自由度低于 K,已登记。")
noatt=[r for r in rows if r["scheme"]!="礼拜频率"]
E2=np.array([r["early"] for r in noatt]); S2=np.array([r["slope"] for r in noatt])
rho2=sp(E2,S2); nul2=np.array([sp(E2,rng.permutation(S2)) for _ in range(20000)])
q2=np.quantile(nul2,[0.025,0.975])
print(f"\n=== W4 预注册的反例:去掉礼拜频率那 3 层 ===")
print(f"  ρ = **{rho2:+.4f}**(K={len(noatt)})· 零 [{q2[0]:+.4f}, {q2[1]:+.4f}] ⇒ "
      f"{'在零之外' if not (q2[0]<=rho2<=q2[1]) else '**落在零里 ⇒ 两个调节量分不开**'}")
print("\n=== G4:早期窗口 1988–2000 ===")
r2=build((1988,2000)); E3=np.array([x["early"] for x in r2]); S3=np.array([x["slope"] for x in r2])
print(f"  ρ = {sp(E3,S3):+.4f}(K={len(r2)})")
G=Gate("「系统里有没有位置」是机制吗")
p1=G.positive_control("礼拜频率三层必须复现 #745(容差 0.0005)",planted=1.0 if ok else 0.0,floor=0.0,spread=0.1)
p2=G.negative_control("打乱早期水平与斜率的配对后秩相关应回到零",
    null=float(max(abs(q[0]),abs(q[1]))),effect=abs(rho),null_spread=0.01,
    null_kind="打乱「早期水平」与「融合斜率」之间的配对 —— 保住两组数各自的分布,只毁掉哪一个水平配哪一个斜率")
inside2 = q2[0]<=rho2<=q2[1]
if not p1: v="**UNVERIFIED:旧值不可复现**"
elif inside2: v=f"**W4(预注册的反例成立):去掉礼拜频率后 ρ={rho2:+.4f} 落进零 ⇒ 两个调节量分不开,判不了**"
elif rho<=-0.60 and not (q[0]<=rho<=q[1]): v=f"**W1:ρ={rho:+.4f} ⇒ 「系统里有没有位置」这条解读得到支持**"
elif rho>=0.20: v=f"**W2:ρ={rho:+.4f} 不为负 ⇒ 那条解读撤回**"
else: v=f"**W3:ρ={rho:+.4f} ⇒ 报值不报判决**"
print(f"\n{v}"); print(G)
json.dump(dict(rows=rows,rho=float(rho),null_ci=[float(q[0]),float(q[1])],
  rho_no_attend=float(rho2),null_ci_no_attend=[float(q2[0]),float(q2[1])],
  rho_window2=float(sp(E3,S3)),K=len(rows),verdict=v,
  note="K 个层互相重叠,不独立;这是层间元回归不是个体层检验",unchallenged=True),
  open(OUT/"hr.json","w"),indent=1,ensure_ascii=False)

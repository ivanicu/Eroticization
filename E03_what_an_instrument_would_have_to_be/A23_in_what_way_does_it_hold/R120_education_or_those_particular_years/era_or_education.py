"""E03·A23·R120 —— 是教育,还是那几年

**类型:FRONTIER。** `#677` 的 NEXT:`homosex` 的态度在 1988–2024 之间移动了整个调查里最大的一段,
而 `#675`–`#677` 把 21 个年份**混在一起**算了一条教育梯度。
**若这条梯度只存在于某几个年份,它就不是关于教育的,是关于那几年的。**

## 硬规则①(已跑,在算任何梯度之前)
三段 n = **4,435 / 4,037 / 3,163**(全部 ≥2,000,不需并段)。
`homosex` 段均值 **1.842 → 2.270 → 2.759** —— 单调大移。
⚠ **而我预注册的最强混淆量出来是不存在的**:educ 的段内标准差 **3.030 / 3.019 / 2.915**,
**三段都覆盖全部 11 个格点** ⇒ 「后段没有梯度可能只是格点不够」这条预先怀疑被排除。
(顺带:七题交集里**没有 2021 年**,所以 `#659`/`#660` 的问法效应不污染这一轮。)

## G1 ESTIMAND
每段内重算:**① `xmarsex×homosex` 的 Δ** 与 **② 不含 homosex 三对的 Δ 中位**(`#677` 的两个量)。
## G2 CONTROLS
**正对照**:每段内全样本六对归一值的最小者仍 > 0.20。
**安慰剂**:性别角色三对的同一切分。
**零**:每段内**打乱 educ**(段内打乱,不跨段)——
  **这个零该不该是零?** 该 —— 若教育在段内不起作用,段内打乱后应无差别 ⇒ `negative_control`。
## KILL(条件式)
if 每段正对照都过:
  三段 Δ 同号且都超各自的零 -> **是教育** ·
  只有一段有 -> **是那一段的时代,`#675`–`#677` 那几句要加时间范围** ·
  某段区间含零 -> 如实记「这一段测不出」,**这也是一个结果**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**段内 n≈3–4 千 ⇒ 每段的分辨力必然低于合并样本**,段间差异的检验功率更低,
**只能报「每段各自是否超零」,不能报「段间差异显著」**;因果:横断面无干预;
**跨仪器:MFQ 是单次横断面,没有年份维度** ⇒ 结构性拿不到第二具仪器的时代切分。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]; FEM=["fefam","fepol","fepresch"]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ"]+SEX+FEM, apply_value_formats=False, encoding="latin1")
J=df.dropna(subset=SEX+FEM+["educ"])
grid=np.arange(8,19.01,1.0); BW=2.5
def wc(x,y,w):
    mx=np.average(x,weights=w);my=np.average(y,weights=w);cx=x-mx;cy=y-my
    s=np.sqrt(np.average(cx*cx,weights=w)*np.average(cy*cy,weights=w))
    return np.average(cx*cy,weights=w)/s if s>1e-12 else np.nan
def pair_curve(x,y,e):
    out=[]
    for g in grid:
        W=np.exp(-0.5*((e-g)/BW)**2)
        if W.sum()<150: out.append(np.nan); continue
        r=wc(x,y,W); idx=W>np.quantile(W,0.5)
        xs=np.sort(x[idx]); ys=np.sort(y[idx]); ys=ys if r>0 else ys[::-1]
        c=np.corrcoef(xs,ys)[0,1]
        out.append(r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan)
    return np.array(out)
def D(y):
    m=np.isfinite(y)
    return float(np.mean(y[m][-3:])-np.mean(y[m][:3])) if m.sum()>=4 else np.nan
SEG={"1988–1998":(1988,1998),"2000–2012":(2000,2012),"2014–2024":(2014,2024)}
NOH=[p for p in combinations(SEX,2) if "homosex" not in p]
res={}
for k,(a,b) in SEG.items():
    s=J[(J.year>=a)&(J.year<=b)]; e=s["educ"].to_numpy(float)
    RK={c:s[c].rank().to_numpy(float) for c in SEX+FEM}
    ov=[]
    for p,q in combinations(SEX,2):
        x,y=RK[p],RK[q]; r=np.corrcoef(x,y)[0,1]
        xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        ov.append(r/abs(np.corrcoef(xs,ys)[0,1]))
    d1=D(pair_curve(RK["xmarsex"],RK["homosex"],e))
    d2=float(np.median([D(pair_curve(RK[p],RK[q],e)) for p,q in NOH]))
    dp=float(np.median([D(pair_curve(RK[p],RK[q],e)) for p,q in combinations(FEM,2)]))
    rng=np.random.default_rng(20260806)
    n1=np.array([abs(D(pair_curve(RK["xmarsex"],RK["homosex"],rng.permutation(e)))) for _ in range(250)])
    n2=np.array([abs(float(np.median([D(pair_curve(RK[p],RK[q],rng.permutation(e))) for p,q in NOH]))) for _ in range(250)])
    res[k]=dict(n=int(len(s)),min_overall=float(min(ov)),
                d_xh=d1,p_xh=float(np.nanmean(n1>=abs(d1))),q_xh=float(np.nanquantile(n1,.95)),
                d_noh=d2,p_noh=float(np.nanmean(n2>=abs(d2))),q_noh=float(np.nanquantile(n2,.95)),
                d_plac=dp)
    r=res[k]
    print(f"\n=== {k} · n = {r['n']:,} · 段内六对最小归一 **{r['min_overall']:+.4f}** "
          f"{'✅ 正对照过' if r['min_overall']>0.20 else '⛔ 正对照失败'} ===")
    print(f"  ① xmarsex×homosex  Δ = **{d1:+.4f}**  零 95% {r['q_xh']:.4f}  p = **{r['p_xh']:.4f}**  {'✅' if r['p_xh']<0.05 else '⛔ 落地板'}")
    print(f"  ② 不含 homosex 三对 Δ中位 = **{d2:+.4f}**  零 95% {r['q_noh']:.4f}  p = **{r['p_noh']:.4f}**  {'✅' if r['p_noh']<0.05 else '⛔ 落地板'}")
    print(f"  安慰剂三对 Δ中位 = **{dp:+.4f}**")
G=Gate("是教育,还是那几年")
p1=G.positive_control("三段内六对最小归一值都 >0.20",
                      planted=float(min(v["min_overall"] for v in res.values())),floor=0.20,spread=0.01)
xh=[v["p_xh"]<0.05 for v in res.values()]; sg=[np.sign(v["d_xh"]) for v in res.values()]
if p1:
    if all(xh) and len(set(sg))==1: verdict=f"**是教育:三段 xmarsex×homosex 的 Δ 同号且都超各自的零({[f'{v[chr(100)+chr(95)+chr(120)+chr(104)]:+.4f}' for v in res.values()]})**"
    elif sum(xh)==1: verdict=f"**是那一段的时代:只有 {[k for k,v in res.items() if v['p_xh']<0.05][0]} 超零 ⇒ `#675`–`#677` 那几句要加时间范围**"
    else: verdict=f"**部分段测不出:{sum(xh)}/3 段超零,同号 {len(set(sg))==1} —— 逐段如实登,这也是一个结果**"
else: verdict="UNVERIFIED —— 正对照失败"
print(f"\n{verdict}"); print(G)
json.dump(dict(res=res,verdict=verdict,unchallenged=True),open(OUT/"era_or_education.json","w"),indent=1,ensure_ascii=False)

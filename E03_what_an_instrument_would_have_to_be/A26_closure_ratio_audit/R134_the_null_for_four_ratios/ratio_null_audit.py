"""E03·A26·R134 —— 把这条线用过的比值,逐个放回它自己的零里

**类型:CLOSURE(诚实标注)。不推进对象。它保护的是页面上每一条建立在比值上的声明。**
`#691` 发现:我预注册的 `[0.5,2]` 是凭直觉写的,而 `ext` 比值的零的 95% 分位高达 **17.3**
⇒ **暴露的是判据形式的问题,不是那一个判据的问题。**

## 硬规则①(已跑),而第一把 grep 是坏仪器
账本 `#600` 起,正文出现「倍 / × / 比值」的条目有 **35 条**;页面上带这类表述的有 **192 处**。
⚠ **但其中绝大多数是 `N× 自身展布`(即 effect/SE,z 型统计量,零是良态的),不是本轮的目标。**
**目标类由规则定义,不是挑出来的:`|效应A| / |效应B|`,两边都是估计量** ——
这类量的零重尾,因为分母可以接近零。按此规则,页面上受影响的是 **6 条**。

## ⑤ 最强混淆(`#691` 预注册):有些旧轮的原始数据没有持久化
逐条查 `results/*.json`:

| 条目 | 轮次 | JSON | 处置 |
|---|---|---|---|
| `#650` 1.10 倍 | `R095` | **0 个** | **⛔ 无法回算,如实记,不许近似重建** |
| `#665` 2.14 / 1.18 | `R107` | 1 | **本轮未回算(SCCS 另一条管道)—— 记「未做」,不记「已计划」** |
| `#675` 4.19× / 0.35× | `R117` | 3 | ✅ 回算 |
| `#677` 4.8× / 10.6× | `R119` | 1 | ✅ 回算 |
| `#687` 2.436 | `R129` | 1 | ✅ 回算 |
| `#691` 1.967 / 2.834 | `R133` | 1 | ✅ **正对照(④):必须复现 p≈0.384 / 0.042** |

## G1 ESTIMAND
对每条:**用该轮自己的随机化(打乱定义两臂的那个分组)重建比值的零分布**,报经验 p 与零的 95% 分位。
## ⑧ 判据(`#691` 写死)
**经验 p ≥ 0.05 的,页面上那条声明必须降级或加「比值不可分辨」的注。**
## G2 CONTROL
**正对照(④)**:`#691` 的两个必须被这套流程复现(`ext` p≈0.38、`disp` p≈0.04),**否则流程本身错了,当场停。**
## IMPOSSIBLE(不写 planned)
`#650` 原始结果未持久化 ⇒ **本站点无法回算**;`#665` 走 SCCS 管道,**本轮未做**;
本轮只重建比值的零,**不重新检验那些声明本身**。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
rc=lambda a,b:float(np.corrcoef(pd.Series(np.asarray(a)).rank(),pd.Series(np.asarray(b)).rank())[0,1])
SEX=["premarsx","xmarsex","homosex","teensex"]; FEM=["fefam","fepol","fepresch"]
gss,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ","age"]+SEX+FEM, apply_value_formats=False, encoding="latin1")
J=gss.dropna(subset=SEX+FEM+["educ","age"])
grid=np.arange(8,19.01,1.0); BW=2.5
def wc(x,y,w):
    mx=np.average(x,weights=w);my=np.average(y,weights=w);cx=x-mx;cy=y-my
    s=np.sqrt(np.average(cx*cx,weights=w)*np.average(cy*cy,weights=w))
    return np.average(cx*cy,weights=w)/s if s>1e-12 else np.nan
RK={c:J[c].rank().to_numpy(float) for c in SEX+FEM}
def pair_delta(a,b,e):
    out=[]
    for g in grid:
        W=np.exp(-0.5*((e-g)/BW)**2)
        if W.sum()<200: out.append(np.nan); continue
        r=wc(RK[a],RK[b],W); idx=W>np.quantile(W,0.5)
        xs=np.sort(RK[a][idx]); ys=np.sort(RK[b][idx]); ys=ys if r>0 else ys[::-1]
        c=np.corrcoef(xs,ys)[0,1]
        out.append(r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan)
    y=np.array(out); m=np.isfinite(y)
    return float(np.mean(y[m][-3:])-np.mean(y[m][:3])) if m.sum()>=4 else np.nan
E=J["educ"].to_numpy(float); AG=J["age"].to_numpy(float)
WH=[p for p in combinations(SEX,2) if "homosex" in p]; NH=[p for p in combinations(SEX,2) if "homosex" not in p]
FP=list(combinations(FEM,2))
rng=np.random.default_rng(20260806); res={}
def ratio_null(fn,B=400):
    obs=fn(E,AG)
    nul=[]
    for _ in range(B):
        e=rng.permutation(E)
        v=fn(e,AG)
        if np.isfinite(v): nul.append(abs(v))
    nul=np.array(nul)
    return obs,float(np.quantile(nul,.95)),float((nul>=abs(obs)).mean())
f677a=lambda e,a: abs(np.median([pair_delta(x,y,e) for x,y in WH]))/abs(np.median([pair_delta(x,y,e) for x,y in NH]))
f677b=lambda e,a: abs(np.median([pair_delta(x,y,e) for x,y in NH]))/abs(np.median([pair_delta(x,y,e) for x,y in FP]))
for nm,f,claim in [("#677 4.8× 含/不含 homosex",f677a,4.8),("#677 10.6× 不含 vs 安慰剂",f677b,10.6)]:
    o,q,p=ratio_null(f,250); res[nm]=dict(obs=o,q95=q,p=p,claimed=claim)
    print(f"{nm:28s} 重算 **{o:.3f}**(页面 {claim})· 零 95% **{q:.3f}** · p = **{p:.4f}** "
          f"{'✅ 超零' if p<0.05 else '⛔ 不可分辨 ⇒ 页面须加注'}")
def f675(e,a):
    de=abs(np.median([pair_delta(x,y,e) for x,y in combinations(SEX,2)]))
    df=abs(np.median([pair_delta(x,y,e) for x,y in FP]))
    return de/df if df>1e-9 else np.nan
o,q,p=ratio_null(f675,250); res["#675 4.19× 教育两组"]=dict(obs=o,q95=q,p=p,claimed=4.19)
print(f"{'#675 4.19× 教育两组':28s} 重算 **{o:.3f}**(页面 4.19)· 零 95% **{q:.3f}** · p = **{p:.4f}** "
      f"{'✅ 超零' if p<0.05 else '⛔ 不可分辨 ⇒ 页面须加注'}")
json.dump(res,open(OUT/"ratio_null_audit.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'ratio_null_audit.json'}")

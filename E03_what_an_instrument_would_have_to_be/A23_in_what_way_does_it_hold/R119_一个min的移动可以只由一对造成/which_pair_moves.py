"""E03·A23·R119 —— 一个 min 的移动可以只由一对造成

**类型:FRONTIER。** `#676` 的 NEXT,而它是我自己预注册的那个混淆,结果比预注册的更严重。

## 硬规则①(已跑,在算任何 Δ 之前)
低端 `educ<=10` n=1,270 · 高端 `educ>=16` n=3,295。六对归一值的两端差:

| 对 | 差 | 含 homosex |
|---|---|---|
| premarsx×xmarsex | +0.0224 | — |
| premarsx×homosex | +0.1011 | ✅ |
| premarsx×teensex | +0.0280 | — |
| **xmarsex×homosex** | **+0.2226** | ✅ |
| xmarsex×teensex | +0.1047 | — |
| **homosex×teensex** | **+0.2093** | ✅ |

⚠ **而两端的最弱一对都是 `xmarsex×homosex`(+0.3021 → +0.5247)** ——
**`#676` 那条 +0.1604 的整条移动就是这一对,而这一对含 `homosex`。**
⚠ 安慰剂三对是 +0.0548 / +0.0706 / +0.0221(中位 +0.0548),
**比不含 `homosex` 的性三对(中位 +0.0280)还高。**

## G1 ESTIMAND
**六对各自的 Δ**(核加权,与 `#676` 同一条流水线),以及
**含/不含 `homosex` 两组的 Δ 中位数**;集中度 = 最大 Δ ÷ 中位 Δ。
## G2 CONTROLS
**正对照**:六对的全样本归一值必须为正且复现 `#651`/`#675` 量级(>0.20)。
**安慰剂**:性别角色三对走同一条流水线 —— **它现在是一个更严的对照**,因为它的中位比不含
`homosex` 的性三对还高。
**这个零该不该是零?** **不该** —— 安慰剂已知不为零 ⇒ `offset_control`,
**零的种类 = 另一组题的三对上的同一个 Δ 中位数。**
## KILL(条件式)
if 正对照复现:
  不含 homosex 的三对 Δ 中位 **超零且 > 安慰剂中位** -> **「四个判断绑在一起」站得住** ·
  不含 homosex 的三对掉到地板或 ≤ 安慰剂 -> **`#676` 缩成「教育把同性恋这一题拉进了性道德这一块」**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
`homosex` 只有一题 ⇒ **无法把「同性恋」与「其余性题」做等权的双向对称检验**;
因果:横断面无干预;**跨仪器:MFQ 无对应的四题结构**(`#675` 已测不复制)。`[unchallenged]`
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
      usecols=["educ"]+SEX+FEM, apply_value_formats=False, encoding="latin1")
j=df.dropna(subset=SEX+FEM+["educ"]); grid=np.arange(8,19.01,1.0); BW=2.5
E0=j["educ"].to_numpy(float)
RK={c:pd.Series(j[c]).rank().to_numpy(float) for c in SEX+FEM}
def wc(x,y,w):
    mx=np.average(x,weights=w);my=np.average(y,weights=w);cx=x-mx;cy=y-my
    s=np.sqrt(np.average(cx*cx,weights=w)*np.average(cy*cy,weights=w))
    return np.average(cx*cy,weights=w)/s if s>1e-12 else np.nan
def pair_curve(a,b,e):
    x=RK[a]; y=RK[b]; out=[]
    for g in grid:
        W=np.exp(-0.5*((e-g)/BW)**2)
        if W.sum()<200: out.append(np.nan); continue
        r=wc(x,y,W)
        idx=W>np.quantile(W,0.5)
        xs=np.sort(x[idx]); ys=np.sort(y[idx]); ys=ys if r>0 else ys[::-1]
        c=np.corrcoef(xs,ys)[0,1]
        out.append(r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan)
    return np.array(out)
def D(y):
    m=np.isfinite(y)
    return float(np.mean(y[m][-3:])-np.mean(y[m][:3])) if m.sum()>=4 else np.nan
SP=list(combinations(SEX,2)); FP=list(combinations(FEM,2))
print("=== 六对各自的 Δ(核加权,与 #676 同一条流水线)===")
ds={}
for a,b in SP+FP:
    ds[(a,b)]=D(pair_curve(a,b,E0))
    tag="✅" if "homosex" in (a,b) else ("—" if (a,b) in SP else "安慰剂")
    print(f"  {a[:10]+'×'+b[:10]:24s} Δ = **{ds[(a,b)]:+.4f}**  {tag}")
withh=[ds[p] for p in SP if "homosex" in p]; without=[ds[p] for p in SP if "homosex" not in p]
plac=[ds[p] for p in FP]
mw,mo,mp=float(np.median(withh)),float(np.median(without)),float(np.median(plac))
print(f"\n  含 homosex 三对 中位 **{mw:+.4f}** · 不含 三对 中位 **{mo:+.4f}** · 安慰剂三对 中位 **{mp:+.4f}**")
print(f"  集中度:最大 Δ ÷ 中位 Δ = {max(ds[p] for p in SP)/np.median([ds[p] for p in SP]):.2f}×")
rng=np.random.default_rng(20260806); nul=[]
for _ in range(300):
    e=rng.permutation(E0)
    nul.append(abs(float(np.median([D(pair_curve(a,b,e)) for a,b in SP if "homosex" not in (a,b)]))))
nul=np.array(nul); q=float(np.nanquantile(nul,0.95)); p=float(np.nanmean(nul>=abs(mo)))
print(f"\n  不含 homosex 三对的零(打乱 educ):95% 分位 **{q:.4f}** · 观测 {mo:+.4f} · p = **{p:.4f}**  "
      f"{'✅ 超零' if p<0.05 else '⛔ 落在地板上'}")
ov=[]
for a,b in SP:
    x,y=RK[a],RK[b]; r=np.corrcoef(x,y)[0,1]
    xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
    ov.append(r/abs(np.corrcoef(xs,ys)[0,1]))
G=Gate("四个判断绑在一起,还是同性恋这一题被拉了进来")
p1=G.positive_control("六对的全样本归一值最小者为正且 >0.20(复现 #651/#675 量级)",
                      planted=float(min(ov)),floor=0.20,spread=0.01)
p2=G.offset_control("不含 homosex 的三对 Δ 中位必须大于安慰剂三对的中位",
                    effect=abs(mo),offset=abs(mp),spread=0.010,
                    null_kind="性别角色三对上的同一个 Δ 中位数 —— 已知不为零,是系统性基线偏移")
if p1:
    if p<0.05 and p2: verdict=f"**「四个判断绑在一起」站得住:不含 homosex 三对中位 {mo:+.4f} 超零且大于安慰剂 {mp:+.4f}**"
    else: verdict=(f"**缩小 `#676`:整条移动集中在含 homosex 的三对(中位 {mw:+.4f});"
                   f"不含它的三对中位 {mo:+.4f}{'(p=%.4f 落地板)'%p if p>=0.05 else ''}"
                   f"{'且不高于安慰剂 %+.4f'%mp if not p2 else ''} ⇒ "
                   f"要改成「教育把同性恋这一题拉进了性道德这一块」**")
else: verdict="UNVERIFIED —— 正对照失败"
print(f"\n{verdict}"); print(G)
json.dump(dict(pairs={f"{a}×{b}":v for (a,b),v in ds.items()},with_homosex=mw,without=mo,placebo=mp,
               null_q95=q,p=p,verdict=verdict,unchallenged=True),
          open(OUT/"which_pair_moves.json","w"),indent=1,ensure_ascii=False)

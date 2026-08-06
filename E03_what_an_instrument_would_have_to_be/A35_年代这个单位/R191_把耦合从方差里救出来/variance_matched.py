"""E03·A35·R191 —— 把耦合从方差里救出来,还是救不出来

**类型:FRONTIER。这是 `#747`①,而它可能把 `#745` 从「降级」推到「撤回」。**

**心理学的那一句(本轮要判的):`#745` 说只有常去教堂的人把两个判断连了起来;
`#747` 发现那一层恰好也是方差升得最多的一层,于是把它降了级。
本轮问:那条耦合梯度,是不是**只是**各层边际轨迹不同造出来的?**

## ⚠ 字面写法有陷阱,已改
`#747`① 写的是「重采样使各层各年方差相等」。**照字面做是危险的**:
**重采样改的是 `homosex` 自己的分布,而它正是被测相关的两个变量之一** ——
**`#746` 刚栽过的「对结果取条件」(Oldham 1962)。**
⇒ 改用 **`#739` 已验证过的那具装置**:**ρ 全程恒定,只让该层自己的真实逐年边际走切点**,
直接问「**边际轨迹本身能造出多少斜率**」。**这具装置在 `#739` 里通过了正对照(切点固定 ⇒ −0.0020),
而且它在该开火的地方开过火(三个小对的伪影份额 60–86%)。**

## ⚠ 按 `#746`① 先算 MDE(装置自己的种子展布),再解释份额
**跑之前先量:同一格换 5 个种子,模拟斜率抖多少。抖动大过份额的差,份额就没有分辨力。**

## 硬规则①
三层可用年 18 / 21 / 25(n≥150);`homosex` 的 sd 轨迹 1.321→1.131 · 1.199→1.237 · **0.862→1.373**。
⚠ **口径必须与 `#745` 一致:耦合是天花板归一的,所以模拟也算归一斜率**(`#721` 的教训)。

## G1 ESTIMAND
每层的 **`artefact_share` = 模拟(ρ 恒定、真实边际走)的归一斜率 ÷ 实测归一斜率**。
## W1–W4(双边)
| 世界 | 判据 | 读法 |
|---|---|---|
| **W1 耦合救回来** | 三层份额都 ≤ 0.35,**且模拟不复现那个梯度** | `#745` 的梯度恢复 |
| **W2 就是边际** | **常去的份额 ≥ 0.65,或模拟复现了梯度的次序** | **`#745` 从降级推到撤回** |
| **W3 之间** | 报份额,不报判决 | |
| **W4 判不了** | 模拟的种子展布 ≥ 份额之差 | 这具装置在这个 n 上分不开 |

⚠ **W2 的正结果会让我撤回两轮前写上页的一条结论 —— 这正是本轮设计成能出它的理由。**

## G2 CONTROLS
**④ 正对照**(两条,都能失败):
(a) **切点固定不动时,模拟的 36 年变化必须 ≈ 0**(阈 |0.02|,`#739` 实测 −0.0020);
(b) **实测的三层归一斜率必须复现 `#745`**(+0.00287 / +0.00178 / −0.00044,容差 0.0005)。
**零** = `negative_control`,**零的种类 = ρ 设为该层早期实测值并全程恒定 ——
保住该层的样本量、档数与两个题各自的真实边际轨迹,只把「潜相关随年代变化」抽掉。**
## G3:三层 × {实测, 模拟} 六格全报。G4:ρ 取早期 / 全期 / 末期三种设定。
## ⑤ 停止条件(**双边**,跑之前写死)
- **正对照 (a) 或 (b) 不过 ⇒ UNVERIFIED 并停。**
- **模拟的种子展布 ≥ 三层份额的极差 ⇒ 判 W4,并停。**
- 之后依 **W2 → W1 → W3** 判。
## IMPOSSIBLE(不写 planned)
模拟是**双变量正态 + 固定切点**;真实作答未必如此 ⇒ **它给的是「边际能解释多少」,不是全部混淆的清单**;
仍是**重复横断面**;**换不了仪器**。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr, norm
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
A,Bc="premarsx","homosex"
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year","attend",A,Bc],encoding="latin1")
R=g.dropna(subset=[A,Bc,"attend"]).copy()
R["rel"]=pd.cut(R.attend,[-1,1,4,8],labels=["几乎不去","偶尔","常去"]); R=R.dropna(subset=["rel"]).copy()
LV=["几乎不去","偶尔","常去"]; FN,FY=150,5
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def cuts(sub,c):
    p=sub[c].value_counts(normalize=True).sort_index()
    cum=np.cumsum([p.get(v,0.0) for v in sorted(p.index)])[:-1]
    return [float(norm.ppf(min(max(x,1e-4),1-1e-4))) for x in cum]
def obs_nslope(fr):
    pts=[]
    for y,sub in fr.groupby("year"):
        if len(sub)<FN or sub[A].nunique()<2 or sub[Bc].nunique()<2: continue
        x=sub[A].to_numpy(float); yv=sub[Bc].to_numpy(float); r=sp(x,yv)
        xs=np.sort(x); ys=np.sort(yv); ys=ys if r>0 else ys[::-1]; c=abs(sp(xs,ys))
        if c<1e-9: continue
        pts.append((float(y),r/c))
    return (float(np.polyfit([p[0] for p in pts],[p[1] for p in pts],1)[0]) if len(pts)>=FY else np.nan), [p[0] for p in pts]
def sim_nslope(fr,yrs,rho,rng,move=True,N=4000,REP=40):
    y0=yrs[0]; ca0=cuts(fr[fr.year==y0],A); cb0=cuts(fr[fr.year==y0],Bc)
    out=[]
    for y in yrs:
        sub=fr[fr.year==y]
        ca=cuts(sub,A) if move else ca0; cb=cuts(sub,Bc) if move else cb0
        vs=[]
        for _ in range(REP):
            L=rng.multivariate_normal([0,0],[[1,rho],[rho,1]],size=N)
            a=np.digitize(L[:,0],ca); b=np.digitize(L[:,1],cb)
            r=sp(a,b); xs=np.sort(a.astype(float)); ys=np.sort(b.astype(float)); ys=ys if r>0 else ys[::-1]
            c=abs(sp(xs,ys))
            if c>1e-9: vs.append(r/c)
        if vs: out.append((float(y),float(np.median(vs))))
    return float(np.polyfit([p[0] for p in out],[p[1] for p in out],1)[0])
LED={"几乎不去":-0.00044,"偶尔":0.00178,"常去":0.00287}
print("=== ④ 正对照 (b):实测三层归一斜率必须复现 `#745`(容差 0.0005)===")
obs={}; YR={}
okb=True
for lv in LV:
    obs[lv],YR[lv]=obs_nslope(R[R.rel==lv]); d=abs(obs[lv]-LED[lv]); okb&= d<=0.0005
    print(f"  {lv:8s} 实测 {obs[lv]:+.5f} · 账本 {LED[lv]:+.5f} · 差 {d:.5f} {'✅' if d<=0.0005 else '⛔'}")
if not okb: print("\n⛔ ⑤ 触发,停"); sys.exit(0)
rng=np.random.default_rng(20260806)
early={lv:float(np.median([sp(*(lambda s:(s[A],s[Bc]))(R[(R.rel==lv)&(R.year==y)][[A,Bc]].dropna()))
        for y in YR[lv][:4]])) for lv in LV}
print("\n=== ④ 正对照 (a):切点固定不动时,模拟 36 年变化必须 ≈ 0(阈 |0.02|)===")
fix=sim_nslope(R[R.rel=="常去"],YR["常去"],early["常去"],rng,move=False)*36
oka=abs(fix)<0.02
print(f"  常去,切点固定:36 年 **{fix:+.4f}** {'✅' if oka else '⛔ ⑤ 触发'}")
if not oka: print("⛔ 停"); sys.exit(0)
print("\n=== ⚠ 按 `#746`① 先算这具装置自己的 MDE(5 个种子的展布)===")
sp5=[sim_nslope(R[R.rel=="常去"],YR["常去"],early["常去"],np.random.default_rng(500+s)) for s in range(5)]
spread=max(sp5)-min(sp5)
print(f"  常去,5 个种子的模拟斜率:"+" ".join(f"{x:+.5f}" for x in sp5)+f"  ⇒ **展布 {spread:.5f}**")
print("\n=== G3:三层 × {实测, 模拟} ===")
print(f"{'层':10s}{'早期 ρ':>9s}{'实测归一斜率':>14s}{'模拟(ρ恒定)':>14s}{'伪影份额':>10s}")
res={}
for lv in LV:
    s=sim_nslope(R[R.rel==lv],YR[lv],early[lv],rng)
    share=s/obs[lv] if abs(obs[lv])>1e-6 else np.nan
    res[lv]=dict(early_rho=early[lv],obs=obs[lv],sim=s,share=float(share))
    print(f"{lv:10s}{early[lv]:>9.3f}{obs[lv]:>+14.5f}{s:>+14.5f}{share:>10.2f}")
shares=[res[lv]["share"] for lv in LV if np.isfinite(res[lv]["share"])]
rng_sh=max(shares)-min(shares)
print(f"\n  份额极差 {rng_sh:.2f} vs 装置展布(换算成份额)约 {spread/abs(obs['常去']):.2f}")
print("\n=== G4:ρ 的三种设定(常去)===")
for lab,rho in (("早期",early["常去"]),("全期",0.5*(early["常去"]+0.6)),("末期",0.63)):
    s=sim_nslope(R[R.rel=="常去"],YR["常去"],rho,rng)
    print(f"  ρ={rho:.3f}({lab}) ⇒ 模拟 {s:+.5f} · 份额 {s/obs['常去']:.2f}")
G=Gate("把耦合从方差里救出来还是救不出来")
p1=G.positive_control("切点固定时模拟≈0,且实测复现 #745",planted=1.0 if (oka and okb) else 0.0,floor=0.0,spread=0.1)
p2=G.negative_control("ρ 恒定的模拟不该复现实测的耦合上升",null=abs(res["常去"]["sim"]),effect=abs(obs["常去"]),
    null_spread=float(spread),
    null_kind="ρ 设为该层早期实测值并全程恒定 —— 保住该层样本量、档数与两个题各自的真实边际轨迹,只把「潜相关随年代变化」抽掉")
sim_order=[res[lv]["sim"] for lv in LV]; obs_order=[obs[lv] for lv in LV]
reproduces = (np.argmax(sim_order)==np.argmax(obs_order)) and (sim_order[2]>sim_order[0])
if not (oka and okb): v="**UNVERIFIED:正对照没过**"
elif spread>=rng_sh: v=f"**W4:装置自己的种子展布 {spread:.5f} ≥ 份额极差 ⇒ 这具装置在这个 n 上分不开,判不了**"
elif res["常去"]["share"]>=0.65 or reproduces:
    v=(f"**W2:常去的伪影份额 {res['常去']['share']:.2f}"
       f"{',且模拟复现了梯度的次序' if reproduces else ''} ⇒ `#745` 从降级推到撤回**")
elif all(s<=0.35 for s in shares) and not reproduces: v=f"**W1:三层份额都 ≤0.35 且模拟不复现梯度 ⇒ `#745` 的梯度恢复**"
else: v=f"**W3:份额 "+" / ".join(f"{res[lv]['share']:.2f}" for lv in LV)+" ⇒ 报份额,不报判决**"
print(f"\n{v}"); print(G)
json.dump(dict(strata=res,fixed_cut_control=fix,seed_spread=float(spread),
  reproduces_order=bool(reproduces),verdict=v,unchallenged=True),open(OUT/"vm.json","w"),indent=1,ensure_ascii=False)

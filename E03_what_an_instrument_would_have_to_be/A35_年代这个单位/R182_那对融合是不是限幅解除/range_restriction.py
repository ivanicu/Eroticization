"""E03·A35·R182 —— 那对「融合」是不是限幅解除

**类型:FRONTIER。这是 `#738`① —— 而它的最强混淆在 `#738` 落地时就写下了。
同时兑现 `#737`①:把 GSS 的真实边际搬进合成世界当切点。**

**心理学的那一句(本轮要判的):`premarsx × homosex` 从 0.288 走到 0.632。
但同一段时间里 `homosex` 从「几乎所有人都说 always wrong」走到中段 ——
一个题从地板挪开,它跟别的题的相关**本来就会**上升。那 0.27 里有多少是这件事?**

## 硬规则①(逐年真实边际,先印后用)
两个题的**逐年四档分布**直接从 GSS 读出来,**不是我造的切点** ——
这既是本轮的输入,也是 `#737`① 要的那件事。

## G1 ESTIMAND
**在潜相关 ρ 恒定不变的合成世界里**,只让两个题的**切点按 GSS 的真实逐年边际移动**,
测得到的 Spearman 会涨多少 ⇒ **`artefact_share = 合成的 36 年上涨 ÷ 实测的 +0.2703`。**
## W1 / W2 / W3(三分,双边)
| 世界 | 合成上涨 | 读法 |
|---|---|---|
| **W1 基本是限幅** | **≥ +0.20**(≥74%) | **「融合」撤回** —— 那是 `homosex` 挪离地板的机械后果 |
| **W2 基本是真的** | **≤ +0.08**(≤30%) | 融合站得住,并报出这部分伪影份额 |
| **W3 一半一半** | 0.08–0.20 | **报份额,不报判决** |
## G2 CONTROLS
**④ 正对照**:**切点固定不动**时,合成的斜率必须 ≈ 0(|36 年变化| < 0.02)——
**否则合成本身在造趋势,后面的份额都不可信。**
⚠ **且它能失败**:若固定切点也造出上涨,本轮直接停。
**零** = `negative_control`,**零的种类 = ρ 设为该对的 1988 年实测值并全程恒定 ——
保住样本量、档数、以及两个题各自的真实边际轨迹,只把「潜相关随年代变化」这件事抽掉。**
**SHAM**:对 `premarsx × xmarsex`(实测 −0.0388,几乎没动)跑同一流程 ——
**若合成也给它造出上涨,说明这套模拟对任何一对都造趋势。**
## G3:六对全跑。G4:ρ 取 1988 值 / 全期均值 / 2024 值 三种设定各跑一次。
## ⑤ 停止条件(**双边**,跑之前写死)
- **④ 固定切点的合成斜率 |36 年变化| ≥ 0.02 ⇒ UNVERIFIED 并停。**
- **合成上涨 ≥ +0.20 ⇒ W1,页上「融合」那一句撤回。**
- **≤ +0.08 ⇒ W2,融合成立并报份额。**
- **落在中间 ⇒ W3,报份额不报判决。**
- ⚠ **三档都写了,没有哑口。**
## IMPOSSIBLE(不写 planned)
合成用的是**双变量正态 + 固定切点**;真实的作答过程未必如此
⇒ **它给的是「限幅能解释多少」的一个估计,不是全部混淆的清单。**
**换不了仪器**:本轮问的是 GSS 这对题自己的历史。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr, norm
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
SEX=["premarsx","xmarsex","homosex","teensex"]
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+SEX,encoding="latin1")
J=g.dropna(subset=SEX); yrs=sorted(J.year.unique()); Y=np.array(yrs,float)
def cuts_from(fr,c):
    """把某年某题的真实四档比例翻成潜正态上的切点(`#737`① 要的那件事)。"""
    p=fr[c].value_counts(normalize=True).sort_index()
    cum=np.cumsum([p.get(v,0.0) for v in sorted(p.index)])[:-1]
    return [float(norm.ppf(min(max(x,1e-4),1-1e-4))) for x in cum]
print("=== 硬规则①:逐年真实边际(四档比例)· 只印首尾两年 ===")
for c in ("premarsx","homosex"):
    for y in (yrs[0],yrs[-1]):
        p=J[J.year==y][c].value_counts(normalize=True).sort_index()
        print(f"  {c:10s} {int(y)}: "+" ".join(f"{v:.0f}档 {p.get(v,0):.3f}" for v in sorted(p.index))
              +f"   切点 {[round(x,3) for x in cuts_from(J[J.year==y],c)]}")
obs={}
for a,b in itertools.combinations(SEX,2):
    v=[sp(*(lambda m:(m[a],m[b]))(J[J.year==y][[a,b]].dropna())) for y in yrs]
    obs[(a,b)]=dict(first=v[0],last=v[-1],slope36=float(np.polyfit(Y,v,1)[0]*36))
rng=np.random.default_rng(20260806); N=4000; REP=60
def simulate(a,b,rho,move_cuts=True):
    """ρ 恒定;切点按真实逐年边际移动(或固定在 1988)。"""
    ca0=cuts_from(J[J.year==yrs[0]],a); cb0=cuts_from(J[J.year==yrs[0]],b)
    out=[]
    for y in yrs:
        ca=cuts_from(J[J.year==y],a) if move_cuts else ca0
        cb=cuts_from(J[J.year==y],b) if move_cuts else cb0
        vs=[]
        for _ in range(REP):
            L=rng.multivariate_normal([0,0],[[1,rho],[rho,1]],size=N)
            vs.append(sp(np.digitize(L[:,0],ca),np.digitize(L[:,1],cb)))
        out.append(float(np.median(vs)))
    return float(np.polyfit(Y,out,1)[0]*36), out
print("\n=== ④ 正对照:切点固定不动时,合成斜率必须 ≈ 0(|36年| < 0.02)===")
s_fix,_=simulate("premarsx","homosex",obs[("premarsx","homosex")]["first"],move_cuts=False)
pc=abs(s_fix)<0.02
print(f"  premarsx × homosex,切点固定:36 年变化 **{s_fix:+.4f}** {'✅' if pc else '⛔ ⑤ 触发'}")
if not pc:
    print("⛔ 停:合成本身在造趋势"); sys.exit(0)
print("\n=== G3 六对全跑:ρ 恒定 = 1988 实测值,切点按真实边际移动 ===")
print(f"{'对':26s}{'实测36年':>10s}{'合成36年':>10s}{'伪影份额':>10s}")
res={}
for (a,b),d in obs.items():
    s_sim,_=simulate(a,b,d["first"],move_cuts=True)
    share=s_sim/d["slope36"] if abs(d["slope36"])>1e-6 else np.nan
    res[f"{a}×{b}"]=dict(obs36=d["slope36"],sim36=s_sim,share=float(share))
    tag=" ★" if (a,b)==("premarsx","homosex") else (" ← SHAM" if (a,b)==("premarsx","xmarsex") else "")
    print(f"{a+'×'+b:26s}{d['slope36']:>+10.4f}{s_sim:>+10.4f}{share:>10.2f}{tag}")
key=res["premarsx×homosex"]
print("\n=== G4:ρ 的三种设定 ===")
for lab,rho in (("1988 实测",obs[("premarsx","homosex")]["first"]),
                ("全期均值",float(np.mean([obs[("premarsx","homosex")]["first"],obs[("premarsx","homosex")]["last"]]))),
                ("2024 实测",obs[("premarsx","homosex")]["last"])):
    s,_=simulate("premarsx","homosex",rho,move_cuts=True)
    print(f"  ρ={rho:.3f}({lab:6s}) ⇒ 合成 36 年 **{s:+.4f}** · 份额 {s/key['obs36']:.2f}")
G=Gate("那对融合是不是限幅解除")
p1=G.positive_control("切点固定时合成斜率必须≈0",planted=float(0.02-abs(s_fix)),floor=0.0,spread=0.001)
p2=G.negative_control("ρ 恒定的合成不该复现出实测的上涨",null=abs(key["sim36"]),effect=abs(key["obs36"]),
    null_spread=0.005,null_kind="ρ 设为该对 1988 年实测值并全程恒定 —— 保住样本量、档数与两个题各自的真实边际轨迹,只把「潜相关随年代变化」抽掉")
sim=key["sim36"]
if not p1: v="**UNVERIFIED:合成本身在造趋势**"
elif sim>=0.20: v=f"**W1:限幅解除就能造出 {sim:+.4f}(占实测 {key['share']:.0%})⇒ 「融合」撤回**"
elif sim<=0.08: v=f"**W2:限幅只造出 {sim:+.4f}(占实测 {key['share']:.0%})⇒ 融合站得住,伪影份额已报**"
else: v=f"**W3:限幅造出 {sim:+.4f},占实测 {key['share']:.0%} ⇒ 报份额,不报判决**"
print(f"\n{v}"); print(G)
json.dump(dict(observed={f"{a}×{b}":obs[(a,b)] for a,b in obs},simulated=res,
  fixed_cut_control=s_fix,verdict=v,unchallenged=True),open(OUT/"rr.json","w"),indent=1,ensure_ascii=False)

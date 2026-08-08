"""E03·A35·R185 —— 「出轨没加入」是关于人的,还是它根本没有空间

**类型:FRONTIER。这是 `#741`① —— 而在谈「自愿 vs 承诺」之前,先排除最便宜的那个解释。**

**心理学的那一句(本轮要判的):同性恋并进来的时候,婚前性和未成年性都跟着动了,
唯独出轨没有。这是因为美国人**不肯**把出轨和它们放在一起,
还是因为出轨这一题几乎所有人都还在说「永远是错的」,**没有地方可动**?**

## ⚠ 假设与它的反例都在 `#741`① 里写下了
假设「加入的是**自愿且不伤害第三方**的那一族」——
**而 `teensex`(14–16 岁)显然不是「自愿」的干净例子,它却加入了(+0.00381)**
⇒ **这个假设很可能是错的,所以本轮不测它,先测那个更便宜的对手。**

## 硬规则①(已跑)
「always wrong」的占比 1988 → 2024:
`premarsx` 0.263 → 0.169 · **`xmarsex` 0.793 → 0.687** · `homosex` 0.768 → **0.349** · `teensex` 0.694 → 0.484。
⇒ **`xmarsex` 到 2024 年仍有 68.7% 说「永远是错的」。**

## G1 ESTIMAND
三对 `homosex × {premarsx, xmarsex, teensex}` 的**世代内合并斜率**,
**生相关**(与 `#741` 同一条路径)与**天花板归一**各一次;
外加每对的**逐年天花板轨迹**(有没有空间)。
## W_a / W_b / W_c(双边)
| 世界 | 判据 | 读法 |
|---|---|---|
| **W_a 是关于人的** | 归一后 `xmarsex×homosex` 的斜率 ≤ **0.25 ×** `premarsx×homosex` 的 | 出轨确实没被并进来 |
| **W_b 是没有空间** | 归一后 ≥ **0.60 ×** | **「出轨没加入」这句话撤回** —— 那是天花板说的,不是人说的 |
| **W_c 之间** | 0.25–0.60 | 报份额,不报判决 |

⚠ **W_b 的正结果会削掉我上一轮刚写上页的那半句 —— 这正是本轮设计成能出它的理由。**

## G2 CONTROLS
**④ 正对照**:三对的**生**斜率必须复现 `#741`(+0.00672 / +0.00075 / +0.00381,容差 0.0005)。
**零** = `negative_control`,**零的种类 = 在每个世代内部打乱受访者的年份标签 ——
保住世代构成、每年 n 与作答分布,只毁掉「同一世代里谁属于哪一年」**(与 `#740`/`#741` 同一具零),
**对归一后的 `xmarsex×homosex` 跑,B=200。**
## G3:三对 × {生, 归一} = 6 格全报。G4:逐年天花板轨迹全报。
## ⑤ 停止条件(**双边**,跑之前写死)
- **三对的生斜率复现不到 0.0005 ⇒ UNVERIFIED 并停。**
- 依 W_b → W_a → W_c 的顺序判(**先看能不能推翻我自己刚写的那句**)。
## IMPOSSIBLE(不写 planned)
天花板归一**修的是边际不对称,不是「这一题还没开始变」** ——
⚠ **若 `xmarsex` 的分布几乎没动,归一也救不回来,那时「没有空间」与「不肯并」在这具仪器上不可分辨,
必须如实记「判不了」,不许挑一个说。** 仍是**重复横断面**;**换不了仪器**。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year","cohort"]+SEX,encoding="latin1")
g=g.dropna(subset=["cohort"]).copy()
g["gen"]=pd.cut(g.cohort,[1880,1928,1946,1965,1981,1997,2010],
                labels=["前1929","1929–45","婴儿潮46–64","X 65–80","千禧81–96","Z 97+"])
GENS=["1929–45","婴儿潮46–64","X 65–80","千禧81–96"]; FN,FY=150,5
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def cells(a,b,frame,norm=False):
    num=[];wt=[];ceils=[]
    for gn in GENS:
        fr=frame[frame.gen==gn][["year",a,b]].dropna(); pts=[]
        for y,sub in fr.groupby("year"):
            if len(sub)<FN or sub[a].nunique()<2 or sub[b].nunique()<2: continue
            x=sub[a].to_numpy(float); yv=sub[b].to_numpy(float); r=sp(x,yv)
            xs=np.sort(x); ys=np.sort(yv); ys=ys if r>0 else ys[::-1]; c=abs(sp(xs,ys))
            if c<1e-9: continue
            ceils.append((float(y),c))
            pts.append((float(y),(r/c) if norm else r,len(sub)))
        if len(pts)>=FY:
            num.append(float(np.polyfit([p[0] for p in pts],[p[1] for p in pts],1)[0])); wt.append(sum(p[2] for p in pts))
    return (float(np.average(num,weights=wt)) if num else np.nan), ceils
PAIRS=[("premarsx","homosex"),("xmarsex","homosex"),("teensex","homosex")]
LED={("premarsx","homosex"):0.00672,("xmarsex","homosex"):0.00075,("teensex","homosex"):0.00381}
print("=== ④ 正对照:三对的**生**斜率必须复现 `#741`(容差 0.0005)===")
raw={}; ce={}
for p in PAIRS:
    raw[p],ce[p]=cells(*p,g,norm=False)
    print(f"  {p[0]:10s}×{p[1]:10s} 实测 **{raw[p]:+.5f}** · 账本 {LED[p]:+.5f} · 差 {abs(raw[p]-LED[p]):.5f} "
          f"{'✅' if abs(raw[p]-LED[p])<=0.0005 else '⛔'}")
maxd=max(abs(raw[p]-LED[p]) for p in PAIRS)
if maxd>0.0005:
    print("\n⛔ ⑤ 触发,停"); json.dump(dict(stop="旧值不可复现",raw={f"{a}×{b}":raw[(a,b)] for a,b in raw}),
        open(OUT/"x.json","w"),indent=1,ensure_ascii=False); sys.exit(0)
print("\n=== G3:三对 × {生, 归一} 六格全报 ===")
nrm={}
print(f"{'对':24s}{'生斜率/年':>12s}{'归一斜率/年':>13s}{'天花板 首→末':>18s}")
for p in PAIRS:
    nrm[p],_=cells(*p,g,norm=True)
    cs=sorted(ce[p]); c0=np.median([c for y,c in cs if y<=1995]); c1=np.median([c for y,c in cs if y>=2015])
    print(f"{p[0]+'×'+p[1]:24s}{raw[p]:>+12.5f}{nrm[p]:>+13.5f}{f'{c0:.3f} → {c1:.3f}':>18s}")
key=nrm[("premarsx","homosex")]; xm=nrm[("xmarsex","homosex")]
share=xm/key if abs(key)>1e-9 else np.nan
print(f"\n  归一后 `xmarsex×homosex` 占 `premarsx×homosex` 的 **{share:.0%}**")
rng=np.random.default_rng(20260806); nul=[]
for _ in range(200):
    P=g.copy(); P["year"]=P.groupby("gen",observed=True)["year"].transform(lambda s: rng.permutation(s.to_numpy()))
    v,_=cells("xmarsex","homosex",P,norm=True)
    if np.isfinite(v): nul.append(v)
q=np.quantile(nul,[0.025,0.975])
print(f"\n=== 零(世代内打乱年份,B={len(nul)},对归一后的 `xmarsex×homosex`)===")
print(f"  95% 区间 [{q[0]:+.5f}, {q[1]:+.5f}] · 实测 {xm:+.5f} ⇒ "
      f"{'在零之外' if not (q[0]<=xm<=q[1]) else '**落在零里 —— 它确实没动**'}")
G=Gate("出轨没加入是关于人还是没有空间")
p1=G.positive_control("三对的生斜率必须复现 #741(容差 0.0005)",planted=float(0.0005-maxd),floor=0.0,spread=0.00002)
p2=G.negative_control("世代内打乱年份后归一的 `xmarsex×homosex` 应回到零",
    null=float(max(abs(q[0]),abs(q[1]))),effect=abs(key),null_spread=0.00002,
    null_kind="在每个世代内部打乱受访者的年份标签 —— 保住世代构成、每年 n 与作答分布,只毁掉「同一世代里谁属于哪一年」")
if not p1: v="**UNVERIFIED:旧值不可复现**"
elif share>=0.60: v=f"**W_b:归一后 `xmarsex×homosex` 占 {share:.0%} ⇒ 「出轨没加入」撤回 —— 那是天花板说的**"
elif share<=0.25: v=f"**W_a:归一后仍只占 {share:.0%} ⇒ 出轨确实没被并进来,天花板不解释它**"
else: v=f"**W_c:归一后占 {share:.0%} ⇒ 报份额,不报判决**"
print(f"\n{v}"); print(G)
json.dump(dict(raw={f"{a}×{b}":raw[(a,b)] for a,b in raw},norm={f"{a}×{b}":nrm[(a,b)] for a,b in nrm},
  share=float(share),null_ci=[float(q[0]),float(q[1])],verdict=v,unchallenged=True),
  open(OUT/"x.json","w"),indent=1,ensure_ascii=False)

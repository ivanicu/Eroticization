"""E03·A35·R190 —— 变宽容和变成一个立场,发生在不同的人身上

**类型:FRONTIER。**

**心理学的那一句(本轮要判的):`#745` 说耦合只发生在常去教堂的人身上。
那他们是不是也「动得最多」?若是,耦合就只是「动得多」的副产品;
若不是 —— 若动得最少的那一群恰恰耦合最强 —— 那么变宽容与变成一个立场,
是两件不同的事,而且发生在不同的人身上。**

## ⚠ 按 `#746`① 先算 MDE,再跑(而不是跑完才发现零比效应大)
逐层逐年 n 中位 234–292,题内 sd 1.14–1.41 ⇒
**水平斜率零的 95% 分位 = 36 年 0.0568(常去)/ 0.0845(几乎不去)。**
**实测的水平变化是它的 6–26 倍 ⇒ 这一问能判。**(而 `#746` 那一问的零比效应大 13 倍,所以判不了。)

## 硬规则①(已跑)
`homosex` 1988 → 2024:**几乎不去 1.788 → 3.258 · 偶尔 1.538 → 3.022 · 常去 1.431 → 1.973**。
逐层 n:6,030 / 6,667 / 8,381;≥150 的调查年 18 / 21 / 25。

## G1 ESTIMAND
每层的 **`homosex` 水平对年份的斜率**(以及 `premarsx` 的,作对照),
与 `#745` 已测的**耦合斜率**并排。
## W1–W4(双边)
| 世界 | 判据 | 读法 |
|---|---|---|
| **W1 耦合是「动得多」的副产品** | 水平斜率与耦合斜率**同序**(常去最大) | 耦合不独立,页上那句要降级 |
| **W2 两件事,不同的人** | **反序**(常去水平最小、耦合最大) | **变宽容与变成一个立场分开了** |
| **W3 水平齐平** | 三层水平斜率彼此在 MDE 内 | 只有耦合在分人 |
| **W4 机械** | 变异度上升的那一层耦合更强 | 耦合差是方差差造出来的 |

⚠ **W1 与 W4 的正结果都削 `#745` —— 这正是本轮设计成能出它们的理由。**
⚠ **而 W4 有一个可先验的方向**:`homosex` 从 1.79 走到 3.26 的那一层**穿过了量表中段**,
方差应当上升、相关能力应当变强 ⇒ **机械解释预测「几乎不去」耦合更强。而 `#745` 实测它是零。**
**⇒ 机械解释的预测与观测相反,但这一条必须量出来,不许只靠推理。**

## G2 CONTROLS
**④ 正对照**:`#745` 的三层耦合斜率必须复现(**+0.00287 / +0.00178 / −0.00044**,容差 0.0005)。
**零** = `negative_control`,**零的种类 = 在每个礼拜频率层内部打乱受访者的年份标签 ——
保住该层构成、每年 n 与作答分布,只毁掉「同一层里谁属于哪一年」。**
**W4 控制(同一迭代内)**:逐层逐年的 `homosex` **标准差轨迹** —— 直接量,不推理。
## G3:2 个题 × 3 层 = 6 格全报 + 方差轨迹。G4:二分切法(去不去)。
## ⑤ 停止条件(**双边**,跑之前写死)
- **耦合三层复现不到 0.0005 ⇒ UNVERIFIED 并停。**
- **常去的水平斜率 ≥ 几乎不去的 ⇒ 判 W1**(耦合是副产品);
- **常去的水平斜率 ≤ 几乎不去的一半,且两者之差 > MDE ⇒ 判 W2**;
- **三层水平斜率两两之差都 < MDE ⇒ 判 W3**;
- **方差上升最多的那一层耦合也最强 ⇒ 判 W4,`#745` 降级。**
## IMPOSSIBLE(不写 planned)
仍是**重复横断面** ⇒ 「这一层的美国人」,不是「同一个人」;
`attend` 是**自报频率**;**换不了仪器**;
⚠ **本轮不解释为什么** —— `#746` 刚因为我给的解释没被支持而撤回过一次,
**本轮只报「是两件事」这个事实,不配解释。**`[unchallenged]`
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
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year","attend",A,Bc],encoding="latin1")
R=g.dropna(subset=[A,Bc,"attend"]).copy()
R["rel"]=pd.cut(R.attend,[-1,1,4,8],labels=["几乎不去","偶尔","常去"]); R=R.dropna(subset=["rel"]).copy()
LV=["几乎不去","偶尔","常去"]; FN,FY=150,5
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def series(fr,col):
    out=[]
    for y,sub in fr.groupby("year"):
        if len(sub)<FN: continue
        out.append((float(y),float(sub[col].mean()),float(sub[col].std()),len(sub)))
    return out
def nslope(fr):
    pts=[]
    for y,sub in fr.groupby("year"):
        if len(sub)<FN or sub[A].nunique()<2 or sub[Bc].nunique()<2: continue
        x=sub[A].to_numpy(float); yv=sub[Bc].to_numpy(float); r=sp(x,yv)
        xs=np.sort(x); ys=np.sort(yv); ys=ys if r>0 else ys[::-1]; c=abs(sp(xs,ys))
        if c<1e-9: continue
        pts.append((float(y),r/c))
    return (float(np.polyfit([p[0] for p in pts],[p[1] for p in pts],1)[0]) if len(pts)>=FY else np.nan)
print("=== ④ 正对照:`#745` 的三层耦合斜率必须复现(容差 0.0005)===")
LED={"常去":0.00287,"偶尔":0.00178,"几乎不去":-0.00044}
cpl={lv:nslope(R[R.rel==lv]) for lv in LV}; ok=True
for lv in LV:
    d=abs(cpl[lv]-LED[lv]); ok&= d<=0.0005
    print(f"  {lv:8s} 实测 {cpl[lv]:+.5f} · 账本 {LED[lv]:+.5f} · 差 {d:.5f} {'✅' if d<=0.0005 else '⛔'}")
if not ok:
    print("\n⛔ ⑤ 触发,停"); sys.exit(0)
print("\n=== G3:水平斜率(2 题 × 3 层)与耦合并排 ===")
print(f"{'层':10s}{'homosex 斜率/年':>16s}{'36 年':>9s}{'premarsx 斜率/年':>17s}{'36 年':>9s}{'耦合斜率/年':>13s}")
lvl={}
for lv in LV:
    fr=R[R.rel==lv]
    sh=series(fr,Bc); spm=series(fr,A)
    a=float(np.polyfit([p[0] for p in sh],[p[1] for p in sh],1)[0])
    b=float(np.polyfit([p[0] for p in spm],[p[1] for p in spm],1)[0])
    lvl[lv]=(a,b)
    print(f"{lv:10s}{a:>+16.5f}{a*36:>+9.3f}{b:>+17.5f}{b*36:>+9.3f}{cpl[lv]:>+13.5f}")
print("\n=== W4 控制:逐层 `homosex` 的标准差轨迹(直接量,不推理)===")
print(f"{'层':10s}{'sd 首':>8s}{'sd 末':>8s}{'sd 斜率/年':>12s}")
sds={}
for lv in LV:
    sh=series(R[R.rel==lv],Bc)
    s0=np.median([p[2] for p in sh[:4]]); s1=np.median([p[2] for p in sh[-4:]])
    ss=float(np.polyfit([p[0] for p in sh],[p[2] for p in sh],1)[0]); sds[lv]=ss
    print(f"{lv:10s}{s0:>8.3f}{s1:>8.3f}{ss:>+12.5f}")
rng=np.random.default_rng(20260806)
print("\n=== 零(层内打乱年份,B=300)· 水平斜率 ===")
nul={}
for lv in LV:
    fr=R[R.rel==lv]; v=[]
    for _ in range(300):
        P=fr.copy(); P["year"]=rng.permutation(P.year.to_numpy())
        sh=series(P,Bc)
        if len(sh)>=FY: v.append(float(np.polyfit([p[0] for p in sh],[p[1] for p in sh],1)[0]))
    q=np.quantile(np.abs(v),0.95); nul[lv]=float(q)
    print(f"  {lv:10s} |零| 95% 分位 {q:.5f}/年 = 36 年 {q*36:.4f} · 实测 {lvl[lv][0]*36:+.3f} ⇒ "
          f"**{abs(lvl[lv][0])/q:.1f} 倍**")
lo,hi=lvl["几乎不去"][0],lvl["常去"][0]
mde=max(nul.values())
G=Gate("变宽容和变成一个立场发生在不同的人身上")
p1=G.positive_control("`#745` 的三层耦合斜率必须复现(容差 0.0005)",planted=1.0 if ok else 0.0,floor=0.0,spread=0.1)
p2=G.negative_control("层内打乱年份后水平斜率应回到零",null=float(max(nul.values())),
    effect=abs(hi),null_spread=0.00005,
    null_kind="在每个礼拜频率层内部打乱受访者的年份标签 —— 保住该层构成、每年 n 与作答分布,只毁掉「同一层里谁属于哪一年」")
w4 = max(sds,key=lambda k:sds[k])==max(cpl,key=lambda k:cpl[k])
if not p1: v="**UNVERIFIED:旧值不可复现**"
elif w4: v=f"**W4:方差上升最多的层({max(sds,key=lambda k:sds[k])})也是耦合最强的层 ⇒ 耦合差可能是方差差造的,`#745` 降级**"
elif hi>=lo: v=f"**W1:常去的水平斜率 {hi:+.5f} ≥ 几乎不去的 {lo:+.5f} ⇒ 耦合是「动得多」的副产品**"
elif hi<=lo/2 and abs(hi-lo)>mde: v=(f"**W2:常去的水平只走 {hi*36:+.3f},不到几乎不去 {lo*36:+.3f} 的一半,"
    f"而它的耦合最强 ⇒ 变宽容与变成一个立场是两件事,而且发生在不同的人身上**")
elif all(abs(lvl[a][0]-lvl[b][0])<mde for a in LV for b in LV): v="**W3:三层水平斜率彼此在 MDE 内 ⇒ 只有耦合在分人**"
else: v=f"**报展布不报判决:水平 {lo*36:+.3f} / {lvl['偶尔'][0]*36:+.3f} / {hi*36:+.3f}**"
print(f"\n{v}"); print(G)
json.dump(dict(coupling=cpl,level={k:dict(homosex=v0[0],premarsx=v0[1]) for k,v0 in lvl.items()},
  sd_slope=sds,null_level={k:nul[k] for k in nul},verdict=v,unchallenged=True),
  open(OUT/"lc.json","w"),indent=1,ensure_ascii=False)

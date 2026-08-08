"""E03·A23·R122 —— 编码极性把上一轮整个翻过来

**类型:FRONTIER。** `#679` 的 NEXT 本来是「把分化当一个量测出来」,
而硬规则①第一步(打印分母)就发现分母是**负的且在变得更负** ——
去查码本,`#679` 的两条结论**都是编码极性的产物**。

## 从码本读出的极性(door ①:看对象,不看我对对象的记忆)
| 组 | 码 | 宽容/容忍方向 | 处理 |
|---|---|---|---|
| 性四题 | 1=always wrong … 4=not wrong at all | **高** | 保持 |
| `spk*` (5) | 1=allowed, 2=not | **低** | **反转** |
| `col*` (5) | 4=allowed, 5=not | **低** | **反转** |
| ⚠ `colcom` | **4=yes FIRED, 5=not fired** | **高** | **不反转 —— 同组内单独一题反向** |
| `lib*` (5) | 1=remove, 2=not remove | **高** | 保持 |
| `suicide*` (4) · `ab*` (7) | 1=yes, 2=no | **低** | **反转** |
| `nat*` (5) | 1=too little … 3=too much | **无宽容轴** | **排除,并写明理由** |

**⇒ 三套极性约定 + 一处组内翻转,而 `#679` 一处都没查。**

## 这一轮撤回 `#679` 的两条
**① 「教育把性与堕胎、安乐死推开」方向反了。** 性题宽容=高、堕胎题宽容=低 ⇒
**负相关是「一致」,更负 = 绑得更紧。** 对齐后重测方向。
**② 「分裂是按题组不是按对象」不成立。** `lib*` 保持而 `spk*`/`col*` 反转,
**这正是那三组符号分裂的来源 —— 是编码极性,不是题组效应。**
⚠ **而 `#661` 已经点名过这个缺陷(「1 在不同题里意思不同」),一窗之内重犯。**

## G1 ESTIMAND(对齐之后)
全部题对齐到「高 = 更宽容/更容忍」,重算 `#679` 的两个量:
**① 逐题 Δ(31→26 题,`nat*` 排除)· ② 预注册的 `*homo` 组内对照。**
## G2 CONTROLS
**正对照**:对齐后性四题内部最弱一环仍 +0.5069(不受对齐影响,因为四题同极性)。
**极性正对照(本轮新增)**:对齐后,**同一题组内五个对象的相关应当同号** ——
若仍异号,说明对齐没做对,**当场记 UNVERIFIED**。
**零**:打乱 educ。
## KILL(条件式)
if 极性正对照过(组内同号) :
  `*homo` 组内对照三组同号且超零 -> **机制是对象** ·
  仍不同号 -> **机制测不出,而 `#679` 的「题组说」已被本轮独立否掉**
else UNVERIFIED —— 对齐本身没做对
## IMPOSSIBLE(不写 planned)
`nat*` 五题**没有宽容轴** ⇒ 结构性无法纳入对齐分析(只能报 |ρ|,本轮不报);
跨仪器:MFQ 无容忍题组;因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]
B={"容忍·言论":["spkath","spkrac","spkcom","spkmil","spkhomo"],
   "容忍·任教":["colath","colrac","colcom","colmil","colhomo"],
   "容忍·藏书":["libath","librac","libcom","libmil","libhomo"],
   "自杀":["suicide1","suicide2","suicide3","suicide4"],
   "堕胎":["abdefect","abnomore","abhlth","abpoor","abrape","absingle","abany"]}
FLIP=lambda c: (c.startswith("spk") or (c.startswith("col") and c!="colcom")
                or c.startswith("suicide") or c.startswith("ab"))
ALL=[c for v in B.values() for c in v]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["educ"]+SEX+ALL, apply_value_formats=False, encoding="latin1")
for c in ALL:
    if FLIP(c): df[c]=-df[c]
grid=np.arange(8,19.01,1.0); BW=2.5
def wc(x,y,w):
    mx=np.average(x,weights=w);my=np.average(y,weights=w);cx=x-mx;cy=y-my
    s=np.sqrt(np.average(cx*cx,weights=w)*np.average(cy*cy,weights=w))
    return np.average(cx*cy,weights=w)/s if s>1e-12 else np.nan
def norm_at(x,y,W):
    r=wc(x,y,W)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    idx=W>=np.quantile(W,0.5)
    xs=np.sort(x[idx]); ys=np.sort(y[idx]); ys=ys if r>0 else ys[::-1]
    c=np.corrcoef(xs,ys)[0,1]
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
CACHE={}
def prep(c):
    if c not in CACHE:
        j=df.dropna(subset=SEX+["educ",c])
        CACHE[c]=(j[c].rank().to_numpy(float),[j[s].rank().to_numpy(float) for s in SEX],
                  j["educ"].to_numpy(float),len(j))
    return CACHE[c]
def item_delta(X,S,e):
    out=[]
    for g in grid:
        W=np.exp(-0.5*((e-g)/BW)**2)
        if W.sum()<150: out.append(np.nan); continue
        v=[u for u in (norm_at(X,s,W) for s in S) if np.isfinite(u)]
        out.append(float(np.median(v)) if v else np.nan)
    y=np.array(out); m=np.isfinite(y)
    return float(np.mean(y[m][-3:])-np.mean(y[m][:3])) if m.sum()>=4 else np.nan
print("=== 极性正对照:对齐后同一题组内五个对象的相关应当同号 ===")
ok=True
for bn,items in B.items():
    lv=[]
    for c in items:
        X,S,e,n=prep(c); W=np.ones(len(X))
        lv.append(float(np.median([u for u in (norm_at(X,s,W) for s in S) if np.isfinite(u)])))
    same=min(lv)*max(lv)>0; ok&=same
    print(f"  {bn:10s} 全样本耦合 {[f'{v:+.3f}' for v in lv]}  {'✅ 同号' if same else '⛔ 异号'}")
print(f"\n=== 逐题 Δ(对齐后,26 题;nat* 已因无宽容轴排除)===")
rows=[]
for bn,items in B.items():
    for c in items:
        X,S,e,n=prep(c); d=item_delta(X,S,e)
        rng=np.random.default_rng(20260806)
        nul=np.array([abs(item_delta(X,S,rng.permutation(e))) for _ in range(120)])
        p=float(np.nanmean(nul>=abs(d)))
        rows.append(dict(item=c,batt=bn,n=n,delta=d,p=p,homo=c.endswith("homo")))
        print(f"  {c:10s} {bn:10s} n={n:>6,} Δ={d:>+8.4f} p={p:.4f} {'✅homo' if c.endswith('homo') else ''}")
ps=sorted(r["p"] for r in rows); C=len(rows)
surv=[r["item"] for r in rows if r["p"]<=0.05*(ps.index(r["p"])+1)/C]
print(f"\n多重性 BH(q=0.05, C={C}):存活 {len(surv)}/{C}")
cons={}
for bn,pre in [("容忍·言论","spk"),("容忍·任教","col"),("容忍·藏书","lib")]:
    h=[r["delta"] for r in rows if r["item"]==pre+"homo"][0]
    o=[r["delta"] for r in rows if r["item"].startswith(pre) and not r["homo"]]
    cons[bn]=h-float(np.median(o))
    print(f"  {bn}  *homo {h:+.4f}  其余中位 {np.median(o):+.4f}  组内对照 **{cons[bn]:+.4f}**")
sg=[np.sign(v) for v in cons.values()]
print(f"\n  三组组内对照 {[f'{v:+.4f}' for v in cons.values()]} —— {'**同号**' if len(set(sg))==1 else '**仍不同号**'}")
G=Gate("编码极性把上一轮整个翻过来")
p1=G.positive_control("极性对齐正对照:每个题组内五个对象的耦合同号",planted=1.0 if ok else 0.0,floor=0.5,spread=0.01)
if p1:
    v=(f"**机制是对象:三组组内对照 {[f'{x:+.4f}' for x in cons.values()]} 同号**" if len(set(sg))==1
       else f"**机制仍测不出,但 `#679` 的「题组说」已被独立否掉:组内对照 {[f'{x:+.4f}' for x in cons.values()]}**")
else: v="UNVERIFIED —— 对齐本身没做对(组内仍异号)"
print(f"\n{v}"); print(G)
json.dump(dict(rows=rows,contrasts=cons,bh=surv,polarity_ok=bool(ok),verdict=v,unchallenged=True),
          open(OUT/"polarity_audit_and_redo.json","w"),indent=1,ensure_ascii=False)

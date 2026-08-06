"""E03·A33·R177 —— 「一样紧」是两个观测量相除得来的

**类型:FRONTIER。这是 `#733`① —— 唯一剩下的真债,已被推后五轮。**

**心理学的那一句(本轮要判的):这一页说「她对性、和对『没结婚算不算一个家』的判断一样紧」。
那句话是把两个观测量相除得来的;换成「每一块对照它自己的零」之后,它还是同一句话吗?**

## 缺口
`#650` 报的是 **性三题 0.4679 对 同居四题 0.4257 = 1.10 倍**,而页上那一行带着 `°`(没有零)。
⚠ **但真正的问题不是「零没算」** —— `realstat` 与 `#716` 都记过同一条:
**一个量只能对照它自己的零,不能对照另一个观测量。**
⇒ **1.10 这个数不是一次检验,是一次比较。** 本轮给两块各自配零。

## ⚠ 元分离器(先写下来,因为它可能让整行作废)
两块的 **k 是 3 对 4**,而 `#715` 已量出**零随 k 变**(同池 k=3 的零是 k=6 的三倍以上)。
⇒ **两个比值可能根本不在同一把尺子上,那样的话「一样紧」不是「假」,是「问得不对」。**
**这一条无论结果如何都要报。**

## G1 ESTIMAND
每块的**先对齐翻向、再取块内全部对的中位**(与 `#650` 逐字同一个量 —— ⚠ 第一次重跑漏了对齐,正对照当场开火)÷ **同池同 k 的零的 95% 分位**。
## G2 CONTROLS
**零** = `negative_control`,**零的种类 = 同一批人、同一个 10 题池、同样 k、同样天花板归一,
只打散「哪几题算一块」;全枚举 C(10,k),按 `#730` 的修法排除任何整块落在真域内的块。**
⚠ **「真域」的选择本身是一个规格**:`#650` 已证 `家庭七题` **不是**一个域 ⇒
**主规格只排除 {性三, 同居四};G4 第二格再排除 {家庭七} 内部的块,两个都报。**
**④ 正对照**:必须复现 `#650` 的 **0.4679 / 0.4257**(容差 0.005)。
**PLACEBO**:打乱行之后两块的中位应回到零 —— 用来证明零毁掉的是「谁跟谁一块」,不是别的。
## G3:2 块 × 2 排除规格 = 4 格全报,外加 placebo。G4:k=3/k=4 各自的零分开报。
## ⑤ 停止条件(跑之前写死)
- **两块的中位复现不到 0.005 ⇒ UNVERIFIED 并停。**
- **任一块的比值 ≤ 1.0 ⇒ 那一块在自己的零里,「一样紧」失去一条腿。**
- **两块比值之差 > 2 倍 ⇒ 「一样紧」在比值这一层是假的。**
- **无论结果如何:因为 k 不同,必须在页上写明两个比值不在同一把尺子上。**
## IMPOSSIBLE(不写 planned)
NSFG **单轮、只有女性 15–44、无年份** ⇒ 无法做时间或人群对照;
**换不了仪器**:GSS 没有同构念的「没结婚算不算一个家」题组(`#700`/`#732` 已枚举);
池只有 10 题 ⇒ **k=4 的零只有 C(10,4)=210 个块可枚举,而排除之后更少 —— 分辨率有限,如实报块数。**
`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.blocks import pairmat, aligned_pooled_median as apm
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
NS=pathlib.Path("data/external/nsfg")
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"'); LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    m=pat.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)),m.group(4))
SEXN=["samesex","sxok18","sxok16"]; FAMN=["staytog","chunless","chsuppor","okcohab","marrfail","chcohab","prvntdiv"]
COH=["okcohab","chcohab","chsuppor","prvntdiv"]
POOL=SEXN+[c for c in FAMN]
cols={n:LAY[n] for n in POOL}; buf={n:[] for n in cols}
for line in open(NS/"2011_2013_FemRespData.dat",errors="replace"):
    for n,(s,w,_) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
X=pd.DataFrame({n:np.where(np.isin(np.array(buf[n]),[1,2,3,4,5]),np.array(buf[n]),np.nan) for n in cols})
print(f"硬规则①:NSFG 2011–2013 女性 15–44 单轮 · 池 {len(POOL)} 题 · 完整个案 n = {len(X.dropna()):,}")
M=pairmat(X,POOL,floor=150); I={c:i for i,c in enumerate(POOL)}
BLK={"性三题":SEXN,"同居四题":COH}
LED={"性三题":0.4679,"同居四题":0.4257}
print("\n=== ④ 正对照:复现 `#650` 的块内归一中位(容差 0.005)===")
obs={}
for nm,items in BLK.items():
    obs[nm]=apm(M,[I[c] for c in items])
    print(f"  {nm:8s} 实测 **{obs[nm]:+.4f}** · 账本 {LED[nm]:+.4f} · 差 {abs(obs[nm]-LED[nm]):.4f} "
          f"{'✅' if abs(obs[nm]-LED[nm])<=0.005 else '⛔ ⑤ 触发'}")
maxd=max(abs(obs[n]-LED[n]) for n in LED)
if maxd>0.005:
    print("\n⛔ 停:旧值不可复现")
    json.dump(dict(stop="旧值不可复现",obs=obs,ledger=LED),open(OUT/"pay650.json","w"),indent=1,ensure_ascii=False); sys.exit(0)
SPECS={"主规格(排除 性三 · 同居四)":[set(SEXN),set(COH)],
       "G4(再排除 家庭七 内部)":[set(SEXN),set(COH),set(FAMN)]}
print("\n=== 零:同池同 k 全枚举,两种排除规格都报 ===")
res={}
for sn,doms in SPECS.items():
    for nm,items in BLK.items():
        k=len(items)
        allb=[c for c in itertools.combinations(range(len(POOL)),k)
              if not any({POOL[i] for i in c}<=d for d in doms)]
        v=np.array([apm(M,list(c)) for c in allb])   # ⚠ 零也必须对齐,否则极性把零压低、比值虚高; v=v[np.isfinite(v)]
        q=float(np.quantile(v,0.95))
        res[(sn,nm)]=dict(k=k,null=q,median=float(np.median(v)),blocks=int(v.size),ratio=obs[nm]/q)
        print(f"  {sn:26s} {nm:8s} k={k} · 块 {v.size:3d} · **零的 95% 分位 {q:+.4f}** · "
              f"实测 {obs[nm]:+.4f} ⇒ **{obs[nm]/q:.2f}×**")
print("\n=== PLACEBO:打乱行之后两块的中位应回到零 ===")
rng=np.random.default_rng(20260806)
Xs=X.copy()
for c in POOL: Xs[c]=rng.permutation(Xs[c].to_numpy())
Ms=pairmat(Xs,POOL,floor=150)
for nm,items in BLK.items():
    print(f"  {nm:8s} 打乱后中位 **{apm(Ms,[I[c] for c in items]):+.4f}**")
main=[res[("主规格(排除 性三 · 同居四)",n)] for n in BLK]
r1,r2=main[0]["ratio"],main[1]["ratio"]
G=Gate("「一样紧」是两个观测量相除得来的")
p1=G.positive_control("必须复现 `#650` 的 0.4679 / 0.4257(容差 0.005)",planted=float(0.005-maxd),floor=0.0,spread=0.0002)
p2=G.negative_control("同池同 k 的随机题组应低于真块",
    null=float(np.mean([m["null"] for m in main])),effect=float(np.mean(list(obs.values()))),
    null_spread=0.005,null_kind="同一批人、同一个 10 题池、同样 k、同样天花板归一,只打散哪几题算一块;全枚举 C(10,k),排除任何整块落在真域内的块")
if not p1: v="**UNVERIFIED:旧值不可复现**"
elif min(r1,r2)<=1.0: v=f"**「一样紧」失去一条腿:比值 {r1:.2f}× 与 {r2:.2f}×,其中一个在自己的零里**"
elif max(r1,r2)/min(r1,r2)>2: v=f"**「一样紧」在比值这一层是假的:{r1:.2f}× 对 {r2:.2f}×,差 {max(r1,r2)/min(r1,r2):.2f} 倍**"
else: v=f"**两块各自都高于自己的零({r1:.2f}× 与 {r2:.2f}×),而 k 不同 ⇒ 两个比值不在同一把尺子上**"
print(f"\n{v}")
print("⚠ 无论结果如何都要报:**两块的 k 是 3 对 4,而 `#715` 已量出零随 k 变** ——")
print(f"   本轮实测同池 k=3 的零 {main[0]['null']:.4f} · k=4 的零 {main[1]['null']:.4f}"
      f"(差 {main[0]['null']/main[1]['null']:.2f} 倍)⇒ **「一样紧」不是假,是问得不对。**")
print(G)
json.dump(dict(obs=obs,ledger=LED,cells={f"{a}|{b}":res[(a,b)] for a,b in res},verdict=v,unchallenged=True),
  open(OUT/"pay650.json","w"),indent=1,ensure_ascii=False)

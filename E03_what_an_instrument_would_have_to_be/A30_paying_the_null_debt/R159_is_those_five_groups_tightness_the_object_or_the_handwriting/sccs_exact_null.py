"""E03·A30·R159 —— 民族志那五组的紧,是「同一手段施于四个对象」,还是编码者的笔迹

**类型:FRONTIER。**

**心理学的那一句:一个社会「怎么管孩子」在民族志里被读成一件事 —— 那是这个社会真的成套,
还是同一个编码团队把整批育儿变量涂成了一片?**

## 为什么现在做(`#716`①)
`#716` 量出 SCCS 五组的**生**数字最大(0.651–0.896)而**比值最小**(2.54–3.50×),
因为**它们的零也最大(0.2560)** —— 任取四个 Barry 1977 的育儿变量本来就相当一致。
⚠ **而那个零上一轮还在抽样(2,500/4,840),本轮全枚举。**
⚠ **且必须说清这个零测的是什么:池内 20 个变量全出自同一编码团队 ⇒
它测的正是 `#528` 的笔迹 —— 「同一支笔写出来的东西本来多一致」。不能只报一个倍数。**

## 池的结构给了一个免费的分离器
池 = **5 手段 × 4 对象**(以身作则 429–432 · 讲课 437–440 · 体罚 453–456 · 放任 465–468 · 疼爱 469–472)。
- **真块** = 同一手段 × 四个对象(5 个)—— `#653` 判「团」的那五个。
- **对照块** = **同一对象 × 四种手段**(固定对象,5 选 4 ⇒ 每对象 5 块 × 4 对象 = 20 块)。
- **全零** = 池内任取 4 个(4,840 块,排除 5 个真块)。

## W1 / W2
| | 真块(手段内) | 同对象块(手段间) | 读法 |
|---|---|---|---|
| **W1 对象是接缝** | 高 | **明显更低** | 「同一手段施于谁」才是一件事 |
| **W2 是笔迹** | 高 | **一样高** | **整批育儿变量都被涂成一片,`#653` 那五个「团」测的是编码者** |

**W2 的正结果我不高兴** —— 页上五条民族志的「团」会从「社会成套」降级成「笔迹」。

## G1 ESTIMAND
最弱一环(归一),与 `#653` 逐字同一条路径。
## G2 CONTROLS
**④ 正对照**:复现 `#653` 五组(0.895 / 0.822 / 0.786 / 0.785 / 0.651),容差 0.001。
**零**:`negative_control`,**零的种类 = 同一编码团队(Barry 1977)的同一个 20 变量池、同样 k=4、同样对齐,
只把「哪四个算一组」打散 —— 所以它测的是笔迹,不是随机。**
## ⑤ 停止条件(跑之前写死)
- 正对照复现不到 0.001 ⇒ 停。
- **同对象块的中位 ≥ 真块中位的 80% ⇒ 判 W2(笔迹)**,页上五条要加限定。
- **≤ 50% ⇒ 判 W1(对象是接缝)。** 落在 50–80% 之间 ⇒ **记「判不了」,报区间。**
## IMPOSSIBLE(不写 planned)
只有一个编码团队 ⇒ **换不了笔**;`#639` 那种跨团队复制只在体罚一个手段上有。
⚠ **本轮结构性地只有一具仪器,而这不是偷懒:换不了仪器** —— 问的就是 SCCS 这一个池**自己的内部结构**
(5 手段 × 4 对象是它独有的编码设计),**别的档案没有同一套「手段 × 对象」网格**(`#700` 已枚举)。
在别处重问这一句,需要的不是另一份数据,是**另一次田野**。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,s=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float)); return sp(x,y[::-1] if s<0 else y)
def weakest(F,items,floor=30):
    A=F[items].dropna()
    if len(A)<floor: return np.nan
    m=A.mean(axis=1); A=A.copy()
    for i in items:
        if sp(A[i],m)<0: A[i]=-A[i]
    v=[]
    for a,b in itertools.combinations(items,2):
        mm=A[[a,b]].dropna()
        if len(mm)<floor or mm[a].nunique()<2 or mm[b].nunique()<2: return np.nan
        r=sp(mm[a],mm[b])
        if not np.isfinite(r) or r==0: return np.nan
        c=rmax(mm[a],mm[b],1 if r>0 else -1)
        if not np.isfinite(c) or abs(c)<1e-9: return np.nan
        v.append(r/abs(c))
    return float(np.min(v))
S=pathlib.Path("data/external/dplace/repo/datasets/SCCS")
W=pd.read_csv(S/"data.csv").pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
TECH={"以身作则":429,"讲课":437,"体罚":453,"放任":465,"疼爱":469}
GRP={k:[f"SCCS{b+i}" for i in range(4)] for k,b in TECH.items()}
POOL=[v for L in GRP.values() for v in L]
print("=== 硬规则①:池 = 5 手段 × 4 对象,逐变量 n ===")
for k,L in GRP.items():
    print(f"  {k:5s} {L}  联合 n={len(W[L].dropna()):3d}  逐变量 n=" + "/".join(str(int(W[c].notna().sum())) for c in L))
LED={"以身作则":0.895,"体罚":0.822,"讲课":0.786,"疼爱":0.785,"放任":0.651}
print("\n=== ④ 正对照(容差 0.001)===")
tru={}
for k in LED:
    tru[k]=weakest(W,GRP[k]); print(f"  {k:5s} 实测 {tru[k]:+.4f} 账本 {LED[k]:+.3f} 差 {abs(tru[k]-LED[k]):.4f} {'✅' if abs(tru[k]-LED[k])<=0.001 else '⛔'}")
maxd=max(abs(tru[k]-LED[k]) for k in LED)
if maxd>0.001:
    print(f"\n⛔ ⑤ 触发:最大差 {maxd:.4f} ⇒ 旧值不可复现,停"); sys.exit(0)
# 同对象块:固定对象序号 i,5 种手段里取 4
same_t=[]
for i in range(4):
    col=[GRP[k][i] for k in TECH]
    for c in itertools.combinations(col,4): same_t.append((i,c,weakest(W,list(c))))
st=np.array([x[2] for x in same_t if np.isfinite(x[2])])
truv=np.array(list(tru.values()))
print(f"\n=== 对照:同一对象 × 四种手段({len(st)}/{len(same_t)} 块可算)===")
for i in range(4):
    r=[x[2] for x in same_t if x[0]==i and np.isfinite(x[2])]
    print(f"  对象 {i}  中位 {np.median(r):+.4f}  范围 [{min(r):+.4f}, {max(r):+.4f}]")
print(f"  **同对象块中位 {np.median(st):+.4f}** vs **真块(同手段)中位 {np.median(truv):+.4f}** "
      f"⇒ 比 **{np.median(st)/np.median(truv)*100:.0f}%**")
# 全零:全枚举
truset={frozenset(v) for v in GRP.values()}
allb=[c for c in itertools.combinations(POOL,4) if frozenset(c) not in truset]
v=np.array([weakest(W,list(c)) for c in allb]); n_bad=int(np.sum(~np.isfinite(v))); v=v[np.isfinite(v)]
q95=float(np.quantile(v,0.95)); q99=float(np.quantile(v,0.99))
print(f"\n=== 全零:**全枚举 {len(allb):,} 块**(排除 5 真块;{n_bad} 块因地板/常数不可算)===")
print(f"  95% 分位 **{q95:+.4f}** · 99% 分位 {q99:+.4f} · 中位 {np.median(v):+.4f} · 最大 {v.max():+.4f}")
print(f"  ⚠ 上一轮抽 2,500 给的是 +0.2560 —— 精确值 **{q95:+.4f}**,差 {abs(q95-0.2560):.4f}")
print("\n五组对照精确零:")
for k in sorted(tru,key=lambda x:-tru[x]):
    print(f"  {k:5s} {tru[k]:+.4f}  **{tru[k]/q95:5.2f}×**  过 q99 {'✅' if tru[k]>q99 else '⛔'}")
ratio=float(np.median(st)/np.median(truv))
from lib.gates import Gate
G=Gate("那五组的紧是对象还是笔迹")
p1=G.positive_control("必须复现 #653 五组(最大差 <0.001)",planted=float(0.001-maxd),floor=0.0,spread=0.00005)
p2=G.negative_control("同池随机四变量应低于真块",null=q95,effect=float(np.median(truv)),null_spread=0.005,
  null_kind="同一编码团队 Barry 1977 的同一个 20 变量池、同样 k=4、同样对齐,只打散「哪四个算一组」—— 所以它测的是笔迹")
if not p1: verdict="**判不了:旧值不可复现**"
elif ratio>=0.80: verdict=f"**W2 笔迹:同对象块中位是真块的 {ratio*100:.0f}%,整批育儿变量被涂成一片**"
elif ratio<=0.50: verdict=f"**W1 对象是接缝:同对象块只有真块的 {ratio*100:.0f}%**"
else: verdict=f"**判不了({ratio*100:.0f}% 落在 50–80% 之间)—— 报区间,不报判决**"
print(f"\n{verdict}"); print(G)
json.dump(dict(true=tru,ledger=LED,max_diff=maxd,same_target_median=float(np.median(st)),
  true_median=float(np.median(truv)),ratio=ratio,null_q95=q95,null_q99=q99,n_blocks=int(v.size),
  sampled_prev=0.2560,verdict=verdict,unchallenged=True),open(OUT/"sccs_exact.json","w"),indent=1,ensure_ascii=False)

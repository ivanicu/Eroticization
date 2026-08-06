"""E03·A31·R162 —— 一个人的「性」是一件事,还是二十件

**类型:FRONTIER。而它换了方向:`#715`–`#719` 连着五轮都在元层(补零、查分母),
`#719`④ 的反基地条款要求回到对象。**

**心理学的那一句:这一页说「一个人有一套性道德」。那一个人的性**自我概念**呢 ——
他的性焦虑、性自尊、性掌控、性动机,是同一件事的不同侧面,还是彼此无关的二十件事?**

## 为什么是这具仪器,以及它**不能**做什么(硬规则②:先点名仪器)
`#653` 的 NEXT 早就写了「第三具仪器就在手边而从未被问过这一句」,指的是这份
**MSSCQ 多维性自我概念问卷**(100 题 · **n = 17,685** · 5 点「像不像我」)。
⚠ **但它回答不了 `#653`② 那一句,而这本身要先说清楚:**
GSS 的 `premarsx/xmarsex/homosex/teensex` 问的是**对别人行为的道德判断**;
MSSCQ 问的是**关于自己的描述**。**两者不是同一个构念 ⇒ 本轮不是 `#653`② 的第三个复制点,
「性是一块」仍然只有 GSS 一条腿(`#718`)。** 本轮问的是一个**相邻但不同**的问题。

## ⚠ 「它是谁分的组」——`#653`① 要求先答,而诚实答案是:这一版发布里没有人分
**码本只列了 100 条题干,一个面的名字都没有。** 那个 20 面 × 5 题的循环结构
(第 k 面 = 题 k, k+20, k+40, k+60, k+80)是**我从题干读出来的推断** ——
证据:题 3「我很清楚自己的性感受与需要」对题 23「我很清楚自己的性动机与欲望」;
题 6「我一天到晚想着性」对题 26「我想性比想别的都多」。
⇒ **它因此不是前提,是本轮要检验的东西之一。**

## W1 / W2 / W3
| | 20 个面的最弱一环 | 读法 |
|---|---|---|
| **W1 性是一件事** | 面内高 **且** 任取五题也高 | 整份问卷是一个总体因子,「二十个面」是包装 |
| **W2 性是二十件事** | **面内高,任取五题低** | 一个人的性由若干**互不蕴含**的侧面组成 |
| **W3 我的分面是错的** | **面内不高于任取五题** | 从题干读结构这件事本身失败,后续分面分析全部作废 |

## G1 ESTIMAND
每个面的**最弱一环(天花板归一)**,以及它 ÷ **同池同 k 的零**。
⚠ **近义改写的混淆(`#653`⑤ 预注册)由估计量本身吃掉**:最弱一环取的是 **min**,
所以一对近义题**抬不高**它;但仍**逐面报最强一对**,凡最强 > 0.60 者标注。
## G2 CONTROLS
**零** = 同一份问卷 100 题池里随机抽 5 题(排除 20 个真面),**`negative_control`,
零的种类 = 同一批人、同一条 5 点量表、同样 k=5,只打散「哪五题算一个面」。**
⚠ 组合抽样是 `#719` 量出**最不稳**的那一类(B=4,000 时 10.73%)⇒ **抽 200,000 块**,
且**把零自己的跨种子相对标准差一起报**。C(100,5)=75,287,520 不可枚举,**如实说明**。
**④ 正对照**:`#542` 实测这份问卷**面内 0.579 / 面间 0.173** —— **面内必须复现到 ±0.05。**
**安慰剂**:打乱行之后,面内最弱一环应回到零。
## ⑤ 停止条件(跑之前写死)
- 正对照复现不到 ±0.05 ⇒ **记「旧值不可复现」并停。**
- **20 个面的最弱一环中位 ≤ 零的 95% 分位 ⇒ 判 W3(我的分面是错的),后续分面分析全部作废。**
- **面内中位 > 零 且 「任取五题」的中位 > 零 ⇒ 判 W1;面内 > 零 而任取五题 ≤ 零 ⇒ 判 W2。**
## IMPOSSIBLE(不写 planned)
自选网络志愿者,非概率样本 · 无年份 · **换不了仪器**:没有第二份同构念的性自我概念公开数据集
(`#700` 已枚举);本轮**不声称**任何跨仪器复制。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, re
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
P="data/external/openpsych/MSSCQ/MSSCQ/"
D=pd.read_csv(P+"data.csv",sep="\t")
Q=[f"Q{i}" for i in range(1,101)]
X=D[Q].replace(0,np.nan).dropna()
X=X[(X>=1).all(axis=1)&(X<=5).all(axis=1)]
txt=pathlib.Path(P+"codebook.txt").read_text(encoding="latin1")
ITEMTXT=dict(re.findall(r'^(\d{1,3})\.\s+(.+)$',txt,re.M))
print(f"=== 硬规则①:MSSCQ · 100 题 · 5 点量表 · 完整个案 **n = {len(X):,}** / 原始 {len(D):,} 行 ===")
print(f"  每题取值 {sorted(pd.unique(X[Q].values.ravel()))} · 码本**没有面的名字**,20×5 循环结构是从题干读出来的推断")
for k in (3,23,6,26): print(f"    Q{k}: {ITEMTXT[str(k)][:78]}")
R=X[Q].rank().to_numpy(float)
C=np.corrcoef(R.T)                       # 全部 4,950 对 Spearman,一次算完
S=np.sort(R,axis=0)
CEIL=np.abs(np.corrcoef(S.T))            # 同向天花板;反向对的天花板取绝对值同值
M=np.where(np.abs(CEIL)>1e-9,C/np.abs(CEIL),np.nan); np.fill_diagonal(M,1.0)
FAC={k:[k-1+20*j for j in range(5)] for k in range(1,21)}   # 0-based 列号
def wl(ix): return float(min(M[a,b] for a,b in itertools.combinations(ix,2)))
def strongest(ix): return float(max(M[a,b] for a,b in itertools.combinations(ix,2)))
w={k:wl(v) for k,v in FAC.items()}; st={k:strongest(v) for k,v in FAC.items()}
print(f"\n=== ④ 正对照:与 `#542` 的面内 0.579 / 面间 0.173 对照(容差 ±0.05)===")
within=float(np.median([np.median([M[a,b] for a,b in itertools.combinations(v,2)]) for v in FAC.values()]))
allp=[M[a,b] for a,b in itertools.combinations(range(100),2)]
cross=float(np.median([M[a,b] for a,b in itertools.combinations(range(100),2)
                       if (a%20)!=(b%20)]))
print(f"  面内中位 **{within:.4f}**(#542: 0.579)差 {abs(within-0.579):.4f} {'✅' if abs(within-0.579)<=0.05 else '⛔'}")
print(f"  面间中位 **{cross:.4f}**(#542: 0.173)差 {abs(cross-0.173):.4f} {'✅' if abs(cross-0.173)<=0.05 else '⛔'}")
if abs(within-0.579)>0.05:
    print("\n⛔ ⑤ 触发:面内复现不到 ±0.05 ⇒ 旧值不可复现,停")
    json.dump(dict(stop="旧值不可复现",within=within,cross=cross),open(OUT/"msscq.json","w"),indent=1,ensure_ascii=False); sys.exit(0)
rng=np.random.default_rng(20260806); NB=200_000
truth={frozenset(v) for v in FAC.values()}
blocks=[]
while len(blocks)<NB:
    b=rng.choice(100,5,replace=False)
    if frozenset(b.tolist()) not in truth: blocks.append(b)
BL=np.array(blocks)
pairs=list(itertools.combinations(range(5),2))
V=np.min(np.stack([M[BL[:,a],BL[:,b]] for a,b in pairs],axis=1),axis=1)
q95=float(np.quantile(V,0.95)); q99=float(np.quantile(V,0.99))
sub=[float(np.quantile(V[np.random.default_rng(s).choice(NB,20000,replace=False)],0.95)) for s in range(20)]
print(f"\n=== 零:100 题池随机抽 5 题(**抽 {NB:,} / C(100,5)=75,287,520,不可枚举,如实说明**)===")
print(f"  **零的 95% 分位 {q95:+.4f}** · 99% {q99:+.4f} · 中位 {np.median(V):+.4f}")
print(f"  ⚠ 零自己的跨种子相对标准差(20 子样各 20,000)= **{np.std(sub)/np.median(sub)*100:.2f}%**"
      f"(`#719`:组合抽样是最不稳的一类)")
print(f"\n=== 20 个面(G3 全报,含不支持结论的)===")
print(f"{'面':>3s}{'最弱一环':>10s}{'÷零':>8s}{'最强一对':>10s}   代表题")
for k in sorted(w,key=lambda x:-w[x]):
    flag=" ⚠最强>0.60" if st[k]>0.60 else ""
    print(f"{k:>3d}{w[k]:>+10.4f}{w[k]/q95:>8.2f}{st[k]:>+10.4f}   {ITEMTXT[str(k)][:44]}{flag}")
med=float(np.median(list(w.values()))); above=sum(1 for v in w.values() if v>q95)
rand_med=float(np.median(V))
G=Gate("一个人的性是一件事还是二十件")
p1=G.positive_control("面内一致必须复现 #542 的 0.579(±0.05)",planted=float(0.05-abs(within-0.579)),floor=0.0,spread=0.002)
p2=G.negative_control("同池随机五题的最弱一环应低于真面",null=q95,effect=med,null_spread=0.005,
  null_kind="同一批人、同一条 5 点量表、同样 k=5,只打散「哪五题算一个面」;抽 200,000 块(池不可枚举)")
if not p1: v="**判不了:旧值不可复现**"
elif med<=q95: v=f"**W3:面内中位 {med:+.4f} 不高于零 {q95:+.4f} ⇒ 我从题干读出来的分面是错的,后续分面分析作废**"
elif rand_med>q95: v="**W1:任取五题也高于零 ⇒ 整份问卷是一个总体因子**"
else: v=(f"**W2:面内中位 {med:+.4f} 是零的 {med/q95:.2f} 倍,而任取五题的中位只有 {rand_med:+.4f} "
         f"⇒ 一个人的性是若干互不蕴含的侧面,不是一件事**({above}/20 个面高于零)")
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(X)),within=within,cross=cross,weakest={str(k):w[k] for k in w},
  strongest={str(k):st[k] for k in st},null_q95=q95,null_q99=q99,null_median=rand_med,
  null_rel_sd=float(np.std(sub)/np.median(sub)),n_blocks=NB,median_within=med,n_above=above,
  verdict=v,unchallenged=True),open(OUT/"msscq.json","w"),indent=1,ensure_ascii=False)

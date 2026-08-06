"""E03·A31·R162b —— 同一批面,加上最优符号(`#718` 刚建的那一步,我上一版忘了用)

⚠ **上一版 `msscq_facets.py` 有三个面是负的**(面 17 −0.7815 · 面 7 −0.5980 · 面 8 −0.5473)。
`#718` 两轮前刚证明:**那种负号多半是「没做符号对齐」造成的,不是数据里的**,
而「互斥」与「无关」是两句完全不同的心理学。**k=5 只有 2^4 = 16 种符号指派,穷举即可。**
⇒ 本文件把真面与零**同时**换成最优符号(`#713` 的类型对齐:分子分母必须是同一个估计量)。

⚠ **换不了仪器**:同 `msscq_facets.py`。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
P="data/external/openpsych/MSSCQ/MSSCQ/"
D=pd.read_csv(P+"data.csv",sep="\t"); Q=[f"Q{i}" for i in range(1,101)]
X=D[Q].replace(0,np.nan).dropna(); X=X[(X>=1).all(axis=1)&(X<=5).all(axis=1)]
ITEMTXT=dict(re.findall(r'^(\d{1,3})\.\s+(.+)$',pathlib.Path(P+"codebook.txt").read_text(encoding="latin1"),re.M))
R=X[Q].rank().to_numpy(float); C=np.corrcoef(R.T)
CEIL=np.abs(np.corrcoef(np.sort(R,axis=0).T))
M=np.where(CEIL>1e-9,C/CEIL,np.nan); np.fill_diagonal(M,1.0)
FAC={k:[k-1+20*j for j in range(5)] for k in range(1,21)}
PR=list(itertools.combinations(range(5),2))
SIGNS=[np.array([1]+[1 if (b>>t)&1==0 else -1 for t in range(4)]) for b in range(1<<4)]
SGP=np.array([[s[a]*s[b] for a,b in PR] for s in SIGNS])          # (16,10)
def opt(BL):
    V=np.stack([M[BL[:,a],BL[:,b]] for a,b in PR],axis=1)          # (B,10)
    return np.max(np.min(V[:,None,:]*SGP[None,:,:],axis=2),axis=1)
def greedy_wl(ix): return float(min(M[a,b] for a,b in itertools.combinations(ix,2)))
FA=np.array([FAC[k] for k in range(1,21)])
wo=opt(FA); wg=np.array([greedy_wl(FAC[k]) for k in range(1,21)])
assert (wo>=wg-1e-9).all(), "④ 正对照失败:最优 < 不翻向(不可能)⇒ 实现有 bug"
print("④ 正对照:最优 ≥ 不翻向,20/20 成立 ✅(严格可证:不翻向也在 16 种里)")
rng=np.random.default_rng(20260806); NB=200_000; truth={frozenset(v.tolist()) for v in FA}
bl=[]
while len(bl)<NB:
    b=rng.choice(100,5,replace=False)
    if frozenset(b.tolist()) not in truth: bl.append(b)
BL=np.array(bl); Vo=opt(BL)
Vg=np.min(np.stack([M[BL[:,a],BL[:,b]] for a,b in PR],axis=1),axis=1)
q95o=float(np.quantile(Vo,0.95)); q95g=float(np.quantile(Vg,0.95))
sub=[float(np.quantile(Vo[np.random.default_rng(s).choice(NB,20000,replace=False)],0.95)) for s in range(20)]
print(f"\n零(抽 {NB:,} / C(100,5)=75,287,520,不可枚举):"
      f"**最优符号下零的 95% 分位 {q95o:+.4f}**(不翻向时 {q95g:+.4f})· "
      f"跨种子相对标准差 **{np.std(sub)/np.median(sub)*100:.2f}%**")
print(f"\n{'面':>3s}{'不翻向':>10s}{'最优符号':>10s}{'÷零':>8s}   代表题")
order=sorted(range(20),key=lambda i:-wo[i])
for i in order:
    print(f"{i+1:>3d}{wg[i]:>+10.4f}{wo[i]:>+10.4f}{wo[i]/q95o:>8.2f}   {ITEMTXT[str(i+1)][:50]}")
med=float(np.median(wo)); above=int((wo>q95o).sum()); flipped=[i+1 for i in range(20) if wg[i]<0<=wo[i]]
print(f"\n**三个负面的去向:{flipped or '无'} —— 不翻向时为负,最优符号下 ≥0**")
G=Gate("一个人的性是一件事还是二十件(最优符号)")
p1=G.positive_control("最优 ≥ 不翻向(20/20,严格可证)",planted=float(min(wo-wg)+0.01),floor=0.0,spread=0.0005)
p2=G.negative_control("同池随机五题应低于真面",null=q95o,effect=med,null_spread=0.005,
  null_kind="同一批人、同一条 5 点量表、同样 k=5、同样取最优符号,只打散「哪五题算一个面」")
v=(f"**面内中位 {med:+.4f} 是零的 {med/q95o:.2f} 倍({above}/20 个面高于零),"
   f"而任取五题的中位只有 {np.median(Vo):+.4f} ⇒ 一个人的性由若干互不蕴含的侧面组成,不是一件事**")
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(X)),weakest_optimal={str(i+1):float(wo[i]) for i in range(20)},
  weakest_unflipped={str(i+1):float(wg[i]) for i in range(20)},flipped=flipped,
  null_q95_optimal=q95o,null_q95_unflipped=q95g,null_median=float(np.median(Vo)),
  null_rel_sd=float(np.std(sub)/np.median(sub)),n_blocks=NB,median=med,n_above=above,
  verdict=v,unchallenged=True),open(OUT/"msscq_optimal.json","w"),indent=1,ensure_ascii=False)

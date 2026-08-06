"""E03·A30·R158b —— MFQ k=6 的零:全枚举,不抽样

**为什么必须全枚举:** 同一个零,抽 1,500 块给 **0.0162**(`#715`),抽 3,000 块给 **0.0313**(`#716` 首版),
而它是页上每一个比值的**分母**。**尾部分位数在抽样下极不稳。**
块空间只有 C(30,6) = 593,775,**可枚举 ⇒ 就不该抽样**:精确、便宜、连种子这个自由度都不存在。

**加速的关键(而它不是近似):** 30 题与池均值的相关**全为正** ⇒ `align()` 在任何块里都不翻向
⇒ 块的最弱一环退化为 **435 个预算好的归一对相关上的 min**,查表即可。脚本会**先验证这个前提**再用它。
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pyreadstat
from scipy.stats import spearmanr
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
ITEM={"emotionally":"HARM","weak":"HARM","cruel":"HARM","compassion":"HARM","animal":"HARM","kill":"HARM",
 "treated":"FAIRNESS","unfairly":"FAIRNESS","rights":"FAIRNESS","fairly":"FAIRNESS","justice":"FAIRNESS","rich":"FAIRNESS",
 "lovecountry":"INGROUP","betray":"INGROUP","loyalty":"INGROUP","history":"INGROUP","family":"INGROUP","team":"INGROUP",
 "respect":"AUTHORITY","traditions":"AUTHORITY","chaos":"AUTHORITY","kidrespect":"AUTHORITY","sexroles":"AUTHORITY","soldier":"AUTHORITY",
 "decency":"PURITY","disgusting":"PURITY","god":"PURITY","harmlessdg":"PURITY","unnatural":"PURITY","chastity":"PURITY"}
mfq,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
X=mfq[list(ITEM)].dropna(); it=list(ITEM); n=len(it)
pos=all(sp(X[i],X.mean(axis=1))>0 for i in it)
print(f"硬规则①:MFQ n={len(X)} · 30 题 / 5 域 · 单次施测")
print(f"前提验证:30 题与池均值的相关全为正 = **{pos}** ⇒ {'查表合法' if pos else '⛔ 前提不成立,必须走 align 全算'}")
assert pos, "前提不成立"
M=np.eye(n)
for a in range(n):
    for b in range(a+1,n):
        x=X[it[a]].to_numpy(float); y=X[it[b]].to_numpy(float)
        r=sp(x,y); xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        M[a,b]=M[b,a]=r/abs(sp(xs,ys))
truth={frozenset(i for i,d in ITEM.items() if d==D) for D in set(ITEM.values())}
vals=[]
for c in itertools.combinations(range(n),6):
    if frozenset(it[j] for j in c) in truth: continue
    m=1.0
    for a,b in itertools.combinations(c,2):
        if M[a,b]<m: m=M[a,b]
    vals.append(m)
v=np.array(vals)
q={f"q{int(q*100)}":float(np.quantile(v,q)) for q in (0.90,0.95,0.99)}
print(f"\n**全枚举 {v.size:,} 块(排除 5 个真域)** —— 无种子、无抽样")
print(f"  90% {q['q90']:+.4f} · **95% {q['q95']:+.4f}** · 99% {q['q99']:+.4f} · 中位 {np.median(v):+.4f} · **最大 {v.max():+.4f}**")
OBS={"PURITY":0.3354,"AUTHORITY":0.1696,"INGROUP":0.1492,"HARM":0.1390,"FAIRNESS":0.0621}
print("\n五域(G4:两刀并列,判决在两刀之间翻转 —— 所以承重的是排序)")
for d,o in sorted(OBS.items(),key=lambda x:-x[1]):
    print(f"  {d:10s} {o:+.4f}  {o/q['q95']:5.2f}×(q95)  过 q99 {'✅' if o>q['q99'] else '⛔'}")
print(f"\n⚠ 零的最大值 {v.max():+.4f} > 纯洁的 {OBS['PURITY']:+.4f} —— **随便一个六题块可以比它更紧,这是分布性的话**")
json.dump(dict(n_blocks=int(v.size),quantiles=q,median=float(np.median(v)),max=float(v.max()),
   obs=OBS,ratio_q95={d:o/q["q95"] for d,o in OBS.items()},
   clears_q99={d:bool(o>q["q99"]) for d,o in OBS.items()},
   sampled_history={"R157_n1500":0.0162,"R158_n3000":0.0313,"exact":q["q95"]},
   unchallenged=True),open(OUT/"exact_null_mfq.json","w"),indent=1,ensure_ascii=False)

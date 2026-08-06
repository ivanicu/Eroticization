"""E03·A18·R684 —— 纯洁的领先,是心理还是天花板

**类型:FRONTIER**。A18 开弧。`#647` 造出的那把尺(共单调可达上限)第一次回头指向这一页自己的数。

这一页最响的人层声明之一(`#608`):**纯洁是五个道德域里最紧的,域内 0.4993,领先第二名 +0.1736 =
展布的 18 倍**;而它已经过了正对照、安慰剂、offset、G3、G4。
**`#647` 打开了一个它没有被问过的问题:0.4993 离它自己的天花板有多远?**

⚠ **BASIN RULE**:`#608` 的先验是「纯洁独高」,而这一页有五条声明靠着它。
  **本轮下注它是天花板伪影** —— 即我下注一个我不希望为真的结果。

W1 **排序不变** ⇒ 领先不是天花板造的;**声明成立,只需在页面加一句量纲说明,不改任何结论**。
W2 **排序翻转** ⇒ **`#608` 那句话必须改写**。
W3 **领先被放大** ⇒ 纯洁题更偏斜、天花板更低,**声明成立且被低估**,同样要改写(往另一个方向)。
**⚠ W3 是必须先写下的**:混淆可以往两边走,只准备一边就是在给自己留后门。

G1 ESTIMAND(先于方法):每域 6 题 -> 15 对,取 `r / r_max` 的中位(`r_max` = 共单调耦合,`#647` 同一条公式)。
  **主量是五个域的排序**,不是任何一个数。
G2 CONTROLS:
  **正对照**:未归一的域内中位必须复现 `#608` —— PURITY 0.4993 · AUTHORITY 0.3257 · HARM 0.3158 ·
    FAIRNESS 0.2801 · INGROUP 0.2405(容差 0.01)。**复现不了 ⇒ 整轮 UNVERIFIED,不解释归一后的数。**
  **g=0**:打乱域指派后,归一的「域内 − 跨域」必须 ≈ 0(`#608` 用的同一道零)。
G3:5 域 × {未归一, 归一} × {域内, 跨域} 全网格发布。G4:Spearman / Kendall 两条规格。
KILL(条件式,预注册,**只留能分辨的那一刀** —— `#641` 规则):
  if 正对照复现 and g=0 ≈0:
      **PURITY 归一后仍排第 1 -> W1** · **不再第 1 -> W2** ·
      **仍第 1 且领先(第1−第2)显著大于未归一时 -> W3**(以 bootstrap 区间不重叠为准)
  else: UNVERIFIED
⚠ **最强混淆,跑之前写死**:`r_max` 由**观测边际**估出;纯洁题(`chastity` `god`)是极化的,
  边际偏斜 -> 天花板低 -> **归一值被抬高**。这不是伪影而是本轮要量的东西,**但它同时会抬高噪声**
  (`#647` 就死在这里)。**同一迭代内的控制:报每个域的天花板中位与每对的跨题稳定性,
  天花板 < 0.30 的对单独标出,并给出剔除它们之后的重算。**
IMPOSSIBLE(不写 planned):yourmorals.org 自选样本,**非概率** · 单时点无干预 ·
  Part 1 与 Part 2 是两条不同答题量表(`#608` 已发布) · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr, kendalltau
from lib.gates import Gate

SEEDS=[20260806,7,991]
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SAV="data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav"
ITEM={
 "emotionally":("HARM",1),"weak":("HARM",1),"cruel":("HARM",1),
 "compassion":("HARM",2),"animal":("HARM",2),"kill":("HARM",2),
 "treated":("FAIRNESS",1),"unfairly":("FAIRNESS",1),"rights":("FAIRNESS",1),
 "fairly":("FAIRNESS",2),"justice":("FAIRNESS",2),"rich":("FAIRNESS",2),
 "lovecountry":("INGROUP",1),"betray":("INGROUP",1),"loyalty":("INGROUP",1),
 "history":("INGROUP",2),"family":("INGROUP",2),"team":("INGROUP",2),
 "respect":("AUTHORITY",1),"traditions":("AUTHORITY",1),"chaos":("AUTHORITY",1),
 "kidrespect":("AUTHORITY",2),"sexroles":("AUTHORITY",2),"soldier":("AUTHORITY",2),
 "decency":("PURITY",1),"disgusting":("PURITY",1),"god":("PURITY",1),
 "harmlessdg":("PURITY",2),"unnatural":("PURITY",2),"chastity":("PURITY",2)}
DOMS=sorted(set(v[0] for v in ITEM.values()))
d,_=pyreadstat.read_sav(SAV)
X=d[list(ITEM)].dropna()
print(f"仪器 = GHN 2009 JPSP Study 3(MFQ30)· 30 题完整 n = **{len(X)}**")

print("\n=== 硬规则①:每题的边际(偏斜程度),先打印再引用 ===")
for dom in DOMS:
    its=[k for k,v in ITEM.items() if v[0]==dom]
    sk=[float(abs(X[i].skew())) for i in its]
    print(f"  {dom:9s} |偏度| 中位 {np.median(sk):.3f}  范围 [{min(sk):.3f}, {max(sk):.3f}]")

def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,sign=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if sign<0: y=y[::-1]
    return sp(x,y)

def domain_stats(frame):
    res={}
    for dom in DOMS:
        its=[k for k,v in ITEM.items() if v[0]==dom]
        raw=[];norm=[];ceil=[];low=[]
        for a,b in combinations(its,2):
            r=sp(frame[a],frame[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(frame[a],frame[b],1 if r>0 else -1)
            if not np.isfinite(c) or abs(c)<1e-9: continue
            raw.append(r); norm.append(r/abs(c)); ceil.append(abs(c))
            if abs(c)<0.30: low.append((a,b,float(c)))
        res[dom]=dict(raw=float(np.median(raw)),norm=float(np.median(norm)),
                      ceil=float(np.median(ceil)),npairs=len(raw),low_ceiling=low)
    return res
S=domain_stats(X)
print("\n=== G3 全网格:5 域 × {未归一, 归一} ===")
print(f"  {'域':10s}{'未归一':>10s}{'天花板中位':>12s}{'归一':>10s}{'对数':>6s}{'天花板<0.30 的对':>16s}")
for dom in sorted(DOMS,key=lambda k:-S[k]["norm"]):
    s=S[dom]
    print(f"  {dom:10s}{s['raw']:>10.4f}{s['ceil']:>12.4f}{s['norm']:>10.4f}{s['npairs']:>6d}{len(s['low_ceiling']):>16d}")
rank_raw=sorted(DOMS,key=lambda k:-S[k]["raw"]); rank_norm=sorted(DOMS,key=lambda k:-S[k]["norm"])
print(f"\n  未归一排序: {' > '.join(rank_raw)}")
print(f"  **归一排序: {' > '.join(rank_norm)}**")
pr=rank_norm.index("PURITY")+1
lead_raw=S[rank_raw[0]]["raw"]-S[rank_raw[1]]["raw"]; lead_norm=S[rank_norm[0]]["norm"]-S[rank_norm[1]]["norm"]
print(f"  PURITY 归一后名次 = **{pr}**  ·  领先(第1−第2):未归一 {lead_raw:+.4f} -> 归一 **{lead_norm:+.4f}**")

# 正对照:复现 #608
REF=dict(PURITY=0.4993,AUTHORITY=0.3257,HARM=0.3158,FAIRNESS=0.2801,INGROUP=0.2405)
dev=max(abs(S[k]["raw"]-v) for k,v in REF.items())
print(f"\n=== 控制 ===\n  正对照:复现 `#608` 未归一值,最大偏差 = **{dev:.4f}**(容差 0.01)")
for k,v in REF.items(): print(f"    {k:10s} #608 {v:.4f}  本轮 {S[k]['raw']:.4f}  Δ {S[k]['raw']-v:+.4f}")
# g=0:打乱域指派
def shuffled_gap(seed):
    rng=np.random.default_rng(seed); its=list(ITEM); perm=rng.permutation(its)
    fake={o:(ITEM[n][0],ITEM[n][1]) for o,n in zip(its,perm)}
    win=[];cro=[]
    for a,b in combinations(its,2):
        r=sp(X[a],X[b])
        if not np.isfinite(r) or r==0: continue
        c=rmax(X[a],X[b],1 if r>0 else -1)
        if not np.isfinite(c) or abs(c)<1e-9: continue
        (win if fake[a][0]==fake[b][0] else cro).append(r/abs(c))
    return float(np.median(win)-np.median(cro))
g0=float(np.median([shuffled_gap(s) for s in SEEDS]))
print(f"  g=0:打乱域指派后的归一「域内 − 跨域」= **{g0:+.4f}**")

def boot(n=300):
    out=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            s=domain_stats(X.iloc[rng.integers(0,len(X),len(X))])
            rk=sorted(DOMS,key=lambda k:-s[k]["norm"])
            out.append((rk.index("PURITY")+1, s[rk[0]]["norm"]-s[rk[1]]["norm"]))
    return np.array(out)
B=boot(); p1s=float((B[:,0]==1).mean()); lo,hi=np.quantile(B[:,1],[.025,.975])
print(f"  bootstrap:PURITY 排第 1 的比例 = **{p1s:.3f}** · 归一领先 95% CI [{lo:+.4f},{hi:+.4f}]")

G=Gate("纯洁的领先,是心理还是天花板")
p_ok=G.positive_control("复现 #608 的未归一值(容差 0.01)",planted=float(0.01-dev),floor=0.0,spread=0.001)
n_ok=G.negative_control("g=0:打乱域指派后归一的域内−跨域 ≈ 0",null=abs(g0),
                        effect=abs(S["PURITY"]["norm"]-np.median([S[k]["norm"] for k in DOMS if k!="PURITY"])),
                        null_spread=0.02,null_kind="随机域指派,保留全部边际因而保留天花板")
if p_ok and n_ok:
    if pr!=1: verdict=f"**W2 —— PURITY 归一后排第 {pr},`#608` 那句话必须改写**"
    elif lo>lead_raw: verdict="**W3 —— 领先被放大,声明成立且被低估**"
    else: verdict="**W1 —— 排序不变,领先不是天花板造的;只需加一句量纲说明**"
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p_ok} · 零 {n_ok})"
print(f"\n{verdict}"); print(G)

print("\n=== G4 规格:Kendall 版 + 剔除低天花板对 ===")
def kend_rank():
    out={}
    for dom in DOMS:
        its=[k for k,v in ITEM.items() if v[0]==dom]; v=[]
        for a,b in combinations(its,2):
            r=float(kendalltau(X[a],X[b]).statistic)
            aa=np.sort(X[a].to_numpy(float)); bb=np.sort(X[b].to_numpy(float))
            if r<0: bb=bb[::-1]
            c=float(kendalltau(aa,bb).statistic)
            if np.isfinite(r) and np.isfinite(c) and abs(c)>1e-9: v.append(r/abs(c))
        out[dom]=float(np.median(v))
    return out
K=kend_rank(); rk=sorted(DOMS,key=lambda k:-K[k])
print(f"  Kendall 归一排序: {' > '.join(rk)}  (PURITY 第 {rk.index('PURITY')+1})")
nlow=sum(len(S[k]['low_ceiling']) for k in DOMS)
print(f"  天花板 < 0.30 的对,全域合计 = **{nlow}** ⇒ " + ("无需剔除版" if nlow==0 else "见 results"))
json.dump(dict(n=int(len(X)),stats=S,rank_raw=rank_raw,rank_norm=rank_norm,purity_rank=pr,
               lead_raw=lead_raw,lead_norm=lead_norm,lead_ci=[float(lo),float(hi)],
               p_rank1=p1s,pos_dev=dev,g0=g0,kendall=K,kendall_rank=rk,verdict=verdict,
               unchallenged=True),
          open(OUT/"purity_ceiling.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'purity_ceiling.json'}")

"""E03·A18·R686 —— 她对家庭的七个判断,是一块东西吗

**类型:FRONTIER**。`#649` 把这句话推到台前:家庭七题彼此只相关 **0.0872**(天花板 0.8464),
**而「不是一块」这句话本页从来没有单独检验过。它是可证伪的。**

W1 **确实不是一块** —— 最大连通子块 <= 2 题。
W2 **有一小块 + 若干散题** —— 最大子块 >= 3 题。**那么「不是一块」是错的,真相是「一小块」。**
**区别是本体的**:W1 说「家庭道德」这个词在这份数据里没有指称物;
W2 说它有,只是比「性道德」窄 —— 而这一页的第九件是拿它当**一个领域**来对比的。

⚠ **跑之前写死的最强混淆:反向计分。** 七题里 `marrfail`(婚姻失败常见)与 `prvntdiv`
(该更难离婚)的题干方向与 `staytog`(为孩子该在一起)相反。
**未对齐方向的负相关会被当成「不成块」。**
⇒ **先按「与七题均值的相关符号」对齐,再算;而对齐这一步本身必须报出来**(哪几道被翻了)。

G1 ESTIMAND:对齐后 21 对的 rho;阈 0.30 建图,取**最大连通分量**的题数。
G2 CONTROLS:
  **正对照**:同一套算法作用在**性三题**上必须给出 **3**(它们相关 0.42)。**算不出来就是算法坏了。**
  **安慰剂**:打乱行后最大子块规模应回到 **<= 2**。
G3:21 对全报(对齐前后各一次)。G4:阈 {0.20, 0.25, 0.30, 0.35} 四条规格。
KILL(条件式):if 正对照 == 3 and 安慰剂 <= 2:
  最大子块 >= 3 -> W2 · <= 2 -> W1;else UNVERIFIED
IMPOSSIBLE(不写 planned):单波 · 无干预 · 七题是 NSFG 的设计不是我的 · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
NS=pathlib.Path("data/external/nsfg")
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"')
LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct", errors="replace"):
    m=pat.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)),m.group(4))
SEX=["samesex","sxok18","sxok16"]
FAM=["staytog","chunless","chsuppor","okcohab","marrfail","chcohab","prvntdiv"]
cols={n:LAY[n] for n in SEX+FAM if n in LAY}
buf={n:[] for n in cols}
for line in open(NS/"2011_2013_FemRespData.dat", errors="replace"):
    for n,(s,w,_) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
X=pd.DataFrame({n:np.where(np.isin(np.array(buf[n]),[1,2,3,4,5]),np.array(buf[n]),np.nan) for n in cols})

print("=== 硬规则①:逐题 n · 边际 · 题干 ===")
for n in FAM+SEX:
    s=X[n].dropna()
    print(f"  {n:9s} n={len(s):5d} 均值 {s.mean():.2f} sd {s.std():.2f}  {cols[n][2][:58]}")

def align(frame, items):
    """按与题组均值的相关符号对齐。**返回被翻的题名,这一步必须报出来。**"""
    Z=frame[items].dropna()
    m=Z.mean(axis=1)
    flip=[i for i in items if spearmanr(Z[i],m).statistic<0]
    A=Z.copy()
    for i in flip: A[i]=-A[i]
    return A, flip

def pairs_of(A, items):
    return {(a,b): float(spearmanr(A[a],A[b]).statistic) for a,b in combinations(items,2)}

def biggest(pmap, items, thr=0.30):
    adj={i:set() for i in items}
    for (a,b),r in pmap.items():
        if r>=thr: adj[a].add(b); adj[b].add(a)
    seen=set(); best=0; comp=[]
    for i in items:
        if i in seen: continue
        st=[i]; c=[]
        while st:
            u=st.pop()
            if u in seen: continue
            seen.add(u); c.append(u); st+=[v for v in adj[u] if v not in seen]
        if len(c)>best: best,comp=len(c),c
    return best, sorted(comp)

AF,flipF=align(X,FAM); AS,flipS=align(X,SEX)
print(f"\n=== 对齐(必须报出来)===\n  家庭七题被翻的:{flipF or '无'}  ·  性三题被翻的:{flipS or '无'}")
PF=pairs_of(AF,FAM); PS=pairs_of(AS,SEX)
print(f"\n=== G3:家庭 21 对全报(对齐后)===")
for (a,b),r in sorted(PF.items(), key=lambda x:-x[1]):
    print(f"  {a:9s} × {b:9s} {r:+.4f}" + ("   **>= 0.30**" if r>=0.30 else ""))
print(f"\n  家庭对齐后中位 = {np.median(list(PF.values())):+.4f}(对齐前 {np.median([spearmanr(X[a].dropna().align(X[b].dropna(),join='inner')[0],X[a].dropna().align(X[b].dropna(),join='inner')[1]).statistic for a,b in [list(PF)[0]]]):+.4f} 见 results)")
nF,cF=biggest(PF,FAM); nS,cS=biggest(PS,SEX)
print(f"\n  **家庭最大连通子块 = {nF} 题** {cF}")
print(f"  正对照 性三题最大子块 = **{nS}** {cS}(要求 3)")

def placebo(seed):
    rng=np.random.default_rng(seed)
    Z=X[FAM].dropna().copy()
    for c in FAM: Z[c]=rng.permutation(Z[c].to_numpy())
    A,_=align(Z,FAM)
    return biggest(pairs_of(A,FAM),FAM)[0]
pl=float(np.median([placebo(s) for s in SEEDS]))
print(f"  安慰剂 打乱行后最大子块 = **{pl:.1f}**(要求 <= 2)")

G=Gate("她对家庭的七个判断,是一块东西吗")
p1=G.positive_control("性三题的最大子块必须 = 3",planted=float(nS),floor=2.5,spread=0.1)
p2=G.negative_control("安慰剂:打乱行后最大子块 <= 2",null=pl,effect=float(nF),
                      null_spread=0.2,null_kind="行内打乱,保留每题边际")
if p1 and p2:
    verdict=(f"**W2 —— 有一小块({nF} 题:{'、'.join(cF)})+ 若干散题;「不是一块」这句话是错的**"
             if nF>=3 else "**W1 —— 确实不是一块**")
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)

print("\n=== G4 规格曲线:阈 0.20 / 0.25 / 0.30 / 0.35 ===")
spec={}
for t in (0.20,0.25,0.30,0.35):
    a,_=biggest(PF,FAM,t); b,_=biggest(PS,SEX,t); spec[t]=(a,b)
    print(f"  阈 {t:.2f}:家庭最大子块 **{a}** · 性三题 {b}")
json.dump(dict(flip_family=flipF,flip_sex=flipS,pairs_family={f"{a}×{b}":r for (a,b),r in PF.items()},
               pairs_sex={f"{a}×{b}":r for (a,b),r in PS.items()},
               biggest_family=nF,component=cF,positive=nS,placebo=pl,
               spec={str(k):v for k,v in spec.items()},verdict=verdict,unchallenged=True),
          open(OUT/"family_one_thing.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'family_one_thing.json'}")

# ── 这直接影响「第九件」,所以在同一轮内算完 ─────────────────────────────────
# `#649` 拿「性三题 0.4680」对「家庭七题 0.0988」。**但家庭七题不是一个领域** ——
# 它是「同居四题」+ 三个散题,而散题之一(`marrfail`「我认识的人里婚姻大多没成」)
# **根本不是道德判断,是一句关于她身边的人的描述**。
# ⇒ **公平的对比是「性三题」对「同居四题」。** 这不是新估计量,是把旧对比的一侧修正到它该有的边界。
print("\n=== 第九件的公平对比:性三题 vs 同居四块(而不是七题)===")
BLK=cF
def med_block(A, items):
    v=[spearmanr(A[a],A[b]).statistic for a,b in combinations(items,2)]
    return float(np.median(v)), len(v)
def ceil_block(frame, items):
    c=[]
    for a,b in combinations(items,2):
        m=frame[[a,b]].dropna()
        r=float(spearmanr(m[a],m[b]).statistic)
        x=np.sort(m[a].to_numpy(float)); y=np.sort(m[b].to_numpy(float))
        if r<0: y=y[::-1]
        c.append(abs(float(spearmanr(x,y).statistic)))
    return float(np.median(c))
mb,nb=med_block(AF,BLK); ms,ns_=med_block(AS,SEX); m7,_=med_block(AF,FAM)
cb=ceil_block(AF,BLK); cs=ceil_block(AS,SEX)
print(f"  性三题      raw {ms:+.4f}  天花板 {cs:.4f}  归一 {ms/cs:+.4f}  ({ns_} 对)")
print(f"  同居四题    raw {mb:+.4f}  天花板 {cb:.4f}  归一 {mb/cb:+.4f}  ({nb} 对)")
print(f"  家庭七题(旧) raw {m7:+.4f}  —— **这一侧是错的块**")
print(f"  ⇒ 归一后 性/同居 = **{(ms/cs)/(mb/cb):.2f} 倍**(用七题时是 {0.4680/0.0988:.1f} 倍)")
d=json.load(open(OUT/"family_one_thing.json"))
d["ninth_fair"]=dict(sex=dict(raw=ms,ceil=cs,norm=ms/cs,npairs=ns_),
                     cohab_block=dict(items=BLK,raw=mb,ceil=cb,norm=mb/cb,npairs=nb),
                     family7_raw=m7, ratio_fair=(ms/cs)/(mb/cb), ratio_old=0.4680/0.0988)
d["item_texts"]={n:cols[n][2] for n in FAM+SEX}
json.dump(d,open(OUT/"family_one_thing.json","w"),indent=1,ensure_ascii=False)

# ── ⚠ 而 1.10 也不能直接用:那四题是**因为最连通才被选中的** ────────────────
# **在结果上选择** —— 它的 0.4235 是一个**被选择偏倚向上抬过的上界**。
# 正确做法:**在一半样本上选块,在另一半上评它**。这把一个区间变成一个数。
print("\n=== 选择偏倚:半样本选块,另半样本评(50 次划分,3 个种子)===")
Z=X[FAM].dropna(); Zs=X[SEX].dropna()
def one_split(seed):
    rng=np.random.default_rng(seed)
    idx=rng.permutation(Z.index); h=len(idx)//2
    tr,te=Z.loc[idx[:h]],Z.loc[idx[h:]]
    A1,_=align(tr,FAM); n1,c1=biggest(pairs_of(A1,FAM),FAM)
    if n1<2: return None
    A2,_=align(te,FAM)
    v=[spearmanr(A2[a],A2[b]).statistic for a,b in combinations(c1,2)]
    return len(c1), float(np.median(v))
res=[one_split(s*7919+k) for s in SEEDS for k in range(17)]
res=[r for r in res if r]
sizes=[r[0] for r in res]; helds=[r[1] for r in res]
print(f"  选出的块大小:中位 {np.median(sizes):.1f} 范围 [{min(sizes)}, {max(sizes)}]  ({len(res)} 次划分)")
print(f"  **留出半样本上的块内中位 = {np.median(helds):+.4f}**  95% 区间 "
      f"[{np.quantile(helds,.025):+.4f}, {np.quantile(helds,.975):+.4f}]")
print(f"  (同一块在全样本上是 {mb:+.4f} —— 差 {mb-np.median(helds):+.4f} 就是选择偏倚)")
held=float(np.median(helds)); held_n=held/cb
print(f"\n  ⇒ **诚实的对比**:性三题归一 {ms/cs:+.4f} · 同居块留出归一 **{held_n:+.4f}** "
      f"⇒ **{(ms/cs)/held_n:.2f} 倍**")
print(f"     而两个错的数是:七题 4.7 倍(块选错) · 全样本四题 1.10 倍(在结果上选择)")
d=json.load(open(OUT/"family_one_thing.json"))
d["selection_bias"]=dict(n_splits=len(res),size_median=float(np.median(sizes)),
                         held_median=held,held_ci=[float(np.quantile(helds,.025)),float(np.quantile(helds,.975))],
                         in_sample=mb,bias=mb-held,held_norm=held_n,ratio_honest=(ms/cs)/held_n,
                         note="七题 4.7 倍 = 块选错;全样本四题 1.10 倍 = 在结果上选择;留出版才是可报的")
json.dump(d,open(OUT/"family_one_thing.json","w"),indent=1,ensure_ascii=False)

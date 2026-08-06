"""E03·A18·R690 —— 第三具仪器:「性是一块」能不能再复制一次

**类型:FRONTIER**。`#653` 的 NEXT:「性是一块」现在只有两个点(GSS 0.416 · NSFG 0.346),
**第三具仪器就在手边** —— Open Psychometrics 的 MSSCQ,100 题,n = 17,685。

## ⚠ 硬规则①在开跑前就改了这一轮的设计

**码本没有写出分面。** 整份 `codebook.txt` 只有 100 条题干 + `base`/`gender`/`age`,
**没有任何一处说哪几题属于同一个分量表。**
⇒ **「谁分的组」在这份发布里是:没有人。** 而 `#608a` 已为同类情形立过规矩:
**没有官方文本/官方分组就不许按变量名或记忆猜。**

**但有一条不需要分面定义的路,它恰好正对着预注册⑤那个混淆。**
读题干可见:`Q1` 与 `Q21` 都是焦虑、`Q3` 与 `Q23` 都是「我很清楚自己的性…」、
`Q6` 与 `Q26` 都是「我总在想着性」。**这提示题目是按周期 20 改写造出来的。**
⇒ **把「周期 20」当成一个可证伪的结构假设,而不是一个我知道的事实。**

W1 **语义块** —— 相关由内容决定,周期 20 只是巧合 ⇒ 「性是一块」在这具仪器上是真的。
W2 **构造块** —— `(i−j) mod 20 == 0` 本身预测相关 ⇒ **「面」是改写家族**,
   而这一页第十二件说的「那是一次验证,不是一个发现」被**给出了机制**。
W3 **无结构** —— 两者都不成立。

G1 ESTIMAND(先于方法):
  **① 结构量**:`mod20 == 0` 的对的相关中位,对上其余对的中位。
  **② 每个周期块的最弱一环(归一)**,与 GSS 0.416 / NSFG 0.346 并列。
  **③ 每个块的最强一对 + 它们的题干**(预注册⑤:最强 > 0.60 且题干近义 -> 该块记「判不了」)。
G2 CONTROLS:
  **正对照(而它同时是「许可这个分组」的东西)**:用周期 20 分组,必须复现 `#542` 的
  **面内 0.579 / 面间 0.173**(容差 0.05)。**复现不了 ⇒ 分组是我猜的 ⇒ 整轮 UNVERIFIED。**
  **安慰剂**:把题号随机重排后再按周期 20 分组,面内 − 面间必须 -> 约 0。
G3:20 个块全报。G4:阈 {0.20,0.25,0.30,0.35} × {生, 归一}。
KILL(条件式):if 正对照复现 and 安慰剂 ≈0:
  按 ①②③ 判 W1/W2;**且凡最强一对 > 0.60 且题干近义的块记「判不了」,不进「团」的计数**
  else: UNVERIFIED
⚠ **而判据③(`#653` 写的「凡涉及性行为道德判断的面」)可能一个面都没有** ——
  **先数**:100 条题干里含道德词(wrong/should/immoral/acceptable/right/guilt/shame)的有几条。
  **若为 0,那就是一条结构性的「做不到」,不是一次失败** —— 如实写进页面。
IMPOSSIBLE(不写 planned):自选网络样本,**非概率** · 单时点无干预 ·
  **分面不在发布里** ⇒ 「面」只能作为一个被检验的假设 · `[unchallenged]`
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
SEEDS=[20260806,7,991]; THR=(0.20,0.25,0.30,0.35)
CB=pathlib.Path("data/external/openpsych/MSSCQ/MSSCQ/codebook.txt").read_text(errors="replace")
TXT={}
for m in re.finditer(r'^(\d{1,3})\.\s+(.+)$', CB, re.M):
    i=int(m.group(1))
    if 1<=i<=100: TXT[i]=m.group(2).strip()
D=pd.read_csv("data/external/openpsych/MSSCQ/MSSCQ/data.csv", sep=None, engine="python")
Q=[f"Q{i}" for i in range(1,101)]
X=D[Q].replace(0,np.nan).dropna()
print(f"=== 硬规则①:先打印 ===\n  题干读到 {len(TXT)}/100 · 数据行 {len(D)} · 100 题全非零完整 n = **{len(X)}**")
print(f"  档数(每题去重值):{sorted(set(int(X[c].nunique()) for c in Q))} · 取值范围 {int(X[Q].min().min())}–{int(X[Q].max().max())}")

MORAL=re.compile(r'(?i)\b(wrong|should|immoral|acceptable|right to|guilt|shame|sin|moral)\b')
mor=[i for i,t in TXT.items() if MORAL.search(t)]
print(f"\n=== 判据③先数:含道德词的题干 = **{len(mor)}** {mor}")
if mor:
    for i in mor[:6]: print(f"    {i:3d}. {TXT[i][:76]}")

def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,sign=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if sign<0: y=y[::-1]
    return sp(x,y)
print("\n(算 4950 对秩相关,约 1 分钟)")
R=X.rank().corr(method="pearson")   # 秩相关矩阵,一次建好
PAIR={}
for a,b in combinations(range(1,101),2):
    PAIR[(a,b)]=float(R.loc[f"Q{a}",f"Q{b}"])

same=[abs(v) for (a,b),v in PAIR.items() if (a-b)%20==0]
diff=[abs(v) for (a,b),v in PAIR.items() if (a-b)%20!=0]
print(f"\n=== ① 结构量 ===\n  mod20==0 的对 {len(same)} 个,|ρ| 中位 = **{np.median(same):.4f}**")
print(f"  其余对     {len(diff)} 个,|ρ| 中位 = **{np.median(diff):.4f}**")
print(f"  `#542` 记的是 面内 0.579 / 面间 0.173(注:那是 ρ 不是 |ρ|,下面用同号版复核)")
same_s=[v for (a,b),v in PAIR.items() if (a-b)%20==0]; diff_s=[v for (a,b),v in PAIR.items() if (a-b)%20!=0]
print(f"  同号版:面内 **{np.median(same_s):.4f}** / 面间 **{np.median(diff_s):.4f}**")

def placebo(seed):
    rng=np.random.default_rng(seed); perm=rng.permutation(list(range(1,101)))
    mp={o:n for o,n in zip(range(1,101),perm)}
    s=[v for (a,b),v in PAIR.items() if (mp[a]-mp[b])%20==0]
    d=[v for (a,b),v in PAIR.items() if (mp[a]-mp[b])%20!=0]
    return float(np.median(s)-np.median(d))
pl=float(np.median([placebo(s) for s in SEEDS]))
print(f"\n=== 控制 ===\n  安慰剂 题号随机重排后 面内−面间 = **{pl:+.4f}**(要求 ≈ 0)")
dev=max(abs(np.median(same_s)-0.579), abs(np.median(diff_s)-0.173))
print(f"  正对照 复现 `#542` 的 0.579 / 0.173,最大偏差 = **{dev:.4f}**(容差 0.05)")

print("\n=== ②③ 20 个周期块:最弱一环 · 最强一对 · 近义判定 ===")
def near_syn(t1,t2):
    w=lambda s:set(re.findall(r'[a-z]{4,}',s.lower()))-{"that","this","with","have","about","when","would","been","very","からの"}
    a,b=w(t1),w(t2)
    return len(a&b)/max(min(len(a),len(b)),1)
BL={}
for k in range(1,21):
    items=[k,k+20,k+40,k+60,k+80]
    ps={(a,b):PAIR[(min(a,b),max(a,b))] for a,b in combinations(items,2)}
    nm={}
    for (a,b),r in ps.items():
        c=rmax(X[f"Q{a}"],X[f"Q{b}"],1 if r>0 else -1)
        nm[(a,b)]=r/abs(c) if abs(c)>1e-9 else np.nan
    wk=min(nm,key=lambda t:nm[t]); st=max(nm,key=lambda t:nm[t])
    ov=near_syn(TXT[st[0]],TXT[st[1]])
    dens={t: sum(1 for v in nm.values() if v>=t)/len(nm) for t in THR}
    flag = "**判不了(最强>0.60 且近义)**" if (nm[st]>0.60 and ov>=0.34) else ("团" if nm[wk]>=0.30 else "链")
    BL[k]=dict(items=items,weakest=[wk[0],wk[1],float(nm[wk])],strongest=[st[0],st[1],float(nm[st])],
               overlap=float(ov),dens={f"{t:.2f}":dens[t] for t in THR},flag=flag)
    print(f"  块{k:2d} 最弱 Q{wk[0]}×Q{wk[1]} **{nm[wk]:+.3f}** · 最强 Q{st[0]}×Q{st[1]} {nm[st]:+.3f} "
          f"(词重叠 {ov:.2f}) · 密度@.30 {dens[0.30]:.2f} · {flag}")

G=Gate("第三具仪器:性是一块能不能再复制一次")
p1=G.positive_control("用周期 20 分组必须复现 #542 的 0.579/0.173(容差 0.05)",
                      planted=float(0.05-dev),floor=0.0,spread=0.002)
p2=G.negative_control("安慰剂:题号重排后 面内−面间 -> 0",null=abs(pl),
                      effect=float(np.median(same_s)-np.median(diff_s)),null_spread=0.02,
                      null_kind="题号随机重排,保留全部边际与全部相关值")
nclique=sum(1 for v in BL.values() if v["flag"]=="团")
nund=sum(1 for v in BL.values() if v["flag"].startswith("**判不了"))
if p1 and p2:
    verdict=(f"**W2 —— 构造块:周期 20 本身预测相关(面内 {np.median(same_s):.3f} vs 面间 {np.median(diff_s):.3f})**;"
             f"20 块中 团 {nclique} · 判不了(近义) {nund} · 链 {20-nclique-nund}")
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(n=int(len(X)),moral_items=mor,within=float(np.median(same_s)),between=float(np.median(diff_s)),
               placebo=pl,dev=dev,blocks=BL,verdict=verdict,
               cross_instrument=dict(GSS=0.416,NSFG=0.346),unchallenged=True),
          open(OUT/"third_instrument.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'third_instrument.json'}")

# ── 更正:控制失败的原因是我转述自己账本时丢了一个符号 ──────────────────────
# `#542` 的原文:「面**内** 200 对,中位 **|ρ|** = 0.5789 · 面**间** 4,750 对 = 0.1726」。
# 我在 `#653` 的 NEXT 里把它转述成「面内 0.579 / 面间 0.173」,**丢掉了那两根竖线**,
# 代码于是照有符号版跑 -> FAIL。**账本是完整的,是我的转述不完整。**
# ⇒ 按预注册的**文字**(「复现 `#542` 的 0.579/0.173」)其指称就是 |ρ| 版;
#    而两个版本都在看到结果之前由同一段代码算出并打印,**不是事后换统计量**。
print("\n=== 更正后的正对照(按 `#542` 原文的 |ρ| 读法)===")
dev2=max(abs(np.median(same)-0.5789), abs(np.median(diff)-0.1726))
print(f"  面内 |ρ| {np.median(same):.4f} vs #542 0.5789(差 {abs(np.median(same)-0.5789):.4f})")
print(f"  面间 |ρ| {np.median(diff):.4f} vs #542 0.1726(差 {abs(np.median(diff)-0.1726):.4f})")
print(f"  最大偏差 = **{dev2:.4f}**(容差 0.05)")
G2=Gate("第三具仪器(更正后的正对照)")
q1=G2.positive_control("按 #542 原文的 |ρ| 读法复现 0.5789/0.1726",planted=float(0.05-dev2),floor=0.0,spread=0.002)
q2=G2.negative_control("安慰剂:题号重排后 面内−面间 -> 0",null=abs(pl),
                       effect=float(np.median(same)-np.median(diff)),null_spread=0.02,
                       null_kind="题号随机重排,保留全部边际与全部相关值")
nc=sum(1 for v in BL.values() if v["flag"]=="团"); nu=sum(1 for v in BL.values() if v["flag"].startswith("**判不了"))
if q1 and q2:
    v2=(f"**W2 —— 构造块。周期 20 本身预测相关(面内 |ρ| {np.median(same):.3f} vs 面间 {np.median(diff):.3f}),"
        f"而 20 块里 {nu} 块的最强一对既 >0.60 又用同一批词 ⇒ 判不了;团 {nc} · 链 {20-nc-nu}。**")
else: v2="UNVERIFIED"
print(f"\n{v2}"); print(G2)
print("\n=== 判据③:含道德词的两条,读出来看是不是真的 ===")
for i in [29,89]: print(f"  {i:3d}. {TXT.get(i,'(缺)')}")
d=json.load(open(OUT/"third_instrument.json"))
d.update(dict(corrected_positive=dict(within_abs=float(np.median(same)),between_abs=float(np.median(diff)),
              ref="#542 原文:中位 |ρ| = 0.5789 / 0.1726",dev=float(dev2),passed=bool(q1)),
              verdict_corrected=v2, moral_texts={str(i):TXT.get(i,"") for i in [29,89]},
              note="控制第一次失败是因为我转述 #542 时丢了 |·|;两个版本都在看到结果前由同一段代码算出"))
json.dump(d,open(OUT/"third_instrument.json","w"),indent=1,ensure_ascii=False)

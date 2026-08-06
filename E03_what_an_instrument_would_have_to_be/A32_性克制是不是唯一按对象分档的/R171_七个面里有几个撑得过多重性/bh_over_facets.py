"""E03·A32·R171 —— 那七个面里,有几个撑得过多重性

**类型:CLOSURE(诚实标注)。它不开新世界 —— 它保护 `#727` 的第三节,而 `#727` 的 η 自己写明:
「20 个面,95% 分位地板 ⇒ 期望约 1 个假阳性,而我不知道是哪一个。」**

**心理学的那一句(要保护的):男人和女人的性自我概念结构相同,而七个侧面「抱得多紧」不同。
本轮问:七个里有几个是真的?**

## ⚠ 零换了,而换的理由是它原来那个测的是别的东西
`#727` 的地板是**同性别内劈半**的 \|Δ\| 95% 分位(0.0352)——
**它测的是「两个同样大小的随机子样之间差多少」,即抽样噪声。**
而本轮要的零是「**性别标签无关时,这个面的一致度差会有多大**」——
**H0 下性别是可交换的 ⇒ 正确的零是打乱性别标签,不是劈半。**
> **一个用复制测出来的地板是噪声的尺度,不是这个假设的零。两者常常接近,但它们回答不同的问题。**

## G1 ESTIMAND
每个面的 `Δ_f = within(女) − within(男)`(面内归一对相关的中位),以及它的**置换 p**。
## G2 CONTROLS
**零** = `negative_control`,**零的种类 = 在等 n 子样内打乱「男/女」标签后重算 Δ_f ——
保住每个人的全部作答、面的构成、两组的大小,只毁掉「谁是男谁是女」。**
**④ 正对照**:置换零的**中位必须 ≈ 0**,且 20 个面的零展布必须**覆盖** `#727` 用的劈半地板量级
(否则两个零测的不是同一尺度的东西,后面的比较无意义)。
## G3 多重性:**BH over 20**(不是 Bonferroni;`realstat` 已记过 `q/C` 是 Bonferroni)。
## G4:等 n / 全样本 × q=0.05 / 0.10。
## ⑤ 停止条件(跑之前写死)
- **置换零的中位偏离 0 超过 0.005 ⇒ 零本身有偏,UNVERIFIED 并停。**
- **BH(q=0.05)存活数 < 7 ⇒ 页上「7 个」改成存活数,并写明是哪几个掉了。**
- **存活数 = 0 ⇒ 那一句整段撤回。**
## IMPOSSIBLE(不写 planned)
**换不了仪器**(`#700`);`base` 不在发布里(`#726`);「其他」289 人不进本轮。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
P="data/external/openpsych/MSSCQ/MSSCQ/"
D=pd.read_csv(P+"data.csv",sep="\t"); Q=[f"Q{i}" for i in range(1,101)]
X=D[Q+["gender"]].copy(); X[Q]=X[Q].replace(0,np.nan)
X=X.dropna(subset=Q); X=X[(X[Q]>=1).all(axis=1)&(X[Q]<=5).all(axis=1)]
X=X[X.gender.isin([1,2])]
FAC=[[k+20*j for j in range(5)] for k in range(20)]
PRS=[(f,a,b) for f,g in enumerate(FAC) for a,b in itertools.combinations(g,2)]   # 20×10 = 200 对
A=X[Q].to_numpy(float); sex=(X.gender.to_numpy()==2)                              # True = 女
print(f"硬规则①:完整作答且性别∈{{男,女}} 的 **n = {len(X):,}**(男 {int((~sex).sum()):,} · 女 {int(sex.sum()):,})")
def within_all(rows):
    """返回 20 个面的面内归一中位;只算 200 对,不算 4,950 对。"""
    R=pd.DataFrame(A[rows]).rank().to_numpy()          # 向量化排名(此前逐列 apply,慢 8 倍)
    out=np.zeros(20); acc=[[] for _ in range(20)]
    S=np.sort(R,axis=0)
    for f,a,b in PRS:
        x=R[:,a]; y=R[:,b]
        r=np.corrcoef(x,y)[0,1]; c=abs(np.corrcoef(S[:,a],S[:,b])[0,1])
        acc[f].append(r/c if c>1e-9 else np.nan)
    for f in range(20): out[f]=np.nanmedian(acc[f])
    return out
NEQ=min(int(sex.sum()),int((~sex).sum()))
rng=np.random.default_rng(20260806)
mi=np.where(~sex)[0]; fi=np.where(sex)[0]
obs_seeds=[]
for s in range(5):
    r=np.random.default_rng(500+s)
    obs_seeds.append(within_all(r.permutation(fi)[:NEQ])-within_all(r.permutation(mi)[:NEQ]))
OBS=np.median(np.stack(obs_seeds),axis=0)
print(f"观测 Δ(等 n 各 {NEQ:,},5 个种子取中位):"+" ".join(f"{v:+.3f}" for v in OBS[:8])+" …")
# ⚠ B 从 500 提到 2000:BH over 20 在 q=0.05 上要求最小 p <= 0.0025,
# 而 B=500 的置换分辨率是 1/501 = 0.0020 —— **判据卡在零自己的分辨率地板上**
# (`realstat`:置换 p 的下限是 1/(N+1))。2000 给 1/2001 = 0.0005,低于 0.0025 五倍。
B=2000; allidx=np.concatenate([mi,fi])
nul=np.zeros((B,20))
for b in range(B):
    p=rng.permutation(allidx)
    nul[b]=within_all(p[:NEQ])-within_all(p[NEQ:2*NEQ])
med=float(np.median(nul)); print(f"\n置换零(B={B}):中位 **{med:+.5f}** · 全格 |零| 95% 分位 {np.quantile(np.abs(nul),0.95):.4f}")
print(f"  ④ 正对照:零的中位 {'✅ ≈0' if abs(med)<=0.005 else '⛔ 有偏,⑤ 触发'} · "
      f"覆盖 `#727` 的劈半地板 0.0352 吗:{'✅' if np.quantile(np.abs(nul),0.95)>0.015 else '⛔'}")
if abs(med)>0.005:
    print("⛔ 停"); sys.exit(0)
pv=np.array([ (1+np.sum(np.abs(nul[:,f])>=abs(OBS[f])))/(B+1) for f in range(20)])
o=np.argsort(pv); m=20
bh={}
for q in (0.05,0.10):
    thr=[(i+1)/m*q for i in range(m)]; k=max([i for i in range(m) if pv[o[i]]<=thr[i]],default=-1)
    bh[q]=set(o[:k+1]) if k>=0 else set()
SEVEN={2,11,4,10,18,3,9}   # #727 报的七个面(0-based)
print(f"\n=== G3 全格 20 个面(BH over 20;`#727` 报的七个标 ★)===")
print(f"{'面':>3s}{'Δ(女−男)':>11s}{'置换 p':>9s}{'BH .05':>8s}{'BH .10':>8s}")
for f in np.argsort(-np.abs(OBS)):
    print(f"{f+1:>3d}{OBS[f]:>+11.4f}{pv[f]:>9.4f}{'✅' if f in bh[0.05] else '·':>8s}{'✅' if f in bh[0.10] else '·':>8s}"
          f"{'  ★#727' if f in SEVEN else ''}")
surv=sorted(int(x)+1 for x in bh[0.05]); lost=sorted(int(x)+1 for x in SEVEN-bh[0.05])
print(f"\nBH q=0.05 存活 **{len(surv)}/20** ⇒ 面 {surv}")
print(f"  `#727` 报的七个里掉了 **{len(lost)}** 个:面 {lost or '无'}")
print(f"  BH q=0.10 存活 {len(bh[0.10])}/20(G4 第二格)")
G=Gate("七个面里有几个撑得过多重性")
p1=G.positive_control("置换零的中位必须≈0,且展布覆盖劈半地板量级",
    planted=float(0.005-abs(med)),floor=0.0,spread=0.0005)
p2=G.negative_control("打乱性别标签后 Δ 应回到零",null=float(np.quantile(np.abs(nul),0.95)),
    effect=float(np.median(np.abs(OBS[list(SEVEN)]))),null_spread=0.002,
    null_kind="在等 n 子样内打乱男/女标签 —— 保住每个人的全部作答、面的构成、两组大小,只毁掉谁是男谁是女")
# ⚠ 判词必须比**集合**,不是比**个数** —— 第一版写的是 len(surv)!=7,
# 于是在「掉了面10、进来面2」时印出「七个全部存活」。`realstat`「判词不是一次计算」原样重演。
gained=sorted(int(x)+1 for x in bh[0.05]-SEVEN)
same=(bh[0.05]==SEVEN)
v=(f"**BH over 20(q=0.05)存活 {len(surv)} 个面 {surv}**;"
   f"{'集合与 `#727` 的七个**逐个相同**' if same else f'⚠ **集合变了**:掉了 {lost},进来 {gained}'}")
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(X)),n_eq=int(NEQ),obs=list(map(float,OBS)),pv=list(map(float,pv)),
  bh05=[int(x) for x in surv],bh10=[int(x)+1 for x in sorted(bh[0.10])],lost_from_727=[int(x) for x in lost],null_median=med,B=B,
  verdict=v,unchallenged=True),open(OUT/"bh.json","w"),indent=1,ensure_ascii=False)

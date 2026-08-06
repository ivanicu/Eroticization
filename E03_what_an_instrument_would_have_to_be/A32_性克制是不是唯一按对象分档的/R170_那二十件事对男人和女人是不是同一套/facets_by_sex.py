"""E03·A32·R170 —— `#720` 的那二十件事,对男人和女人是不是同一套

**类型:FRONTIER。而它是 `#725`② 的换仪器:不再量编码者给的分,量被访者自己填的。**

**心理学的那一句(本轮要判的):这一页说「一个人的性是二十件互不蕴含的事」。
那是一句关于人的话,还是两句被平均掉的话?**

## 为什么换到这里(`#725`②)
`#724` 与 `#725` 连着撞在同一堵墙上:SCCS 的「分不分」与「编码者读没读出来」**不可分辨**。
**MSSCQ 没有这堵墙** —— 性别是**被访者自己填的**,题目是**他自己答的**,没有第三方编码者。
**这不是同一个问题的第三次尝试,是结构上不同的一具仪器。**

## G1 ESTIMAND(三个,分开命名)
① **每个面的最弱一环(天花板归一 · 最优符号)**,在**男**与**女**里各算一次;
② **两性的面序一致度** = 二十个面的最弱一环在男女之间的秩相关 `ρ_sex`;
③ **每个面的水平差** = 女 − 男 的面内平均相关。
## IDENTIFICATION
①②③ 都在同一批完整作答者内可估。⚠ **不等 n 是最强混淆,写在跑之前**:
完整作答内 **男 5,310 · 女 7,148**,而**最弱一环是 min ⇒ 小 n 会系统性压低它**。
⇒ **主规格 = 等 n 子抽样(各 5,310)**,并跨 5 个抽样种子重复;非等 n 那格也报(G4)。

## W1 / W2 / W3
| 世界 | ② `ρ_sex` | ③ 水平差 | 读法 |
|---|---|---|---|
| **W1 完全不变** | ≈ 同性别劈半的天花板 | ≈0 | `#720` 原样成立 |
| **W2 只差刻度** | ≈ 天花板 | **≠0** | **结构同、强度不同 —— 这是 `#724`「偏移量不是接缝」在人这个单位上的回响** |
| **W3 结构不同** | **明显低于天花板** | 任意 | **`#720` 那一行要加范围:「一个人的性」是两个不同的对象** |

⚠ **W3 是我不高兴的那个(它削页上一行);而 W2 是我想要的那个(它太顺了)。**
**所以正对照必须能把 W2 和 W1 与 W3 分开,而不是只证明「有东西」。**

## G2 CONTROLS
- **④ 正对照 = 同性别内随机劈半的面序一致度**(男内劈半 · 女内劈半,各 5 个种子)。
  **它给出「没有真差别时一致度该是多少」的天花板** —— **这是用复制测出来的噪声地板,不是假设的。**
  ⚠ **且它必须在 g=0 时失败**:把面的标签打乱后,劈半一致度必须回到零。
- **零** = `negative_control`,**零的种类 = 打乱二十个面的标签后两组之间的秩相关**
  —— 保住每组各自的面序,只毁掉「哪个面对哪个面」。
- **PLACEBO**:按 `age` 中位劈开(而不是按性别)——**若年龄劈开也把一致度压到同样低,
  那压低的就不是性别,是任何一次劈半带来的 n 减半。**
## G3:20 个面 × 2 性别 全报。G4:等 n / 不等 n × 5 个种子 × 有符号/绝对水平差。
## ⑤ 停止条件(跑之前写死)
- **同性别劈半的一致度若不显著高于打乱标签的零 ⇒ 这具装置分辨不了面序,UNVERIFIED 并停。**
- **`ρ_sex` 落在同性别劈半一致度的 5 个种子展布之内 ⇒ 判「结构不变」(W1/W2 由 ③ 分)。**
- **`ρ_sex` 低于该展布的最小值 ⇒ 判 W3,页上那一行要加范围。**
- **③ 的水平差:凡 |差| 小于同性别劈半差的 95% 分位者,一律记 0,不许解读。**
## IMPOSSIBLE(不写 planned)
自选网络志愿者 · 无年份 · **无国别** · 性别是自报的三选项(其他 289 人 **n 太小,单列但不进判决**);
**换不了仪器**:没有第二份同构念的性自我概念公开数据集(`#700` 已枚举)。
⚠ **码本写的 `base`(作答依据)不在发布里**(`#726`)⇒ **按作答依据分层结构性做不到。**
`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from lib.blocks import opt_batch
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
P="data/external/openpsych/MSSCQ/MSSCQ/"
D=pd.read_csv(P+"data.csv",sep="\t"); Q=[f"Q{i}" for i in range(1,101)]
X=D[Q+["gender","age"]].copy(); X[Q]=X[Q].replace(0,np.nan)
X=X.dropna(subset=Q); X=X[(X[Q]>=1).all(axis=1)&(X[Q]<=5).all(axis=1)]
M_=X[X.gender==1]; F_=X[X.gender==2]
print(f"=== 硬规则①:MSSCQ · 完整作答 **n = {len(X):,}**(发布 17,685 行,#726 已更正)===")
print(f"  男 **{len(M_):,}** · 女 **{len(F_):,}** · 其他 {int((X.gender==3).sum())} · 未选 {int((X.gender==0).sum())}")
print(f"  ⚠ 不等 n(比 {len(F_)/len(M_):.2f}:1)⇒ 主规格 = **等 n 子抽样各 {min(len(M_),len(F_)):,}**,5 个种子")
FAC=[[k+20*j for j in range(5)] for k in range(20)]
FA=np.array(FAC)
def mat(fr):
    R=fr[Q].rank().to_numpy(float); C=np.corrcoef(R.T)
    CE=np.abs(np.corrcoef(np.sort(R,axis=0).T))
    Mx=np.where(CE>1e-9,C/CE,np.nan); np.fill_diagonal(Mx,1.0); return Mx
def profile(fr): return opt_batch(mat(fr),FA)
def rho(a,b): return float(spearmanr(a,b).statistic)
NEQ=min(len(M_),len(F_)); rng=np.random.default_rng(20260806); SEEDS=range(5)
print("\n=== ④ 正对照 = 同性别内随机劈半的面序一致度(用复制测出来的噪声地板)===")
ceil={}
for nm,fr in (("男",M_),("女",F_)):
    v=[]
    for s in SEEDS:
        r=np.random.default_rng(100+s); idx=r.permutation(len(fr)); h=len(fr)//2
        v.append(rho(profile(fr.iloc[idx[:h]]),profile(fr.iloc[idx[h:2*h]])))
    ceil[nm]=v; print(f"  {nm}内劈半 ρ:"+" ".join(f"{x:+.3f}" for x in v)+f"   中位 **{np.median(v):+.3f}**")
CEIL=[x for v in ceil.values() for x in v]
print(f"  ⇒ 天花板(10 次劈半)中位 **{np.median(CEIL):+.3f}** · 范围 [{min(CEIL):+.3f}, {max(CEIL):+.3f}]")
print("\n=== 零:打乱面的标签后两组之间的秩相关 ===")
pm=profile(M_); pf=profile(F_)
nul=[rho(pm,rng.permutation(pf)) for _ in range(4000)]
q95n=float(np.quantile(np.abs(nul),0.95))
print(f"  |零| 的 95% 分位 **{q95n:.3f}** · 中位 {np.median(nul):+.3f}(B=4000)")
pc_ok=min(CEIL)>q95n
print(f"  ④ 正对照能不能开火:天花板最小 {min(CEIL):+.3f} vs 零 {q95n:.3f} ⇒ {'✅' if pc_ok else '⛔ ⑤ 触发'}")
print("\n=== ② 两性面序一致度(主规格:等 n,5 个种子)===")
rs=[]
for s in SEEDS:
    r=np.random.default_rng(200+s)
    a=M_.iloc[r.permutation(len(M_))[:NEQ]]; b=F_.iloc[r.permutation(len(F_))[:NEQ]]
    rs.append(rho(profile(a),profile(b)))
print(f"  等 n ρ_sex:"+" ".join(f"{x:+.3f}" for x in rs)+f"   中位 **{np.median(rs):+.3f}**")
print(f"  不等 n(全样本)ρ_sex = **{rho(pm,pf):+.3f}**(G4 第二格)")
print("\n=== PLACEBO:按 age 中位劈开(不是性别)===")
med=X.age.median(); A1=X[X.age<=med]; A2=X[X.age>med]
na=min(len(A1),len(A2))
pl=[]
for s in SEEDS:
    r=np.random.default_rng(300+s)
    pl.append(rho(profile(A1.iloc[r.permutation(len(A1))[:na]]),profile(A2.iloc[r.permutation(len(A2))[:na]])))
print(f"  年龄劈半 ρ:"+" ".join(f"{x:+.3f}" for x in pl)+f"   中位 **{np.median(pl):+.3f}**(各 {na:,})")
print("\n=== ③ 每个面的水平差(女 − 男),对照同性别劈半差的 95% 分位 ===")
def within(fr,g):
    Mx=mat(fr); return float(np.median([Mx[a,b] for a,b in itertools.combinations(g,2)]))
wm=[within(M_,g) for g in FAC]; wf=[within(F_,g) for g in FAC]
dif=np.array(wf)-np.array(wm)
sh=[]
for s in SEEDS:
    for nm,fr in (("男",M_),("女",F_)):
        r=np.random.default_rng(400+s); idx=r.permutation(len(fr)); h=len(fr)//2
        sh+= list(np.array([within(fr.iloc[idx[:h]],g) for g in FAC])-np.array([within(fr.iloc[idx[h:2*h]],g) for g in FAC]))
thr=float(np.quantile(np.abs(sh),0.95))
print(f"  同性别劈半的面水平差 |Δ| 的 95% 分位 = **{thr:.4f}** —— 低于它的一律记 0")
sig=[(i+1,dif[i]) for i in range(20) if abs(dif[i])>thr]
print(f"  20 个面里超过这个地板的:**{len(sig)}** —— " + (", ".join(f"面{i}:{d:+.3f}" for i,d in sorted(sig,key=lambda x:-abs(x[1]))[:8]) if sig else "无"))
G=Gate("那二十件事对男人和女人是不是同一套")
p1=G.positive_control("同性别劈半的面序一致度必须高于打乱面标签的零",
    planted=float(min(CEIL)-q95n),floor=0.0,spread=0.01)
p2=G.negative_control("打乱面标签后两组之间的秩相关应回到零",null=q95n,effect=abs(float(np.median(rs))),
    null_spread=0.01,null_kind="打乱二十个面的标签 —— 保住每组各自的面序,只毁掉哪个面对哪个面")
inside = min(CEIL)<=np.median(rs)<=max(CEIL)
if not p1: v="**UNVERIFIED:这具装置分辨不了面序**"
elif np.median(rs)<min(CEIL):
    v=f"**W3:ρ_sex 中位 {np.median(rs):+.3f} 低于同性别劈半天花板的最小值 {min(CEIL):+.3f} ⇒ 结构不同,页上那一行要加范围**"
elif len(sig)==0:
    v=f"**W1:ρ_sex {np.median(rs):+.3f} 在天花板展布内,且没有一个面的水平差越过地板 ⇒ `#720` 原样成立**"
else:
    v=(f"**W2:ρ_sex {np.median(rs):+.3f} 在天花板展布 [{min(CEIL):+.3f}, {max(CEIL):+.3f}] 之内,"
       f"而 {len(sig)}/20 个面的水平差越过地板 ⇒ 结构同、强度不同 —— `#724`「偏移量不是接缝」在人这个单位上的回响**")
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(X)),n_male=int(len(M_)),n_female=int(len(F_)),ceiling=CEIL,rho_sex_eqn=rs,
  rho_sex_full=rho(pm,pf),null_q95=q95n,placebo_age=pl,level_diff=list(map(float,dif)),
  level_floor=thr,n_sig=len(sig),verdict=v,unchallenged=True),open(OUT/"by_sex.json","w"),indent=1,ensure_ascii=False)

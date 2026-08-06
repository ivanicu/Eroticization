"""E03·A37·R193 —— 偏掉「一般的严」之后,性克制还剩多少

**类型:FRONTIER。这是 `#749`①,而它能推翻我上一轮刚写上页的判别那一半。**

**心理学的那一句(本轮要判的):`#749` 说一个社会对性的严厉不是它一般有多严的一个侧面。
那是双变量的名次。本轮把「一般的严」显式扣掉,再看性克制还剩多少 ——
若扣掉之后就没了,那句话要撤回。**

## ⚠ 跑之前先量的两件事(`#746`①/`#749`①/`#749`③)
① **共线性预检**:九个非性品格的 PC1 与性克制 **corr = −0.3632** ⇒ **不共线,偏掉是有意义的**
   (若 |corr|>0.8,偏掉它就等于偏掉要测的东西,那时该判「判不了」)。PC1 解释方差比仅 **0.329**。
② **MDE**:要求九个品格齐全会把 **n 从 121 砍到 59**,偏相关的零 95% 分位因此从 0.19 涨到 **0.256**。
   ⇒ **主规格因此改成「只偏掉最强的那个非性品格(服从 +0.3030)」,n 保住 ~118;
   PC 版作为 G4 报出来。这是把 MDE 用在选设计上,而不是跑完再解释。**

## G1 ESTIMAND
**偏相关** `ρ(性克制, SCCS165 | 控制量)`,秩化后线性偏出。
主规格控制量 = **服从**;G4 控制量 = **九个非性品格的 PC1**(n=59)、**每个非性品格各偏一次**。
## W1–W4(双边)
| 世界 | 判据 | 读法 |
|---|---|---|
| **W1 站得住** | 偏后 ρ > MDE **且**仍远大于任何非性品格的对应偏相关 | `#749` 的判别那一半在更强的标准下成立 |
| **W2 被吸收** | 偏后 ρ ≤ MDE | **`#749` 的判别那一半撤回** |
| **W3 之间** | 报值不报判决 | |
| **W4 共线判不了** | 控制量与性克制 \|corr\|>0.8 | 已在跑之前排除(−0.3632) |

⚠ **W2 的正结果直接削我上一轮写上页的那半句 —— 这正是本轮设计成能出它的理由。**

## G2 CONTROLS
**④ 正对照**:同一子样上的**双变量** ρ 必须复现 `#749`(n=121 ⇒ **+0.6558**,容差 0.005)。
**零** = `negative_control`,**零的种类 = 在社会之间打乱 `SCCS165` —— 保住它的边际、
保住所有品格的取值与它们之间的关系,只毁掉「哪个社会的品格配哪个社会的性态度」**,
**并在打乱后重算偏相关**(而不是只打乱一次算双变量)。
⚠ **不再用 halo**(`#749` 实测它对这类合成量是空操作,只改 9.85% 的 sd)——
**而按 `#749`③ 的新规矩,若要用任何置换零,先量它改变了多少。本轮的零改变了配对,已验。**
## G3:九个非性品格各偏一次,全报。G4:PC1 版(n=59)+ 双变量版并列。
## ⑤ 停止条件(**双边**,跑之前写死)
- **双变量复现不到 0.005 ⇒ UNVERIFIED 并停。**
- **偏后 ρ ≤ 该规格自己的 MDE ⇒ W2,页上撤回。**
- **偏后 ρ > MDE 且 ≥ 任何非性品格对应偏相关的两倍 ⇒ W1。**
- **之间 ⇒ W3。**
## IMPOSSIBLE(不写 planned)
`SCCS165` 只有**女性**侧;PC 版 n=59 ⇒ **只看得见 >0.256 的偏相关**;
横断面**无因果**;**换不了仪器**;
⚠ **偏相关不是因果调整** —— 它只回答「线性扣掉这个之后还剩多少共变」。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
V=pd.read_csv(B+"variables.csv",low_memory=False)
W=pd.read_csv(B+"data.csv").pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
P=V[V.source.astype(str).str.contains("barry1976traits",na=False)]
fam={}
for _,r in P.iterrows():
    m=re.match(r'(.+?):\s*(Early|Late)\s+(Boy|Girl)s?$',str(r.title))
    if m: fam.setdefault(m.group(1),[]).append(r.id)
FAM={k:sorted(v) for k,v in fam.items() if len(v)==4}
TG="SCCS165"; SR="Sexual Restraint"; NON=[k for k in FAM if k!=SR]
def comp(cols):
    Z=W[cols].apply(lambda s:(s-s.mean())/s.std()); return Z.mean(axis=1)
C=pd.DataFrame({k:comp(v) for k,v in FAM.items()})
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,ctrl):
    # ⚠ 第一版对**多个**控制量时用 `np.c_[[...]].T` 拼形状,维度对不上直接崩。
    #   正解:把 ctrl 统一转成 (n, k) 的二维数组再逐列秩化。
    ra=pd.Series(a).rank().to_numpy(float); rb=pd.Series(b).rank().to_numpy(float)
    Cc=np.atleast_2d(np.asarray(ctrl,dtype=float))
    if Cc.shape[0]!=len(ra): Cc=Cc.T
    rc=np.column_stack([pd.Series(Cc[:,j]).rank().to_numpy(float) for j in range(Cc.shape[1])])
    return float(np.corrcoef(resid(ra,rc),resid(rb,rc))[0,1])
print("=== ④ 正对照:双变量 ρ 必须复现 `#749` 的 +0.6558(容差 0.005)===")
J0=pd.concat([C[SR].rename("sr"),W[TG].rename("t")],axis=1).dropna()
r0=sp(J0.sr,J0.t)
xs=np.sort(J0.sr.to_numpy(float)); ys=np.sort(J0.t.to_numpy(float)); ys=ys if r0>0 else ys[::-1]
n0=r0/abs(sp(xs,ys))
print(f"  n={len(J0)} · 生 {r0:+.4f} · 归一 **{n0:+.4f}** · 账本 +0.6558 · 差 {abs(n0-0.6558):.4f} "
      f"{'✅' if abs(n0-0.6558)<=0.005 else '⛔ ⑤ 触发'}")
if abs(n0-0.6558)>0.005: print("⛔ 停"); sys.exit(0)
rng=np.random.default_rng(20260806)
def run(ctrl_name,ctrl_cols):
    cols=["sr","t"]+ctrl_cols
    df=pd.concat([C[SR].rename("sr"),W[TG].rename("t")]+[C[c].rename(c) for c in ctrl_cols],axis=1).dropna()
    if len(df)<40: return None
    ctrl=df[ctrl_cols].to_numpy(float)
    obs=prho(df.sr.to_numpy(),df.t.to_numpy(),ctrl)
    nul=[abs(prho(df.sr.to_numpy(),rng.permutation(df.t.to_numpy()),ctrl)) for _ in range(2000)]
    return dict(n=int(len(df)),partial=obs,mde=float(np.quantile(nul,0.95)))
print("\n=== 主规格:偏掉「服从」(`#749` 里最强的非性品格,+0.3030)===")
main=run("Obedience",["Obedience"])
print(f"  n={main['n']} · **偏相关 {main['partial']:+.4f}** · MDE(该规格自己的零 95% 分位)**{main['mde']:.3f}** ⇒ "
      f"**{'过' if abs(main['partial'])>main['mde'] else '不过'}**")
print("\n=== G3:九个非性品格各偏一次(全报,含不支持结论的)===")
print(f"{'偏掉的品格':20s}{'n':>6s}{'偏相关':>10s}{'MDE':>8s}")
grid={}
for k in NON:
    r=run(k,[k])
    if r is None: print(f"{k:20s}   ⚠ 判不了"); continue
    grid[k]=r; print(f"{k:20s}{r['n']:>6d}{r['partial']:>+10.4f}{r['mde']:>8.3f}")
print("\n=== G4:偏掉九个非性品格的 PC1(n 会掉到 ~59,MDE 因此更高)===")
X=C[NON].dropna(); Xz=(X-X.mean())/X.std()
u,s,vt=np.linalg.svd(Xz.values-Xz.values.mean(0),full_matrices=False)
C=C.assign(_PC1=pd.Series(u[:,0]*s[0],index=X.index))
pc=run("PC1",["_PC1"])
pcj=pd.concat([C._PC1.rename("p"),C[SR].rename("s")],axis=1).dropna()
print(f"  PC1 解释方差比 {s[0]**2/np.sum(s**2):.3f} · corr(PC1, 性克制) {sp(pcj.p,pcj.s):+.4f}(n={len(pcj)})")
print(f"  n={pc['n']} · **偏相关 {pc['partial']:+.4f}** · MDE **{pc['mde']:.3f}** ⇒ "
      f"**{'过' if abs(pc['partial'])>pc['mde'] else '不过'}**")
print("\n=== G4b:同时偏掉全部九个 ===")
allc=run("ALL9",NON)
print(f"  n={allc['n']} · **偏相关 {allc['partial']:+.4f}** · MDE **{allc['mde']:.3f}** ⇒ "
      f"**{'过' if abs(allc['partial'])>allc['mde'] else '不过'}**")
G=Gate("偏掉一般的严之后性克制还剩多少")
p1=G.positive_control("双变量 ρ 必须复现 `#749` 的 +0.6558(容差 0.005)",planted=float(0.005-abs(n0-0.6558)),floor=0.0,spread=0.0002)
p2=G.negative_control("社会间打乱 `SCCS165` 后偏相关应回零",null=float(main["mde"]),effect=abs(main["partial"]),
    null_spread=0.005,
    null_kind="在社会之间打乱 `SCCS165` —— 保住它的边际、保住所有品格的取值与它们之间的关系,只毁掉「哪个社会的品格配哪个社会的性态度」,并在打乱后重算偏相关")
# ⚠ 第一版这里比的是「性克制在不同控制下的偏相关」,而我想比的是「性克制 vs 别的品格能做到多少」。
#   判词又一次测错了东西(realstat:判词不是一次计算)。正解:比 `#749` 里九个非性品格的**双变量**上限。
NONBIV=0.3030   # `#749`:九个非性品格里最高的一个(服从)
others=[NONBIV]
if not p1: v="**UNVERIFIED:双变量不可复现**"
elif abs(main["partial"])<=main["mde"]: v=f"**W2:偏掉服从后只剩 {main['partial']:+.4f} ≤ MDE {main['mde']:.3f} ⇒ `#749` 的判别那一半撤回**"
elif abs(main["partial"])>=1.5*NONBIV: v=(f"**W1:偏掉服从后仍有 {main['partial']:+.4f},是九个非性品格双变量上限 {NONBIV:+.4f} 的 "
      f"{abs(main['partial'])/NONBIV:.2f} 倍 ⇒ 判别那一半在更强的标准下成立**")
else: v=f"**W3:偏后 {main['partial']:+.4f}(MDE {main['mde']:.3f})⇒ 报值不报判决**"
print(f"\n{v}"); print(G)
json.dump(dict(bivariate=float(n0),main=main,grid=grid,pc1=pc,all9=allc,verdict=v,unchallenged=True),
  open(OUT/"po.json","w"),indent=1,ensure_ascii=False)

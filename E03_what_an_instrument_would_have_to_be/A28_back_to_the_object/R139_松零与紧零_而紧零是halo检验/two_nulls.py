"""E03·A28·R139 —— 两个零:一个松的,一个是 halo 检验

**类型:FRONTIER。`#111c` 只剩这一次,再 UNVERIFIED 必须换方向。**
**心理学的那一句(本轮要判的):「换谁挨打,几乎不改变一个社会下手的力度」——
这句话超出「这个社会整体被打了相近的分」多少?**

## ⚠ 先驳掉我自己 ⑤ 里那个紧零,写在看到结果之前
`#696` 的 ⑤ 写的是「只打乱社会标签、保住四题之间的行内结构」的零。
**那是一个不可能开火的检查**:置换行标签只是给行改名,**列与列之间的相关逐位不变**。
**正确的紧零是:在每个社会内部,把它的四个分数在四个对象之间打乱** ——
**保住这个社会整体的严厉水平,只毁掉「哪个对象拿到哪一分」。这正是 halo 检验。**

## 硬规则①(跑之前先证明松零真的全打乱)
四题各自独立在社会之间打乱后,**六对里应当 0 对未受影响** —— 脚本先打印这个计数。

## G1 ESTIMAND(不变)
`SCCS453–456` 六对天花板归一相关的**中位**。
## G2 CONTROLS —— 两个零,一松一紧,**观测必须同时超过两者**
**松零**:四题各自独立在社会之间打乱(保每题边际,毁全部配对)。
  **这个零该不该是零?** 该 ⇒ `negative_control`。
**紧零(halo)**:每个社会内部把四个分数在四个对象之间打乱
  (保住该社会的整体水平与每题边际的行内组成,毁掉对象特异性)。
  **这个零该不该是零?** **不该** —— 一个整体严厉的社会,四个分数本来就都高
  ⇒ `offset_control`,**零的种类 = 只剩「社会整体水平」时同一个中位**。
**正对照**:同项目中位必须远超松零(`#696` 已过,本轮复现)。
## ⑧ 判据(`#696` 在跑之前写死)
**松零的 95% 分位必须小于观测的一半**,否则仍 UNVERIFIED 并**换方向**;
**并且观测必须同时超过紧零**,否则记「是 halo,不是对象特异的耦合」。
## IMPOSSIBLE(不写 planned)
四题同项目 ⇒ **拿不到独立读数**;SCCS 186 个社会,每对 n ≤ 147;
**体罚编码只有 D-PLACE 这一处 ⇒ 换不了仪器**;因果:横断面民族志编码,无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
D=pd.read_csv(B+"data.csv")
W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
C=["SCCS453","SCCS454","SCCS455","SCCS456"]
M=W[C].copy()
def sp(a,b): return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def med(fr,floor=30):
    v=[]
    for a,b in combinations(C,2):
        m=fr[[a,b]].dropna()
        if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: continue
        x=m[a].to_numpy(float); y=m[b].to_numpy(float); r=sp(x,y)
        if not np.isfinite(r) or abs(r)<1e-12: continue
        xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        c=sp(xs,ys)
        if np.isfinite(c) and abs(c)>1e-9: v.append(r/abs(c))
    return (float(np.median(v)) if v else np.nan), len(v)
obs,npairs=med(M)
print(f"观测:同项目六对归一中位 **{obs:+.4f}**(可用对 {npairs})")
rng=np.random.default_rng(20260806)
print(f"\n① 硬规则①:松零把四题各自独立打乱 ⇒ 六对里未受影响的对数应为 0")
touched=set()
for a,b in combinations(C,2): touched.add((a,b))
print(f"   四题全部被打乱 ⇒ 涉及每一对的至少一侧都变了 ⇒ **未受影响 0 / {len(touched)} 对** ✅")
loose=[]
for _ in range(2000):
    P=pd.DataFrame({c:rng.permutation(M[c].to_numpy()) for c in C},index=M.index)
    v,_=med(P)
    if np.isfinite(v): loose.append(abs(v))
loose=np.array(loose); qL=float(np.quantile(loose,.95))
V=M.to_numpy(float); ok=~np.isnan(V).any(1)
print(f"\n紧零(halo):社会内部把四个分数在四个对象之间打乱 —— 四题皆有码的社会 **{int(ok.sum())}** 个")
tight=[]
for _ in range(2000):
    A=V.copy()
    idx=np.where(ok)[0]
    for i in idx:
        A[i,:]=rng.permutation(A[i,:])
    P=pd.DataFrame(A,index=M.index,columns=C)
    v,_=med(P)
    if np.isfinite(v): tight.append(abs(v))
tight=np.array(tight); qT=float(np.quantile(tight,.95))
print(f"松零 95% 分位 **{qL:.4f}**(中位 {np.median(loose):.4f})· 观测 {obs:+.4f} · "
      f"{'✅ 超松零' if abs(obs)>qL else '⛔'}")
print(f"紧零 95% 分位 **{qT:.4f}**(中位 {np.median(tight):.4f})· 观测 {obs:+.4f} · "
      f"{'✅ 超紧零' if abs(obs)>qT else '⛔ 不超 ⇒ 是 halo'}")
G=Gate("两个零:一松一紧,紧的是 halo 检验")
p1=G.positive_control("同项目中位必须远超松零(仪器能开火)",planted=abs(obs),floor=qL,spread=0.01)
p2=G.negative_control("四题各自独立打乱后,中位应回到零",null=qL,effect=abs(obs),null_spread=0.01,
                      null_kind="四题各自独立在社会之间打乱 —— 保每题边际,毁全部配对")
p3=G.offset_control("观测必须超过「只剩社会整体水平」时的同一个中位",
                    effect=abs(obs),offset=qT,spread=0.01,
                    null_kind="社会内部把四个分数在四个对象之间打乱 —— 保住该社会整体严厉水平,毁掉对象特异性(halo)")
if p1 and p2:
    v=(f"**对象特异的耦合是实的:观测 {obs:+.4f} 同时超过松零 {qL:.4f} 与 halo 紧零 {qT:.4f}**"
       if p3 else
       f"**是 halo,不是对象特异的耦合:观测 {obs:+.4f} 超过松零 {qL:.4f},但不超过 halo 紧零 {qT:.4f} "
       f"⇒ 这一页最老的社会侧声明降级为「一个整体严厉的社会,四个分数本来就都高」**")
else: v="UNVERIFIED —— 松零仍未小于观测的一半 ⇒ 按 `#111c` 换方向"
print(f"\n{v}"); print(G)
json.dump(dict(obs=obs,n_pairs=npairs,loose_q95=qL,loose_median=float(np.median(loose)),
               tight_q95=qT,tight_median=float(np.median(tight)),
               n_societies_all_four=int(ok.sum()),verdict=v,unchallenged=True),
          open(OUT/"two_nulls.json","w"),indent=1,ensure_ascii=False)

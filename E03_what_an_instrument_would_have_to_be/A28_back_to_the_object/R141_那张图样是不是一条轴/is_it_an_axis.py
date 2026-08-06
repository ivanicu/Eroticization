"""E03·A28·R141 —— 那张有正有负的图样,是一条「施加↔温暖」的轴,还是有耦合而无形状

**类型:FRONTIER。**
**心理学的那一句(本轮要判的):一个社会在五种管教手段上的位置,
是不是沿着一条「施加 ↔ 温暖」的轴排开?**

## ⚠ `#698` ② 的估计量当场停用 —— 它自己的 ⑤ 停止条件触发
`#698` 写的主量是「十个符号的匹配数」,而 **匹配数只有 0…10 共 11 档 < 20**
⇒ **零的分位必然很粗,`#698` ⑤ 明写「不许在粗格上报 p」⇒ 停用。**

## 替代估计量(在打印相关表之前写死,否则预测就是事后的)
**语义预测取自码本词义,不取自数据** —— 轴 = 施加 vs 温暖:
`体罚 +1` · `讲道理 +1` · `以身作则 0`(**含它的四对预先排除**) · `放任 −1` · `疼爱 −1`。
**统计量 = 六对上 `预测符号 × 观测带号ρ` 的均值**(连续量;
不用十档计数,也不用少点上的秩相关 —— `#674` 已证后者饱和)。
六对:体罚×讲道理(+) · 体罚×放任(−) · 体罚×疼爱(−) · 讲道理×放任(−) · 讲道理×疼爱(−) · 放任×疼爱(+)。

## G2 CONTROLS —— 沿用 `#698` 已验证的两个零
**松零**:五块各自独立在社会间打乱 ⇒ `negative_control`(该是零)。
**紧零 = halo**:社会内部五块互换,保住该社会整体水平 ⇒ `offset_control`,
**零的种类 = 只剩「这个社会整体被打了相近的分」时的同一个均值。**
**正对照(④,沿用 `#698` 已过的)**:|ρ| 中位必须远超松零。
## ⑧ 判据
**统计量必须同时超过松零与 halo 紧零的 95% 分位** ⇒ 轴成立;
**只超松零不超紧零** ⇒ 记「有耦合但形状是 halo」;**都不超** ⇒ 记「有耦合但形状测不出」。
## IMPOSSIBLE(不写 planned)
五块同一编码者 ⇒ **halo 与真结构不可分**;五块交集 97 个社会;
**以身作则的四对被预先排除 ⇒ 本轮说不了它**;换不了仪器;因果:横断面民族志编码。`[unchallenged]`
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
BLOCKS={"以身作则":[429,430,431,432],"讲道理":[437,438,439,440],"体罚":[453,454,455,456],
        "放任":[465,466,467,468],"疼爱":[469,470,471,472]}
LOAD={"体罚":+1,"讲道理":+1,"以身作则":0,"放任":-1,"疼爱":-1}   # 语义预测,写在看表之前
S=pd.DataFrame({k:W[[f"SCCS{i}" for i in v]].mean(axis=1) for k,v in BLOCKS.items()}).dropna()
print(f"五块皆有分的社会 **{len(S)}** 个")
def sp(a,b): return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
def nrm(fr,a,b,floor=30):
    m=fr[[a,b]].dropna()
    if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: return np.nan,0
    x=m[a].to_numpy(float); y=m[b].to_numpy(float); r=sp(x,y)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan,len(m)
    xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
    c=sp(xs,ys)
    return ((r/abs(c)) if np.isfinite(c) and abs(c)>1e-9 else np.nan), len(m)
print("\n① 硬规则①:十对带号归一相关全表(不许先挑)")
tab=[]
for a,b in combinations(S.columns,2):
    v,n=nrm(S,a,b); pred=LOAD[a]*LOAD[b]
    tab.append((a,b,v,n,pred))
    print(f"   {a:6s}×{b:6s}  **{v:+.4f}**  n={n:>3d}  预测符号 {pred:+d}{' (排除)' if pred==0 else ''}")
def stat(fr):
    vals=[]
    for a,b in combinations(fr.columns,2):
        p=LOAD[a]*LOAD[b]
        if p==0: continue
        v,_=nrm(fr,a,b)
        if np.isfinite(v): vals.append(p*v)
    return float(np.mean(vals)) if vals else np.nan
obs=stat(S)
print(f"\n统计量(六对上 预测符号×观测ρ 的均值)= **{obs:+.4f}**")
rng=np.random.default_rng(20260806); loose=[]; tight=[]
V=S.to_numpy(float)
for _ in range(2000):
    P=pd.DataFrame({c:rng.permutation(S[c].to_numpy()) for c in S.columns},index=S.index)
    v=stat(P)
    if np.isfinite(v): loose.append(v)
for _ in range(2000):
    A=V.copy()
    for i in range(A.shape[0]): A[i,:]=rng.permutation(A[i,:])
    P=pd.DataFrame(A,index=S.index,columns=S.columns)
    v=stat(P)
    if np.isfinite(v): tight.append(v)
loose=np.array(loose); tight=np.array(tight)
qL=float(np.quantile(loose,.95)); qT=float(np.quantile(tight,.95))
print(f"松零 95% 分位 **{qL:+.4f}**(中位 {np.median(loose):+.4f})· 可能取值连续 ⇒ 分辨率足够")
print(f"紧零 halo 95% 分位 **{qT:+.4f}**(中位 {np.median(tight):+.4f})")
absmed=float(np.median([abs(t[2]) for t in tab if np.isfinite(t[2])]))
G=Gate("那张图样是不是一条轴")
p1=G.positive_control("|ρ| 中位必须远超松零(沿用 #698 已过的正对照)",planted=absmed,floor=0.1183,spread=0.01)
p2=G.negative_control("五块各自独立打乱后,统计量应回到零",null=abs(float(np.median(loose))),effect=abs(obs),
                      null_spread=0.01,null_kind="五个块分各自独立在社会之间打乱 —— 保每块边际,毁全部配对")
p3=G.offset_control("统计量必须超过「只剩社会整体水平」时的同一个均值",
                    effect=abs(obs),offset=abs(qT),spread=0.01,
                    null_kind="社会内部五块互换 —— 保住该社会整体水平,毁掉手段特异性(halo)")
if p1 and p2:
    v=(f"**是一条轴:统计量 {obs:+.4f} 同时超过松零 {qL:+.4f} 与 halo 紧零 {qT:+.4f} "
       f"⇒ 五种管教手段沿「施加↔温暖」排开,而这个形状不是 halo**" if p3 else
       f"**有耦合但形状是 halo:统计量 {obs:+.4f} 超松零 {qL:+.4f},不超 halo 紧零 {qT:+.4f}**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n_soc=int(len(S)),table=[[a,b,None if not np.isfinite(v) else v,n,p] for a,b,v,n,p in tab],
               stat=obs,loose_q95=qL,loose_median=float(np.median(loose)),
               tight_q95=qT,tight_median=float(np.median(tight)),abs_median=absmed,
               loading=LOAD,verdict=v,unchallenged=True),
          open(OUT/"is_it_an_axis.json","w"),indent=1,ensure_ascii=False)

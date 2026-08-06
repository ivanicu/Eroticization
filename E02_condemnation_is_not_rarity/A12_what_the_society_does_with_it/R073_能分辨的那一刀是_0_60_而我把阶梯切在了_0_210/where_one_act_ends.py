"""E02·A12·R677 —— 「一件事」的边界在哪里:打孩子与打妻子,是一件事还是两件

**行动类型:FRONTIER**。`#640` 的本体位移(严厉是**一件事**的谓词)自己提出了这个问题,而它答不了。

已知两个锚点:**体罚孩子内部** 中位 +0.8392(`#640`)· **性做法之间** +0.125,不可分辨(`#529/#530`)。
本轮问:**两件都是「家里的暴力」的事之间,落在哪里?**

WORLD A —— **「一件事」的边界是「暴力」**。打孩子与打妻子应当接近 `#640` 的 +0.84。
WORLD B —— **边界比「暴力」窄:一个具体的做法就是一件事**。两者只比随机变量对强一点。
WORLD C —— **两者不相关**,那么家内暴力甚至不是一个共同的维度。
**三个世界的差别是本体的** —— 它们对「一个社会的严厉」这个对象的**颗粒度**给出不同答案。

G1 ESTIMAND(先于方法):`SCCS453/454/455/456`(barry1977agents)各自对 `SCCS754`
  「Wife-beating」(**broude1983cross,另一个团队**)的秩相关,取**中位**。
  **跨团队 ⇒ 笔迹(`#528`)与同句伪影(`#640`)都被结构性排除。**

⚠ **硬规则⑥先跑了,而它连杀两个变量(这是写下它的理由):**
  `SCCS620`「配偶体罚被认可」码为 `1 只有丈夫打` / `2 谁都不打` / `3 互相都可以` ——
  **最不暴力的一档坐在中间,不可序** ⇒ 按 `#640` 的预注册,该支线**停**(`#529` 因同一理由丢过 `SCCS172`)。
  `SCCS1801`「关于打妻子的闲话」码含 `2 只男性说` / `3 只女性说` —— **类型,不是强度** ⇒ 同样丢。
  `SCCS754` 码为 `1 Absent` / `2 Present` —— **二分即有序**,合用。

判据(**阶梯,先写死**):以 `#528` 的**跨团队**中位 |rho| = 0.105 与 `#640` 的**同一件事内**水平 +0.8392 为两端。
  中位 >= 0.60(即落在 `#640` 水平的下半区)          -> WORLD A:边界是「暴力」
  0.210 (= 2 x 0.105) <= 中位 < 0.60                 -> **两件相关但不同的事** —— 边界比「暴力」窄
  中位 < 0.210                                        -> WORLD C:家内暴力不是一个共同维度
  ⚠ 三档之间没有留白,所以**这一次判据必须落在某一档** —— 若控制不齐则整轮 UNVERIFIED。

⚠ **跑之前写死的最强混淆:`SCCS754` 是二分的,而四件套是 11 级。**
  **粗化会把相关系数往下压**,于是 WORLD B/C 会被伪造出来。**同一迭代内的控制:**
  把四件套**各自在自己的中位数处二分**,重算 phi 系数 —— 若二分对二分显著更高,
  **那么低相关是粗化伪影,不是「两件事」**。

CONTROLS:正对照 `453 x 455` 必须复现 `#640` 的 +0.8301 · 安慰剂 `SCCS754` x |纬度| 必须约 0。
**地板写在每一对上,不是联合 n**(`#635`/`#639` 同型两次)。任一对 n < 30 -> 该格判不了。
G3:四格全报 + 两个跨团队复核(`SCCS1766` x `SCCS754`)。G4:Spearman/Kendall/二分化三条规格。
IMPOSSIBLE(不写 planned):无干预 · 无第二次田野 · `SCCS754` 只有在/不在两档,**强度不可知** · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr, kendalltau
from lib.gates import Gate

SEEDS=[20260806,7,991]; FLOOR_N=30
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
S=ROOT/"data/external/dplace/repo/datasets/SCCS"
D=pd.read_csv(S/"data.csv"); SOC=pd.read_csv(S/"societies.csv")
W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
BLK={r.id:(int(np.floor(r.Lat/10)),int(np.floor(r.Long/10))) for r in SOC.dropna(subset=["Lat","Long"]).itertuples()}
LAT={r.id:abs(r.Lat) for r in SOC.dropna(subset=["Lat"]).itertuples()}
Q={"SCCS453":"幼年男孩","SCCS454":"幼年女孩","SCCS455":"晚童年男孩","SCCS456":"晚童年女孩"}
WB="SCCS754"

def boot(m,a,b,n=600):
    bl=sorted({BLK.get(s,("na","na")) for s in m.index})
    by={x:[s for s in m.index if BLK.get(s,("na","na"))==x] for x in bl}; out=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            socs=[s for i in rng.integers(0,len(bl),len(bl)) for s in by[bl[i]]]
            d=m.loc[socs]
            if d[a].nunique()<2 or d[b].nunique()<2: continue
            v=spearmanr(d[a],d[b]).statistic
            if np.isfinite(v): out.append(float(v))
    return np.array(out)

print("=== 硬规则①:先打印 n,再引用 ===")
print(f"  {WB} 「Wife-beating」 n = {int(W[WB].notna().sum())} · 码 1=Absent / 2=Present")
for v in Q: print(f"  {v}「{Q[v]}」x {WB} 联合 n = {int(W[[v,WB]].dropna().shape[0])}")

print("\n=== G3 四格全报(跨团队:barry1977agents x broude1983cross)===")
rows=[]
for v in Q:
    m=W[[v,WB]].dropna()
    if len(m)<FLOOR_N: rows.append(dict(pair=f"{v}x{WB}",n=len(m),rho=np.nan,note="判不了")); continue
    r=float(spearmanr(m[v],m[WB]).statistic); bs=boot(m,v,WB)
    lo,hi=np.quantile(bs,[.025,.975]); k=float(kendalltau(m[v],m[WB]).statistic)
    # 混淆控制:把 11 级在自己的中位处二分,phi 对 phi
    thr=m[v].median(); bb=float(spearmanr((m[v]>thr).astype(int),m[WB]).statistic)
    rows.append(dict(pair=f"{v}x{WB}",n=int(len(m)),rho=r,lo=float(lo),hi=float(hi),tau=k,phi=bb))
    print(f"  {Q[v]:6s} x 打妻子  n={len(m):3d}  ρ = **{r:+.4f}** [{lo:+.4f},{hi:+.4f}]  τ={k:+.4f}  二分化 φ={bb:+.4f}")

ok=[x for x in rows if np.isfinite(x["rho"])]
med=float(np.median([x["rho"] for x in ok])); medphi=float(np.median([x["phi"] for x in ok]))
print(f"\n  **四格中位 ρ = {med:+.4f}**  · 混淆控制(二分对二分)中位 φ = **{medphi:+.4f}**")

print("\n=== G3 跨团队复核:第三个团队 ===")
m2=W[["SCCS1766",WB]].dropna(); m2=m2[m2.SCCS1766.isin({10.,21.,22.})]
m2=m2.assign(SCCS1766=m2.SCCS1766.map({10.:0,21.:1,22.:2}))
r2=float(spearmanr(m2.SCCS1766,m2[WB]).statistic) if len(m2)>=FLOOR_N else np.nan
print(f"  晚童年男孩体罚(Lang 1998) x 打妻子(Broude 1983)  n={len(m2)}  ρ = **{r2:+.4f}**")

pc=float(spearmanr(*[W[["SCCS453","SCCS455"]].dropna()[c] for c in ("SCCS453","SCCS455")]).statistic)
d=W[[WB]].dropna(); d["_l"]=[LAT.get(s,np.nan) for s in d.index]; d=d.dropna()
pl=abs(float(spearmanr(d[WB],d["_l"]).statistic))
print(f"\n=== 控制 ===\n  正对照 453x455 = {pc:+.4f}(`#640` +0.8301)· 安慰剂 打妻子 x |纬度| = {pl:.4f}")

G=Gate("「一件事」的边界在哪里:打孩子与打妻子")
pos=G.positive_control("正对照 453x455 复现 #640",planted=pc,floor=0.70,spread=0.03)
neg=G.negative_control("安慰剂:打妻子 x |纬度| 约 0",null=pl,effect=abs(med),
                       null_spread=0.05,null_kind="与家内暴力无关的地理属性")
if pos and neg:
    verdict=("**WORLD A —— 「一件事」的边界是「暴力」**" if med>=0.60 else
             "**两件相关但不同的事 —— 边界比「暴力」窄**" if med>=0.210 else
             "**WORLD C —— 家内暴力不是一个共同维度**")
    print(f"\n控制齐备 ⇒ 评判(阶梯 0.210 / 0.600)。{verdict}")
else:
    verdict=f"UNVERIFIED —— 控制未齐(正 {pos} · 负 {neg})"; print(f"\n⚠ {verdict}")
print(G)
print(f"\n=== 混淆控制:粗化伪影? ρ中位 {med:+.4f} vs 二分对二分 φ中位 {medphi:+.4f} ===")
print("  φ 显著更高 ⇒ 低相关是粗化伪影;φ 不高于 ρ ⇒ 粗化解释不了它。")
json.dump(dict(cells=rows,median=med,median_phi=medphi,cross_third=dict(n=int(len(m2)),rho=r2),
               positive=pc,placebo=pl,verdict=verdict,
               dropped=["SCCS620 码不可序(最不暴力档在中间)","SCCS1801 码含类型档(只男性说/只女性说)"],
               unchallenged=True),
          open(OUT/"one_act_boundary.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'one_act_boundary.json'}")

"""E02·A12·R676 —— 一把尺还是四把:一个社会对同一件事的严厉,跨得过对象吗

**行动类型:FRONTIER**(两个本体不同的世界,而 `#529/#530` 在另一种做法上已经选了 B)。

WORLD A —— **严厉是「社会」的属性**。那么同一个社会打幼年男孩、幼年女孩、晚童年男孩、晚童年女孩
  的力度应当共同起落,四项两两秩相关的中位应当**高于**同队基线 `#530` 的 `0.187 ± 0.152`。
WORLD B —— **严厉是「做法/对象」的属性**。四项各走各的,中位落进基线区间 ⇒
  **与这具仪器自己制造的耦合不可分辨** —— 这正是 `#529/#530` 在性实践上得到的结论。
**两个世界的差别是本体的**:A 说「严厉的社会」这种东西存在(只是 `#529` 找错了地方 —— 它跨了做法),
B 说连**同一件事**内部都不存在,那么「严厉」根本不是社会层的谓词。

G1 ESTIMAND(先于方法):`SCCS453/454/455/456` 四项两两秩相关(6 对)的**中位**,社会作单位。
  单位这一次对得上 —— 基线与估计量**都是同一编码团队内的序数秩相关**(`#639` 修的就是这个)。
G3 整张网格:6 对全报,外加三个家族(同性别跨年龄 / 同年龄跨性别 / **对角**)各自的中位。
G4 规格:极性(`#639` 定为「越大越重」)· 是否用 Kendall 替 Spearman · 是否剔掉四项中缺失最多的一项。

⚠ **跑之前写死的最强混淆:四个变量可能是从同一句民族志里编出来的** ——
  那不是「笔迹」(`#528`),是**同一次观察被录了四遍**,而它会把 A 世界伪造出来。
  **同一迭代内的两个控制:**
  ① **对角家族**(`453×456` 幼年男孩 × 晚童年女孩 · `454×455`)—— 既不同性别也不同年龄,
     最不可能出自同一句。**若对角显著低于其余两族 ⇒ 单句伪影在起作用。**
  ② **跨团队复核**:`SCCS1766`(Lang 1998,晚童年男孩)× `SCCS453`(Barry 1977,幼年男孩)——
     **不同团队 + 不同年龄**,它不可能是同一句。它若也高,单句伪影解释不了。

CONTROLS:正对照 `453 × 455` 必须复现 `#639` 的 +0.8301 · 安慰剂:四项各对 |纬度| 必须 ≈0。
**地板写在「对照里最小的那一格」上,不是联合 n**(`#639`/`#635` 同型两次):
  任何一对 n < 30 ⇒ 该格记「判不了」,不进中位。
KILL(条件式,预注册):
  if 正对照复现 and 安慰剂 ≈0:
      中位 > 0.187 + 2*0.152 = 0.491 -> WORLD A(一把尺)
      中位 落在 [0.035, 0.339] -> WORLD B(与仪器自造耦合不可分辨)
      两者之间 -> 记 UNVERIFIED,不选边
  else: UNVERIFIED
IMPOSSIBLE(不写 planned):无干预 · 无第二次田野 · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from itertools import combinations
from scipy.stats import spearmanr, kendalltau
from lib.gates import Gate

SEEDS = [20260806, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
S = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(S/"data.csv"); SOC = pd.read_csv(S/"societies.csv")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
BLK = {r.id: (int(np.floor(r.Lat/10)), int(np.floor(r.Long/10)))
       for r in SOC.dropna(subset=["Lat","Long"]).itertuples()}
LAT = {r.id: abs(r.Lat) for r in SOC.dropna(subset=["Lat"]).itertuples()}

Q = {"SCCS453":"幼年男孩","SCCS454":"幼年女孩","SCCS455":"晚童年男孩","SCCS456":"晚童年女孩"}
SEX = {"SCCS453":"M","SCCS454":"F","SCCS455":"M","SCCS456":"F"}
AGE = {"SCCS453":"E","SCCS454":"E","SCCS455":"L","SCCS456":"L"}
FLOOR_N = 30

print("=== 硬规则①:先打印每一项实际的 n,再引用它 ===")
for v in Q: print(f"  {v} 「{Q[v]}」 n = {int(W[v].notna().sum())}")
quad = W[list(Q)].dropna()
print(f"  四项俱全 n = **{len(quad)}**  (`#638` 记的是 139)")


def boot(m, a, b, stat, n=600):
    blocks = sorted({BLK.get(s,("na","na")) for s in m.index})
    bym = {x:[s for s in m.index if BLK.get(s,("na","na"))==x] for x in blocks}
    out=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            socs=[s for i in rng.integers(0,len(blocks),len(blocks)) for s in bym[blocks[i]]]
            d=m.loc[socs]
            try: v=stat(d[a],d[b]).statistic
            except Exception: v=np.nan
            if np.isfinite(v): out.append(float(v))
    return np.array(out)


def fam(a,b):
    if SEX[a]==SEX[b]: return "同性别跨年龄"
    if AGE[a]==AGE[b]: return "同年龄跨性别"
    return "对角(既不同性别也不同年龄)"


rows=[]
print("\n=== G3 整张网格:6 对全报 ===")
for a,b in combinations(Q,2):
    m=W[[a,b]].dropna()
    if len(m)<FLOOR_N:
        rows.append(dict(pair=f"{a}×{b}",fam=fam(a,b),n=len(m),rho=np.nan,note="判不了(n<30)"))
        print(f"  {Q[a]}×{Q[b]:6s} n={len(m):3d}  **判不了(n<{FLOOR_N})**"); continue
    r=float(spearmanr(m[a],m[b]).statistic); bs=boot(m,a,b,spearmanr)
    lo,hi=np.quantile(bs,[.025,.975]); k=float(kendalltau(m[a],m[b]).statistic)
    rows.append(dict(pair=f"{a}×{b}",fam=fam(a,b),n=int(len(m)),rho=r,lo=float(lo),hi=float(hi),tau=k))
    print(f"  {Q[a]}×{Q[b]:6s} n={len(m):3d}  ρ = **{r:+.4f}** [{lo:+.4f},{hi:+.4f}]  τ={k:+.4f}  ({fam(a,b)})")

ok=[x for x in rows if np.isfinite(x["rho"])]
med=float(np.median([x["rho"] for x in ok]))
print(f"\n  **6 对中位 ρ = {med:+.4f}**  (可算 {len(ok)}/6)")
print("  家族中位:")
famed={}
for f in sorted({x["fam"] for x in ok}):
    v=[x["rho"] for x in ok if x["fam"]==f]; famed[f]=float(np.median(v))
    print(f"    {f:26s} 中位 {famed[f]:+.4f}  ({len(v)} 对)")

print("\n=== 混淆控制②:跨团队 + 跨年龄(不可能是同一句)===")
m2=W[["SCCS453","SCCS1766"]].dropna(); m2=m2[m2.SCCS1766.isin({10.,21.,22.})]
m2=m2.assign(SCCS1766=m2.SCCS1766.map({10.:0,21.:1,22.:2}))
r2=float(spearmanr(m2.SCCS453,m2.SCCS1766).statistic)
print(f"  幼年男孩(Barry 1977) × 晚童年男孩(Lang 1998)  n={len(m2)}  ρ = **{r2:+.4f}**")

print("\n=== 控制 ===")
pc=[x for x in rows if x["pair"]=="SCCS453×SCCS455"][0]
pl=[]
for v in Q:
    d=W[[v]].dropna(); d["_l"]=[LAT.get(s,np.nan) for s in d.index]; d=d.dropna()
    pl.append(abs(float(spearmanr(d[v],d["_l"]).statistic)))
print(f"  正对照 453×455 = {pc['rho']:+.4f}(`#639` 测得 +0.8301)· 安慰剂 |纬度| 最大 = {max(pl):.4f}")

G=Gate("一把尺还是四把:严厉跨得过对象吗")
pos=G.positive_control("正对照 453×455 必须复现 #639 的 +0.8301", planted=pc["rho"], floor=0.70, spread=0.03)
neg=G.negative_control("安慰剂:四项 × |纬度| 必须 ≈0", null=max(pl), effect=abs(med),
                       null_spread=0.05, null_kind="与体罚无关的地理属性")
LO,HI=0.187-0.152,0.187+2*0.152
if pos and neg:
    verdict=("**WORLD A —— 一把尺**" if med>HI else
             "**WORLD B —— 与这具仪器自己制造的耦合不可分辨**" if LO<=med<=0.187+0.152 else
             f"**UNVERIFIED —— 中位 {med:+.4f} 落在两条判据之间,不选边**")
    print(f"\n控制齐备 ⇒ 评判(阈 {HI:.3f})。{verdict}")
else:
    verdict=f"UNVERIFIED —— 控制未齐(正 {pos} · 负 {neg})"; print(f"\n⚠ {verdict}")
print(G)

diag=famed.get("对角(既不同性别也不同年龄)",np.nan)
oth=[v for k,v in famed.items() if not k.startswith("对角")]
print(f"\n=== 混淆控制①:对角 {diag:+.4f} vs 其余两族中位 {np.median(oth):+.4f} ===")
print("  对角显著低 ⇒ 单句伪影在起作用;对角同样高 ⇒ 解释不了。")
json.dump(dict(quad_n=int(len(quad)),pairs=rows,median=med,family=famed,
               cross_team=dict(n=int(len(m2)),rho=r2),verdict=verdict,
               positive=pc["rho"],placebo_max=max(pl),unchallenged=True),
          open(OUT/"one_ruler.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'one_ruler.json'}")

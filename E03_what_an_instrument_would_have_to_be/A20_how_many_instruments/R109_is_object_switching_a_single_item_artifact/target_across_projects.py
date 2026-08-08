"""E03·A20·R109 —— 「换对象几乎不改变」是不是单项目伪影

**类型:FRONTIER**。`#664` 查出单仪器声明 8/11,其中「换对象 +0.845」两侧都出自 `barry1977agents`。

⚠ **BASIN**:连续五次下注反对自己喜欢的结果、四次输(`#665` 已把这条校准记下)。
   **本轮仍下注不受欢迎的那一面:W2(它是单项目伪影)** —— 下注是纪律,不是预测。

## ⑤ 的预言应验了一半,而穷举的 bug 差点让我看错另一半

**`lang1998conan` 只编了一个对象**(`SCCS1766` 晚童年男孩;`SCCS620` 是配偶不是儿童)
⇒ **体罚上的跨项目「换对象」构不成。**

⚠ **而我第一次穷举把每个变量的 n 都印成 186 —— 因为我数的是 `data.csv` 的行数,不是非缺失值。**
**`#617` 记过这个缺陷一模一样**,而这一轮的穷举正是为了「先查再说话」而设的。
按真实覆盖重算:75 个命中里 **21 个行数 ≥150 而真实 <60**。
`rohner1981parental` 按性别分列的温暖变量真实 n = **5 / 6**,联合 **4** ⇒ **不可用。**

## 而真实覆盖也翻出了一条真的路
**`barry1976traits` 的 `SCCS322–325`「Obedience:幼年/晚童年 男/女」,真实 n = 160–162** ——
**一个完整的四对象切分,出自与 `barry1977agents` 不同的出版物。**
⚠ **但两者第一作者同为 Barry ⇒ 跨出版物,不是跨独立研究组。这一点必须先测,不能假定。**

## G1 ESTIMAND(先于方法)
**主量 = `SCCS322–325`(服从,Barry 1976)四对象两两秩相关的中位**,与 `#640` 的体罚 **+0.8392** 并列。
## G2 CONTROLS
**正对照**:`SCCS453–456`(体罚,Barry 1977)必须复现 `#640` 的 **+0.8392**。
**安慰剂**:服从四对象 × |纬度|,必须 ≈0。**这个零该不该是零?** 该 ⇒ `negative_control`。
**独立性检查(不能假定)**:`ρ(服从, 体罚)` 的跨构念中位 —— **若 ≥0.90,两套码是一套,判不了。**
## G3/G4:两个四件套各 6 对全报;{原始, 天花板归一} 两条规格。
## KILL(条件式)
if 正对照复现 and 安慰剂≈0 and 独立性检查 <0.90:
  服从四对象中位 **≥0.40** -> **W1:「换对象几乎不改变」不是单项目伪影(但仅跨出版物)**
  **<0.40** -> **W2:`#640` 必须缩小成「在 Barry 1977 这一个项目内」**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**同一第一作者 ⇒ 跨出版物而非跨独立团队**,这是这一轮拿不到的东西 ·
**体罚上的跨项目换对象构不成**(Lang 只有一个对象),可证伪形式:
  *若出现一个 `source ∉ {barry1977agents}` 且按 ≥2 个对象分列的儿童体罚变量、真实 n ≥ 60,即被推翻* ·
**跨仪器:换不了仪器,只此一具** · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
B="data/external/dplace/repo/datasets/SCCS/"
D=pd.read_csv(B+"data.csv"); SOC=pd.read_csv(B+"societies.csv")
W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
BLK={r.id:(int(np.floor(r.Lat/10)),int(np.floor(r.Long/10))) for r in SOC.dropna(subset=["Lat","Long"]).itertuples()}
LAT={r.id:abs(r.Lat) for r in SOC.dropna(subset=["Lat"]).itertuples()}
OBEY=["SCCS322","SCCS323","SCCS324","SCCS325"]; CORP=["SCCS453","SCCS454","SCCS455","SCCS456"]
print("=== 硬规则①:真实非缺失 n,以及每一对的联合 n ===")
for nm,Q in [("服从 Barry 1976",OBEY),("体罚 Barry 1977",CORP)]:
    print(f"  {nm}: 各自 {[int(W[c].notna().sum()) for c in Q]}")
    print(f"      6 对联合 n {[int(W[[a,b]].dropna().shape[0]) for a,b in combinations(Q,2)]}")

def rmax(a,b,s=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if s<0: y=y[::-1]
    return float(spearmanr(x,y).statistic)
def quad(Q,floor=30):
    out=[]
    for a,b in combinations(Q,2):
        m=W[[a,b]].dropna()
        if len(m)<floor: out.append((a,b,np.nan,np.nan,len(m))); continue
        r=float(spearmanr(m[a],m[b]).statistic)
        c=rmax(m[a],m[b],1 if r>0 else -1)
        out.append((a,b,r,r/abs(c) if abs(c)>1e-9 else np.nan,len(m)))
    return out
res={}
for nm,Q in [("服从 Barry 1976",OBEY),("体罚 Barry 1977",CORP)]:
    o=quad(Q); v=[x[2] for x in o if np.isfinite(x[2])]; nvals=[x[3] for x in o if np.isfinite(x[3])]
    res[nm]=dict(pairs=[(a,b,float(r),float(n) if np.isfinite(n) else None,int(k)) for a,b,r,n,k in o],
                 med=float(np.median(v)),med_norm=float(np.median(nvals)))
    print(f"\n=== {nm}:6 对全报 ===")
    for a,b,r,n,k in o: print(f"  {a}×{b} n={k:3d} ρ={r:+.4f} 归一={n:+.4f}")
    print(f"  **中位 {np.median(v):+.4f} · 归一中位 {np.median(nvals):+.4f}**")

# 独立性检查
cross=[]
for a in OBEY:
    for b in CORP:
        m=W[[a,b]].dropna()
        if len(m)>=30: cross.append(abs(float(spearmanr(m[a],m[b]).statistic)))
indep=float(np.median(cross))
print(f"\n=== 独立性检查(不能假定):服从 × 体罚 的跨构念 |ρ| 中位 = **{indep:.4f}** (16 对)")
print(f"  {'⚠ ≥0.90 ⇒ 两套码是一套,判不了' if indep>=0.90 else '**<0.90 ⇒ 是两套不同的码**'}")
d=W[[OBEY[0]]].dropna(); d["_l"]=[LAT.get(x,np.nan) for x in d.index]; d=d.dropna()
pl=abs(float(spearmanr(d[OBEY[0]],d["_l"]).statistic))
print(f"\n=== 控制 ===\n  正对照 体罚中位 {res['体罚 Barry 1977']['med']:+.4f}(`#640`: +0.8392)")
print(f"  安慰剂 服从 × |纬度| = {pl:.4f}")
G=Gate("「换对象几乎不改变」是不是单项目伪影")
p1=G.positive_control("体罚四对象必须复现 #640 的 +0.8392(容差 0.02)",
                      planted=float(0.02-abs(res['体罚 Barry 1977']['med']-0.8392)),floor=0.0,spread=0.002)
p2=G.negative_control("安慰剂:服从 × |纬度|",null=pl,effect=abs(res['服从 Barry 1976']['med']),
                      null_spread=0.05,null_kind="与管教无关的地理属性")
om=res['服从 Barry 1976']['med']
if p1 and p2 and indep<0.90:
    verdict=(f"**W1 —— 不是单项目伪影:服从四对象中位 {om:+.4f}(跨出版物,同一第一作者)**" if abs(om)>=0.40
             else f"**W2 —— `#640` 必须缩小成「在 Barry 1977 这一个项目内」:服从四对象只有 {om:+.4f}**")
elif indep>=0.90: verdict=f"**判不了 —— 两套码 |ρ|={indep:.4f} ≈ 一套**"
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(quads=res,independence=indep,placebo=pl,verdict=verdict,
               lang_single_target=True,coverage_bug="穷举第一版用行数(186)而非非缺失值,#617 同型;75 个命中里 21 个假覆盖",
               unchallenged=True),open(OUT/"target_across_projects.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'target_across_projects.json'}")

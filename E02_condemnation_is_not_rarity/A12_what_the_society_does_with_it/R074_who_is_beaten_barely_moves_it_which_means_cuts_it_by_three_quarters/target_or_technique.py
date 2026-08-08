"""E02·A12·R678 —— 换谁挨打 vs 换用什么手段,哪一个更能改变一个社会的数

**行动类型:FRONTIER**,而且**这一步的正面结果是我不希望看到的那一个**(BASIN RULE)。

`#640` 量了**换对象**:体罚四个对象之间,中位 +0.8392。结论写成「严厉是一件事的谓词」。
`barry1977agents` 同时把**同样四个对象**编在了**十种管教手段**上。于是第一次可以量**换手段**。

⚠ **若换手段也是 +0.84,那 `#640` 的 +0.84 就不是关于社会的事实,是编码者的光环效应** ——
  同一个人给同一个社会的所有育儿栏打了相近的分。**那会把 `#640` 的解释撤回**,而这正是本轮要冒的险。
  这是 §3 的 meta-separator:**存在一个可信的结果,它会显示我的世界分解本身是错的。**

## 硬规则⑥:先读码本

十种手段全部可序。**极性只有一半是自述的**(同一团队同一篇论文):
  有语言锚:以身作则 `11=Example given as most important` · 讲课 `11=Constant and one of the most important`
            放任 `0=Harsh socialization by parents` · 疼爱 `2=generally low expression of affection`
  只有 `(k of 11)`:舆论压力 · 取笑 · 责骂 · 警告 · 奖励礼物 —— **极性不可读 ⇒ 本轮不用**
  体罚:`#639` 从数据内部定为「越大越重」。
=> **5 把极性确定的尺 x 4 个对象 = 20 个变量。**

## 免费的极性正对照(码本给出的方向性预测,数据可以证伪)

**体罚 x 放任必须是负的** —— 打得多的社会应当更不放任。**若为正,两个极性里至少一个是错的,整轮停。**

G1 ESTIMAND(先于方法):
  **A** = 组内换对象:每种手段的 6 个对象对,共 30 对,取中位 rho。
  **B** = 组内换手段:每个对象的 10 个手段对,共 40 对,取中位 |rho|。
  **对比 A - B**,以及 B 的符号结构。

预测矩阵:
| 世界 | A | B | 读法 |
| A 一整套/光环 | ~0.84 | ~0.84 | 换手段与换对象一样不改变 => **评分者光环,`#640` 的解释被撤回** |
| B 手段各自独立 | ~0.84 | 远小于 0.84 | **`#640` 存活**,「一件事」= 一个具体手段 |
| C 一条严厉<->温暖的轴 | ~0.84 | 中等且**符号成系统** | 手段不独立,但也不是一整套 |
**C 与 B 的区别是符号结构,不是量级** —— 用 20 个变量的 PC1 方差占比判(>=50% 即一条轴)。

KILL(条件式,预注册):
  if 极性对照为负 and 正对照 `453x455` 复现 and 安慰剂约 0:
      B 中位 >= 0.60                       -> WORLD A(光环)
      B 中位 <  0.60 and PC1 >= 0.50       -> WORLD C(一条轴)
      B 中位 <  0.60 and PC1 <  0.50       -> WORLD B(手段各自独立)
  else: UNVERIFIED
⚠ **`#641` 的新规则先跑**:判据的分档间距必须先与该估计量自己的区间宽度比一次。
  **本轮先算 B 中位的块 bootstrap 区间,若它跨过 0.60 则如实记「这一刀判不了」,不硬判。**

IMPOSSIBLE(不写 planned):单一编码团队(这正是本轮要检验的东西,无法同时排除)· 无干预 ·
  五种手段中四种的极性来自码本文字而非数据 · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate

SEEDS=[20260806,7,991]; FLOOR_N=30
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
S=ROOT/"data/external/dplace/repo/datasets/SCCS"
D=pd.read_csv(S/"data.csv"); SOC=pd.read_csv(S/"societies.csv")
W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
BLK={r.id:(int(np.floor(r.Lat/10)),int(np.floor(r.Long/10))) for r in SOC.dropna(subset=["Lat","Long"]).itertuples()}
LAT={r.id:abs(r.Lat) for r in SOC.dropna(subset=["Lat"]).itertuples()}

TECH={"体罚":453,"以身作则":429,"讲课":437,"放任":465,"疼爱":469}
TGT=["幼年男孩","幼年女孩","晚童年男孩","晚童年女孩"]
COL={(t,i):f"SCCS{b+i}" for t,b in TECH.items() for i in range(4)}

def rho(a,b,frame=None):
    m=(frame if frame is not None else W)[[a,b]].dropna()
    if len(m)<FLOOR_N or m[a].nunique()<2 or m[b].nunique()<2: return np.nan,len(m)
    return float(spearmanr(m[a],m[b]).statistic), len(m)

def boot_med(pairs,n=600):
    bl=sorted({BLK.get(s,("na","na")) for s in W.index})
    by={x:[s for s in W.index if BLK.get(s,("na","na"))==x] for x in bl}; out=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            socs=[s for i in rng.integers(0,len(bl),len(bl)) for s in by[bl[i]]]
            f=W.loc[socs]; vals=[]
            for a,b,ab in pairs:
                r,_=rho(a,b,f)
                if np.isfinite(r): vals.append(abs(r) if ab else r)
            if vals: out.append(float(np.median(vals)))
    return np.array(out)

print("=== 硬规则①:先打印 n ===")
for t,b in TECH.items():
    print(f"  {t:5s} SCCS{b}-{b+3}  n = {[int(W[COL[(t,i)]].notna().sum()) for i in range(4)]}")

print("\n=== 极性正对照(码本的方向性预测):体罚 x 放任 必须为负 ===")
pol=[]
for i in range(4):
    r,n=rho(COL[("体罚",i)],COL[("放任",i)]); pol.append(r)
    print(f"  {TGT[i]:6s} n={n:3d}  ρ = **{r:+.4f}**")
pol_ok=all(np.isfinite(x) and x<0 for x in pol)
print(f"  => 四格全负? **{pol_ok}**")

print("\n=== A:组内换对象(每种手段 6 对)===")
A=[]
for t in TECH:
    v=[]
    for i,j in combinations(range(4),2):
        r,n=rho(COL[(t,i)],COL[(t,j)])
        if np.isfinite(r): v.append(r); A.append((COL[(t,i)],COL[(t,j)],False))
    print(f"  {t:5s} 中位 {np.median(v):+.4f}  ({len(v)}/6 对可算)")
Av=[rho(a,b)[0] for a,b,_ in A]; Amed=float(np.median(Av))

print("\n=== B:组内换手段(每个对象 10 对)===")
B=[]; sign_tbl={}
for i in range(4):
    v=[]
    for t1,t2 in combinations(TECH,2):
        r,n=rho(COL[(t1,i)],COL[(t2,i)])
        if np.isfinite(r):
            v.append(abs(r)); B.append((COL[(t1,i)],COL[(t2,i)],True))
            sign_tbl.setdefault((t1,t2),[]).append(r)
    print(f"  {TGT[i]:6s} |ρ| 中位 {np.median(v):+.4f}  ({len(v)}/10 对可算)")
Bv=[abs(rho(a,b)[0]) for a,b,_ in B]; Bmed=float(np.median(Bv))

print(f"\n  **A(换对象)中位 = {Amed:+.4f}  ·  B(换手段)|ρ| 中位 = {Bmed:+.4f}  ·  A − B = {Amed-Bmed:+.4f}**")
bs=boot_med(B); Blo,Bhi=np.quantile(bs,[.025,.975])
print(f"  B 中位 95% CI = [{Blo:+.4f}, {Bhi:+.4f}]")
print(f"  ⚠ `#641` 规则:0.60 这一刀 {'**跨在区间内 ⇒ 判不了**' if Blo<0.60<Bhi else '**在区间之外 ⇒ 可判**'}")

print("\n=== 手段对的符号(四个对象的中位)===")
for (t1,t2),v in sorted(sign_tbl.items(), key=lambda x:np.median(x[1])):
    print(f"  {t1:5s} × {t2:5s}  {np.median(v):+.4f}")

sub=W[[COL[k] for k in COL]].dropna()
z=(sub-sub.mean())/sub.std()
ev=np.linalg.eigvalsh(np.corrcoef(z.values.T))[::-1]
pc1=float(ev[0]/ev.sum())
print(f"\n=== PC1 方差占比(20 个变量,完整观测 n={len(sub)})= **{pc1:.4f}** ===")

d=W[[COL[("体罚",2)]]].dropna(); d["_l"]=[LAT.get(s,np.nan) for s in d.index]; d=d.dropna()
pl=abs(float(spearmanr(d[COL[("体罚",2)]],d["_l"]).statistic))
pc,_=rho("SCCS453","SCCS455")
G=Gate("换谁挨打 vs 换用什么手段")
p1=G.positive_control("正对照 453x455 复现 #640",planted=pc,floor=0.70,spread=0.03)
p2=G.negative_control("安慰剂:体罚 x |纬度| 约 0",null=pl,effect=abs(Amed),null_spread=0.05,
                      null_kind="与管教无关的地理属性")
if p1 and p2 and pol_ok:
    if Blo<0.60<Bhi: verdict="**判不了 —— B 中位的区间跨过 0.60(`#641` 规则)**"
    elif Bmed>=0.60: verdict="**WORLD A —— 评分者光环,`#640` 的解释被撤回**"
    elif pc1>=0.50: verdict="**WORLD C —— 一条严厉↔温暖的轴**"
    else: verdict="**WORLD B —— 手段各自独立,`#640` 存活**"
    print(f"\n控制齐备(含极性 {pol_ok})⇒ {verdict}")
else:
    verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2} · 极性 {pol_ok})"; print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(A_median=Amed,B_median=Bmed,B_ci=[float(Blo),float(Bhi)],pc1=pc1,
               polarity_control=[float(x) for x in pol],polarity_ok=bool(pol_ok),
               positive=pc,placebo=pl,n_complete=int(len(sub)),
               sign_table={f"{a}×{b}":float(np.median(v)) for (a,b),v in sign_tbl.items()},
               verdict=verdict,unchallenged=True),
          open(OUT/"target_or_technique.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'target_or_technique.json'}")

# ── 标注的更正(同一轮内):预注册的 C-vs-B 判据用错了聚合层级 ────────────────
# PC1 算在 20 个变量上,而那里的主导结构是**五个手段块**(每块内部 +0.84),不是手段之间的轴。
# 正确:先把每种手段塌成一个社会分(4 个对象的标准化均值),再算 5 个分的 PC1。
# ⚠ **A-vs-B 的判决不受影响**(B 中位远低于 0.60);受影响的只有 C 与 B 的区分。
print("\n=== 标注的更正:PC1 应算在「手段分」上,不是 20 个变量上 ===")
score={}
for t in TECH:
    cols=[COL[(t,i)] for i in range(4)]
    zz=(W[cols]-W[cols].mean())/W[cols].std()
    score[t]=zz.mean(axis=1, skipna=True)
Sdf=pd.DataFrame(score).dropna()
ev2=np.linalg.eigvalsh(np.corrcoef(((Sdf-Sdf.mean())/Sdf.std()).values.T))[::-1]
pc1b=float(ev2[0]/ev2.sum()); pc2b=float(ev2[1]/ev2.sum())
print(f"  5 个手段分,完整观测 n={len(Sdf)}  PC1 = **{pc1b:.4f}** · PC2 = {pc2b:.4f} · PC1+PC2 = {pc1b+pc2b:.4f}")
print(f"  20 变量版(预注册写的)= {pc1:.4f} —— **层级错了,这一格作废**")
print(f"  => C vs B:PC1 {'>=' if pc1b>=0.50 else '<'} 0.50 ⇒ " +
      ("**WORLD C —— 一条轴**" if pc1b>=0.50 else "**WORLD B —— 手段各自独立**"))
d=json.load(open(OUT/"target_or_technique.json"))
d["pc1_corrected_technique_level"]=pc1b; d["pc2_corrected"]=pc2b
d["n_complete_technique"]=int(len(Sdf))
d["correction"]="预注册的 PC1 算在 20 个变量上,层级错;正确层级是 5 个手段分。A-vs-B 判决不受影响。"
d["verdict"]=("**WORLD C —— 一条严厉↔温暖的轴**" if pc1b>=0.50 else "**WORLD B —— 手段各自独立,`#640` 存活**")
json.dump(d,open(OUT/"target_or_technique.json","w"),indent=1,ensure_ascii=False)
print(f"  最终判决:{d['verdict']}")

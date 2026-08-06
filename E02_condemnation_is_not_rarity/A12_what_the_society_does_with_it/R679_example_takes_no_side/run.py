"""E02·A12·R679 —— 以身作则为什么不站队:它可能根本不是一种管教手段

**类型:FRONTIER**。`#642` 里唯一没有被任何判据覆盖、而又最反直觉的一格。

`#642` 实测:体罚/讲课/放任/疼爱四者两两 |rho| 在 **0.31–0.47**;
而**以身作则对这四种全在 +0.0175 ~ +0.1393**。
**而它自己四个对象之间是 +0.9215 —— 五种手段里最高。所以「不站队」不是「测不准」。**

两个世界:
  **W1 它是一种手段,只是与别的手段正交** —— 一个社会可以又打又讲又榜样,互不牵连。
  **W2 它量到的根本不是管教手段,是一种社会结构** —— 孩子看得见成年人在做什么。
  **区别是本体的**:W2 说这一栏被放错了类别,那么 `#642` 的「五种手段」本身就是错的分组。

G1 ESTIMAND(两个,都预注册,多重性 = 2 个族,两个都报):
  **① E vs O**:E = 以身作则 x 其余四种的 |rho| 中位(4 对象 x 4 = 16 对);
               O = 其余四种**彼此之间**的 |rho| 中位(4 对象 x 6 = 24 对)。对比 **O − E**。
  **② S**:以身作则 x 四个**结构**变量的 |rho| 中位,对上 **体罚 x 同样四个**的 |rho| 中位。
     结构变量(极性全部可读,同一来源同样四个对象):
       `SCCS425+i` 正规学校教育  `1=Informal training, with minimal guidance` -> `6=Formal schooling typical`
       `SCCS409+i` 教育中非父母参与 · `SCCS361+i` 养育中非父母参与 · `SCCS377+i` 权威中非父母参与
       三者同码:`1=Exclusively parental` -> `7=Exclusively non-parental`
     **若以身作则贴着结构而体罚不贴 => W2。**

⚠ **跑之前写死的最强混淆:量程。** 实测 以身作则 `SCCS429` 唯一值 **6** 档、sd **1.386**;
  体罚 `SCCS453` 唯一值 **10** 档、sd **2.106**。**窄量程会压低相关**,足以伪造 E 低。
  **同一迭代内的控制:把每个变量都按分位裁成 6 档再算一遍两个中位。**

CONTROLS:**正对照** 以身作则自己四对象之间必须复现 `#642` 的 +0.9215 ——
  **这一格证明「不站队」不是「测不准」,是本轮最重要的控制。**
  **安慰剂** 以身作则 x |纬度| 约 0。**地板** 每一对 n < 30 记判不了。
KILL(条件式,预注册,并遵守 `#641` 规则:先看区间再判):
  if 正对照复现 and 安慰剂约 0:
      **O − E 的块 bootstrap 区间含零 -> 判不了,如实记**;不含零且 O > E -> 以身作则确实是异类
      **S 的差的区间含零 -> W1 与 W2 判不了**;不含零且 以身作则 > 体罚 -> **W2**
  else: UNVERIFIED
IMPOSSIBLE(不写 planned):单一编码团队 · 无干预 · 「结构」四变量本身也是这个团队编的 · `[unchallenged]`
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

EX="以身作则"; TECH={EX:429,"体罚":453,"讲课":437,"放任":465,"疼爱":469}
STRUCT={"正规学校":425,"非父母教育":409,"非父母养育":361,"非父母权威":377}
TGT=["幼年男孩","幼年女孩","晚童年男孩","晚童年女孩"]
OTH=[t for t in TECH if t!=EX]

def rho(a,b,f=None):
    m=(f if f is not None else W)[[a,b]].dropna()
    if len(m)<FLOOR_N or m[a].nunique()<2 or m[b].nunique()<2: return np.nan,len(m)
    return float(spearmanr(m[a],m[b]).statistic), len(m)

def rebin(s,k=6):
    try: return pd.qcut(s.rank(method="first"),k,labels=False,duplicates="drop")
    except Exception: return s

WB=W.copy()
for t,b in list(TECH.items())+list(STRUCT.items()):
    for i in range(4):
        c=f"SCCS{b+i}"
        if c in WB.columns: WB[c]=rebin(WB[c])

def med(pairs,f=None):
    v=[abs(rho(a,b,f)[0]) for a,b in pairs]; v=[x for x in v if np.isfinite(x)]
    return (float(np.median(v)) if v else np.nan), len(v)

E_pairs=[(f"SCCS{TECH[EX]+i}",f"SCCS{TECH[t]+i}") for i in range(4) for t in OTH]
O_pairs=[(f"SCCS{TECH[a]+i}",f"SCCS{TECH[b]+i}") for i in range(4) for a,b in combinations(OTH,2)]
Sx_pairs=[(f"SCCS{TECH[EX]+i}",f"SCCS{STRUCT[s]+i}") for i in range(4) for s in STRUCT]
Sc_pairs=[(f"SCCS{TECH['体罚']+i}",f"SCCS{STRUCT[s]+i}") for i in range(4) for s in STRUCT]

print("=== 硬规则①:以身作则的 16 个 rho 全报(不是中位)===")
for i in range(4):
    row=[]
    for t in OTH:
        r,n=rho(f"SCCS{TECH[EX]+i}",f"SCCS{TECH[t]+i}"); row.append(f"{t} {r:+.4f}(n={n})")
    print(f"  {TGT[i]:6s} " + " · ".join(row))

Em,En=med(E_pairs); Om,On=med(O_pairs); Sxm,Sxn=med(Sx_pairs); Scm,Scn=med(Sc_pairs)
print(f"\n=== 估计量① ===\n  E 以身作则×其余四种 |ρ| 中位 = **{Em:+.4f}** ({En}/16)")
print(f"  O 其余四种彼此    |ρ| 中位 = **{Om:+.4f}** ({On}/24)   **O − E = {Om-Em:+.4f}**")
print(f"\n=== 估计量② ===\n  以身作则 × 结构 |ρ| 中位 = **{Sxm:+.4f}** ({Sxn}/16)")
print(f"  体罚   × 结构 |ρ| 中位 = **{Scm:+.4f}** ({Scn}/16)   **差 = {Sxm-Scm:+.4f}**")

def boot_diff(p1,p2,n=600):
    bl=sorted({BLK.get(s,("na","na")) for s in W.index})
    by={x:[s for s in W.index if BLK.get(s,("na","na"))==x] for x in bl}; out=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            socs=[s for i in rng.integers(0,len(bl),len(bl)) for s in by[bl[i]]]
            f=W.loc[socs]; a,_=med(p1,f); b,_=med(p2,f)
            if np.isfinite(a) and np.isfinite(b): out.append(b-a)
    return np.array(out)

b1=boot_diff(E_pairs,O_pairs); l1,h1=np.quantile(b1,[.025,.975])
b2=boot_diff(Sc_pairs,Sx_pairs); l2,h2=np.quantile(b2,[.025,.975])
print(f"\n=== `#641` 规则:先看区间 ===")
print(f"  O − E = {Om-Em:+.4f}  95% CI [{l1:+.4f},{h1:+.4f}]  -> {'**含零 ⇒ 判不了**' if l1<0<h1 else '**不含零 ⇒ 可判**'}")
print(f"  S 差  = {Sxm-Scm:+.4f}  95% CI [{l2:+.4f},{h2:+.4f}]  -> {'**含零 ⇒ 判不了**' if l2<0<h2 else '**不含零 ⇒ 可判**'}")

print("\n=== 混淆控制:全部裁成 6 档分位后重算 ===")
Em2,_=med(E_pairs,WB); Om2,_=med(O_pairs,WB); Sx2,_=med(Sx_pairs,WB); Sc2,_=med(Sc_pairs,WB)
print(f"  裁档后 O − E = **{Om2-Em2:+.4f}**(原 {Om-Em:+.4f})· S 差 = **{Sx2-Sc2:+.4f}**(原 {Sxm-Scm:+.4f})")
print("  差若基本不变 ⇒ 量程解释不了它。")

pc=float(np.median([rho(f"SCCS{TECH[EX]+i}",f"SCCS{TECH[EX]+j}")[0] for i,j in combinations(range(4),2)]))
d=W[[f"SCCS{TECH[EX]}"]].dropna(); d["_l"]=[LAT.get(s,np.nan) for s in d.index]; d=d.dropna()
pl=abs(float(spearmanr(d[f"SCCS{TECH[EX]}"],d["_l"]).statistic))
G=Gate("以身作则为什么不站队")
p1=G.positive_control("以身作则自己四对象之间复现 #642 的 +0.9215(证明不站队≠测不准)",
                      planted=pc,floor=0.80,spread=0.03)
p2=G.negative_control("安慰剂:以身作则 x |纬度| 约 0",null=pl,effect=abs(Om-Em),
                      null_spread=0.05,null_kind="与管教无关的地理属性")
if p1 and p2:
    v1="**判不了**" if l1<0<h1 else ("**以身作则确实是异类**" if Om>Em else "**方向相反**")
    v2="**W1 与 W2 判不了**" if l2<0<h2 else ("**W2 —— 它量的是结构,不是手段**" if Sxm>Scm else "**W1 —— 它是一种手段,只是正交**")
    verdict=f"① {v1} · ② {v2}"
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(E=Em,O=Om,O_minus_E=Om-Em,ci_OE=[float(l1),float(h1)],
               S_example=Sxm,S_corporal=Scm,S_diff=Sxm-Scm,ci_S=[float(l2),float(h2)],
               rebinned=dict(O_minus_E=Om2-Em2,S_diff=Sx2-Sc2),
               positive=pc,placebo=pl,verdict=verdict,unchallenged=True),
          open(OUT/"example_no_side.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'example_no_side.json'}")

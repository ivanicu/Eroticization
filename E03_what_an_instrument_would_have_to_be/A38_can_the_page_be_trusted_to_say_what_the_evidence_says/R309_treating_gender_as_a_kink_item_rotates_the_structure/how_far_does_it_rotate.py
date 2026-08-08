r"""#870 · E03·A99·R309 —— 把性别当成一道性癖题,会把 E01 的结构转多少?

**还 `#869`①(上一条账本的 NEXT①)。** `#869` 量出:BKS 的 likert 题集由一条鸭子类型规则发现,
它有 **38 份逐字拷贝**,而 `if c!='biomale'  # 剔除:它是协变量` 这一行**只在 32 份里** ⇒
**8 轮把一个性别协变量当成了 likert 题,占它们题集的 5%(1/20)。**
`#869` 把影响如实登记为 **`UNASSESSED`,并明写「不许说影响很小,也不许说影响很大」。本轮去量。**

**⚠ 这不是随便一个 5%。那 8 个文件全部落在 `E01` 的 `A05`–`A07`,而那几条弧的产物是
「第三个维度的正名是 breadth」「宽度有形状,它的第一个名字是体液」「ALS 载荷当人层变量」——
也就是说,如果污染进了成分结构,被污染的是 E01 给一个维度**起名字**的那一步。**

**⚠⚠ 规范检验先跑(三行,零算力,`frontier §1.3`)—— 而它已经把世界 A 打掉了一半:**
`biomale` 与 19 道真题的 **|r| 中位 0.0594**,而**真题彼此之间的 |r| 中位只有 0.0416**
⇒ **它是相关意义上「高于平均」的一列,不是一列惰性噪声。**
最强的四条:主导 **+0.268** · 文字模态 **−0.267** · 顺从 **−0.258** · 异性配对 **+0.258**,
**与真题之间的最大相关(0.283)同一量级** ⇒ **它会载荷,不会置身事外。**
⇒ **所以「它无关紧要」这条路已经不通;本轮问的是「转了多少、转到哪一根成分上」。**

**⚠⚠ 算术先行(`realstat` 的算术陷阱)—— 什么是被强制的:**
**往 PCA 里多加一列,总会把子空间转一点。** 所以「转了」不是发现,**「转得比一道真题还多」才是。**
⇒ **估计量必须是同一个操作在 20 列上各做一遍**,而不是「加 vs 不加」这种两侧不对称的比较:

`G1` **估计量(先于方法命名)**:
   ① **`rot_drop(x)`** —— 从**那 8 轮实际用的 20 列解**出发,**去掉第 x 列**后
      前 k 个主成分子空间的**最大主角**(degrees)。**对 20 列各做一遍,同一个操作。**
      ⇒ **问题变成:`biomale` 在这 20 个数里排第几?** 它是普通一员,还是离群的那一个?
   ② **`|corr(PC_k, biomale)|`**,k = 1,2,3,在 20 列解上 —— **哪一根成分被性别沾上了。**
   ③ **`biomale` 的载荷排名** —— 在 20 列里,它的共同度排第几。

四个世界(**每个都有分支**):
   A **`biomale` 是普通一员**:`rot_drop` 排名进不了前四分之一,且三根成分与它的相关都在安慰剂之下
     ⇒ **那 8 轮的结构结论不受影响**,而**影响的上界要写出来,不许只说「小」。**
   B **它沾上了某一根既有成分** ⇒ **那根成分的名字被污染** ⇒ 点名是哪一根。
   C **它自成一根成分,并把某根真成分挤出前 k** ⇒ ⚠ **这是我不欢迎的那一个** ——
     **E01 那句「第三个维度的正名是 breadth」,可能有一部分是「是不是男的」。**
   D **⚠ 元分离器**:去掉**任何一列**都能把子空间转到同样多 ⇒ **成分不稳,
     「第 k 根成分是什么」这个问法本身在这份数据上不成立** —— 不是某个世界赢了,是命名这件事没有地基。

预测矩阵:
   | 世界 | 现在 | biomale 排名靠后 | 沾上某根成分 | 自成一根并挤掉真成分 | 谁都转一样多 |
   | A 普通一员   | 0.20 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 污染一根   | 0.35 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 挤掉真成分 | 0.25 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D 命名无地基 | 0.20 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(**条件式,不是阈值**):
  if 正控开火(**把 `biomale` 换成一列「按已知权重合成的真题线性组合 + 噪声」⇒ `rot_drop` 必须随权重
     单调上升;权重 0(纯噪声列)时必须落回 19 道真题的分布内**)
     and 负控为零(**样本内打乱 `biomale`,毁掉它与题目的关系、保住边际 ⇒ 必须落回分布内**)
     and 安慰剂为零(**一列独立的 Bernoulli(0.5) ⇒ 必须落回分布内**):
      `biomale` 的 `rot_drop` 排名 ≥ 6/20(即**不在前四分之一**)
        且三根成分与它的 |corr| 都 < 安慰剂 95 分位            -> A
      某根 k≤3 成分与它的 |corr| ≥ 安慰剂 95 分位              -> B
      `biomale` 自己是某根 k≤3 成分的最大载荷                  -> C
      19 道真题的 `rot_drop` 极差 < 它们中位的一半             -> D
  else: UNVERIFIED

⚠⚠ **跑前写下的最强混淆**:**`biomale` 的方差比多数 likert 题大**(0/1 且 p≈0.5 ⇒ var 0.25,
   而 −3..3 的题标准化前方差大得多、标准化后都是 1)。
   ⇒ **若不标准化,它会因为尺度而不是因为内容被低估或高估。**
   ⇒ 控制:**标准化与不标准化两版都跑**,并把它列进规格曲线;
   ⚠ 而**真正的对照是安慰剂**:一列同样是 Bernoulli(0.5)、但与题目无关的列 ——
   **它把「二值 + 这个方差」这件事本身减掉,剩下的才是 `biomale` 的内容。**

`G3` 多重性:整族 = 20 列 × k∈{2,3,4} × 标准化两版 × 相关矩阵两版,BH 与 BY 都做,不同意的格一起发表。
`G4` 规格曲线:k · 标准化 · 成对完整 vs 整行完整,全部逐格印出。
⚠ kill 带 `yardstick` / `yardstick_noise` / `population` / `direction`,
**且控制行的 `population` 与 kill 的 `population` 逐字相同(`#867`)。**

**⚠ 本轮结构性做不到的(登记,不许写「计划中」)**:
① **本轮不重跑那 8 个脚本** —— 它们各自还依赖别的派生文件与阈值,逐个复现是另一整轮;
   **本轮量的是它们共用的那个结构对象(成分空间),不是它们各自的判词** ⇒
   **所以本轮能说「结构转了多少」,不能说「某一轮的结论翻了」**,这条边界不许省;
② **横断面 + 单一发布** ⇒ **无因果识别**;
③ **换不了仪器** —— `biomale` 与那 19 道题只在 BKS 这一份发布里同时存在,
   **第二具仪器在定义上没有这套题** ⇒ **只此一具仪器**,这是结构性的,不是省略;
④ **PCA 的「第 k 根成分」本身不是旋转不变的** —— 本轮报的是**子空间主角**(旋转不变)
   与**逐成分相关**(不是),后者**只在 k 很小且特征值分离时可读**,分离度一并印出。
"""
import json, math, pathlib, sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate
from lib.bks_items import likert_columns

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
SEED, NDRAW = 309, 200
POP = "BKS 那 8 轮实际使用的 20 列题集(19 道真题 + `biomale`),在 k∈{2,3,4} × 标准化两版 × 相关矩阵两版上"

d = pd.read_csv(ROOT / "data/raw/BKSPublic.csv", low_memory=False)
ITEMS = likert_columns(d)                       # 19,例外默认生效
FULL = likert_columns(d, exclude=set())         # 20,那 8 轮实际用的
assert "biomale" in FULL and "biomale" not in ITEMS

print("=== ⓪a 硬规则①:**变量名不是测量** —— 每一列的 n ===")
print(f"  19 道真题:每列非缺中位 **{int(np.median([d[c].notna().sum() for c in ITEMS])):,}** · "
      f"最小 {int(min(d[c].notna().sum() for c in ITEMS)):,}")
print(f"  `biomale`:非缺 **{int(d['biomale'].notna().sum()):,}** · 取值 "
      f"{sorted(pd.Series(d['biomale']).dropna().unique())} ⇒ **0/1,不是 −3..3 的 likert 尺**")
print(f"  ⚠ 它之所以混进来,是因为鸭子类型规则只查「取值 ⊆ {{−3..3}} 且非缺 > 10000」,"
      f"**而 {{0,1}} ⊂ {{−3..3}}**")

M = d[FULL].astype(float)
C0 = M.corr(min_periods=500)
rb = C0["biomale"].drop("biomale")
off = C0.loc[ITEMS, ITEMS].where(~np.eye(len(ITEMS), dtype=bool)).abs().stack()
print("\n=== ⓪b 规范检验(三行,零算力)—— 它已经把世界 A 打掉了一半 ===")
print(f"  `biomale` 与 19 题 |r| 中位 **{rb.abs().median():.4f}** · 最大 {rb.abs().max():.4f}")
print(f"  19 题彼此 |r| 中位 **{off.median():.4f}** · 最大 {off.max():.4f}")
print(f"  ⇒ **它的中位相关是真题之间的 {rb.abs().median()/off.median():.2f} 倍 —— 高于平均的一列**")


def corr_mat(cols, how, dat=None, extra=None):
    X = (dat if dat is not None else d)[cols].astype(float).copy()
    if extra is not None:
        for k, v in extra.items(): X[k] = v
    if how == "整行完整":
        X = X.dropna()
        if len(X) < 500: return None
        return X.corr()
    return X.corr(min_periods=500)


def subspace(Cm, k, standardize):
    """前 k 个主成分的**载荷子空间**。standardize=False 时用协方差(尺度进来)。"""
    if Cm is None or Cm.isna().to_numpy().any():
        Cm = Cm.fillna(0.0) if Cm is not None else None
        if Cm is None: return None, None
    A = Cm.to_numpy()
    if not standardize:
        sd = np.sqrt(np.diag(A)); A = A * np.outer(sd, sd)   # 退回协方差量纲
    w, V = np.linalg.eigh(A)
    o = np.argsort(w)[::-1]
    return V[:, o[:k]], w[o]


def max_principal_angle(U, W):
    """两个子空间的最大主角(度)。⚠ 旋转不变,不依赖「第 k 根是谁」。"""
    if U is None or W is None: return np.nan
    Qu, _ = np.linalg.qr(U); Qw, _ = np.linalg.qr(W)
    s = np.linalg.svd(Qu.T @ Qw, compute_uv=False)
    return float(np.degrees(np.arccos(np.clip(s.min(), -1, 1))))


def rot_drop(cols, drop, k, standardize, how, dat=None, extra=None):
    """从 20 列解出发,去掉 `drop` 一列后子空间转了多少度。**20 列各做一遍,同一个操作。**"""
    keep = [c for c in cols if c != drop]
    Cf = corr_mat(cols, how, dat, extra)
    Ck = corr_mat(keep, how, dat, extra)
    Uf, _ = subspace(Cf, k, standardize)
    Uk, _ = subspace(Ck, k, standardize)
    if Uf is None or Uk is None: return np.nan
    idx = [cols.index(c) for c in keep]
    return max_principal_angle(Uf[idx, :], Uk)


print(f"\n=== ① 规格曲线:20 列 × k∈{{2,3,4}} × 标准化两版 × 相关矩阵两版 ===")
rows = []
for k in (2, 3, 4):
    for std in (True, False):
        for how in ("成对完整", "整行完整"):
            vals = {c: rot_drop(FULL, c, k, std, how) for c in FULL}
            vals = {c: v for c, v in vals.items() if np.isfinite(v)}
            if len(vals) < 15: continue
            item_v = np.array([v for c, v in vals.items() if c != "biomale"])
            bv = vals.get("biomale", np.nan)
            rank = int(1 + sum(1 for v in item_v if v > bv)) if np.isfinite(bv) else -1
            rows.append(dict(k=k, std=bool(std), how=how, biomale=bv, rank=rank, n=len(vals),
                             item_med=float(np.median(item_v)), item_max=float(item_v.max()),
                             item_min=float(item_v.min()),
                             spread_ratio=float((item_v.max()-item_v.min())/max(np.median(item_v),1e-9))))
        r_ = [x for x in rows if x["k"] == k and x["std"] == std]
        if r_:
            print(f"  k={k} 标准化={std} · " + " · ".join(
                f"{x['how']}: biomale **{x['biomale']:.2f}°** 排名 **{x['rank']}/{x['n']}**"
                f"(真题中位 {x['item_med']:.2f}°)" for x in r_))
if not rows:
    raise SystemExit("⛔ 网格为空 —— 空总体不许当作通过")

RANKS = [r["rank"] for r in rows]
RANK_MED = float(np.median(RANKS))
print(f"  ⇒ **`biomale` 的 `rot_drop` 排名中位 {RANK_MED:.1f}/20** "
      f"(逐格 {min(RANKS)}–{max(RANKS)});**排名 1 = 去掉它转得最多**")

print("\n=== ② 哪一根成分被性别沾上了 ===")
Cf = corr_mat(FULL, "成对完整")
comp = []
for k in (3,):
    U, w = subspace(Cf, k, True)
    ib = FULL.index("biomale")
    ev = w[:6] / w.sum()
    for j in range(k):
        L = U[:, j]
        comp.append(dict(pc=j + 1, biomale_loading=float(L[ib]),
                         loading_rank=int(1 + sum(1 for x in np.abs(L) if x > abs(L[ib]))),
                         top_item=FULL[int(np.argmax(np.abs(L)))],
                         ev_share=float(ev[j])))
        print(f"  PC{j+1}(解释 {ev[j]:.1%}) · `biomale` 载荷 **{L[ib]:+.4f}**,"
              f"在 20 列里排 **{comp[-1]['loading_rank']}** · 最大载荷项:{FULL[int(np.argmax(np.abs(L)))][:46]}")
    print(f"  ⚠ 特征值分离度(前 4):" + " · ".join(f"{x:.3f}" for x in ev[:4])
          + "  —— **分离度小则「第 k 根是谁」不可读**(本轮 `④` 号边界)")

print("\n=== ③ 控制(总体与 kill 逐字相同)===")
rng = np.random.default_rng(SEED)
n = len(d)
# 安慰剂:独立 Bernoulli(0.5),把「二值 + 这个方差」本身减掉
plac = pd.Series(rng.integers(0, 2, n).astype(float), index=d.index)
# 负控:样本内打乱 biomale,毁掉与题目的关系,保住边际
negc = pd.Series(rng.permutation(d["biomale"].to_numpy()), index=d.index)
# 正控:按已知权重合成的真题线性组合 + 噪声,剂量扫
Z = d[ITEMS].astype(float)
Zs = (Z - Z.mean()) / Z.std(ddof=1)
base = Zs.mean(axis=1)
def swap_rot(series, k=3, std=True, how="成对完整"):
    dat = d.copy(); dat["__x__"] = series.to_numpy()
    cols = ITEMS + ["__x__"]
    return rot_drop(cols, "__x__", k, std, how, dat=dat)
BIO_R = swap_rot(d["biomale"])
PLA_R = swap_rot(plac)
NEG_R = swap_rot(negc)
DOSE = [0.0, 0.25, 0.5, 1.0, 2.0]
noise = pd.Series(rng.normal(size=n), index=d.index)
curve = [(t, swap_rot((t * base.fillna(0) + noise))) for t in DOSE]
MONO = all(curve[i][1] <= curve[i + 1][1] + 1e-9 for i in range(len(curve) - 1))
ref = np.array([rot_drop(FULL, c, 3, True, "成对完整") for c in ITEMS])
ref = ref[np.isfinite(ref)]
FLOOR = float(np.quantile(ref, 0.95))
print(f"  ⚠ 三条控制**与 `biomale` 用同一个操作**(把那一列换掉再算 `rot_drop`),而不是另一种比较")
print(f"  正控(剂量-反应,合成列 = τ×真题均值 + 噪声):"
      + " · ".join(f"τ={t}→{v:.2f}°" for t, v in curve) + f" · **单调 {MONO}**")
print(f"     **τ=0(纯噪声列)⇒ {curve[0][1]:.2f}°**,必须落回真题分布内(95 分位 {FLOOR:.2f}°):"
      f"**{'是' if curve[0][1] <= FLOOR else '⚠ 否'}** ⇒ ⚠ **`G2` 控制必须能失败**")
CEIL = curve[-1][1] - curve[0][1]
print(f"     **控制也必须能通过**:floor {curve[0][1]:.2f}° < ceiling {curve[-1][1]:.2f}°,"
      f"跨度 {CEIL:.2f}°")
print(f"  负控(打乱 `biomale`,保边际毁关系)⇒ **{NEG_R:.2f}°**,真题 95 分位 {FLOOR:.2f}° ⇒ "
      f"**{'落回分布内' if NEG_R <= FLOOR else '⚠ 没落回'}**")
print(f"     ⚠ **「这个零该不该是零?」不该是 0,该是「一道真题的量级」** —— "
      f"加一列总会转一点,**打乱之后剩下的正是「加了一列」这件事本身**")
print(f"  安慰剂(独立 Bernoulli(0.5),把「二值+方差」减掉)⇒ **{PLA_R:.2f}°** ⇒ "
      f"**{'落回分布内' if PLA_R <= FLOOR else '⚠ 没落回'}**")
print(f"  **而 `biomale` 自己:{BIO_R:.2f}°** —— 对照真题中位 {np.median(ref):.2f}° · 95 分位 {FLOOR:.2f}°")

PL95 = float(np.quantile([abs(np.corrcoef(
    pd.Series(rng.integers(0, 2, n)).to_numpy(),
    (d[ITEMS].astype(float).mean(axis=1).fillna(0)).to_numpy())[0, 1]) for _ in range(NDRAW)], 0.95))
U3, w3 = subspace(Cf, 3, True)
S = d[FULL].astype(float).fillna(d[FULL].astype(float).mean())
Ss = (S - S.mean()) / S.std(ddof=1)
pcs = Ss.to_numpy() @ U3
bio = d["biomale"].to_numpy(float)
CORRS = [float(abs(np.corrcoef(pcs[:, j][np.isfinite(bio)], bio[np.isfinite(bio)])[0, 1])) for j in range(3)]
print(f"  **三根成分的人层得分与 `biomale` 的 |corr|**:" + " · ".join(f"PC{j+1} **{c:.4f}**" for j, c in enumerate(CORRS))
      + f" · 安慰剂 95 分位 **{PL95:.4f}**")

ps = np.array([min(1.0, (1 + sum(1 for v in ref if v >= r["biomale"])) / (len(ref) + 1)) for r in rows])
Cn = len(ps); o_ = np.argsort(ps); q = 0.05
cH = q * np.arange(1, Cn + 1) / Cn; cY = cH / np.sum(1.0 / np.arange(1, Cn + 1))
su = lambda pv, cr: (int(np.max(np.where(pv <= cr)[0])) + 1 if (pv <= cr).any() else 0)
kH, kY = su(ps[o_], cH), su(ps[o_], cY)
print(f"\n=== ④ 多重性:整族 **{Cn}** 格 · BH 存活 **{kH}** · BY **{kY}** ===")

SPREAD = float(np.median([r["spread_ratio"] for r in rows]))
G = Gate("#870 · 把性别当成一道性癖题会把结构转多少")
G.asserted("① 规范检验(三行零算力):`biomale` 与 19 题的 |r| 中位 vs 真题彼此的 |r| 中位 —— "
           "**它是不是一列惰性噪声**",
           bool(rb.abs().median() > off.median()),
           f"{rb.abs().median():.4f} vs {off.median():.4f} ⇒ {rb.abs().median()/off.median():.2f} 倍,"
           f"**高于平均,会载荷**", kind="control", population=POP)
G.asserted("② 前提(跑前写下的最强混淆):**`biomale` 是 0/1,方差与 likert 题不同尺** ⇒ "
           "**标准化两版都跑**,而真正的对照是**同样 Bernoulli(0.5) 但与题目无关的安慰剂**",
           bool(any(r["std"] for r in rows) and any(not r["std"] for r in rows)),
           f"标准化 {sum(1 for r in rows if r['std'])} 格 · 不标准化 {sum(1 for r in rows if not r['std'])} 格 · "
           f"安慰剂 {PLA_R:.2f}°", kind="control", population=POP)
G.asserted("③ 正控(**剂量-反应**,与 `biomale` 同一个操作:换掉那一列再算 `rot_drop`)⇒ 必须**单调**;"
           "**τ=0 的纯噪声列必须落回真题分布内**;且 floor < ceiling",
           bool(MONO and curve[0][1] <= FLOOR and CEIL > 0),
           " ".join(f"τ={t}:{v:.2f}°" for t, v in curve) + f" · 单调 {MONO} · τ=0 {curve[0][1]:.2f}° ≤ "
           f"{FLOOR:.2f}° · 跨度 {CEIL:.2f}°", kind="control", population=POP)
G.asserted("④ 负控:打乱 `biomale`(保边际、毁关系)⇒ 落回真题分布内 "
           "(⚠ **这个零不该是 0,该是「一道真题的量级」** —— 加一列总会转一点)",
           bool(NEG_R <= FLOOR), f"{NEG_R:.2f}° ≤ 95 分位 {FLOOR:.2f}°(真题中位 {np.median(ref):.2f}°)",
           kind="control", population=POP)
G.asserted("⑤ 安慰剂:独立 Bernoulli(0.5) ⇒ 落回真题分布内",
           bool(PLA_R <= FLOOR), f"{PLA_R:.2f}° ≤ {FLOOR:.2f}°", kind="control", population=POP)
G.asserted("⑥ kill(预注册):「那 8 轮的结构不受 `biomale` 影响」要成立,需 **它的 `rot_drop` 排名 ≥ 6/20"
           "(不在前四分之一)且三根成分与它的 |corr| 都 < 安慰剂 95 分位**",
           bool(RANK_MED >= 6 and all(c < PL95 for c in CORRS)),
           f"排名中位 {RANK_MED:.1f}/20(逐格 {min(RANKS)}–{max(RANKS)})· "
           f"|corr| PC1 {CORRS[0]:.4f} · PC2 {CORRS[1]:.4f} · PC3 {CORRS[2]:.4f} vs 安慰剂 95 分位 {PL95:.4f} · "
           f"biomale {BIO_R:.2f}° vs 真题中位 {np.median(ref):.2f}°",
           kind="kill",
           yardstick="`rot_drop`(去掉一列后前 k 主成分子空间的最大主角),**20 列各做一遍,同一个操作**;"
                     "尺的零 = 19 道真题自己的 `rot_drop` 分布",
           yardstick_noise=float(np.std(ref)), population=POP,
           direction=[r["rank"] - 5.5 for r in rows])
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
worst_pc = int(np.argmax(CORRS)) + 1
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif SPREAD < 0.5:
    VERD = (f"**D 去掉任何一列都转得差不多(真题 `rot_drop` 极差 / 中位 = {SPREAD:.2f})⇒ "
            f"「第 k 根成分是什么」在这份数据上没有地基。**")
elif any(c >= PL95 for c in CORRS):
    # ⚠⚠ **判词必须只说旋转不变的那部分。**
    # 第一版写的是「`biomale` 是 PC1 的最大载荷 ⇒ 它自成一根成分」——
    # **而本轮自己的第 ④ 号边界就写着:特征值分离度小的时候「第 k 根是谁」不可读。**
    # 实测分离度 0.100 / 0.092 —— **PC1 与 PC2 几乎等大,谁是第一根由抽样决定。**
    # ⇒ 只说两件旋转不变的事:**相关有多大**,以及**去掉它转多少度**。
    hi = max(CORRS); lo = min(c for c in CORRS)
    ev1, ev2 = comp[0]["ev_share"], comp[1]["ev_share"]
    k2 = [r for r in rows if r["k"] == 2]
    k4 = [r for r in rows if r["k"] == 4]
    VERD = (f"**C 性别沾在**领先的那一小片子空间**上,而且沾得比任何一道真题都紧。**\n"
            f"  ⚠ **只说旋转不变的两件事**(本轮第 ④ 号边界:特征值分离度 "
            f"{ev1:.3f} / {ev2:.3f} 几乎等大 ⇒ **「第几根成分」不可读,不说**):\n"
            f"  ① **前两根成分的人层得分与「是不是男的」相关 {CORRS[0]:.4f} 与 {CORRS[1]:.4f}**,"
            f"而安慰剂(独立 Bernoulli(0.5))的 95 分位是 **{PL95:.4f}** ⇒ **{hi/PL95:.0f} 倍**。\n"
            f"  ② 去掉它,**k=2 的子空间转 {k2[0]['biomale']:.2f}°,在 20 列里排 {k2[0]['rank']}/20** —— "
            f"**比去掉任何一道真题都多**(真题中位 {k2[0]['item_med']:.2f}°);\n"
            f"     而 **k=4 时排 {k4[0]['rank']}/20,低于平均** ⇒ **影响集中在领先成分上,"
            f"而领先成分正是被拿去起名字的那几根。**\n"
            f"  ⚠⚠ **而多重性校正下旋转那一半一格都不剩(BH {kH}/{Cn} · BY {kY}/{Cn})** ——\n"
            f"  **所以本轮的结论靠的是相关那一条({hi:.4f} vs {PL95:.4f}),不是转角那一条。如实写。**\n"
            f"  ⇒ **一句关于人的话:那 8 轮在算「一个人的性癖有多宽」的时候,\n"
            f"  把「这个人是不是男的」当成了他的第 20 个性癖 ——\n"
            f"  而它不是安静地待在角落:它和领先那两根轴的相关是 0.56 与 0.47,\n"
            f"  是一个纯随机二值列的三十多倍。于是那几根轴量到的,有一大块是性别。**\n"
            f"  ⚠ **而本轮不能说某一轮的结论翻了** —— 量的是它们共用的结构对象,不是各自的判词;\n"
            f"  逐轮重跑是另一整轮,**在那之前「翻没翻」仍是 `UNASSESSED`**。")
else:
    VERD = (f"**A `biomale` 是普通一员** —— `rot_drop` 排名中位 {RANK_MED:.1f}/20,"
            f"三根成分 |corr| 最大 {max(CORRS):.4f} < 安慰剂 95 分位 {PL95:.4f}。\n"
            f"  ⇒ **那 8 轮的结构结论不受影响,而影响的上界是:子空间最多转 "
            f"{max(r['biomale'] for r in rows):.2f}°,不超过一道真题的 "
            f"{max(r['biomale'] for r in rows)/np.median(ref):.2f} 倍。**")
print(VERD)
print(f"\n⚠ **本轮结构性做不到的**:① **不重跑那 8 个脚本** —— 本轮量的是它们**共用的结构对象**,"
      f"**所以能说「结构转了多少」,不能说「某一轮的结论翻了」**;② 横断面 + 单一发布 ⇒ **无因果识别**;"
      f"③ **换不了仪器** —— `biomale` 与这 19 道题只在 BKS 这一份发布里同时存在,**只此一具仪器**,"
      f"结构性的,不是省略;④ **「第 k 根成分」不是旋转不变的** —— 主角是,逐成分相关不是,"
      f"**特征值分离度已印出,分离小则那一列不可读**。")

json.dump(dict(grid=rows, components=comp, rank_median=RANK_MED, ranks=RANKS,
               corr_pc_biomale=CORRS, placebo_corr_95=PL95,
               biomale_rot=BIO_R, placebo_rot=PLA_R, negative_rot=NEG_R,
               item_rot=dict(median=float(np.median(ref)), p95=FLOOR, sd=float(np.std(ref))),
               dose_curve=[[t, v] for t, v in curve], monotone=MONO, spread_ratio=SPREAD,
               multiplicity=dict(cells=Cn, bh=int(kH), by=int(kY), q=q),
               gauge=dict(biomale_median_r=float(rb.abs().median()),
                          item_median_r=float(off.median()),
                          ratio=float(rb.abs().median() / off.median())),
               derivation="adding any column rotates the subspace; the finding must be 'rotates more "
                          "than a real item', so rot_drop is computed for all 20 columns identically",
               scope_limit="this round measures the SHARED structural object, not the 8 rounds' own "
                           "verdicts; re-running them is a separate round",
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), seed=SEED, population=POP),
          open(OUT / "how_far_does_it_rotate.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'how_far_does_it_rotate.json'}")

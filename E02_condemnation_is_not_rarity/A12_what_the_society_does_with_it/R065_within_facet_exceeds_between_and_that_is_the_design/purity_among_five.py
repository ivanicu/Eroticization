"""E02·A219·R589 — 我差点在盘上躺着一份数据的情况下去提议采数据

`#543` 的 NEXT 本来是「写一份采数据的可行性表」。**在写之前先清点自己的磁盘,结果计划作废。**
行动类型:**FRONTIER**(主检验分离世界)+ 一条 PRODUCTION 记录(第三例结构性阻断)。

**A · 清点的结果(P4 的先前技艺闸,指向我自己的硬盘):**
`data/external/dataverse/` 里五个未打开的包,其中 `10.7910_DVN_SJTRBI` 是
**Graham · Haidt · Nosek 2009 JPSP** —— 道德基础问卷,`Study_3` 有 **8,193 人 × 234 列**。
⚠ 而它**同样只发量表分,不发条目文本**(`HARM_REL` `FAIRNESS_REL` `INGROUP_REL`
`AUTHORITY_REL` `PURITY_REL`,以及它们的 `_STA` / `_AVG` / `SACRED_*` 版本)。
⇒ **这是同一个阻断的第三例**:RWAS 无题目文本(`#541c`)· MSSCQ 无非性条目(`#541b`)· 此处无条目。
**文件在,语义不在 —— 三份外部心理学数据集,三次。**

**B · 所以本轮不做需要分类的事,只做量表分本身能答的那个问题。**
⚠ **我不主张 `PURITY` 是「性的那个」** —— 那需要条目文本,而它不在包里(`#541c` 的教训)。
本轮只用**文件里自带的名字**问一个不需要任何分类的问题:

G1 ESTIMAND(先于方法):对五个基础中的每一个 `f`,
   **`sep(f) = 1 − median|ρ(f, 其他四个)|`** —— 它与其余四个的中位相关有多低。
   **主量 = `sep` 的排名**;`PURITY` 排第几,是一个**读得出来的事实**,不是一个分类。

WORLDS:
  W-PURITY-APART  `PURITY` 的 `sep` 最高 ⇒ 那个基础与其余道德最分得开
  W-PURITY-TYPICAL 排在中间 ⇒ 它不比别的更独立
  W-PURITY-CENTRAL 最低 ⇒ 它反而是最居中的
⚠ BASIN:`W-PURITY-APART` 与页面上「性是更紧、更独立的领域」呼应,**所以是我想要的**,
   **不是**下注方向。本轮下注 `W-PURITY-TYPICAL`。

CONTROLS(G2):
  正对照 同一基础的**两种算法**(`_REL` 与 `_AVG`)必须高度相关 —— 它给出该仪器的上限;
  安慰剂 每个基础 × **随机整数标签** ≈ 0(该是零 ⇒ negative_control);
  关键零 **把五列随机配对**(打乱受访者行)后重算 sep 排名,300 次 —— 排名是否只是噪声;
  规格曲线 三套算法(`_REL` · `_AVG` · `_STA`)各跑一遍,**全格公布**。
KILL(条件式):if 正对照通过 and 安慰剂 ≈ 0:三套算法下 `PURITY` 的排名一致 -> 按三分判
   else UNVERIFIED(排名不稳 = 不可判,不是「居中」)
IMPOSSIBLE:**无条目文本 ⇒ 不能断言任何基础的内容是什么** · 自选网络样本(YourMorals)⇒ 无人群外推 ·
   横断面非因果 · 基础的定义来自量表作者 · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import rankdata
from lib.gates import Gate
SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
df, meta = pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
F = ["HARM", "FAIRNESS", "INGROUP", "AUTHORITY", "PURITY"]
SUITES = {"_REL": "宗教/相关性算法", "_AVG": "平均分", "_STA": "标准化分"}
print(f"=== 硬规则 1:{len(df)} 行 · {len(df.columns)} 列 ===")
for suf in SUITES:
    have = [f + suf for f in F if f + suf in df.columns]
    print(f"  {suf}: {len(have)}/5 列存在  n(全五列非缺失) = "
          f"{int(df[have].notna().all(1).sum()) if len(have)==5 else 'NA'}")
def rho(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 300: return np.nan
    return float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1])
rows, ranks = [], {}
for suf, label in SUITES.items():
    cols = [f + suf for f in F]
    if not all(c in df.columns for c in cols): print(f"  {suf} 缺列,跳过"); continue
    V = {f: df[f + suf].values.astype(float) for f in F}
    sep = {}
    for f in F:
        rs = [abs(rho(V[f], V[g])) for g in F if g != f]
        rs = [r for r in rs if np.isfinite(r)]
        sep[f] = 1 - float(np.median(rs))
        rows.append(dict(suite=suf, foundation=f, sep=sep[f], median_r=float(np.median(rs)),
                         n=int(np.isfinite(V[f]).sum()), k_pairs=len(rs),
                         inclusion=[f"{label}", f"与其余四个基础配对 {len(rs)} 对",
                                    f"每对 n>=300", "秩相关,取绝对值"]))
    order = sorted(F, key=lambda x: -sep[x])
    ranks[suf] = order
    print(f"\n  --- {suf} ({label}) --- sep 由高到低(越高=越与其他分得开)")
    for i, f in enumerate(order): print(f"    {i+1}. {f:10s} sep={sep[f]:.4f} (中位|ρ|={1-sep[f]:.4f})")
G = Gate("五个道德基础里,哪一个与其余最分得开?(Graham·Haidt·Nosek 2009,Study 3)")
pcs = [abs(rho(df[f + "_REL"].values.astype(float), df[f + "_AVG"].values.astype(float)))
       for f in F if f + "_REL" in df.columns and f + "_AVG" in df.columns]
allr = [1 - r["sep"] for r in rows]
G.positive_control("正对照:同一基础的两种算法(_REL × _AVG)", planted=float(np.median(pcs)),
                   floor=float(np.median(allr)), spread=1e-9)
rng = np.random.default_rng(SEEDS[0])
tg = rng.integers(0, 7, len(df)).astype(float)
zs = [abs(rho(df[f + "_AVG"].values.astype(float), tg)) for f in F if f + "_AVG" in df.columns]
G.negative_control("安慰剂:基础分 × 随机整数标签", null=float(np.median(zs)),
                   effect=float(np.median(pcs)), null_spread=float(np.std(zs)),
                   null_kind="与问卷无关的随机整数标签")
# 关键零:打乱行 -> sep 排名是否只是噪声
nullpos = []
for sd in SEEDS:
    r2 = np.random.default_rng(sd)
    for _ in range(100):
        V = {f: df[f + "_AVG"].values.astype(float) for f in F if f + "_AVG" in df.columns}
        V = {f: (v[r2.permutation(len(v))] if i else v) for i, (f, v) in enumerate(V.items())}
        sp = {f: 1 - np.median([abs(rho(V[f], V[g])) for g in V if g != f]) for f in V}
        nullpos.append(sorted(V, key=lambda x: -sp[x]).index("PURITY") + 1)
G.negative_control("关键零:打乱行后 PURITY 的名次应随机(均值≈3)",
                   null=abs(float(np.mean(nullpos)) - 3.0), effect=2.0, null_spread=float(np.std(nullpos)),
                   null_kind="受访者行置换,破坏基础之间的配对")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['suite']}|{r['foundation']}": r for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {f"{r['suite']}|{r['foundation']}": r for r in rows})
pr = {s: o.index("PURITY") + 1 for s, o in ranks.items()}
print(f"\n  PURITY 在三套算法下的名次:{pr}  (1 = 最分得开)")
print(f"  打乱行后的名次均值 = {np.mean(nullpos):.2f}(应 ≈3)")
print("\n" + "=" * 76)
if float(np.median(pcs)) > float(np.median(allr)) and np.median(zs) < 0.5 * float(np.median(pcs)):
    vs = set(pr.values())
    if len(vs) > 1:
        world = "UNVERIFIED"; verdict = f"三套算法给出不同名次 {pr} -> **排名不稳,不可判**(不是「居中」)"
    elif pr["_AVG"] == 1:
        world = "W-PURITY-APART"; verdict = f"三套算法一致把 PURITY 排第 1 -> **它与其余道德最分得开**"
    elif pr["_AVG"] == 5:
        world = "W-PURITY-CENTRAL"; verdict = f"三套一致排第 5 -> **它反而最居中**"
    else:
        world = "W-PURITY-TYPICAL"; verdict = f"三套一致排第 {pr['_AVG']} -> **它不比别的更独立**"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:**包里没有条目文本**,所以「PURITY 这个基础的内容是什么」"
          "在本轮无法验证 —— 本轮只报**文件里自带的名字**的名次,"
          "**不主张它是「性的那个」**,那需要条目,而条目不在包里。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rows=rows, ranks={k: v for k, v in ranks.items()}, purity_rank=pr,
               null_rank_mean=float(np.mean(nullpos)), world=world, verdict=verdict, seeds=SEEDS,
               third_block="第三例:文件在,条目文本不在(RWAS · MSSCQ · MFQ)",
               instrument="Graham·Haidt·Nosek 2009 JPSP Study 3,YourMorals 自选样本",
               impossible=["无条目文本 -> 不能断言任何基础的内容", "自选网络样本无人群外推",
                           "横断面非因果", "基础定义来自量表作者"], unchallenged=True),
          open(OUT / "purity_among_five.json", "w"), indent=1)
print(f"\nwrote {OUT/'purity_among_five.json'}")

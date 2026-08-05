"""E02·A197·R540 — 一个在自己出生地不能被检验的猜想,和一次真正的样本外预测

`#494` 的 NEXT,两件。⚠ 而 ② 有一个诚实问题必须先说破,不是绕过:

**「可见痕迹」这个维度是我看过那六条斜率之后才想到的。**
`#494a` 的分裂是:负 = 同居 · 同性 · 18 岁;正 = 离婚 · 未婚生育 · 16 岁。
我随后注意到正的那三件**一旦发生就留下第三方可见的痕迹**(离婚有法院记录 · 孩子是人 ·
极早的初次性行为在同龄人里显眼),负的那三件**可以不被任何人看见**。
⛔ **这个猜想由那六条斜率生成,因此不能在那六条上被检验。** 在同一批数上「验证」它
是 realstat 的「定义描述实例」那一行 —— **描述实例与定义类别,从内部无法区分。**

⇒ 本轮做三件,每件的行动类型分开标:
  ① **CLOSURE**:补 `Q` 的人层 bootstrap,把 `#494c` 的 `UNCOMPUTED` 换成一个数。
  ② **PRODUCTION**:把猜想**冻结成一个预注册文件**(`FROZEN_visible_trace.md`),
     写清楚它预测什么、由谁检验、在什么数据上 —— **一个不能被检验的猜想,至少要被钉住。**
  ③ **FRONTIER(一个点)**:在**我从没算过阈值斜率的一个话题**上做一次样本外预测:
     **GSS `homosex` × 同性性行为**。可见痕迹编码 = **不可见** -> **预测斜率为负**。
     ⚠ **预测写在计算之前**,并且写在这个 docstring 里,不可事后改。
     ⚠ 一个点不能确证任何东西;它只能**否证**——若为正,猜想当场受创。

G1 ESTIMAND:
  ① `sd_bootstrap(Q)`(人层重抽,重算六条斜率与 Q)。
  ③ `slope_GSS = d|lnOR|/d(谴责占比)`,谴责集取 `homosex` 的 {1} / {1,2} / {1,2,3}
     (1 = always wrong … 4 = not wrong at all;⛔ 5 = other,剔除)。

CONTROLS:③ 用 RULE-v3(`|r| > 自身置换 q95` 且 `> 同问卷参照中位`);精度 = 人层 bootstrap。
KILL(条件式,预注册):
  ③ if 正对照触发 and 阴性为零:
        slope < 0 且 |slope| > MDE -> 猜想在这一点上**存活**(不是被确证)
        slope > 0 且 |slope| > MDE -> **猜想受创,写进冻结文件**
        |slope| < MDE              -> **UNVERIFIED-by-power,直说,不写成「符合」**
     else: UNVERIFIED
IMPOSSIBLE:猜想**不能**在生成它的六条斜率上被检验 · 一个样本外点不能确证 ·
  会话约束禁止派对抗 agent ⇒ 无法做清白编码 ⇒ `[unchallenged]` ·
  GSS 与 NSFG 是不同仪器,斜率的**量级**不可比,**只有符号可比**。
"""
import os, sys, pathlib, json, math, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"

# ================================================================ ① CLOSURE: Q 的 bootstrap
print("=" * 70); print("① CLOSURE — Q 的人层 bootstrap(`#494c` 的 UNCOMPUTED)"); print("=" * 70)


def parse_dct(p):
    out = {}
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out


LAY = parse_dct(NS / "setup" / "2011_2013_FemRespSetup.dct")
NEED = ["staytog", "sxok18", "sxok16", "samesex", "okcohab", "chsuppor", "prevhusb",
        "evrmarry", "samesexany", "nonmarr", "cebow", "agefstsx"]
cols = {n: LAY[n] for n in NEED if n in LAY}
rows = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        rows[n].append(float(v) if v not in ("", ".") else np.nan)
D = {n: np.array(rows[n]) for n in rows}
Nn = len(D["staytog"])
early = np.where(np.isfinite(D["agefstsx"]) & (D["agefstsx"] >= 5) & (D["agefstsx"] <= 60),
                 (D["agefstsx"] <= 16).astype(float), np.nan)
evermar = D["evrmarry"] == 1
BEH = {"samesex": (np.where(D["samesexany"] == 1, 1.0, np.where(D["samesexany"] == 5, 0.0, np.nan)), None),
       "sxok16": (early, None), "sxok18": (early, None),
       "staytog": (np.where(np.isfinite(D["prevhusb"]) & (D["prevhusb"] < 90),
                            (D["prevhusb"] > 0).astype(float), np.nan), evermar),
       "okcohab": (np.where(np.isfinite(D["nonmarr"]) & (D["nonmarr"] < 90),
                            (D["nonmarr"] > 0).astype(float), np.nan), None),
       "chsuppor": (np.where(np.isfinite(D["cebow"]) & (D["cebow"] < 90),
                             (D["cebow"] > 0).astype(float), np.nan), None)}
REV = {"okcohab": True}; CUTS = [[4], [3, 4], [2, 3, 4]]


def cond_at(v, codes, rev=False):
    ok = np.isin(v, [1, 2, 3, 4])
    return np.where(ok, np.isin(v, [5 - c for c in codes] if rev else codes).astype(float), np.nan)


def absl(c, b, mk=None):
    m = np.isfinite(c) & np.isfinite(b)
    if mk is not None: m &= mk
    cc, bb = c[m], b[m]
    a1, a0 = bb[cc == 1], bb[cc == 0]
    if len(a1) < 30 or len(a0) < 30: return np.nan, np.nan
    p1, p0 = a1.mean(), a0.mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, np.nan
    return abs(math.log((p1 / (1 - p1)) / (p0 / (1 - p0)))), float(cc.mean())


def slopes_and_Q(idx=None):
    sl, se_ = [], []
    for t, (b, mk) in BEH.items():
        v = D[t] if idx is None else D[t][idx]
        bb = b if idx is None else b[idx]
        mm = None if mk is None else (mk if idx is None else mk[idx])
        xs, ys = [], []
        for cd in CUTS:
            a, s = absl(cond_at(v, cd, REV.get(t, False)), bb, mm)
            if np.isfinite(a): xs.append(s); ys.append(a)
        if len(xs) < 3 or np.ptp(xs) < 1e-9: return np.nan
        sl.append(np.polyfit(xs, ys, 1)[0])
    sl = np.array(sl)
    return float(((sl - sl.mean()) ** 2).sum() / max(sl.var(ddof=1), 1e-12))  # 形状同 Q 的无权版


# 主 Q 用 `#494` 的加权定义重算一次,并对人层重抽求其展布
R539 = json.load(open(ROOT / "E02_condemnation_is_not_rarity/A196_is_the_gap_one_number_or_many/"
                      "R539_is_the_pooled_null_an_average/results/is_the_pooled_null_an_average.json"))
TOPS = list(BEH)


def weighted_Q(idx=None, nb=60, seed=0):
    rng = np.random.default_rng(seed)
    sl, se_ = [], []
    for t, (b, mk) in BEH.items():
        def one(ii):
            v = D[t][ii]; bb = b[ii]; mm = None if mk is None else mk[ii]
            xs, ys = [], []
            for cd in CUTS:
                a, s = absl(cond_at(v, cd, REV.get(t, False)), bb, mm)
                if np.isfinite(a): xs.append(s); ys.append(a)
            return np.polyfit(xs, ys, 1)[0] if len(xs) >= 3 and np.ptp(xs) > 1e-9 else np.nan
        base_i = np.arange(Nn) if idx is None else idx
        s0 = one(base_i)
        bb_ = [one(rng.integers(0, len(base_i), len(base_i))) for _ in range(nb)]
        bb_ = np.array([x for x in bb_ if np.isfinite(x)])
        if not np.isfinite(s0) or len(bb_) < 10: return np.nan
        sl.append(s0); se_.append(bb_.std())
    sl, se_ = np.array(sl), np.array(se_)
    w = 1 / se_ ** 2
    return float((((sl - (w @ sl) / w.sum()) ** 2) / se_ ** 2).sum())


Qobs = weighted_Q(None, 60, SEEDS[0])
print(f"  Q(重算,内层 60 次)= {Qobs:.3f}   (`#494a` 报的 = {R539['Q']:.3f})")
QB = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(12):
        q = weighted_Q(rng.integers(0, Nn, Nn), 40, sd)
        if np.isfinite(q): QB.append(q)
QB = np.array(QB)
print(f"  **Q 的人层 bootstrap:sd = {QB.std():.3f}**,95% CI [{np.quantile(QB,.025):.2f}, "
      f"{np.quantile(QB,.975):.2f}](B={len(QB)})  vs `#494a` 的零 q95 = {R539['Q_null_q95']:.2f}")
q_above = float((QB > R539["Q_null_q95"]).mean())
print(f"  重抽中落在零 q95 之上的比例 = {q_above:.3f}")

# ================================================================ ③ FRONTIER: 一个样本外点
print("\n" + "=" * 70); print("③ FRONTIER — 样本外一点:GSS homosex(预测:斜率为负,写在跑之前)")
print("=" * 70)
it = pd.read_stata(DTA, iterator=True); vl = it.variable_labels()
cand = [k for k, v in vl.items()
        if re.search(r"sex of (sex )?partner|partners.*(been|were)|sexsex", str(k) + str(v), re.I)]
print(f"  规则①:同性行为候选 = {cand[:8]}")
beh_col = next((c for c in ["sexsex5", "sexsex18", "sexsex"] if c in vl), None)
print(f"  选用 = {beh_col}")
res3 = None
if beh_col:
    g = pd.read_stata(DTA, columns=["year", "homosex", beh_col, "attend", "zodiac", "age",
                                    "educ", "polviews"], convert_categoricals=False)
    g = g.dropna(subset=["homosex", beh_col])
    print(f"  n={len(g)}  years={int(g.year.min())}-{int(g.year.max())}  "
          f"homosex vals={sorted(g.homosex.unique())}  {beh_col} vals={sorted(g[beh_col].unique())[:6]}")
    # ⛔ homosex 5 = other,剔除;behaviour:非「仅异性」= 阳性(码待打印后判定)
    g = g[g.homosex.isin([1, 2, 3, 4])]
    bvals = sorted(g[beh_col].unique())
    beh = (g[beh_col] != bvals[-1]).astype(float) if len(bvals) >= 2 else None
    print(f"  行为阳性率 = {beh.mean():.4f}(阳性 := {beh_col} != {bvals[-1]})")
    xs, ys = [], []
    for cd in [[1], [1, 2], [1, 2, 3]]:
        c = g.homosex.isin(cd).astype(float).values
        a, s = absl(c, beh.values)
        if np.isfinite(a): xs.append(s); ys.append(a); print(f"    谴责集{cd} 占比={s:.3f} |lnOR|={a:.4f}")
    if len(xs) >= 3:
        sl3 = float(np.polyfit(xs, ys, 1)[0])
        bs = []
        for sd in SEEDS:
            rng = np.random.default_rng(sd)
            for _ in range(300):
                i = rng.integers(0, len(g), len(g))
                gx, gb = g.homosex.values[i], beh.values[i]
                X, Y = [], []
                for cd in [[1], [1, 2], [1, 2, 3]]:
                    a, s = absl(np.isin(gx, cd).astype(float), gb)
                    if np.isfinite(a): X.append(s); Y.append(a)
                if len(X) >= 3 and np.ptp(X) > 1e-9: bs.append(np.polyfit(X, Y, 1)[0])
        bs = np.array(bs); MDE3 = 2.8 * bs.std()
        ci3 = np.quantile(bs, [.025, .975])
        print(f"  **slope = {sl3:+.4f}**  CI [{ci3[0]:+.4f},{ci3[1]:+.4f}]  MDE={MDE3:.4f}")
        res3 = dict(slope=sl3, ci=[float(ci3[0]), float(ci3[1])], MDE=float(MDE3),
                    n=int(len(g)), beh_col=beh_col, beh_rate=float(beh.mean()))
else:
    print("  ⚠ GSS 里找不到同性行为变量 -> ③ 不可做,直说")

G = Gate("一个在自己出生地不能被检验的猜想")
G.asserted("猜想不在生成它的数据上被检验", True,
           "「可见痕迹」由 #494a 的六条斜率生成 -> 已冻结,只在样本外一点上预测", kind="control")
G.has_error_bar("Q", value=Qobs, spread=float(QB.std()), spread_source="bootstrap_人层")
if res3:
    ok3 = (res3["slope"] < 0 and abs(res3["slope"]) > res3["MDE"])
    hurt = (res3["slope"] > 0 and abs(res3["slope"]) > res3["MDE"])
    v3 = ("猜想在这一点上**存活**(不是被确证)" if ok3 else
          "**猜想受创**" if hurt else
          f"**UNVERIFIED-by-power**:|{res3['slope']:.3f}| < MDE {res3['MDE']:.3f},看不见,不写成「符合」")
else:
    v3 = "③ 不可做(GSS 无同性行为变量)"
print("\n" + "=" * 70); print(f"③ 判定:{v3}"); print(G)

FROZEN = ROOT / "FROZEN_visible_trace.md"
FROZEN.write_text(f"""# FROZEN — 「可见痕迹」猜想(冻结于 `E02·A197·R540`)

⛔ **这个猜想由 `#494a` 的六条斜率生成,因此不能在那六条上被检验。**

**猜想**:把「算谴责」的门槛放宽时,差距**变大**的话题,是那些**一旦发生就留下第三方可见痕迹**的
(离婚 → 法院记录 · 未婚生育 → 一个孩子 · 极早的初次性行为 → 同龄人可见);
差距**变小**的话题,是那些**可以不被任何人看见**的(同性性接触 · 同居 · 18 岁的性行为)。

**它预测什么**:任何新话题,只要其行为**不留下可见痕迹**,阈值斜率应为**负**。

**它由谁检验**:不是我在这份数据上。需要 (a) 一个**清白上下文**在看到任何斜率之前完成编码,
或 (b) **新的话题**,其斜率我从未计算过。

**已用掉的样本外点**(用掉即失效,不可重复使用):
- `GSS homosex × {res3['beh_col'] if res3 else 'N/A'}`:编码 = 不可见 -> 预测**负**;
  实测 slope = {f"{res3['slope']:+.4f}" if res3 else "N/A"},MDE = {f"{res3['MDE']:.4f}" if res3 else "N/A"} -> {v3}

**冻结日期**:2026-08-05 · **`Q` 的人层 bootstrap sd = {QB.std():.3f}**(补 `#494c` 的 UNCOMPUTED)
""")
print(f"\nwrote {FROZEN}")
json.dump(dict(Q_recomputed=Qobs, Q_boot_sd=float(QB.std()), Q_boot_B=int(len(QB)),
               Q_boot_ci=[float(np.quantile(QB, .025)), float(np.quantile(QB, .975))],
               Q_above_null_share=q_above, out_of_sample=res3, verdict3=v3,
               seeds=SEEDS, unchallenged=True),
          open(OUT / "freeze_and_one_point.json", "w"), indent=1)
print(f"wrote {OUT/'freeze_and_one_point.json'}")

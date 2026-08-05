"""E02·A197·R541 — 把性别读进去之后,那个样本外点说什么

`#495` 的 NEXT。⛔ `#495b` 作废了上一次的样本外点:`sexsex5` 是**伴侣的性别**,
含义随**受访者**性别翻转,而我忽略了 `sex`,阳性率 0.5368 ≈ 男性占比就是告密者。

⚠ **本轮开头再更正一次我自己在 `#495` 的 NEXT 里写反的映射**:
`sexsex5` 码本:**1 = 仅男性伴侣 · 2 = 男女都有 · 3 = 仅女性伴侣**。
⇒ **男性受访者的同性阳性 = {1,2};女性受访者 = {2,3}。**
(`#495` 的 NEXT 写成了「男性 {2,3},女性 {1,2}」—— 反了。按码本写死,并 assert。)

G1 ESTIMAND(先于方法):`slope = d|lnOR|/d(谴责占比)`,
  谴责集取 `homosex` 的 {1} / {1,2} / {1,2,3}(1 = always wrong … 4 = not wrong at all;⛔ 5 剔除),
  行为 = 过去五年有同性伴侣(**按受访者性别正确编码**)。**分性别各跑一次,不合并。**

⛔ **预测已冻结,不可改**(`FROZEN_visible_trace.md`):
  同性性接触**不留下第三方可见痕迹** ⇒ **预测斜率为负**。

⚠ 上一版缺的那道守卫,本轮加上并写在最前:
  **assert 分性别阳性率 ∈ [0.02, 0.30]** —— 同性接触的合理量级。
  **上一版若有这道 assert,`#495b` 那次自伤就不会发生。**

WORLDS:
  W-TRACE  猜想对   -> 两个性别的斜率都为负且超 MDE
  W-BROKEN 猜想错   -> 为正且超 MDE
  W-BLIND  看不见   -> |slope| < MDE -> **UNVERIFIED-by-power,不写成「符合」**
CONTROLS(RULE-v3):正对照 `homosex` × `attend`(`|r| > 自身置换 q95` 且 `> 同问卷参照中位`);
  安慰剂 `homosex` × `zodiac`;精度 = 人层 bootstrap;先算 MDE 再看点估计。
KILL(条件式,预注册):
  if 正对照与安慰剂在**两个性别**都通过:
      两性别斜率都 < 0 且都超 MDE -> 猜想**存活**(用掉这个点)
      任一性别 > 0 且超 MDE       -> **猜想受创**(用掉这个点)
      否则                        -> **UNVERIFIED-by-power,点不算用掉**
  else: UNVERIFIED
IMPOSSIBLE:一个样本外点不能确证 · 无干预 ⇒ 非因果 · 未派对抗 agent ⇒ `[unchallenged]` ·
  GSS 与 NSFG 量级不可比,**只有符号可比**
"""
import os, sys, pathlib, json, math, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)

lab = pd.read_stata(DTA, columns=["sexsex5", "homosex", "sex"], convert_categoricals=True)
raw = pd.read_stata(DTA, columns=["sexsex5", "homosex", "sex"], convert_categoricals=False)
print("=== 规则①:先读码,不猜 ===")
for c in ["sexsex5", "homosex", "sex"]:
    mp = pd.DataFrame({"code": raw[c], "label": lab[c]}).dropna().drop_duplicates().sort_values("code")
    print(f"  {c}: " + " | ".join(f"{int(r.code)}={r.label}" for r in mp.itertuples()))

g = pd.read_stata(DTA, columns=["year", "sex", "homosex", "sexsex5", "attend", "zodiac",
                                "age", "educ"], convert_categoricals=False)
g = g.dropna(subset=["homosex", "sexsex5", "sex"])
g = g[g.homosex.isin([1, 2, 3, 4])]
# 码本:1=仅男性伴侣 2=男女都有 3=仅女性伴侣;男性(sex=1)同性 = {1,2},女性(sex=2) = {2,3}
g["same"] = np.where(g.sex == 1, g.sexsex5.isin([1, 2]),
                     np.where(g.sex == 2, g.sexsex5.isin([2, 3]), np.nan)).astype(float)
print(f"\nn={len(g)}  years {int(g.year.min())}-{int(g.year.max())}")
rates = {}
for s, nm in [(1, "男性"), (2, "女性")]:
    sub = g[g.sex == s]
    rates[nm] = float(sub.same.mean())
    print(f"  {nm}: n={len(sub):6d}  同性伴侣(过去五年)阳性率 = {sub.same.mean():.4f}")
    assert 0.02 <= sub.same.mean() <= 0.30, \
        f"{nm} 阳性率 {sub.same.mean():.4f} 不在 [0.02,0.30] -> 码又读错了(`#495b` 的守卫)"
print("  ✅ 两个性别的阳性率都落在合理量级 —— 上一版缺的正是这道 assert")


def absl(c, b):
    m = np.isfinite(c) & np.isfinite(b)
    cc, bb = c[m], b[m]
    a1, a0 = bb[cc == 1], bb[cc == 0]
    if len(a1) < 30 or len(a0) < 30: return np.nan, np.nan
    p1, p0 = a1.mean(), a0.mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, np.nan
    return abs(math.log((p1 / (1 - p1)) / (p0 / (1 - p0)))), float(cc.mean())


CUTS = [[1], [1, 2], [1, 2, 3]]


def slope_of(hx, bh):
    xs, ys = [], []
    for cd in CUTS:
        a, s = absl(np.isin(hx, cd).astype(float), bh)
        if np.isfinite(a): xs.append(s); ys.append(a)
    if len(xs) < 3 or np.ptp(xs) < 1e-9: return np.nan, None
    return float(np.polyfit(xs, ys, 1)[0]), list(zip(xs, ys))


print("\n=== 分性别,不合并 ===")
res = {}
for s, nm in [(1, "男性"), (2, "女性")]:
    sub = g[g.sex == s]
    sl, pts = slope_of(sub.homosex.values, sub.same.values)
    bs = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(400):
            i = rng.integers(0, len(sub), len(sub))
            v, _ = slope_of(sub.homosex.values[i], sub.same.values[i])
            if np.isfinite(v): bs.append(v)
    bs = np.array(bs); MDE = 2.8 * bs.std(); ci = np.quantile(bs, [.025, .975])
    res[nm] = dict(n=int(len(sub)), slope=sl, MDE=float(MDE),
                   ci=[float(ci[0]), float(ci[1])], sd=float(bs.std()),
                   points=[[float(a), float(b)] for a, b in (pts or [])])
    print(f"  {nm}: n={len(sub):6d}  **MDE={MDE:.4f}**(先算)  slope={sl:+.4f}  "
          f"CI [{ci[0]:+.4f},{ci[1]:+.4f}]")
    for a, b in (pts or []): print(f"        占比={a:.3f} |lnOR|={b:.4f}")

# ---------------------------------------------------------------- 控制 RULE-v3
G = Gate("把性别读进去之后,那个样本外点说什么?(GSS homosex)")
REFV = ["age", "educ", "zodiac"]
ok_all = []
for s, nm in [(1, "男性"), (2, "女性")]:
    sub = g[g.sex == s]
    c = sub.homosex.isin([1, 2]).astype(float).values
    ref = []
    for cn in REFV:
        m = sub[cn].notna().values & np.isfinite(c)
        if m.sum() > 500 and sub[cn][m].nunique() >= 3:
            ref.append(abs(float(np.corrcoef(c[m], sub[cn].values[m])[0, 1])))
    RM = float(np.median(ref))
    m = sub.attend.notna().values
    r = abs(float(np.corrcoef(c[m], sub.attend.values[m])[0, 1]))
    rg = np.random.default_rng(SEEDS[0])
    q = float(np.quantile([abs(np.corrcoef(c[m][rg.permutation(m.sum())], sub.attend.values[m])[0, 1])
                           for _ in range(300)], .95))
    ok_all.append(G.positive_control(f"正对照-v3[{nm}]:谴责×礼拜出席",
                                     planted=r, floor=max(q, RM), spread=1e-9))
    mz = sub.zodiac.notna().values
    zr = float(np.corrcoef(c[mz], sub.zodiac.values[mz])[0, 1])
    ok_all.append(G.negative_control(f"安慰剂[{nm}]:谴责×星座", null=zr, effect=r,
                                     null_spread=float(np.std(ref)), null_kind="个体层标签置换"))
for nm in res:
    G.has_error_bar(f"斜率[{nm}]", value=res[nm]["slope"], spread=res[nm]["sd"],
                    spread_source="bootstrap_人层")

print("\n" + "=" * 70)
if all(ok_all):
    neg = [nm for nm in res if res[nm]["slope"] < 0 and abs(res[nm]["slope"]) > res[nm]["MDE"]]
    pos = [nm for nm in res if res[nm]["slope"] > 0 and abs(res[nm]["slope"]) > res[nm]["MDE"]]
    blind = [nm for nm in res if abs(res[nm]["slope"]) <= res[nm]["MDE"]]
    if len(neg) == 2:
        verdict, used = "两性别斜率都为负且超 MDE -> **猜想在这一点上存活**", True
    elif pos:
        verdict, used = f"{pos} 斜率为正且超 MDE -> **猜想受创**", True
    else:
        verdict, used = (f"{blind} 的 |slope| <= MDE -> **UNVERIFIED-by-power:看不见,"
                         f"不写成「符合」;这个点不算用掉**", False)
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:一个样本外点不能确证;"
          "而 `sexsex5` 问的是**过去五年的伴侣性别**,不是「曾经有过」——"
          "一个只在更早发生过同性接触的人在这里算阴性。")
else:
    verdict, used = f"UNVERIFIED —— 控制未齐 {ok_all}", False
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rates=rates, per_sex=res, verdict=verdict, point_used=bool(used),
               prediction="negative (frozen in FROZEN_visible_trace.md)",
               seeds=SEEDS, unchallenged=True),
          open(OUT / "out_of_sample_with_sex.json", "w"), indent=1)
print(f"\nwrote {OUT/'out_of_sample_with_sex.json'}")

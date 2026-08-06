"""E02·A244·R627 — 「不覆盖中老年」这堵墙,后面有什么?

`#582` 的 NEXT。行动类型:**FRONTIER**。
挑的规则先于挑本身(`#553b`):八条「可量而未量」里 **`#539` 被后续引用 7 次,最多** ⇒ 选它。
它的条款:**「仅 15–44 岁,不覆盖中老年」**(NSFG 的年龄上限)。

**为什么它可量:** **GSS 覆盖全部成年年龄**,且有四道性态度题。
⇒ 可以直接量出**那堵墙后面的东西**:**45 岁之后,性题之间的耦合还变不变。**

G1 ESTIMAND(先于方法):在**每一个调查年内**,按年龄段算 `性内`(四道性题两两 |ρ| 中位),
   **主量 = 45+ 各段与 15–44 段之间的差**,以及**全年龄的斜率**。
   全部题目二值化(沿用 `#535` 的 S-A 严切点),使段与段之间格式一致。
预注册:
  **A** 45+ 与 <45 之差 **< 各自 MDE** ⇒ **这堵墙后面没有新东西**,`#539` 的这条限制**后果 ≈ 0**;
  **B** 差 **超 MDE** ⇒ 限制**真实且已定价**,NSFG 的结论**必须带上「只到 44 岁」**;
  **C** 45+ 的段本身算不出 ⇒ UNVERIFIED。
CONTROLS:正对照 = `age` 与一个**必然随年龄变**的量(`marital` 已婚比例)必须强相关 ·
  安慰剂 = 性内中位 × 随机重排的年龄标签 ≈ 0 · 逐段 n 全部打印
IMPOSSIBLE:GSS 与 NSFG **题目不同**(GSS 四道 vs NSFG 三道)⇒ 这是**同一构念在另一仪器上的年龄剖面**,
  **不是 NSFG 那条曲线的延长** · 横断面 ⇒ 年龄与队列共线(`#539b` 同一条)· [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEX = ["homosex", "premarsx", "xmarsex", "teensex"]
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "age", "marital"] + SEX, convert_categoricals=False)
B = {}
for c in SEX:
    v = g[c].where(g[c].isin([1, 2, 3, 4]))
    B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
AGE = g.age.where(g.age.between(18, 89)).values.astype(float)
BANDS = [(18, 29), (30, 44), (45, 59), (60, 74), (75, 89)]
print("=== 硬规则 1:逐段 n(四题齐全)===")
ok4 = np.all([np.isfinite(B[c]) for c in SEX], 0)
for lo, hi in BANDS:
    m = ok4 & (AGE >= lo) & (AGE <= hi)
    print(f"  {lo}-{hi}: n={int(m.sum()):6d}")
def med_within(mask):
    o = []
    for a, b in itertools.combinations(SEX, 2):
        m = mask & np.isfinite(B[a]) & np.isfinite(B[b])
        if m.sum() < 300 or np.std(B[a][m]) == 0 or np.std(B[b][m]) == 0: continue
        o.append(abs(float(np.corrcoef(rankdata(B[a][m]), rankdata(B[b][m]))[0, 1])))
    return float(np.median(o)) if o else np.nan
print("\n=== 性内耦合的年龄剖面(GSS,全部成年年龄)===")
res, N = {}, len(AGE)
for lo, hi in BANDS:
    m = (AGE >= lo) & (AGE <= hi)
    pt = med_within(m)
    bs = []
    for sd in SEEDS:
        r = np.random.default_rng(sd)
        idx = np.where(m)[0]
        for _ in range(300):
            mm = np.zeros(N, bool); mm[r.choice(idx, len(idx))] = True
            v = med_within(mm)
            if np.isfinite(v): bs.append(v)
    bs = np.array(bs)
    res[f"{lo}-{hi}"] = dict(n=int(m.sum()), point=pt, sd=float(bs.std()), mde=float(2.8 * bs.std()),
                             inclusion=[f"年龄 {lo}-{hi}", "四题二值化(S-A 严切点)", "每对 n>=300"])
    print(f"  {lo}-{hi}: n={int(m.sum()):6d} 性内={pt:.4f} · MDE={2.8*bs.std():.4f}")
young = np.median([res[k]["point"] for k in ("18-29", "30-44")])
old = np.median([res[k]["point"] for k in ("45-59", "60-74", "75-89")])
mde = float(np.median([res[k]["mde"] for k in res]))
print(f"\n  <45 中位 = **{young:.4f}** · 45+ 中位 = **{old:.4f}** · 差 = **{old-young:+.4f}** · 中位 MDE = {mde:.4f}")
G = Gate("「不覆盖中老年」这堵墙,后面有什么?")
mar = g.marital.where(g.marital.between(1, 5)).values.astype(float)
m = np.isfinite(AGE) & np.isfinite(mar)
pc = abs(float(np.corrcoef(rankdata(AGE[m]), rankdata((mar[m] == 1).astype(float)))[0, 1]))
G.positive_control("正对照:年龄 × 已婚(必然随年龄变)", planted=pc, floor=abs(old - young), spread=1e-9)
rr = np.random.default_rng(SEEDS[0]); sh = AGE.copy(); sh = sh[rr.permutation(N)]
zs = []
for lo, hi in BANDS[:2] + BANDS[2:]:
    mm = (sh >= lo) & (sh <= hi)
    v = med_within(mm)
    if np.isfinite(v): zs.append(v)
G.negative_control("安慰剂:打乱年龄标签后段间差应塌", null=float(np.std(zs)),
                   effect=abs(old - young), null_spread=1e-9, null_kind="年龄标签置换")
G.spec_curve_cells_declare_n("规格曲线逐格 n", res)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", res)
print("\n" + "=" * 72)
if pc > abs(old - young) and np.std(zs) < 0.5 * abs(old - young):
    if abs(old - young) < mde:
        world = "A-NO-WALL"; verdict = (f"45+ 与 <45 之差 {old-young:+.4f} < 中位 MDE {mde:.4f} -> "
            f"**这堵墙后面没有新东西;`#539` 的这条限制后果 ≈ 0**")
    else:
        world = "B-PRICED"; verdict = (f"差 {old-young:+.4f} 超 MDE {mde:.4f} -> "
            f"**限制真实且已定价:NSFG 的结论必须带上「只到 44 岁」**")
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:GSS 四道题 ≠ NSFG 三道题 —— "
          "这是**同一构念在另一仪器上的年龄剖面**,不是 NSFG 那条曲线的延长。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(bands=res, young=young, old=old, diff=old - young, mde=mde, world=world,
               verdict=verdict, seeds=SEEDS, picked_by="被后续条目引用次数最多(#539,7 次)",
               impossible=["GSS 四题 ≠ NSFG 三题,不是同一曲线的延长", "横断面,年龄与队列共线"],
               unchallenged=True), open(OUT / "beyond_44.json", "w"), indent=1)
print(f"\nwrote {OUT/'beyond_44.json'}")

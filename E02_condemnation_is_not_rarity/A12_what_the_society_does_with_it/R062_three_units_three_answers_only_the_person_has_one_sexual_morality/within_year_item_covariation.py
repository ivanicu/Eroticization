"""E02·A213·R579 — 在同一年里,一个人的严厉是一把尺还是几把?

`#533` 的 NEXT,**改方向**:不再用时间差分(它在长窗口退化成趋势,短窗口淹在噪声里)。
行动类型:**FRONTIER**。这是 `#529`(社会)/`#532`(年代)同一个量在**第三个单位:人**上的版本,
而且**年份固定 ⇒ 不受长期趋势污染**,正是前两轮的共同弱点。

G1 ESTIMAND(先于方法):在**每一个调查年内**,
   `ρ_sex(t) = 四道性态度题两两个体层 Spearman 的中位`;
   `ρ_ref(t) = 同一份问卷里**非性**态度题两两的中位`(同年、同受访者、同作答模式)。
   **主量 = `ρ_sex(t) − ρ_ref(t)`**,即**扣掉同问卷共同方法方差之后**,性话题之间还多耦合多少。
   **次量 = 该差随年代的斜率**(严厉是否从多维收拢为一维)。

⚠ 硬规则 2:四题同问卷、同受访者、同一次作答 ⇒ **共同方法方差必然抬高 `ρ_sex`**。
   所以**绝不看 `ρ_sex` 本身**,只看它减去同问卷参照后的余量。参照必须是**态度题**,
   不能是人口学变量 —— 否则扣掉的不是方法方差,是别的东西。

WORLDS:
  W-ONE-RULER  `ρ_sex − ρ_ref` 明显 > 0 ⇒ **性话题之间另有一把共同的尺**,`#529`/`#532` 在人层不成立
  W-SAME-AS-ANY 差 ≈ 0 ⇒ 性话题之间的耦合**不比任意两道态度题更强** ⇒ 三个单位一致
  W-COLLAPSING 差随年代**上升** ⇒ 严厉正在从多维收拢为一维(一个关于历史的结论,不是关于结构的)
⚠ BASIN:`W-SAME-AS-ANY` 会让三个单位漂亮一致,所以**不是**本轮下注方向。本轮下注 `W-ONE-RULER`。

CONTROLS(G2):
  正对照 **堕胎题组**(`abany`/`abdefect`/`abnomore`/`abhlth`/`abpoor`/`abrape`/`absingle`)——
     公认的近重复题组,个体层必须强相关;它给出**这份问卷在个体层能达到的上限**;
  安慰剂 性态度题 × `zodiac`(星座)必须 ≈ 0(该是零 ⇒ negative_control);
  参照   非性态度题两两(`cappun` 死刑 · `grass` 大麻 · `letdie1` 安乐死 · `suicide1` 自杀 ·
     `fepol` 女性从政 · `natfare` 福利支出),**同年同人**;
  ⚠ 堕胎题**不进参照**:它本身与性道德高度相关,放进参照会把要测的东西扣掉。
KILL(条件式):if 正对照 > 参照中位 and 安慰剂 ≈ 0:
     差 > 参照分布的散度 -> W-ONE-RULER;差 ≈ 0 -> W-SAME-AS-ANY;斜率显著 -> W-COLLAPSING
   else UNVERIFIED
IMPOSSIBLE:同问卷同人 ⇒ 方法方差只能**部分**扣除(参照题的方法方差未必与性题相同)·
   一国一仪器 · 观察性非因果 · 未派对抗 agent ⇒ [unchallenged]
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
SEX = {"homosex": [1, 2, 3, 4], "premarsx": [1, 2, 3, 4], "xmarsex": [1, 2, 3, 4], "teensex": [1, 2, 3, 4]}
REF = {"cappun": [1, 2], "grass": [1, 2], "letdie1": [1, 2], "suicide1": [1, 2],
       "fepol": [1, 2], "natfare": [1, 2, 3]}
POS = {"abany": [1, 2], "abdefect": [1, 2], "abnomore": [1, 2], "abhlth": [1, 2],
       "abpoor": [1, 2], "abrape": [1, 2], "absingle": [1, 2]}
cols = ["year", "zodiac"] + list(SEX) + list(REF) + list(POS)
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=cols, convert_categoricals=False)
print("=== 硬规则 1:逐题 n 与年份跨度 ===")
for c in list(SEX) + list(REF) + list(POS):
    d = g[g[c].isin((SEX | REF | POS)[c])]
    print(f"  {c:9s} n={len(d):6d}  {int(d.year.min())}-{int(d.year.max())}  {d.year.nunique()} 年")

def rho(d, a, b, va, vb):
    m = d[a].isin(va) & d[b].isin(vb)
    if m.sum() < 200: return np.nan, int(m.sum())
    return float(np.corrcoef(rankdata(d[a][m]), rankdata(d[b][m]))[0, 1]), int(m.sum())

def med_within(d, spec):
    out = []
    for a, b in itertools.combinations(spec, 2):
        r, n = rho(d, a, b, spec[a], spec[b])
        if np.isfinite(r): out.append(abs(r))
    return (float(np.median(out)), len(out)) if out else (np.nan, 0)

rows = []
print("\n=== 逐年:性题中位 vs 同问卷非性参照中位(差 = 扣掉方法方差后的余量)===")
for y, d in g.groupby("year"):
    s, ks = med_within(d, SEX); r, kr = med_within(d, REF); p, kp = med_within(d, POS)
    if not (np.isfinite(s) and np.isfinite(r)): continue
    rows.append(dict(year=int(y), rho_sex=s, rho_ref=r, rho_pos=p if np.isfinite(p) else None,
                     diff=s - r, k_sex=ks, k_ref=kr, k_pos=kp, n=int(len(d)),
                     inclusion=[f"{int(y)} 年全部受访者 n={len(d)}", f"性题对 {ks} 个 · 参照对 {kr} 个",
                                "每对要求 n>=200", "同年同受访者,年份固定 -> 无趋势污染"]))
    print(f"  {int(y)}  n={len(d):5d}  性题中位={s:.4f}({ks}对)  参照中位={r:.4f}({kr}对)  "
          f"**差={s-r:+.4f}**  正对照(堕胎)={p:.4f}" if np.isfinite(p) else "")

DIFF = float(np.median([x["diff"] for x in rows]))
SEXM = float(np.median([x["rho_sex"] for x in rows]))
REFM = float(np.median([x["rho_ref"] for x in rows]))
POSM = float(np.median([x["rho_pos"] for x in rows if x["rho_pos"] is not None]))
SPREAD = float(np.std([x["rho_ref"] for x in rows]))
yrs = np.array([x["year"] for x in rows]); dfs = np.array([x["diff"] for x in rows])
slope = float(np.polyfit(yrs, dfs, 1)[0])
bs = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(2000):
        i = rng.integers(0, len(yrs), len(yrs))
        if np.ptp(yrs[i]) > 0: bs.append(np.polyfit(yrs[i], dfs[i], 1)[0])
SMDE = 2.8 * np.std(bs)
print(f"\n  **性题中位 {SEXM:.4f} · 参照中位 {REFM:.4f} · 差 {DIFF:+.4f}**  "
      f"(参照的年间散度 {SPREAD:.4f})")
print(f"  正对照(堕胎题组)中位 = {POSM:.4f} —— 这是本问卷个体层的上限")
print(f"  差随年代的斜率 = {slope:+.6f}/年  MDE={SMDE:.6f}  {'超MDE' if abs(slope)>SMDE else '看不见'}")

G = Gate("在同一年里,一个人的严厉是一把尺还是几把?(GSS 个体层)")
G.positive_control("正对照:堕胎题组(问卷个体层上限)", planted=POSM, floor=REFM, spread=1e-9)
zs = []
for a in SEX:
    r, n = rho(g, a, "zodiac", SEX[a], list(range(1, 13)))
    if np.isfinite(r): zs.append(abs(r))
G.negative_control("安慰剂:性态度题 × 星座", null=float(np.median(zs)), effect=POSM,
                   null_spread=float(np.std(zs)), null_kind="无关的个体层标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {str(x["year"]): x for x in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {str(x["year"]): x for x in rows})

print("\n" + "=" * 76)
if POSM > REFM and np.median(zs) < 0.5 * POSM:
    if abs(slope) > SMDE:
        world = "W-COLLAPSING"; verdict = (f"差随年代斜率 {slope:+.6f}/年 超 MDE {SMDE:.6f} -> "
            f"**性话题之间的额外耦合正在{'上升' if slope>0 else '下降'}**")
    elif DIFF > SPREAD:
        world = "W-ONE-RULER"; verdict = (f"差 {DIFF:+.4f} > 参照年间散度 {SPREAD:.4f} -> "
            f"**扣掉同问卷方法方差后,性话题之间仍另有一把共同的尺**")
    else:
        world = "W-SAME-AS-ANY"; verdict = (f"差 {DIFF:+.4f} <= 参照散度 {SPREAD:.4f} -> "
            f"**性话题之间的耦合不比任意两道态度题更强**")
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:参照题(死刑/大麻/安乐死/…)的**方法方差未必与性题相同** ——"
          "若性题彼此在题目措辞上更相似(都用「always wrong…not wrong at all」四级),"
          "那么差里有一部分是**格式相似**,不是道德结构。本轮无法把这两者分开。")
else:
    world, verdict = "UNVERIFIED", f"控制未齐 正对照={POSM:.4f} 参照={REFM:.4f} 安慰剂={np.median(zs):.4f}"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rows=rows, sex_median=SEXM, ref_median=REFM, diff_median=DIFF, pos_median=POSM,
               ref_spread=SPREAD, slope=slope, slope_MDE=float(SMDE), world=world, verdict=verdict,
               placebo=[float(x) for x in zs], seeds=SEEDS,
               instrument="GSS,同一份问卷同一批受访者,年份固定",
               impossible=["方法方差只能部分扣除:参照题格式与性题不同", "一国一仪器",
                           "观察性非因果"], unchallenged=True),
          open(OUT / "within_year_covariation.json", "w"), indent=1)
print(f"\nwrote {OUT/'within_year_covariation.json'}")

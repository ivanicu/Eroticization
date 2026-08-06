"""E02·A236·R619 — 「一个领域一个领域地动」,换一种领域划分还成立吗?

`#574` 的 NEXT。行动类型:**CLOSURE**(保护 `#557` 的结论,不分离新世界)。
`#557e` 自己写着:领域划分由我从题干文本读出,**未经第二人复核** ——
而 `#574` 刚证明:**同一个指派检验,池大小一变,答案会翻。** 所以这一格必须同样验。

预注册(写在跑之前):
  域内与跨域之差的**跨方案极差 < 0.05** -> `#557` 的结论**不带保留**;
  **≥ 0.05** -> 必须带上「领域划分敏感」并写进页面。
⚠ **而这一次先给极差本身做 bootstrap**(`#574d` 的教训:极差自己也有抽样波动,
  上一轮没给它误差棒,只能说「敏感」不能说「影响是 0.038」)。
方案:S1 严格按标题 · S2 把「心理」并入「物质」 · S3 剔除「交通」
方向统一沿用 `#557a`:从 format 库读码 1 的含义,统一成「高风险 = 1」。
IMPOSSIBLE:CLOSURE 不分离世界 · 青少年 · 校内抽样 · 未加权 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate
Y = ROOT / "data/external/yrbs"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
sas = (Y / "2023-SADC-SAS-Input-Program.sas").read_text(errors="replace")
fm = (Y / "2023-SADC-SAS-Formats-Program.sas").read_text(errors="replace")
pos = {m.group(1): (int(m.group(3)) - 1, int(m.group(4)))
       for m in re.finditer(r'^\s*(\w+)\s+(\$\s+)?(\d+)-(\d+)\s*$', sas, re.M)}
S1 = {"性活跃": ["q56", "q59"], "暴力": ["q12", "q16", "q13"], "物质": ["q33", "q42", "q48"],
      "交通": ["q8", "q9"], "心理": ["q26", "q27", "q29"]}
ITEMS = [q for v in S1.values() for q in v]
def code1(q):
    m = re.search(r'value\s+\$H' + q[1:] + r'S(.*?)(?=\nvalue\s|\Z)', fm, re.S | re.I)
    if not m: return None
    o = re.search(r'"1"\s*=\s*"([^"]*)"', m.group(1))
    return o.group(1) if o else None
YES1 = {q for q in ITEMS if (c := code1(q)) and not re.match(r'^(never|none|0 |did not|no )', c.strip(), re.I)}
print(f"=== 方向(从 format 库读)——码1=高风险的题:{sorted(YES1)} ===")
NEED = ["year", "sex"] + ITEMS
rec = {n: [] for n in NEED}
for f in sorted(Y.glob("sadc_2023_state_*.dat")):
    for line in open(f, errors="replace"):
        for n in NEED:
            s, e = pos[n]; v = line[s:e].strip()
            rec[n].append(float(v) if v not in ("", ".") else np.nan)
A = {n: np.array(rec[n]) for n in NEED}
def series(q, sx):
    m = np.isfinite(A[q]) & (A["sex"] == sx)
    risk = (lambda v: (v == 1) if q in YES1 else (v > 1))
    return {y: float(risk(A[q][m & (A["year"].astype(int) == y)]).mean())
            for y in sorted(set(A["year"][m].astype(int)))
            if (m & (A["year"].astype(int) == y)).sum() >= 300}
SER = {(q, sx): series(q, sx) for q in ITEMS for sx in (1, 2)}
def dc(a, b, yrs=None):
    yr = sorted(set(a) & set(b)) if yrs is None else [y for y in yrs if y in a and y in b]
    if len(yr) < 8: return np.nan
    va = np.array([a[y] for y in yr]); vb = np.array([b[y] for y in yr])
    da, db = np.diff(va), np.diff(vb)
    if np.std(da) == 0 or np.std(db) == 0: return np.nan
    return float(np.corrcoef(da, db)[0, 1])
def gap(MAP, boot=None):
    """域内中位 − 跨域中位。boot: 年份的 bootstrap 索引。"""
    inn, out = [], []
    yrs = sorted({y for q in ITEMS for y in SER[(q, 1)]})
    use = [yrs[i] for i in boot] if boot is not None else None
    for a, b in itertools.combinations([q for v in MAP.values() for q in v], 2):
        da = [dc(SER[(a, i)], SER[(b, j)], use) for i, j in ((1, 2), (2, 1))]
        da = [x for x in da if np.isfinite(x)]
        if not da: continue
        r = float(np.mean(da))
        dom = {q: d for d, qs in MAP.items() for q in qs}
        (inn if dom[a] == dom[b] else out).append(r)
    if not inn or not out: return np.nan
    return float(np.median(inn) - np.median(out))
S2 = {k: v for k, v in S1.items() if k not in ("心理", "物质")}; S2["物质"] = S1["物质"] + S1["心理"]
S3 = {k: v for k, v in S1.items() if k != "交通"}
SCH = {"S1 严格按标题": S1, "S2 心理并入物质": S2, "S3 剔除交通": S3}
print("\n=== 三种领域划分 ===")
g = {}
for k, M in SCH.items():
    g[k] = gap(M)
    print(f"  {k:14s} 域内−跨域 = **{g[k]:+.4f}**")
rng_ = max(g.values()) - min(g.values())
print(f"\n  **跨方案极差 = {rng_:.4f}**(预注册门槛 0.05)")
# ⚠ #574d 的教训:给极差本身做 bootstrap
yrs = sorted({y for q in ITEMS for y in SER[(q, 1)]})
bs = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(120):
        idx = sorted(rng.choice(len(yrs), len(yrs), replace=True))
        gg = [gap(M, boot=idx) for M in SCH.values()]
        gg = [x for x in gg if np.isfinite(x)]
        if len(gg) == len(SCH): bs.append(max(gg) - min(gg))
bs = np.array(bs)
print(f"  极差的 bootstrap:中位 {np.median(bs):.4f} · 95% CI [{np.quantile(bs,.025):.4f},{np.quantile(bs,.975):.4f}]")
G = Gate("「一个领域一个领域地动」,换一种领域划分还成立吗?")
cells = {k: dict(n=len(yrs), gap=g[k], inclusion=[k, f"{len(yrs)} 个调查年", "方向由 format 库统一"]) for k in SCH}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
G.positive_control("三种方案都算得出(非退化)", planted=float(len([x for x in g.values() if np.isfinite(x)])),
                   floor=2.5, spread=1e-9)
print("\n" + "=" * 72)
sens = rng_ >= 0.05
verdict = (f"极差 {rng_:.4f} {'≥' if sens else '<'} 0.05 -> "
           + ("**领域划分敏感,`#557` 的结论必须带上这条**" if sens else "**领域划分不承重,`#557` 不带保留**"))
print(f"CLOSURE 结论:{verdict}")
print(f"⚠ 而极差本身的 95% CI 是 [{np.quantile(bs,.025):.4f},{np.quantile(bs,.975):.4f}] —— "
      f"{'它跨过门槛,所以「敏感」这个判决本身也不稳' if np.quantile(bs,.025) < 0.05 < np.quantile(bs,.975) else '它整段在门槛一侧,判决稳'}")
print(G)
json.dump(dict(gaps=g, range=rng_, threshold=0.05, sensitive=bool(sens),
               range_boot_median=float(np.median(bs)),
               range_boot_ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
               yes1=sorted(YES1), n_years=len(yrs), verdict=verdict, action_type="CLOSURE",
               impossible=["CLOSURE 不分离世界", "青少年", "校内抽样", "未加权"],
               unchallenged=True), open(OUT / "domain_sensitivity.json", "w"), indent=1)
print(f"\nwrote {OUT/'domain_sensitivity.json'}")

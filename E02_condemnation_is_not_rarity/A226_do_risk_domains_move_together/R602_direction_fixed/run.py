"""E02·A226·R602 — 方向统一之后,风险领域是一起动还是一件一件?

`#556` 的 NEXT。行动类型:**FRONTIER**。
**编码方向已从 format 库读出(对象,不是记忆)**:
  `$H56S` 1=Yes · `$H26S` 1=Yes  -> **码 1 就是高风险**
  `$H57S`/`$H58S`/`$H59S`/`$H61S` 1=Never had sex · `$H12S`/`$H33S` 1=0 days -> **码 1 是低风险**
⇒ 统一成 **「高风险 = 1」**:`YES1` 组用 `x==1`,其余用 `x>1`(即不属于「从不/0 天」那一档)。
⚠ `#556c`:`#554b` 的 `+0.8849` 是**符号不一致集合上的中位**;本轮按新方向**重报**。

判决逻辑的修(`#556b`):**对绝对值判,并先 `assert 域内 > 0`** —— 域内为负则该分类作废,不进三分。
发布规则(沿用 `#556` NEXT ④):网格至少三格可判**且方向已核对**才上页面。
IMPOSSIBLE:青少年 · 校内抽样 · 未加权 · 领域划分由我从标签读出,未经复核 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
Y = ROOT / "data/external/yrbs"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
sas = (Y / "2023-SADC-SAS-Input-Program.sas").read_text(errors="replace")
fm = (Y / "2023-SADC-SAS-Formats-Program.sas").read_text(errors="replace")
pos = {m.group(1): (int(m.group(3)) - 1, int(m.group(4)))
       for m in re.finditer(r'^\s*(\w+)\s+(\$\s+)?(\d+)-(\d+)\s*$', sas, re.M)}
DOM = {"性活跃": ["q56", "q59"], "暴力": ["q12", "q16", "q13"], "物质": ["q33", "q42", "q48"],
       "交通": ["q8", "q9"], "心理": ["q26", "q27", "q29"]}
ACT_EXTRA = ["q57", "q58", "q61"]   # 仅用于 #554b 的重报,不进领域网格
ITEMS = [q for v in DOM.values() for q in v] + ACT_EXTRA

def code1(q):
    """从 format 库读码 1 的含义 —— 对象,不是记忆。"""
    m = re.search(r'value\s+\$H' + q[1:] + r'S(.*?)(?=\nvalue\s|\Z)', fm, re.S | re.I)
    if not m: return None
    one = re.search(r'"1"\s*=\s*"([^"]*)"', m.group(1))
    return one.group(1) if one else None
print("=== 码 1 的含义(逐题从 format 库读)与方向判定 ===")
YES1 = set()
for q in ITEMS:
    c = code1(q)
    hi = bool(c) and not re.match(r'^(never|none|0 |did not|no )', c.strip(), re.I)
    if hi: YES1.add(q)
    print(f"  {q:4s} 码1 = {str(c)[:34]:36s} -> {'**码1=高风险**' if hi else '码1=低风险(用 x>1)'}")
NEED = ["year", "sex"] + ITEMS
rec = {n: [] for n in NEED}
for f in sorted(Y.glob("sadc_2023_state_*.dat")):
    for line in open(f, errors="replace"):
        for n in NEED:
            s, e = pos[n]; v = line[s:e].strip()
            rec[n].append(float(v) if v not in ("", ".") else np.nan)
A = {n: np.array(rec[n]) for n in NEED}
def risk(q, v):   # 统一成「高风险 = 1」
    return (v == 1) if q in YES1 else (v > 1)
def series(q, sx):
    m = np.isfinite(A[q]) & (A["sex"] == sx)
    return {y: float(risk(q, A[q][m & (A["year"].astype(int) == y)]).mean())
            for y in sorted(set(A["year"][m].astype(int)))
            if (m & (A["year"].astype(int) == y)).sum() >= 300}
keep = {d: [q for q in qs if len(series(q, 1)) >= 8] for d, qs in DOM.items()}
keep = {d: v for d, v in keep.items() if len(v) >= 2}
K = [q for v in keep.values() for q in v]
S = {(q, sx): series(q, sx) for q in K for sx in (1, 2)}
def dc(a, b, perm=None):
    yr = sorted(set(a) & set(b))
    if len(yr) < 8: return np.nan
    va = np.array([a[y] for y in yr]); vb = np.array([b[y] for y in yr])
    if perm is not None: vb = vb[perm.permutation(len(yr))]
    da, db = np.diff(va), np.diff(vb)
    if np.std(da) == 0 or np.std(db) == 0: return np.nan
    return float(np.corrcoef(da, db)[0, 1])
def pm(pairs):
    o = []
    for a, b in pairs:
        v = [dc(S[(a, i)], S[(b, j)]) for i, j in ((1, 2), (2, 1))]
        v = [x for x in v if np.isfinite(x)]
        if v: o.append(float(np.mean(v)))
    return (float(np.median(o)), len(o), o) if o else (np.nan, 0, [])
same = [dc(S[(q, 1)], S[(q, 2)]) for q in K]
SM = float(np.median([x for x in same if np.isfinite(x)]))
inn = [(a, b) for v in keep.values() for a, b in itertools.combinations(v, 2)]
out = [(a, b) for d1, d2 in itertools.combinations(keep, 2) for a in keep[d1] for b in keep[d2]]
IN, kin, inv = pm(inn); OU, kout, ouv = pm(out)
shuf = []
for sd in SEEDS:
    r = np.random.default_rng(sd)
    for _ in range(40):
        for a, b in out:
            x = dc(S[(a, 1)], S[(b, 2)], perm=r)
            if np.isfinite(x): shuf.append(abs(x))
Q95 = float(np.quantile(shuf, .95))
print(f"\n=== 方向统一后 ===\n  同题上限 {SM:+.4f} · 安慰剂 q95 {Q95:.4f}")
print(f"  **域内 {IN:+.4f}**({kin} 对) · **跨域 {OU:+.4f}**({kout} 对) · 比值 {OU/IN if IN else float('nan'):.3f}")
print(f"  域内逐值:{[f'{x:+.3f}' for x in sorted(inv)]}")
# `#554b` 按新方向重报
ACT = ["q56", "q57", "q58", "q59", "q61"]
Sa = {(q, sx): series(q, sx) for q in ACT for sx in (1, 2)}
S.update(Sa)
am, ak, av = pm(list(itertools.combinations(ACT, 2)))
print(f"\n  **`#554b` 重报**:五道性行为题两两,方向统一后中位 = **{am:+.4f}**({ak} 对),"
      f"全部值 {[f'{x:+.3f}' for x in sorted(av)]}")
print("\n" + "=" * 74)
if IN <= 0:
    world, verdict = "UNVERIFIED", f"**域内中位 {IN:+.4f} ≤ 0 -> 这个领域划分作废,不进三分**(`#556b` 的修)"
    print(f"⚠ {verdict}")
elif abs(SM) > Q95:
    world = ("W-ONE-FRONT" if abs(OU) >= 0.7 * abs(IN) else
             ("W-BY-DOMAIN" if abs(OU) < Q95 else "W-PARTIAL"))
    verdict = f"{world}: 域内 {IN:+.4f} · 跨域 {OU:+.4f} · 比值 {abs(OU)/abs(IN):.3f}"
    print(f"控制齐备 ⇒ 评判。**{world}** —— {verdict}")
    print("⚠ 这个 KILL 会怎样失败:**领域的划分是我从标签读出来的**,未经第二人复核。")
else:
    world, verdict = "UNVERIFIED", f"上限 {SM:.4f} 未越过安慰剂 q95 {Q95:.4f}"
    print(f"⚠ {verdict}")
json.dump(dict(yes1=sorted(YES1), domains=keep, ceiling=SM, within=IN, between=OU,
               k_within=kin, k_between=kout, placebo_q95=Q95, within_values=inv,
               activity_remeasured=dict(median=am, k=ak, values=av),
               world=world, verdict=verdict, seeds=SEEDS,
               direction_source="2023-SADC-SAS-Formats-Program.sas 的 $H..S 块,逐题读出",
               impossible=["青少年", "校内抽样", "未加权", "领域划分未经第二人复核"],
               unchallenged=True), open(OUT / "direction_fixed.json", "w"), indent=1)
print(f"\nwrote {OUT/'direction_fixed.json'}")

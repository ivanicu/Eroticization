"""E02·A226·R601 — 青少年的各个风险领域,是一起松开的,还是一件一件?

`#555` 的 NEXT。行动类型:**FRONTIER**。
`#555a` 证明:五道性行为题不是五种行为 ⇒ 用它们回答「一起还是一件一件」是**循环的**。
本轮换成**彼此不下游的五个领域**(由**标签文本**分类,不由题号):

  性活跃  q56 曾有性交 · q59 当前性活跃            (`#555` 已证:同一个量)
  暴力    q12 校内携带武器 · q16 打架 · q13 携枪
  物质    q33 当前吸烟 · q42 当前饮酒 · q48 当前大麻
  交通    q8 系安全带 · q9 坐酒后司机的车
  心理    q26 悲伤无望 · q27 考虑自杀 · q29 自杀未遂

**这是 `#529`(社会:严厉是否跨做法耦合)与 `#532`(年代:态度是否同步)
  在「行为 × 年代」上的形式** —— 而这一次,各项**确实可以独立变动**。

G1 ESTIMAND(先于方法):`ρ_in` = **同一领域内**两题的 Δ 相关中位;
   `ρ_out` = **跨领域**两题的 Δ 相关中位。**主量 = 两者之差。**
   噪声配平与 `#577`/`#599` 同:**男的题 A × 女的题 B**,受访者不相交;上限 = 同题男×女。
WORLDS:
  W-ONE-FRONT   跨领域 ≈ 同领域 ⇒ **所有风险行为一起动**(一条战线)
  W-BY-DOMAIN   跨领域 << 同领域,且跨领域 ≈ 安慰剂 ⇒ **一个领域一个领域地动**
  W-PARTIAL     介于两者 ⇒ 报比值
⚠ BASIN:`W-BY-DOMAIN` 与 `#529`/`#532` 一路呼应,**所以不是下注方向**。本轮下注 `W-ONE-FRONT`。
CONTROLS:上限 同题男Δ×女Δ · 安慰剂 打乱年份 · **年份覆盖 <8 年的题一律剔除(硬规则 1,先打印)**
发布规则(预注册):**网格至少三格可判**才上页面;否则只进账本。
KILL(条件式):if 上限 > 安慰剂 q95: 按三分判 else UNVERIFIED
IMPOSSIBLE:青少年 · 校内抽样 · 未加权 · **领域的划分由我从标签读出,未经第二人复核** · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, re
Y = ROOT / "data/external/yrbs"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
sas = (Y / "2023-SADC-SAS-Input-Program.sas").read_text(errors="replace")
lab = dict(re.findall(r'^\s*(\w+)\s*=\s*"([^"]*)"', sas, re.M))
pos = {m.group(1): (int(m.group(3)) - 1, int(m.group(4)))
       for m in re.finditer(r'^\s*(\w+)\s+(\$\s+)?(\d+)-(\d+)\s*$', sas, re.M)}
DOM = {"性活跃": ["q56", "q59"], "暴力": ["q12", "q16", "q13"], "物质": ["q33", "q42", "q48"],
       "交通": ["q8", "q9"], "心理": ["q26", "q27", "q29"]}
ITEMS = [q for v in DOM.values() for q in v]
NEED = ["year", "sex"] + ITEMS
assert all(n in pos for n in NEED), [n for n in NEED if n not in pos]
rec = {n: [] for n in NEED}
for f in sorted(Y.glob("sadc_2023_state_*.dat")):
    for line in open(f, errors="replace"):
        for n in NEED:
            s, e = pos[n]; v = line[s:e].strip()
            rec[n].append(float(v) if v not in ("", ".") else np.nan)
A = {n: np.array(rec[n]) for n in NEED}
print(f"=== 硬规则 1:{len(A['year']):,} 行,逐题年份与 n(<8 年一律剔除)===")
keep = {}
for d, qs in DOM.items():
    kq = []
    for q in qs:
        ok = np.isfinite(A[q]); ys = sorted(set(A["year"][ok].astype(int)))
        good = len(ys) >= 8
        print(f"  [{d:4s}] {q:4s} n={int(ok.sum()):8,} 年 {len(ys):2d} 个 "
              f"{ys[0]}–{ys[-1]}  {'✅' if good else '⛔剔除'}  {lab.get(q,'')[:28]}")
        if good: kq.append(q)
    if len(kq) >= 2: keep[d] = kq
print(f"⇒ 保留 {len(keep)} 个领域:{ {k: v for k, v in keep.items()} }\n")
K = [q for v in keep.values() for q in v]

def series(q, sx):
    m = np.isfinite(A[q]) & (A["sex"] == sx)
    return {y: float((A[q][m & (A["year"].astype(int) == y)] == 1).mean())
            for y in sorted(set(A["year"][m].astype(int)))
            if (m & (A["year"].astype(int) == y)).sum() >= 300}
S = {(q, sx): series(q, sx) for q in K for sx in (1, 2)}
def dc(a, b, perm=None):
    yr = sorted(set(a) & set(b))
    if len(yr) < 8: return np.nan
    va = np.array([a[y] for y in yr]); vb = np.array([b[y] for y in yr])
    if perm is not None: vb = vb[perm.permutation(len(yr))]
    da, db = np.diff(va), np.diff(vb)
    if np.std(da) == 0 or np.std(db) == 0: return np.nan
    return float(np.corrcoef(da, db)[0, 1])
def pairmed(pairs):
    out = []
    for a, b in pairs:
        v = [dc(S[(a, i)], S[(b, j)]) for i, j in ((1, 2), (2, 1))]
        v = [x for x in v if np.isfinite(x)]
        if v: out.append(float(np.mean(v)))
    return (float(np.median(out)), len(out)) if out else (np.nan, 0)
same = [dc(S[(q, 1)], S[(q, 2)]) for q in K]
SM = float(np.median([x for x in same if np.isfinite(x)]))
inn = [(a, b) for v in keep.values() for a, b in itertools.combinations(v, 2)]
out = [(a, b) for d1, d2 in itertools.combinations(keep, 2) for a in keep[d1] for b in keep[d2]]
IN, kin = pairmed(inn); OUTv, kout = pairmed(out)
shuf = []
for sd in SEEDS:
    r = np.random.default_rng(sd)
    for _ in range(40):
        for a, b in out:
            x = dc(S[(a, 1)], S[(b, 2)], perm=r)
            if np.isfinite(x): shuf.append(abs(x))
Q95 = float(np.quantile(shuf, .95))
print(f"=== 结果 ===\n  同题上限 = **{SM:+.4f}**  安慰剂 q95 = {Q95:.4f}")
print(f"  **领域内 = {IN:+.4f}**({kin} 对) · **跨领域 = {OUTv:+.4f}**({kout} 对) · "
      f"差 = **{IN-OUTv:+.4f}** · 比值 {OUTv/IN if IN else float('nan'):.3f}")
print("\n  逐领域对(全格公布):")
cells = {}
for d1, d2 in itertools.combinations(keep, 2):
    m, n = pairmed([(a, b) for a in keep[d1] for b in keep[d2]])
    cells[f"{d1}×{d2}"] = dict(median=None if not np.isfinite(m) else m, k=n,
                               inclusion=[f"{d1} 与 {d2} 的全部跨域对", "每题年份≥8", "每年每性别 n≥300"])
    print(f"    {d1:4s} × {d2:4s}  {m:+.4f}({n} 对)")
for d in keep:
    m, n = pairmed(list(itertools.combinations(keep[d], 2)))
    cells[f"{d}(域内)"] = dict(median=None if not np.isfinite(m) else m, k=n,
                              inclusion=[f"{d} 域内全部对", "每题年份≥8", "每年每性别 n≥300"])
    print(f"    {d:4s} 域内     {m:+.4f}({n} 对)")
judg = sum(1 for v in cells.values() if v["median"] is not None)
print("\n" + "=" * 74)
if abs(SM) > Q95:
    world = ("W-ONE-FRONT" if OUTv >= 0.7 * IN else
             ("W-BY-DOMAIN" if abs(OUTv) < Q95 else "W-PARTIAL"))
    verdict = f"{world}: 域内 {IN:+.4f} · 跨域 {OUTv:+.4f} · 比值 {OUTv/IN:.3f} · 安慰剂 q95 {Q95:.4f}"
    print(f"控制齐备 ⇒ 评判。**{world}** —— {verdict}")
    print(f"  发布规则(预注册):可判格 {judg} 个 -> {'**上页面**' if judg >= 3 else '**只进账本**'}")
    print("⚠ 这个 KILL 会怎样失败:**领域的划分是我从标签读出来的**;"
          "「打架」与「携带武器」算一个领域是一个判断,而它直接决定了「域内」这个量。")
else:
    world, verdict = "UNVERIFIED", f"上限 {SM:.4f} 未越过安慰剂 q95 {Q95:.4f}"
    print(f"⚠ {verdict}")
json.dump(dict(domains=keep, same_ceiling=SM, within=IN, between=OUTv, k_within=kin, k_between=kout,
               ratio=OUTv / IN if IN else None, placebo_q95=Q95, cells=cells, judgeable=judg,
               world=world, verdict=verdict, seeds=SEEDS,
               prereg="网格至少三格可判才上页面",
               impossible=["青少年不可外推", "校内抽样", "未加权", "领域划分由我从标签读出,未经复核"],
               unchallenged=True), open(OUT / "domain_grid.json", "w"), indent=1)
print(f"\nwrote {OUT/'domain_grid.json'}")

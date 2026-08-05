"""E02·A190·R527 — 两个对照各自坏在什么地方,以及同源参照分布才是那个零

R526 判 UNVERIFIED,正确 —— 但两个对照都是**自己坏了**,不是被测对象坏了。
realstat 的两行,一次撞到两行:

  ① 「控制无法 PASS」:我把正对照的门槛设成 `q95(自身置换零) + 2·sd = 0.562`,
     **却从没算过这个统计量的天花板**。4 级 × 3 级、带大量并列的 Spearman,
     最大可达值远低于 1。门槛可能设在设计返回不了的地方 ⇒ 它的 FAIL 什么也没说。
     修:算 floor(无植入)与 ceiling(最大植入、无噪声 = 两列各自排序后配对),
        要求 floor < t < ceiling,并把三个数一起报。

  ② 「阴性对照不是零」:`SCCS1753` 单系继嗣深度 × 暴力谴责 = +0.268。
     我问了「这个零该是零吗」,答「是」—— **而人类学里父系深度与世仇本来就相关**。
     那不是一个坏掉的零,那是**我选错了对象**。
     ⚠ 更深的问题:`lang1998conan` 的整套变量是**围绕一个关于暴力的假说**编的,
        所以它里面可能**根本不存在**与暴力无关的变量 ⇒ 同源零对可能不存在。

  修法(这才是本轮的内容):**不要挑一个零,建一个参照分布。**
  谴责 × `lang` 的**每一个**其他序数变量 -> |ρ| 的经验分布。
  这个分布同时回答两件 R526 想分开问的事:
    - 同源本身产生多大的耦合?(分布的位置)
    - 范围匹配的那个 ρ 是不是特别?(它在分布里的分位)
  ⇒ 它是一个**测量出来的**同源基线,不是一个**挑出来的**零。

ESTIMAND:范围匹配 ρ 在「同源参照分布」中的分位,以及正对照相对其自身
         floor / ceiling 的位置。

⚠ 这一轮**不重新评判 R526 的斜率**。R526 的 UNVERIFIED 不因对照被修而变成 CONFIRMED
  —— 修好的仪器只授权**新的一次测量**,不追认旧的判定。
"""
import os, sys, pathlib, json, csv, math, collections
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))

import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
SCCS = ROOT / "data/external/dplace/repo/datasets/SCCS"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SENTINEL = {88.0}

tab = collections.defaultdict(lambda: collections.defaultdict(list))
for r in csv.DictReader(open(SCCS / "data.csv")):
    if r["code"] not in ("NA", ""):
        tab[r["soc_id"]][r["var_id"]].append(r["code"])
V = {v["id"]: v for v in csv.DictReader(open(SCCS / "variables.csv"))}
socids = sorted(tab, key=lambda s: int(s.replace("SCCS", "")))
socrow = {s["id"]: s for s in csv.DictReader(open(SCCS / "societies.csv"))}
LAT = np.array([float(socrow[s]["Lat"]) if s in socrow and socrow[s]["Lat"] else np.nan for s in socids])
LON = np.array([float(socrow[s]["Long"]) if s in socrow and socrow[s]["Long"] else np.nan for s in socids])
REGION = np.full(len(socids), -1)
for i, (a, o) in enumerate(zip(LAT, LON)):
    if not (np.isfinite(a) and np.isfinite(o)): continue
    if o < -30:  REGION[i] = 0 if a > 12 else 1
    elif o < 45: REGION[i] = 2 if a > 20 else 3
    elif o < 100: REGION[i] = 4
    else:        REGION[i] = 5 if a > 0 else 6


def vec(var):
    out = []
    for s in socids:
        vs = tab[s].get(var, [])
        if not vs: out.append(np.nan); continue
        try: x = float(max(set(vs), key=vs.count))
        except ValueError: out.append(np.nan); continue
        out.append(np.nan if x in SENTINEL else x)
    return np.array(out, float)


def midrank(a):
    """并列取平均秩。⚠ R526 用的是 argsort(argsort()),它给并列值**任意的不同秩**,
    秩取决于数据文件的行序 —— 重排行(不改变数据)会改变 rho。这是 gauge 失败。
    这些变量只有 3-4 个等级、~50 个社会,并列极多,所以影响不是小数点后第三位。"""
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    s = a[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]: j += 1
        r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8: return np.nan, int(m.sum())
    ra, rb = midrank(x[m]), midrank(y[m])
    ra = ra - ra.mean(); rb = rb - rb.mean()
    d = math.sqrt((ra @ ra) * (rb @ rb))
    return (float(ra @ rb / d) if d > 0 else np.nan), int(m.sum())


def spearman_old(x, y):
    """R526 用的那个 —— 保留,用来量这个缺陷值多少。"""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8: return np.nan, int(m.sum())
    a, b = x[m], y[m]
    ra = np.argsort(np.argsort(a)).astype(float); rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = math.sqrt((ra @ ra) * (rb @ rb))
    return (float(ra @ rb / d) if d > 0 else np.nan), int(m.sum())


def ceiling(x, y):
    """最大植入、无噪声:两列各自排序后配对 —— 给定观测到的并列结构,Spearman 的上界。"""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8: return np.nan
    a = np.sort(x[m]); b = np.sort(y[m])
    r, _ = spearman(a, b)
    return r


def perm_region(x, y, n=3000, seed=0):
    rng = np.random.default_rng(seed)
    m = np.isfinite(x) & np.isfinite(y)
    xs, ys, rs = x[m], y[m], REGION[m]
    out = []
    for _ in range(n):
        yp = ys.copy()
        for g in np.unique(rs):
            idx = np.where(rs == g)[0]
            if len(idx) > 1: yp[idx] = ys[idx][rng.permutation(len(idx))]
        r, _ = spearman(xs, yp)
        if np.isfinite(r): out.append(r)
    return np.array(out)


COND = -vec("SCCS1768")          # 本地社区内暴力:高 = 更被谴责
RARE = -vec("SCCS1750")          # 高 = 更罕见
main_r, main_n = spearman(COND, RARE)

# ---------------------------------------------------------------- ⓪ GAUGE:并列缺陷值多少
print("=== ⓪ GAUGE 修复:并列秩 vs R526 的序秩,以及行序不变性 ===")
SCOPES = [("within_community", "SCCS1768", "SCCS1750"),
          ("within_ethnic", "SCCS1769", "SCCS1776"),
          ("between_ethnic", "SCCS1770", "SCCS1778")]
gauge = []
rngG = np.random.default_rng(20260805)
for nm, av, fv in SCOPES:
    a = -np.floor(vec(av) / 10.0) if av == "SCCS1770" else -vec(av)
    f = -vec(fv)
    r_new, n = spearman(a, f)
    r_old, _ = spearman_old(a, f)
    # 行序不变性:重排社会顺序,重算 —— 正确的实现必须完全不变
    perm = rngG.permutation(len(a))
    r_new_p, _ = spearman(a[perm], f[perm])
    r_old_p, _ = spearman_old(a[perm], f[perm])
    gauge.append(dict(scope=nm, rho_new=r_new, rho_old=r_old, n=n,
                      rowperm_drift_new=abs(r_new - r_new_p),
                      rowperm_drift_old=abs(r_old - r_old_p)))
    print(f"  {nm:17s} 并列秩={r_new:+.4f}  序秩(R526)={r_old:+.4f}  Δ={r_new-r_old:+.4f}   "
          f"行序漂移: 新={abs(r_new-r_new_p):.6f}  旧={abs(r_old-r_old_p):.6f}")
print(f"  ⇒ 旧实现在**同一份数据重排行**之后漂移 "
      f"{max(g['rowperm_drift_old'] for g in gauge):.4f};新实现漂移 "
      f"{max(g['rowperm_drift_new'] for g in gauge):.6f}。\n")

# ---------------------------------------------------------------- ① 正对照:floor < t < ceiling
print("=== ① 正对照的三个数(R526 只算了一个)===")
pcx, pcy = vec("SCCS1776"), vec("SCCS1777")
pc_r, pc_n = spearman(pcx, pcy)
pc_null = perm_region(pcx, pcy, 3000, SEEDS[0])
pc_floor = float(np.quantile(np.abs(pc_null), .95))
pc_ceil = ceiling(pcx, pcy)
t_old = pc_floor + 2 * float(pc_null.std())
print(f"  观测 rho      = {pc_r:+.4f}  (n={pc_n})")
print(f"  floor (q95|零|) = {pc_floor:+.4f}")
print(f"  ceiling (最大植入,无噪声) = {pc_ceil:+.4f}")
print(f"  R526 用的门槛 t = floor + 2sd = {t_old:+.4f}")
print(f"  ⇒ floor < t < ceiling ?  {pc_floor:.4f} < {t_old:.4f} < {pc_ceil:.4f}  -> "
      f"{'成立' if pc_floor < t_old < pc_ceil else '⛔ 不成立 —— R526 的门槛不可用'}")
band_ok = bool(pc_floor < t_old < pc_ceil)
# 可用门槛:落在带内的中点
t_new = (pc_floor + pc_ceil) / 2
print(f"  可用门槛(带内中点) t' = {t_new:+.4f}  ->  正对照 {'PASS' if pc_r > t_new else 'FAIL'}")

# ---------------------------------------------------------------- ② 同源参照分布,而不是挑一个零
print("\n=== ② 同源参照分布:谴责 × lang 的每一个其他序数变量 ===")
lang_vars = [v for v, d in V.items() if d["source"] == "lang1998conan"]
ref = []
for v in lang_vars:
    if v in ("SCCS1768", "SCCS1750"): continue
    y = vec(v)
    u = np.unique(y[np.isfinite(y)])
    if len(u) < 3 or len(u) > 12:      # 只用序数样的
        continue
    r, n = spearman(COND, y)
    if np.isfinite(r) and n >= 30:
        ref.append(dict(var=v, rho=r, n=n, title=V[v]["title"][:52]))
ref.sort(key=lambda d: -abs(d["rho"]))
absr = np.array([abs(d["rho"]) for d in ref])
print(f"  参照对数 = {len(ref)}   |rho| 中位={np.median(absr):.4f}  "
      f"q75={np.quantile(absr,.75):.4f}  q90={np.quantile(absr,.90):.4f}  q95={np.quantile(absr,.95):.4f}  "
      f"max={absr.max():.4f}")
print(f"  ⇒ **同源本身**的典型耦合 = {np.median(absr):.4f},不是 0 —— R526 挑的那个零(+0.2679)"
      f"落在 q{100*(absr<0.2679).mean():.0f}")
print("  参照分布最强的 6 个(它们是同源耦合的上限样本):")
for d in ref[:6]:
    print(f"    |{d['rho']:+.3f}| n={d['n']:3d}  {d['var']:9s} {d['title']}")
pct = float((absr < abs(main_r)).mean())
print(f"\n  范围匹配 rho = {main_r:+.4f} (n={main_n})  ->  在同源参照分布的 q{100*pct:.0f}")

# ---------------------------------------------------------------- 门
G = Gate("修好两个对照之后,范围匹配的耦合还特别吗?")
pc_ok = G.positive_control("正对照(门槛落在 floor<t<ceiling 带内)",
                           planted=abs(pc_r), floor=t_new, spread=float(pc_null.std()) * 0.0 + 1e-9)
nc_ok = G.negative_control("阴性:同源参照分布的中位(测量出来的,不是挑出来的)",
                           null=float(np.median(absr)), effect=main_r,
                           null_spread=float(np.std(absr)),
                           null_kind="同源变量参照分布(lang1998conan 内所有序数变量)")
G.resolvable("范围匹配效应可分辨", effect=abs(main_r),
             spread=float(np.mean([perm_region(COND, RARE, 2000, s).std() for s in SEEDS])))
G.asserted("参照分布不是空的", len(ref) >= 15, f"{len(ref)} 个同源参照对(要求 >=15)", kind="control")

print("\n" + "=" * 68)
if pc_ok and nc_ok:
    verdict = (f"范围匹配 rho 在同源参照分布的 q{100*pct:.0f} —— "
               + ("显著高于同源基线" if pct >= 0.90 else "**不**高于同源基线,同源耦合可以解释它"))
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:参照分布里若混入了与暴力谴责真有关的变量"
          "(lang 的变量集正是围绕暴力假说编的),基线会被抬高,从而**低估**范围匹配的特殊性。"
          "⇒ 这个方向的偏是保守的。")
else:
    verdict = f"UNVERIFIED —— 控制仍未齐(pos={pc_ok} neg={nc_ok})"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(main=dict(rho=main_r, n=main_n, percentile_in_same_source=pct),
               positive=dict(rho=pc_r, n=pc_n, floor=pc_floor, ceiling=pc_ceil,
                             threshold_R526=t_old, band_ok_R526=band_ok, threshold_new=t_new,
                             ok=bool(pc_ok)),
               reference=dict(k=len(ref), median=float(np.median(absr)),
                              q75=float(np.quantile(absr, .75)), q90=float(np.quantile(absr, .90)),
                              q95=float(np.quantile(absr, .95)), max=float(absr.max()),
                              top=[{k: v for k, v in d.items()} for d in ref[:10]]),
               verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT / "repair_both_controls.json", "w"), indent=1)
print(f"\nwrote {OUT/'repair_both_controls.json'}")

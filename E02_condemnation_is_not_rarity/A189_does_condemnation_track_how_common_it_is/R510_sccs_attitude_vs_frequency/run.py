"""E02·A189·R510 — 谴责是不是「有多少人做」的函数?被试 = 186 个社会。

起因(Ivan, 2026-08-05,原话):
  「听黑死亡金属他就觉得很帅,但他如果喜欢吃屎他就会觉得很羞耻,都是小众为什么不同,
   这不才是那个问题吗?」
  「习近平一个人当皇帝这个很小众,这个小众不行,他为什么不羞耻呢?」

E01 的四百轮把「稀有度 → 羞耻」当机制,而 corr=0.758 只是**在性这个领域内**的一个巧合:
被谴责的东西同时又稀有,因为谴责让它稀有,或者两者同为「社会怎么看」的下游。
E01 拿不到「社会怎么看」这一列 —— 它不在那张表里。

SCCS 里它在。Broude & Greene (1976) 对同一批社会**分别**编码了
「这件事被怎么看待」和「这件事有多常见」—— 两列,四对实践。

ESTIMAND(先于方法):
  在有两列的社会上,**谴责等级**与**罕见程度**的 Spearman ρ。
  两个量表都朝「更罕见 / 更被谴责」递增 ⇒ ρ>0 表示两者同行。

三个世界:
  W-A(E01 的操作子):谴责是稀有度的函数        -> ρ 强正
  W-B(Ivan):谴责与稀有度无关                  -> ρ ≈ 0
  W-C(反向):谴责压低了频率                    -> ρ 强正
⚠ 这个设计把 {A,C} 与 {B} 分开,**分不开 A 和 C**。事先写下,不在结果出来后才说。

⚠ 最强混淆,无法在此控制,必须命名:
  态度与频率是**同一批编码者**从**同一批民族志**里编出来的(Broude & Greene 1976)。
  共享方法方差会**朝正方向**抬高 ρ。⇒ 一个正的 ρ 有一个非社会学的解释。
  ⇒ 因此本轮只有**零结果**是干净的;正结果必须降级为 UNVERIFIED-for-mechanism。

⚠ Galton 问题:社会不独立(共享历史/语言/地理)。
  零分布用**区域内置换**(保留区域聚集),而不是自由置换。

预注册条件式 kill:
  仅当 positive_control 触发 且 negative_control 为零时才评判阈值,否则 UNVERIFIED。
"""
import os, sys, pathlib, json, csv, math
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))

import numpy as np
from lib.gates import Gate, check_columns

RG = np.random.default_rng(20260805)
SCCS = ROOT / "data/external/dplace/repo/datasets/SCCS"
OUT = pathlib.Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------------- 载入
rows = list(csv.DictReader(open(SCCS / "data.csv")))
socs = list(csv.DictReader(open(SCCS / "societies.csv")))
vardef = {v["id"]: v["title"] for v in csv.DictReader(open(SCCS / "variables.csv"))}
codes = {}
for c in csv.DictReader(open(SCCS / "codes.csv")):
    codes.setdefault(c["var_id"], {})[c["code"]] = c["name"]

# soc -> var -> code (取第一条;sub_case 重复时用众数)
tab = {}
for r in rows:
    tab.setdefault(r["soc_id"], {}).setdefault(r["var_id"], []).append(r["code"])


def vec(var, socids):
    out = []
    for s in socids:
        vs = [x for x in tab.get(s, {}).get(var, []) if x not in ("NA", "")]
        if not vs:
            out.append(np.nan); continue
        try:
            out.append(float(max(set(vs), key=vs.count)))
        except ValueError:
            out.append(np.nan)
    return np.array(out, float)


socids = sorted(tab.keys(), key=lambda s: int(s.replace("SCCS", "")))
lat = {s["id"]: float(s["Lat"]) for s in socs if s.get("Lat")}
lon = {s["id"]: float(s["Long"]) for s in socs if s.get("Long")}
LAT = np.array([lat.get(s, np.nan) for s in socids])
LON = np.array([lon.get(s, np.nan) for s in socids])
print(f"societies with any data: {len(socids)}   lat/long known: {np.isfinite(LAT).sum()}")

# 粗区域(用于区域内置换 —— Galton)。经纬度分箱,粗但可复算,且**说明它是粗的**。
REGION = np.full(len(socids), -1)
for i, (a, o) in enumerate(zip(LAT, LON)):
    if not (np.isfinite(a) and np.isfinite(o)):
        continue
    if o < -30:   REGION[i] = 0 if a > 12 else 1          # N / S America
    elif o < 45:  REGION[i] = 2 if a > 20 else 3          # Eurasia-W / Africa
    elif o < 100: REGION[i] = 4                            # S/C Asia
    else:         REGION[i] = 5 if a > 0 else 6            # E Asia / Oceania
print("region sizes:", {int(k): int((REGION == k).sum()) for k in sorted(set(REGION))})


# ---------------------------------------------------------------- 统计
def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8:
        return np.nan, int(m.sum())
    a, b = x[m], y[m]
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = math.sqrt((ra @ ra) * (rb @ rb))
    return (float(ra @ rb / d) if d > 0 else np.nan), int(m.sum())


def perm_within_region(x, y, region, n=4000, rng=RG):
    """区域内置换 y —— 保留区域聚集(Galton),打破区域内的配对。"""
    m = np.isfinite(x) & np.isfinite(y)
    xs, ys, rs = x[m], y[m], region[m]
    out = []
    for _ in range(n):
        yp = ys.copy()
        for g in np.unique(rs):
            idx = np.where(rs == g)[0]
            if len(idx) > 1:
                yp[idx] = ys[idx][rng.permutation(len(idx))]
        r, _ = spearman(xs, yp)
        if np.isfinite(r):
            out.append(r)
    return np.array(out)


def perm_free(x, y, n=4000, rng=RG):
    m = np.isfinite(x) & np.isfinite(y)
    xs, ys = x[m], y[m]
    out = []
    for _ in range(n):
        r, _ = spearman(xs, ys[rng.permutation(len(ys))])
        if np.isfinite(r): out.append(r)
    return np.array(out)


# ---------------------------------------------------------------- 变量族
# (名字, 态度变量, 频率变量, 态度需排除的码, 说明)
FAMILY = [
    ("premarital_F", "SCCS165", "SCCS167", set(),
     "婚前性(女):态度 1 期待→6 强烈不许 × 频率 1 普遍→4 罕见"),
    ("premarital_M", "SCCS165", "SCCS166", set(),
     "⚠ 态度只编了女性一侧 —— 男性频率配女性态度,是一个 MISMATCH,列出但不算入主族"),
    ("extramarital", "SCCS169", "SCCS171", set(),
     "⚠ SCCS169 是类别(单/双重标准),不是谴责序数 —— 列出但不算入主族"),
    ("rape",         "SCCS173", "SCCS174", set(),
     "强奸:态度 × 频率"),
    ("homosexual",   "SCCS176", "SCCS177", {2.0},
     "同性:⚠ 态度 code 2='None' 不是谴责等级而是「不存在」—— 必须剔除"),
]
PRIMARY = {"premarital_F", "rape", "homosexual"}   # 预注册的主族(多重性在这三个上)

print("\n=== 编码检查 ===")
for nm, av, fv, drop, note in FAMILY:
    print(f"{nm:14s} {av}={vardef[av][:34]:36s} {fv}={vardef[fv][:30]}")
    if drop: print(f"               drop codes {drop}: " +
                   ", ".join(f"{c}={codes[av].get(str(int(c)),'?')}" for c in drop))

# ---------------------------------------------------------------- 主计算
res = {}
for nm, av, fv, drop, note in FAMILY:
    a = vec(av, socids); f = vec(fv, socids)
    n_drop = 0
    if drop:
        bad = np.isin(a, list(drop)); n_drop = int(bad.sum()); a[bad] = np.nan
    r, n = spearman(a, f)
    res[nm] = dict(rho=r, n=n, n_dropped=n_drop, att=av, freq=fv, note=note,
                   att_title=vardef[av], freq_title=vardef[fv])
    print(f"\n{nm:14s} rho={r:+.4f}  n={n:3d}" + (f"  (dropped {n_drop} 'None')" if n_drop else ""))
    print(f"               {note}")

# ---------------------------------------------------------------- 零分布(族内最大 |rho|)
print("\n=== 族内 max|rho| 零分布(区域内置换,保留 Galton 聚集)===")
NPERM = 4000
maxnull_reg, maxnull_free = [], []
prim = [x for x in FAMILY if x[0] in PRIMARY]
per_null = {}
for nm, av, fv, drop, note in prim:
    a = vec(av, socids); f = vec(fv, socids)
    if drop: a[np.isin(a, list(drop))] = np.nan
    nr = perm_within_region(a, f, REGION, NPERM)
    nf = perm_free(a, f, NPERM)
    per_null[nm] = dict(reg_sd=float(nr.std()), free_sd=float(nf.std()),
                        reg_q95=float(np.quantile(np.abs(nr), .95)),
                        free_q95=float(np.quantile(np.abs(nf), .95)))
    maxnull_reg.append(np.abs(nr)); maxnull_free.append(np.abs(nf))
    print(f"{nm:14s} region-perm sd={nr.std():.4f} q95|.|={np.quantile(np.abs(nr),.95):.4f}   "
          f"free-perm sd={nf.std():.4f} q95|.|={np.quantile(np.abs(nf),.95):.4f}")

L = min(len(x) for x in maxnull_reg)
MAXNULL = np.max(np.vstack([x[:L] for x in maxnull_reg]), axis=0)
fw_q95 = float(np.quantile(MAXNULL, .95))
print(f"\n族内 max|rho| 的 95% 分位 (family-wise) = {fw_q95:.4f}   (family = {sorted(PRIMARY)})")

for nm in PRIMARY:
    r = res[nm]["rho"]
    p = float((MAXNULL >= abs(r)).mean()) if np.isfinite(r) else np.nan
    res[nm]["p_fw"] = p
    print(f"  {nm:14s} |rho|={abs(r):.4f}  p_fw={p:.4f}  " +
          ("超过族零" if abs(r) > fw_q95 else "在族零之内"))

# ---------------------------------------------------------------- 控制
G = Gate("谴责是不是稀有度的函数?(186 个社会)")

# 阳性对照:同一实践、两个性别的频率 —— 必须强相关,否则仪器/数据死了
pcx = vec("SCCS166", socids); pcy = vec("SCCS167", socids)
pc_r, pc_n = spearman(pcx, pcy)
pc_null = perm_within_region(pcx, pcy, REGION, 2000)
print(f"\n阳性对照 男/女婚前性频率  rho={pc_r:+.4f} n={pc_n}  null sd={pc_null.std():.4f}")
pc_ok = G.positive_control("阳性:同一实践两性别的频率必须同行",
                           planted=abs(pc_r), floor=float(np.quantile(np.abs(pc_null), .95)),
                           spread=float(pc_null.std()))

# 阴性对照:谴责 × 与性无关的序数(社会分层量表)。这个零**应该**是零 ⇒ negative_control
ncx = vec("SCCS165", socids); ncy = vec("SCCS158", socids)   # Scale 10 - Social Stratification
nc_r, nc_n = spearman(ncx, ncy)
nc_null = perm_within_region(ncx, ncy, REGION, 2000)
print(f"阴性对照 婚前性态度 × 社会分层  rho={nc_r:+.4f} n={nc_n}  null sd={nc_null.std():.4f}")
main_r = res["premarital_F"]["rho"]
nc_ok = G.negative_control("阴性:谴责 × 社会分层(与性无关)",
                           null=nc_r, effect=main_r,
                           null_spread=float(nc_null.std()),
                           null_kind="区域内置换(保留 Galton 聚集)")

# 分辨率:主效应能否与它自己的零分开
G.resolvable("主效应可分辨", effect=abs(main_r), spread=float(np.std(perm_within_region(
    vec("SCCS165", socids), vec("SCCS167", socids), REGION, 2000))))

# ---------------------------------------------------------------- 规格曲线
print("\n=== 规格曲线(不是一格)===")
spec = []
for nm in sorted(PRIMARY):
    av, fv = res[nm]["att"], res[nm]["freq"]
    drop = dict((x[0], x[3]) for x in FAMILY)[nm]
    a0 = vec(av, socids); f0 = vec(fv, socids)
    if drop: a0[np.isin(a0, list(drop))] = np.nan
    for spec_name, sel in [
        ("all", np.ones(len(socids), bool)),
        ("lat|<30", np.abs(LAT) < 30),
        ("lat|>=30", np.abs(LAT) >= 30),
        ("oldworld", LON > -30),
        ("newworld", LON <= -30),
    ]:
        r, n = spearman(np.where(sel, a0, np.nan), np.where(sel, f0, np.nan))
        spec.append(dict(pair=nm, spec=spec_name, rho=r, n=n))
        print(f"  {nm:14s} {spec_name:10s} rho={r:+.4f} n={n:3d}" if np.isfinite(r)
              else f"  {nm:14s} {spec_name:10s} rho=  n/a  n={n:3d}")

fin = [s for s in spec if np.isfinite(s["rho"])]
if fin:
    signs = [np.sign(s["rho"]) for s in fin]
    dom = max(set(signs), key=signs.count)
    share = signs.count(dom) / len(signs)
    print(f"\nspec_survival: {signs.count(dom)}/{len(signs)} = {share:.0%} 同号 (sign={dom:+.0f});"
          f" 全格公布,含不同号的格")

# ---------------------------------------------------------------- 条件式 kill
print("\n" + "=" * 66)
KILL_THRESH = None
if pc_ok and nc_ok:
    r = abs(res["premarital_F"]["rho"])
    verdict = ("W-B(与稀有度无关)存活" if r <= fw_q95
               else "W-B 被否;{W-A,W-C} 存活且本设计分不开")
    KILL_THRESH = fw_q95
    print(f"控制齐备 ⇒ 评判阈值。|rho|={r:.4f} vs 族零 q95={fw_q95:.4f}")
    print(f"判定:{verdict}")
else:
    verdict = "UNVERIFIED —— 控制未齐,预注册禁止评判阈值"
    print(f"⚠ {verdict}  (pos={pc_ok} neg={nc_ok})")

print(G)

json.dump(dict(res={k: {kk: (None if isinstance(vv, float) and not np.isfinite(vv) else vv)
                        for kk, vv in v.items()} for k, v in res.items()},
               per_null=per_null, fw_q95=fw_q95, spec=[
                   {**s, "rho": (None if not np.isfinite(s["rho"]) else s["rho"])} for s in spec],
               positive=dict(rho=pc_r, n=pc_n, ok=bool(pc_ok)),
               negative=dict(rho=nc_r, n=nc_n, ok=bool(nc_ok)),
               verdict=verdict, kill_threshold=KILL_THRESH,
               n_societies=len(socids), seed=20260805, nperm=NPERM),
          open(OUT / "sccs_attitude_vs_frequency.json", "w"), indent=1)
print(f"\nwrote {OUT/'sccs_attitude_vs_frequency.json'}")

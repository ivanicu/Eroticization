"""E02·A190·R526 — 那个我宣布无法控制的混淆,可以在非性的地方被控制

FRONTIER。SAMPLE:E02 里最承重且未被攻击的断言 —— `#466d`
「规范能压住的行为,谴责与稀有度绑在一起;压不住的,它们脱钩」。
它是 D4、由数据生成,而 `#466c` 已把产生它的 +0.78 降级为 UNVERIFIED-for-mechanism,
理由是**共享编码者**:态度与频率来自同一批人读同一批民族志,可能来自同一句话。
`#466c` 当时写的是「无法在此控制」。**那句话是错的,而 `source` 那一列证明它错。**

GRADIENT。SCCS 的 10 个性变量全部来自 `broude1976cross`;
`lang1998conan` 里有**三对非性的态度×频率**,而且**范围逐对匹配**:

  范围            态度(1 拒斥 → 3 赞赏)   频率(1 罕见 → 4 永久)   规范可控性
  本地社区内      SCCS1768                 SCCS1750                 最高
  本族群内        SCCS1769                 SCCS1776                 中
  跨族群          SCCS1770                 SCCS1778                 最低

  跨族群那一对的可控性最低,因为**对方的行为不受你的规范管**。

ESTIMAND(先于方法):三个范围各自的 Spearman ρ(谴责, 稀有),
以及**它们随范围外移的斜率**。

⚠ GAUGE(frontier §3,三行零算力):
  SCCS1768 高 = 更被**接受**;SCCS1750 高 = 更**常见**。
  谴责 := 反向,稀有 := 反向。ρ(−a,−b) = ρ(a,b) ⇒ 原始 ρ 已等于 ρ(谴责,稀有)。
  脚本仍显式构造朝向后的变量,好让读者能核,而不是相信这句话。

WORLDS(预测矩阵,粗数,形状才是重点):
  W-CODER    耦合是编码者读一句话发两个码产生的  -> 三个范围**持平**
  W-SUPPRESS `#466d`:耦合跟随规范能否压住行为   -> ρ **随范围外移递减**
  W-MIXED    两者都在                            -> 递减但幅度小于自身展布
  | World      | now | 递减 | 持平 |
  | W-CODER    | 0.4 | 0.1  | 0.8  |
  | W-SUPPRESS | 0.3 | 0.8  | 0.1  |
  | W-MIXED    | 0.3 | 0.4  | 0.3  |
  没有平行,最差分支仍移动 ~0.3。

⚠ STRONGEST CONFOUND,先写下:三对的 n 不同(态度 68/62/41)⇒ ρ 精度不同,
  「递减」可能是精度伪影。控制:每个 ρ 报自己的置换展布,斜率对**展布**判,不对零判。
⚠ 哨兵码:lang 大量使用 `88 = 不适用`。**88 必须排除,不能当等级** ——
  这与 `#466` 里 SCCS176 的 `code 2 = None` 是同一个陷阱。

CONTROLS:
  正对照   SCCS1776 频率 × SCCS1777 强度(同现象两面,同编码者)必须强正
  阴性对照 谴责 × SCCS1753 单系继嗣深度(同源、无实质关联)—— 这个零**应该**是零
           ⇒ negative_control。它同时是 W-CODER 需要的对照:
           **同源本身不产生耦合。**
  SHAM     3×3 全格:对角=范围匹配(主),非对角=范围错配。
           若耦合是范围特异的,非对角应弱于对角。

KILL(条件式,预注册,写在跑之前):
  if 正对照触发 and 阴性对照为零:
      斜率显著为负(超自身展布) -> W-CODER 减弱
      斜率在展布之内            -> W-SUPPRESS 在非性域减弱
      else                      -> MIXED
  else: UNVERIFIED

IMPOSSIBLE(本站点结构上做不到,各附所需):
  causally identified        -> 需要对机制的干预
  independently replicated   -> 需要第二个编码团队;本轮**未派对抗 agent**(会话约束)
                                ⇒ 全部行标 [unchallenged]
  temporally resolved        -> SCCS 每社会只有一个焦点年
  cross-dataset              -> EA 不含这些码
"""
import os, sys, pathlib, json, csv, math, collections
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))

import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
SCCS = ROOT / "data/external/dplace/repo/datasets/SCCS"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)

SENTINEL = {88.0}          # lang 的「不适用」

# ---------------------------------------------------------------- 载入
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


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8: return np.nan, int(m.sum())
    a, b = x[m], y[m]
    ra = np.argsort(np.argsort(a)).astype(float); rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = math.sqrt((ra @ ra) * (rb @ rb))
    return (float(ra @ rb / d) if d > 0 else np.nan), int(m.sum())


def perm_region(x, y, n=4000, seed=0):
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


# ---------------------------------------------------------------- 朝向
SCOPES = [
    ("within_community", "SCCS1768", "SCCS1750", 3, 4),   # 态度 1..3, 频率 1..4
    ("within_ethnic",    "SCCS1769", "SCCS1776", 3, 4),
    ("between_ethnic",   "SCCS1770", "SCCS1778", 3, 4),   # ⚠ 1770 用 10/11/20/21/22 两位码
]
print("=== GAUGE:朝向构造(高 = 更被谴责 / 更罕见)===")
COND, RARE = {}, {}
for nm, av, fv, amax, fmax in SCOPES:
    a_raw, f_raw = vec(av), vec(fv)
    # 1770 用两位码 10/11/20/21/22 -> 取十位作为等级
    if av == "SCCS1770": a_raw = np.floor(a_raw / 10.0)
    COND[nm] = -a_raw          # 原码高=更接受 -> 取负 = 更谴责
    RARE[nm] = -f_raw          # 原码高=更常见 -> 取负 = 更罕见
    ua = sorted(set(a_raw[np.isfinite(a_raw)])); uf = sorted(set(f_raw[np.isfinite(f_raw)]))
    print(f"  {nm:17s} {av} levels={ua}  n={np.isfinite(a_raw).sum():3d} | "
          f"{fv} levels={uf} n={np.isfinite(f_raw).sum():3d}")

# ---------------------------------------------------------------- 3x3 全格(对角=主,非对角=SHAM)
print("\n=== G4 规格曲线:3×3 全格(对角 = 范围匹配 = 主;非对角 = SHAM 范围错配)===")
grid = {}
names = [s[0] for s in SCOPES]
print(f"{'cond \\ rare':20s}" + "".join(f"{n:>19s}" for n in names))
for cn in names:
    line = f"{cn:20s}"
    for fn in names:
        r, n = spearman(COND[cn], RARE[fn])
        grid[(cn, fn)] = (r, n)
        line += f"{r:+.3f} (n={n:3d})".rjust(19)
    print(line)

diag = [grid[(n, n)][0] for n in names]
off = [grid[(a, b)][0] for a in names for b in names if a != b]
print(f"\n对角(范围匹配) mean={np.nanmean(diag):+.4f}   "
      f"非对角(SHAM 错配) mean={np.nanmean(off):+.4f}   "
      f"差 = {np.nanmean(diag)-np.nanmean(off):+.4f}")

# ---------------------------------------------------------------- 主:斜率
print("\n=== 主检验:ρ 随范围外移的斜率 ===")
rows = []
for i, nm in enumerate(names):
    r, n = grid[(nm, nm)]
    nulls = [perm_region(COND[nm], RARE[nm], 3000, s) for s in SEEDS]
    sd = float(np.mean([x.std() for x in nulls]))
    q95 = float(np.mean([np.quantile(np.abs(x), .95) for x in nulls]))
    seed_spread = float(np.std([x.std() for x in nulls]))
    rows.append(dict(scope=nm, order=i, rho=r, n=n, null_sd=sd, null_q95=q95,
                     seed_spread=seed_spread))
    print(f"  [{i}] {nm:17s} rho={r:+.4f} n={n:3d}  null sd={sd:.4f} q95={q95:.4f} "
          f"seed_spread={seed_spread:.5f}")

x = np.array([r["order"] for r in rows], float)
y = np.array([r["rho"] for r in rows], float)
slope = float(np.polyfit(x, y, 1)[0]) if np.all(np.isfinite(y)) else np.nan
# 斜率自身展布:每个 ρ 从自己的零里抽,重算斜率
rng = np.random.default_rng(20260805)
sl_null = []
for _ in range(4000):
    yy = np.array([rng.normal(0, r["null_sd"]) for r in rows])
    sl_null.append(np.polyfit(x, yy, 1)[0])
sl_null = np.array(sl_null)
print(f"\n斜率 = {slope:+.4f}   斜率零展布 sd={sl_null.std():.4f}  "
      f"|slope|/sd = {abs(slope)/max(sl_null.std(),1e-12):.2f}x")

# ---------------------------------------------------------------- 控制
G = Gate("耦合是编码者造的,还是规范造的?(lang1998conan,非性)")

pcx, pcy = vec("SCCS1776"), vec("SCCS1777")     # 族内暴力 频率 × 强度
pc_r, pc_n = spearman(pcx, pcy)
pc_null = perm_region(pcx, pcy, 2000, SEEDS[0])
print(f"\n正对照 族内暴力 频率×强度  rho={pc_r:+.4f} n={pc_n} null q95={np.quantile(np.abs(pc_null),.95):.4f}")
pc_ok = G.positive_control("正对照:同现象两面必须同行",
                           planted=abs(pc_r), floor=float(np.quantile(np.abs(pc_null), .95)),
                           spread=float(pc_null.std()))

ncx, ncy = COND["within_community"], vec("SCCS1753")   # 谴责 × 单系继嗣深度
nc_r, nc_n = spearman(ncx, ncy)
nc_null = perm_region(ncx, ncy, 2000, SEEDS[0])
main_r = grid[("within_community", "within_community")][0]
print(f"阴性对照 谴责 × 单系继嗣深度  rho={nc_r:+.4f} n={nc_n} null sd={nc_null.std():.4f}")
nc_ok = G.negative_control("阴性:同源但无实质关联 —— 同源本身不产生耦合",
                           null=nc_r, effect=main_r, null_spread=float(nc_null.std()),
                           null_kind="区域内置换(保留 Galton 聚集)")

G.resolvable("主效应(社区内)可分辨", effect=abs(main_r), spread=rows[0]["null_sd"])
G.has_error_bar("斜率", value=slope, spread=float(sl_null.std()),
                spread_source="每个 ρ 从自身区域内置换零重抽,重算斜率,4000 次")

# ---------------------------------------------------------------- 条件式 KILL
print("\n" + "=" * 68)
if pc_ok and nc_ok:
    if np.isfinite(slope) and abs(slope) > 2 * sl_null.std() and slope < 0:
        verdict = "斜率显著为负 -> W-CODER 减弱;耦合跟随范围而非编码者"
    elif np.isfinite(slope) and abs(slope) <= 2 * sl_null.std():
        verdict = "斜率在自身展布之内 -> 三范围持平 -> W-SUPPRESS 在非性域减弱"
    else:
        verdict = "斜率显著为正 -> 两个世界都没预测这个方向;MIXED/未知"
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:若三个范围的可控性排序本身错了"
          "(即跨族群其实同样受规范支配),斜率的解释就反了 —— 排序是先验的,未被独立评定。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pc_ok} neg={nc_ok}),预注册禁止评判"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(
    scopes=[{k: (None if isinstance(v, float) and not np.isfinite(v) else v)
             for k, v in r.items()} for r in rows],
    grid={f"{a}|{b}": (None if not np.isfinite(grid[(a,b)][0]) else grid[(a,b)][0])
          for a in names for b in names},
    grid_n={f"{a}|{b}": grid[(a,b)][1] for a in names for b in names},
    diag_mean=float(np.nanmean(diag)), off_mean=float(np.nanmean(off)),
    slope=slope, slope_null_sd=float(sl_null.std()),
    positive=dict(rho=pc_r, n=pc_n, ok=bool(pc_ok)),
    negative=dict(rho=nc_r, n=nc_n, ok=bool(nc_ok)),
    verdict=verdict, seeds=SEEDS, source="lang1998conan", sentinel_excluded=sorted(SENTINEL),
    unchallenged=True),
    open(OUT / "violence_scope_gradient.json", "w"), indent=1)
print(f"\nwrote {OUT/'violence_scope_gradient.json'}")

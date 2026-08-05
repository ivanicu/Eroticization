"""E02·A191·R528 — +0.7966 在它自己那批编码者手里,是不是特别的

`#482` 的 NEXT,逐字执行:参照分布建在 `lang` 里,而头条在 `broude1976cross` 里,
**跨源外推未被授权** -> 在 broude 内部建同一个分布。

G1 ESTIMAND(先于方法):`+0.7966`(婚前性谴责 × 稀有)在
**`broude1976cross` 内部所有变量对的 |ρ| 分布**中的分位。
IDENTIFICATION:18 个合格变量(3-12 级,n>=30)-> 18*17/2 = **153 对**,过 `#482` 定的 15 对门槛。

WORLDS:
  W-CODER-BROUDE  这批编码者本来就让自己的变量互相强耦合 -> +0.7966 落在分布**体内**
  W-SPECIAL       耦合是「态度×频率」这个配对特有的        -> 落在**上尾 (q95+)**
  W-DOMAIN        broude 的整套变量都关于性规范,所以什么都相关 -> 也落在体内
  ⚠ **本设计把 W-SPECIAL 与 {W-CODER, W-DOMAIN} 分开,分不开后两者。** 先写下。
  | World          | now | q95+ | 体内 |
  | W-CODER-BROUDE | 0.4 | 0.10 | 0.80 |
  | W-SPECIAL      | 0.4 | 0.85 | 0.10 |
  | W-DOMAIN       | 0.2 | 0.05 | 0.10 |

⚠ STRONGEST CONFOUND,写在跑之前(两个):
  ① 参照分布若**锚在 SCCS165 一个变量上**,量到的是那个变量的边际,不是编码者。
     -> 用**全部 153 对**,不是「谴责 × 其他」。(`#482c` 在 lang 里就是锚着的,这里修掉。)
  ② 153 对里含**同一实践的近重复对**(166×167 男女婚前频率 · 170×171 男女婚外频率),
     它们必然强相关,会抬高上尾 -> **含/不含两版都报**。
  ③ 参照分布还含**被测的那 4 对态度×频率本身** -> 必须剔除,否则用被测对象当基线。

⚠ 预注册规格网格(写在跑之前,`#482` 的 NEXT 要求):
  ordinal filter k ∈ {[3,8], [3,12], [2,20]} × min_n ∈ {20, 30, 40} = 9 格,**全格公布**。

CONTROLS:
  正对照  SCCS166 × SCCS167(同一实践,两性别)—— 这一次算 floor **和** ceiling(`#482d`)
  阴性    参照分布的中位 —— **测量出来的,不是挑出来的**(`#482c`)
  ⚠ 正对照本身也在近重复对里 -> 它同时是「近重复必然强相关」的示范,这是有意的。

KILL(条件式,预注册):
  if 正对照触发 and 阴性为零:
      +0.7966 在**九格全部**达 q95 -> W-CODER-BROUDE 减弱,`#466c` 的降级可撤
      任一格未达            -> 降级保留
  else: UNVERIFIED

预注册探索项(**不评判 kill**,只列出):`broude1983cross` 的
  SCCS743 (Neg) Attitude towards Divorce × SCCS744 (Neg) Frequency of Divorce
  —— 第五对态度×频率,同作者不同出版物,离婚是社会**直接规管**的行为。
  ⚠ "(Neg)" 是 D-PLACE 的反向编码标记 -> **先读码,不假设方向**。

IMPOSSIBLE(附所需):
  causally identified      -> 需要干预
  independently replicated -> 需要第二个编码团队;**本轮未派对抗 agent(会话约束)** ⇒ [unchallenged]
  W-CODER vs W-DOMAIN      -> 需要一个**非性主题**但同编码者的变量集;broude 没有
"""
import os, sys, pathlib, json, csv, math, collections, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))

import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
SCCS = ROOT / "data/external/dplace/repo/datasets/SCCS"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SENTINEL = {88.0}
HEADLINE = ("SCCS165", "SCCS167")          # 婚前性 谴责 × 频率
AF_PAIRS = {("SCCS165", "SCCS167"), ("SCCS165", "SCCS166"),
            ("SCCS169", "SCCS170"), ("SCCS169", "SCCS171"),
            ("SCCS173", "SCCS174"), ("SCCS176", "SCCS177")}
NEAR_DUP = {("SCCS166", "SCCS167"), ("SCCS170", "SCCS171"), ("SCCS163", "SCCS164")}

V = {v["id"]: v for v in csv.DictReader(open(SCCS / "variables.csv"))}
codes = collections.defaultdict(dict)
for c in csv.DictReader(open(SCCS / "codes.csv")): codes[c["var_id"]][c["code"]] = c["name"]
tab = collections.defaultdict(lambda: collections.defaultdict(list))
for r in csv.DictReader(open(SCCS / "data.csv")):
    if r["code"] not in ("NA", ""): tab[r["soc_id"]][r["var_id"]].append(r["code"])
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


def vec(var, drop=()):
    out = []
    for s in socids:
        vs = tab[s].get(var, [])
        if not vs: out.append(np.nan); continue
        try: x = float(max(set(vs), key=vs.count))
        except ValueError: out.append(np.nan); continue
        out.append(np.nan if (x in SENTINEL or x in drop) else x)
    return np.array(out, float)


def midrank(a):                                   # `#482a`
    order = np.argsort(a, kind="mergesort"); r = np.empty(len(a)); s = a[order]; i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]: j += 1
        r[order[i:j + 1]] = (i + j) / 2.0; i = j + 1
    return r


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8: return np.nan, int(m.sum())
    ra, rb = midrank(x[m]), midrank(y[m])
    ra = ra - ra.mean(); rb = rb - rb.mean()
    d = math.sqrt((ra @ ra) * (rb @ rb))
    return (float(ra @ rb / d) if d > 0 else np.nan), int(m.sum())


def ceiling(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8: return np.nan
    r, _ = spearman(np.sort(x[m]), np.sort(y[m])); return r


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


# ---------------------------------------------------------------- 头条 + 探索项
h, hn = spearman(vec(HEADLINE[0]), vec(HEADLINE[1]))
print(f"头条 {HEADLINE[0]}×{HEADLINE[1]} rho={h:+.4f} n={hn}\n")

print("=== 预注册探索项:离婚(broude1983cross)—— 先读码,不假设方向 ===")
for v in ("SCCS743", "SCCS744"):
    print(f"  {v} {V[v]['title'][:46]}")
    for c, nm in sorted(codes[v].items(), key=lambda kv: (kv[0] != 'NA', kv[0])):
        if c != "NA": print(f"      {c} = {nm[:62]}")
dv_r, dv_n = spearman(vec("SCCS743"), vec("SCCS744"))
dv_null = perm_region(vec("SCCS743"), vec("SCCS744"), 3000, SEEDS[0])
print(f"  rho(态度, 频率) = {dv_r:+.4f}  n={dv_n}  零 q95|.|={np.quantile(np.abs(dv_null),.95):.4f}")

# ---------------------------------------------------------------- 规格网格
print("\n=== G4 规格曲线:9 格全公布 ===")
BR = [v for v, d in V.items() if d["source"] == "broude1976cross"]
grid_rows = []
for (klo, khi) in [(3, 8), (3, 12), (2, 20)]:
    for minn in (20, 30, 40):
        ok = []
        for v in BR:
            y = vec(v); u = np.unique(y[np.isfinite(y)])
            if klo <= len(u) <= khi and np.isfinite(y).sum() >= minn: ok.append(v)
        refs, refs_nodup = [], []
        for a, b in itertools.combinations(sorted(ok), 2):
            if (a, b) in AF_PAIRS or (b, a) in AF_PAIRS: continue      # 剔除被测对象
            r, n = spearman(vec(a), vec(b))
            if not np.isfinite(r) or n < minn: continue
            refs.append((abs(r), a, b))
            if (a, b) not in NEAR_DUP and (b, a) not in NEAR_DUP: refs_nodup.append((abs(r), a, b))
        if len(refs) < 15:
            print(f"  k[{klo},{khi}] n>={minn}: 只有 {len(refs)} 对 -> 不可识别,跳过"); continue
        A = np.array([x[0] for x in refs]); Ad = np.array([x[0] for x in refs_nodup])
        pct = float((A < abs(h)).mean()); pctd = float((Ad < abs(h)).mean())
        grid_rows.append(dict(k=[klo, khi], min_n=minn, n_vars=len(ok), n_pairs=len(refs),
                              median=float(np.median(A)), q95=float(np.quantile(A, .95)),
                              maxv=float(A.max()), pct=pct, n_pairs_nodup=len(refs_nodup),
                              pct_nodup=pctd, q95_nodup=float(np.quantile(Ad, .95))))
        print(f"  k[{klo:2d},{khi:2d}] n>={minn}: vars={len(ok):2d} pairs={len(refs):3d} "
              f"med={np.median(A):.3f} q95={np.quantile(A,.95):.3f} max={A.max():.3f} "
              f"-> 头条在 q{100*pct:.0f}   (去近重复 pairs={len(refs_nodup):3d} q{100*pctd:.0f})")

# 主格 = k[3,12], n>=30(`#482` NEXT 里预注册的那个)
main = [g for g in grid_rows if g["k"] == [3, 12] and g["min_n"] == 30][0]
print(f"\n主格 k[3,12] n>=30:参照 {main['n_pairs']} 对,中位 {main['median']:.4f},"
      f"q95 {main['q95']:.4f},max {main['maxv']:.4f}")
print(f"  ⇒ 头条 |{h:+.4f}| 在 q{100*main['pct']:.0f}"
      f"(去近重复 q{100*main['pct_nodup']:.0f})")

# 全分布最强的几对 —— 看上尾装的是什么(meta-separator:问题问对了吗)
ok = [v for v in BR if 3 <= len(np.unique(vec(v)[np.isfinite(vec(v))])) <= 12
      and np.isfinite(vec(v)).sum() >= 30]
allp = []
for a, b in itertools.combinations(sorted(ok), 2):
    if (a, b) in AF_PAIRS or (b, a) in AF_PAIRS: continue
    r, n = spearman(vec(a), vec(b))
    if np.isfinite(r) and n >= 30: allp.append((abs(r), r, a, b, n))
allp.sort(reverse=True)
print("\n  上尾装的是什么(前 8 对):")
for ab, r, a, b, n in allp[:8]:
    print(f"    |{r:+.3f}| n={n:3d}  {a}×{b}  {V[a]['title'][:30]} × {V[b]['title'][:30]}")

# ---------------------------------------------------------------- 控制
G = Gate("+0.7966 在它自己那批编码者手里特别吗?")
pcx, pcy = vec("SCCS166"), vec("SCCS167")
pc_r, pc_n = spearman(pcx, pcy)
pc_null = perm_region(pcx, pcy, 3000, SEEDS[0])
pc_floor = float(np.quantile(np.abs(pc_null), .95)); pc_ceil = ceiling(pcx, pcy)
t = (pc_floor + pc_ceil) / 2
print(f"\n正对照 男/女婚前性频率 rho={pc_r:+.4f} n={pc_n}  floor={pc_floor:.4f} "
      f"ceiling={pc_ceil:.4f} 门槛(带内中点)={t:.4f}")
pc_ok = G.positive_control("正对照:同实践两性别(门槛在 floor<t<ceiling 带内)",
                           planted=abs(pc_r), floor=t, spread=1e-9)
nc_ok = G.negative_control("阴性:参照分布中位(测量,非挑选)",
                           null=main["median"], effect=h,
                           null_spread=float(np.std([g["median"] for g in grid_rows])),
                           null_kind="broude1976cross 内部全变量对参照分布")
G.resolvable("头条可分辨", effect=abs(h),
             spread=float(np.mean([perm_region(vec(*HEADLINE[:1]), vec(HEADLINE[1]), 1500, s).std()
                                   for s in SEEDS])))
G.asserted("参照分布够大", main["n_pairs"] >= 15, f"{main['n_pairs']} 对 (要求>=15)", kind="control")

print("\n" + "=" * 70)
all_q95 = all(g["pct"] >= 0.95 for g in grid_rows) and all(g["pct_nodup"] >= 0.95 for g in grid_rows)
if pc_ok and nc_ok:
    verdict = ("九格全部达 q95 -> W-CODER-BROUDE 减弱,`#466c` 的降级可撤"
               if all_q95 else
               f"并非九格全达 q95(最低 q{100*min(min(g['pct'],g['pct_nodup']) for g in grid_rows):.0f})-> 降级保留")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:参照分布仍无法把 W-CODER 与 W-DOMAIN 分开 ——"
          " broude 的整套变量都关于性规范,所以「同源强耦合」与「同主题强耦合」在此同义。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pc_ok} neg={nc_ok})"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(headline=dict(rho=h, n=hn), grid=grid_rows, main=main,
               divorce=dict(rho=dv_r, n=dv_n, null_q95=float(np.quantile(np.abs(dv_null), .95)),
                            note="预注册探索项,不评判 kill"),
               top_pairs=[dict(rho=r, a=a, b=b, n=n) for _, r, a, b, n in allp[:12]],
               positive=dict(rho=pc_r, n=pc_n, floor=pc_floor, ceiling=pc_ceil,
                             threshold=t, ok=bool(pc_ok)),
               verdict=verdict, all_cells_q95=bool(all_q95), seeds=SEEDS, unchallenged=True),
          open(OUT / "broude_reference_distribution.json", "w"), indent=1)
print(f"\nwrote {OUT/'broude_reference_distribution.json'}")

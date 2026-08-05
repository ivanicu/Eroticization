"""E02·A196·R539 — 我把六条符号矛盾的斜率平均掉了,那个零还算不算零

`#493` 的 NEXT,并且它**指向我不愿意看到的结果**:
`#493b` 我刚把「门槛无关(−0.18,MDE 0.28)」发上了页面,而 `#493c②` 自己写下
分话题斜率**符号相反**(`okcohab` 降、`staytog` 升)。
**若六条斜率彼此可区分,那个合并零就是「平均掉了真实异质性」的产物**,页面那句必须改写。

G1 ESTIMAND(先于方法):
  ① 每话题一条斜率 `s_i = d|lnOR|/d(谴责占比)`,及其人层 bootstrap `se_i`;
     **每条的 MDE = 2.8·se_i,先算再看离散度**(`#485c` 的教训)。
  ② 异质性 `Q = Σ (s_i − s̄_w)² / se_i²`,`s̄_w` = 逆方差加权均值。
     **零不是 χ² 表,而是参数自助**:在「六个话题共用一条斜率 `s̄_w`」下,
     从 `N(s̄_w, se_i²)` 重抽六条,重算 Q,取 4000 次 -> 经验零。

⚠ 结构上限,写在跑之前:AGDGFMT 只有 4 个可用等级 ⇒ 谴责集只能取 {4}/{3,4}/{2,3,4}
  ({1,2,3,4} 会让所有人都是谴责者,退化)⇒ **每话题最多 3 个点,一条 3 点斜率。**
  这不是我偷懒,是量表的结构。`k_points = 3` 是本站点的天花板。

WORLDS:
  W-HOMO   六条是同一条斜率 -> Q 落在零内 -> `#493b` 的合并零成立
  W-HETERO 六条彼此可区分   -> Q 超出零   -> **`#493b` 降级为「平均意义上无关」,页面改写**
  | World     | now | Q 超零 | Q 在零内 |
  | W-HOMO    | 0.4 | 0.10   | 0.85     |
  | W-HETERO  | 0.6 | 0.85   | 0.10     |

⚠ STRONGEST CONFOUND,写在跑之前:三个切点的**行为基率不变、谴责占比变**,
  但每个话题的**占比跨度不同**(`chsuppor` 0.052→0.747,`sxok16` 0.501→0.985)——
  跨度小的话题斜率噪声大 ⇒ `se_i` 差异大 ⇒ **Q 会被少数几个高精度话题主导**。
  控制:同时报**未加权**的斜率方差,以及**逐条 leave-one-topic-out** 的 Q。

CONTROLS(RULE-v3):正对照 = 每题谴责 × 宗教出席;阴性 = 参照分布中位;精度 = 人层 bootstrap。
KILL(条件式,预注册):
  if 正对照全触发 and 阴性为零:
      Q > 零的 q95 -> **W-HETERO,`#493b` 降级并改写页面**
      Q <= q95 且 每条 MDE < 1.40 -> W-HOMO(有内容的零)
      Q <= q95 但 MDE >= 1.40 -> **UNVERIFIED-by-power:直说看不见,不写成「一样」**
  else: UNVERIFIED
IMPOSSIBLE:每话题只有 3 点(量表结构)· 占比与谴责强度不可拆(`#493c①`)·
  无干预 ⇒ 非因果 · 未派对抗 agent(会话约束)⇒ `[unchallenged]`
"""
import os, sys, pathlib, json, math, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
MEANINGFUL = 1.40
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)


def parse_dct(p):
    out = {}
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out


LAY = parse_dct(NS / "setup" / "2011_2013_FemRespSetup.dct")
NEED = ["staytog", "sxok18", "sxok16", "samesex", "okcohab", "chsuppor", "prevhusb",
        "evrmarry", "samesexany", "nonmarr", "cebow", "agefstsx", "attndnow",
        "age_r", "educat", "poverty", "hisp", "religion"]
cols = {n: LAY[n] for n in NEED if n in LAY}
rows = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        rows[n].append(float(v) if v not in ("", ".") else np.nan)
D = {n: np.array(rows[n]) for n in rows}
N = len(D["staytog"])
early = np.where(np.isfinite(D["agefstsx"]) & (D["agefstsx"] >= 5) & (D["agefstsx"] <= 60),
                 (D["agefstsx"] <= 16).astype(float), np.nan)
evermar = D["evrmarry"] == 1
BEH = {
    "samesex": (np.where(D["samesexany"] == 1, 1.0, np.where(D["samesexany"] == 5, 0.0, np.nan)), None),
    "sxok16": (early, None), "sxok18": (early, None),
    "staytog": (np.where(np.isfinite(D["prevhusb"]) & (D["prevhusb"] < 90),
                         (D["prevhusb"] > 0).astype(float), np.nan), evermar),
    "okcohab": (np.where(np.isfinite(D["nonmarr"]) & (D["nonmarr"] < 90),
                         (D["nonmarr"] > 0).astype(float), np.nan), None),
    "chsuppor": (np.where(np.isfinite(D["cebow"]) & (D["cebow"] < 90),
                          (D["cebow"] > 0).astype(float), np.nan), None),
}
REV = {"okcohab": True}
CUTS = [[4], [3, 4], [2, 3, 4]]
print(f"行 {N};话题 {len(BEH)};每话题切点 {len(CUTS)}(量表结构上限)")


def cond_at(v, codes, reverse=False):
    ok = np.isin(v, [1, 2, 3, 4])
    inv = [5 - c for c in codes] if reverse else codes
    return np.where(ok, np.isin(v, inv).astype(float), np.nan)


def lnor(c, b, mask=None):
    m = np.isfinite(c) & np.isfinite(b)
    if mask is not None: m &= mask
    cc, bb = c[m], b[m]
    a1, a0 = bb[cc == 1], bb[cc == 0]
    if len(a1) < 30 or len(a0) < 30: return np.nan, np.nan
    p1, p0 = a1.mean(), a0.mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, np.nan
    return abs(math.log((p1 / (1 - p1)) / (p0 / (1 - p0)))), float(cc.mean())


def topic_slope(topic, idx=None):
    b, mk = BEH[topic]
    v = D[topic] if idx is None else D[topic][idx]
    bb = b if idx is None else b[idx]
    mm = None if mk is None else (mk if idx is None else mk[idx])
    xs, ys = [], []
    for codes in CUTS:
        a, s = lnor(cond_at(v, codes, REV.get(topic, False)), bb, mm)
        if np.isfinite(a): xs.append(s); ys.append(a)
    if len(xs) < 3 or np.ptp(xs) < 1e-9: return np.nan, np.nan
    return float(np.polyfit(xs, ys, 1)[0]), float(np.ptp(xs))


print("\n=== 每话题一条斜率(3 点)+ 自身 MDE ===")
TOP = list(BEH)
slopes, ses, spans = {}, {}, {}
for t in TOP:
    s, span = topic_slope(t)
    bs = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(400):
            v, _ = topic_slope(t, rng.integers(0, N, N))
            if np.isfinite(v): bs.append(v)
    bs = np.array(bs)
    slopes[t], ses[t], spans[t] = s, float(bs.std()), span
    print(f"  {t:9s} slope={s:+.4f}  se={bs.std():.4f}  **MDE={2.8*bs.std():.4f}**  "
          f"占比跨度={span:.3f}  {'⚠ MDE≥有意义量' if 2.8*bs.std() >= MEANINGFUL else ''}")

sl = np.array([slopes[t] for t in TOP]); se = np.array([ses[t] for t in TOP])
w = 1 / se ** 2
sbar = float((w @ sl) / w.sum())
Q = float(((sl - sbar) ** 2 / se ** 2).sum())
print(f"\n逆方差加权均值 s̄_w = {sbar:+.4f};**Q = {Q:.3f}**(df=5)")

rng = np.random.default_rng(SEEDS[0])
Qnull = np.array([(((rng.normal(sbar, se) - sbar) ** 2 / se ** 2).sum()) for _ in range(4000)])
q95 = float(np.quantile(Qnull, .95))
pQ = float((Qnull >= Q).mean())
print(f"  参数自助零(六话题共用一条斜率):q95 = {q95:.3f}  **p = {pQ:.4f}**")

# 混淆控制:未加权方差 + leave-one-topic-out
unw = float(np.var(sl, ddof=1))
print(f"\n混淆控制:未加权斜率方差 = {unw:.4f};leave-one-topic-out 的 Q:")
loo = {}
for i, t in enumerate(TOP):
    m = np.ones(len(TOP), bool); m[i] = False
    ww = 1 / se[m] ** 2; sb = (ww @ sl[m]) / ww.sum()
    loo[t] = float(((sl[m] - sb) ** 2 / se[m] ** 2).sum())
    print(f"  去掉 {t:9s} -> Q={loo[t]:.3f}")

G = Gate("我把六条符号矛盾的斜率平均掉了,那个零还算不算零?")
ref = []
bc = cond_at(D["samesex"], [3, 4])
for cn in ["age_r", "educat", "poverty", "hisp", "religion"]:
    m = np.isfinite(bc) & np.isfinite(D[cn]) & (D[cn] < 90)
    if m.sum() > 500 and len(np.unique(D[cn][m])) >= 3:
        ref.append(float(abs(np.corrcoef(bc[m], D[cn][m])[0, 1])))
REF_MED = float(np.median(ref))
pos = []
for t in TOP:
    c = cond_at(D[t], [3, 4], REV.get(t, False))
    m = np.isfinite(c) & np.isfinite(D["attndnow"]) & (D["attndnow"] < 90)
    r = abs(float(np.corrcoef(c[m], D["attndnow"][m])[0, 1]))
    rg = np.random.default_rng(SEEDS[0])
    q = float(np.quantile([abs(np.corrcoef(c[m][rg.permutation(m.sum())], D["attndnow"][m])[0, 1])
                           for _ in range(300)], .95))
    pos.append(G.positive_control(f"正对照-v3[{t}]", planted=r, floor=max(q, REF_MED), spread=1e-9))
nc = G.negative_control("阴性:参照分布中位(测量,非挑选)", null=REF_MED,
                        effect=float(np.mean(np.abs(sl))), null_spread=float(np.std(ref)),
                        null_kind="同问卷无关变量参照分布")
G.has_error_bar("Q", value=Q, spread=float(Qnull.std()), spread_source="null_零臂")

max_mde = float(2.8 * se.max())
print("\n" + "=" * 70)
if all(pos) and nc:
    if Q > q95:
        verdict = (f"Q={Q:.2f} > 零 q95={q95:.2f}(p={pQ:.4f}) -> "
                   f"**W-HETERO:`#493b` 的合并零降级为「平均意义上无关」,页面须改写**")
    elif max_mde < MEANINGFUL:
        verdict = f"Q={Q:.2f} <= q95={q95:.2f} 且每条 MDE 最大 {max_mde:.2f} < {MEANINGFUL} -> W-HOMO(有内容的零)"
    else:
        verdict = (f"Q={Q:.2f} <= q95={q95:.2f} **但最大 MDE {max_mde:.2f} >= {MEANINGFUL}** -> "
                   f"**UNVERIFIED-by-power:看不见,不是一样**")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:各话题的**占比跨度**差很多 -> `se_i` 差异大 -> "
          "Q 会被少数高精度话题主导(已用未加权方差与 leave-one-out 部分吸收);"
          "而每话题只有 3 个点,是量表的结构上限,不是设计选择。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos} neg={nc})"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(slopes={t: slopes[t] for t in TOP}, ses={t: ses[t] for t in TOP},
               spans={t: spans[t] for t in TOP}, mdes={t: 2.8 * ses[t] for t in TOP},
               sbar_w=sbar, Q=Q, Q_null_q95=q95, p=pQ, unweighted_var=unw, loo=loo,
               max_mde=max_mde, meaningful=MEANINGFUL, verdict=verdict,
               seeds=SEEDS, unchallenged=True),
          open(OUT / "is_the_pooled_null_an_average.json", "w"), indent=1)
print(f"\nwrote {OUT/'is_the_pooled_null_an_average.json'}")

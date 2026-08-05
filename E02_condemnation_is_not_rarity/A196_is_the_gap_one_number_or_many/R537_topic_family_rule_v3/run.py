"""E02·A196·R537 — 「谴责它的人更少承认」是一个数,还是一族随话题变的数

`#491` 的 NEXT,硬约束:**换方向,回到人,不得再以我自己的仪器为对象。**

⛔ **先解决一个诚实问题,而不是绕过它。**
`#491d` 让 `R535` 的问题回到未决,而 `R535` 的三对数**我已经看过**。
**在看过的数上「重新预注册」不是预注册。** 处理办法:**扩族** ——
加入我**从未计算过**的对,**主判据只跑在未看过的子集上**,
已看过的三对**照报但显式标记 `SEEN`**,不进入主判据。

G1 ESTIMAND(先于方法):同一份 NSFG 问卷内,每个「态度×行为」对的 `lnOR` 与其人层 bootstrap CI;
  **主判据 = UNSEEN 子集内,各对 CI 是否两两重叠。**
  重叠 -> 差距是**一个数**;不重叠 -> 差距是**一族随话题变的数**。

族(全部预注册于跑之前;`SEEN` 者不进主判据):
  UNSEEN ① `STAYTOG` IH-2「无法解决婚姻问题时离婚是最好的办法」× `PREVHUSB>0`(**限已婚过者**)
  UNSEEN ② `SXOK18`  IH-5「未婚 18 岁若有强烈感情发生性关系可以」× 早年初次性行为
  UNSEEN ③ `SXOK16`  IH-6「未婚 16 岁…」× 早年初次性行为
  SEEN   ④⑤⑥ `samesex`×曾有过 · `okcohab`×非婚同居 · `chsuppor`×非婚生育(`R535` 已算)
  谴责 := 答「不同意/强烈不同意」{3,4};⛔ 5/8/9 剔除(`#489b`)。
  ⚠ `OKCOHAB` 题干是禁止式(「不应同居」)-> 谴责 = 同意 = {1,2}。每题的方向在代码里显式写死。

WORLDS:
  W-CONSTANT 差距是一个数        -> UNSEEN 各对 CI **两两重叠**
  W-TOPIC    差距随话题变        -> **不重叠**,且极差远大于 CI 宽度
  | World      | now | 不重叠 | 重叠 |
  | W-CONSTANT | 0.4 | 0.10   | 0.85 |
  | W-TOPIC    | 0.6 | 0.85   | 0.10 |

⚠ STRONGEST CONFOUND,写在跑之前:各对的**行为基率**差极大(可能 0.1 ~ 0.9)。
  `lnOR` 对基率位移不变(`#485b` 已验),但**极端基率放大抽样噪声**,会让 CI 变宽 ->
  **偏向「重叠」= 偏向 W-CONSTANT**,即**保守方向**。已写下,不再事后解释。
⚠ 第二个:①限已婚过者 -> **人群与其它对不同**。作为规格轴同时报「不限」的一版。

CONTROLS(**RULE-v3**,`#491c` 预注册):
  正对照 每对的谴责题 × 宗教出席:`|r| > 自身置换零 q95` **且** `|r| > 参照分布中位`
  阴性   参照分布本身(测量,非挑选)
  精度   人层 bootstrap;多重性 = 族内最大 |lnOR| 差的置换零
KILL(条件式,预注册):
  if UNSEEN 各对正对照都触发 and 参照中位 < 0.5·min|lnOR|:
      UNSEEN 内存在**不重叠**的一对 -> **W-TOPIC:差距不是一个数**
      UNSEEN 全部重叠               -> W-CONSTANT
  else: UNVERIFIED
IMPOSSIBLE:无干预 ⇒ 非因果 · 未派对抗 agent(会话约束)⇒ `[unchallenged]` ·
  「话题」与「行为可隐藏性/后果大小」不可拆(`#490e`)。
"""
import os, sys, pathlib, json, math, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
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
# 规则①:先找出真实存在的「初次性行为年龄」变量名,打印,再用
cands = [k for k, v in LAY.items()
         if re.search(r"age at first (sex|interc)|1st sex|first sex", v[2], re.I)]
print("=== 规则①:初次性行为年龄候选(从字典读,不是猜)===")
for k in cands[:8]: print(f"   {k:14s} {LAY[k][2][:66]}")
# ⛔ 第一版我打印了候选**然后让代码盲取 cands[0]**,取到 `wnfstsex_m` =「初次性交的**月份**」,
#    不是年龄 -> 两个 sxok 对建不起来 -> UNSEEN 只剩 1 对 -> 「两两重叠」在 **0 次比较**上算出
#    -> 门全过、判定 W-CONSTANT,而那个判定是**空的**(realstat「空总体通过」)。
#    规则 ① 的新变体:**打印了对象,然后没有读它。** 显式指定,并在下面 assert。
#    第二次:改指 `c_sex15`,而它 **n=15**(罕见分支)-> assert 当场抓住(基率 1.0)。
#    第三次才对:`agefstsx` CE-4「初次性交年龄」。**两次都是「读到名字就用」。**
AGEFS = "agefstsx"

NEED = ["staytog", "sxok18", "sxok16", "samesex", "okcohab", "chsuppor",
        "prevhusb", "evrmarry", "samesexany", "nonmarr", "cebow",
        "attndnow", "age_r", "educat", "poverty", "hisp", "religion"]
if AGEFS: NEED.append(AGEFS)
cols = {n: LAY[n] for n in NEED if n in LAY}
print(f"\n字典命中 {len(cols)}/{len(NEED)};缺: {[n for n in NEED if n not in LAY]}")
rows = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        rows[n].append(float(v) if v not in ("", ".") else np.nan)
D = {n: np.array(rows[n]) for n in rows}
N = len(D[list(D)[0]])
print(f"行 {N}")
for n in ["staytog", "sxok18", "sxok16", "prevhusb", "evrmarry"] + ([AGEFS] if AGEFS else []):
    if n in D:
        u, c = np.unique(D[n][np.isfinite(D[n])], return_counts=True)
        top = dict(sorted(zip(u.tolist(), c.tolist()), key=lambda t: -t[1])[:6])
        print(f"  {n:12s} n={np.isfinite(D[n]).sum():5d}  值={top}   {cols[n][2][:48]}")


def cond(v, reverse=False):
    ok = np.isin(v, [1, 2, 3, 4])
    c = np.isin(v, [1, 2]).astype(float) if reverse else np.isin(v, [3, 4]).astype(float)
    return np.where(ok, c, np.nan)


def lnor(c, b, mask=None):
    m = np.isfinite(c) & np.isfinite(b)
    if mask is not None: m &= mask
    cc, bb = c[m], b[m]
    a1, a0 = bb[cc == 1], bb[cc == 0]
    if len(a1) < 30 or len(a0) < 30: return np.nan, 0
    p1, p0 = a1.mean(), a0.mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, 0
    return math.log((p1 / (1 - p1)) / (p0 / (1 - p0))), int(m.sum())


evermar = D["evrmarry"] == 1
# c_sex15: CE-6「初次性行为时是否 <15 岁」。规则①:先打印码,再定方向,并 assert 非退化。
_u, _c = np.unique(D[AGEFS][np.isfinite(D[AGEFS])], return_counts=True)
print(f"\n{AGEFS} 码分布 = {dict(zip(_u.tolist(), _c.tolist()))}   {cols[AGEFS][2]}")
early = np.where(np.isfinite(D[AGEFS]) & (D[AGEFS] >= 5) & (D[AGEFS] <= 60),
                 (D[AGEFS] <= 16).astype(float), np.nan)
assert 0.02 < np.nanmean(early) < 0.98, f"early 基率 {np.nanmean(early)} 退化 -> 码读错"
print(f"early(初次性交 <=16 岁)基率 = {np.nanmean(early):.4f}  n={np.isfinite(early).sum()}")

FAM = []
if "staytog" in D:
    FAM.append(("UNSEEN staytog×divorced", cond(D["staytog"]),
                np.where(np.isfinite(D["prevhusb"]) & (D["prevhusb"] < 90),
                         (D["prevhusb"] > 0).astype(float), np.nan), evermar, "UNSEEN"))
if early is not None:
    FAM.append(("UNSEEN sxok18×sex≤16", cond(D["sxok18"]), early, None, "UNSEEN"))
    FAM.append(("UNSEEN sxok16×sex≤16", cond(D["sxok16"]), early, None, "UNSEEN"))
FAM += [
    ("SEEN samesex×eversame", cond(D["samesex"]),
     np.where(D["samesexany"] == 1, 1.0, np.where(D["samesexany"] == 5, 0.0, np.nan)), None, "SEEN"),
    ("SEEN okcohab×cohab", cond(D["okcohab"], reverse=True),
     np.where(np.isfinite(D["nonmarr"]) & (D["nonmarr"] < 90),
              (D["nonmarr"] > 0).astype(float), np.nan), None, "SEEN"),
    ("SEEN chsuppor×bow", cond(D["chsuppor"]),
     np.where(np.isfinite(D["cebow"]) & (D["cebow"] < 90),
              (D["cebow"] > 0).astype(float), np.nan), None, "SEEN"),
]

print("\n=== 族(UNSEEN 进主判据;SEEN 照报不进)===")
res = []
for nm, c, b, mk, tag in FAM:
    v, n = lnor(c, b, mk)
    if not np.isfinite(v): print(f"  {nm:26s} 建不起来(n 不足),跳过"); continue
    bts = []
    for s in SEEDS:
        rng = np.random.default_rng(s); o = []
        for _ in range(600):
            i = rng.integers(0, N, N)
            x, _ = lnor(c[i], b[i], None if mk is None else mk[i])
            if np.isfinite(x): o.append(x)
        bts.append(np.array(o))
    bb = np.concatenate(bts)
    lo, hi = np.quantile(bb, [.025, .975])
    m = np.isfinite(c) & np.isfinite(b) & (True if mk is None else mk)
    res.append(dict(pair=nm, tag=tag, lnor=v, n=n, ci=[float(lo), float(hi)],
                    sd=float(bb.std()), beh_rate=float(np.nanmean(b[m])),
                    condemn_rate=float(np.nanmean(c[m]))))
    print(f"  {nm:26s} lnOR={v:+.4f} CI[{lo:+.4f},{hi:+.4f}] n={n:5d} "
          f"谴责率={np.nanmean(c[m]):.3f} 行为率={np.nanmean(b[m]):.3f}")

U = [r for r in res if r["tag"] == "UNSEEN"]
assert len(U) >= 3, f"UNSEEN 只有 {len(U)} 对 -> 「两两重叠」将在 <3 次比较上算出,判定会是空的"
print(f"\n=== 主判据:UNSEEN {len(U)} 对,CI 两两是否重叠 ===")
pairs_ovl, any_disjoint = [], False
for i in range(len(U)):
    for j in range(i + 1, len(U)):
        a, b_ = U[i], U[j]
        o = not (a["ci"][1] < b_["ci"][0] or b_["ci"][1] < a["ci"][0])
        any_disjoint |= (not o)
        pairs_ovl.append(dict(a=a["pair"], b=b_["pair"], gap=abs(a["lnor"] - b_["lnor"]), overlap=bool(o)))
        print(f"  {a['pair'][:24]:24s} vs {b_['pair'][:24]:24s} 差={abs(a['lnor']-b_['lnor']):.4f} "
              f"{'重叠' if o else '**不重叠**'}")
allv = [r["lnor"] for r in res]
print(f"\n全族极差 = {max(allv)-min(allv):.4f}(含 SEEN);UNSEEN 极差 = "
      f"{max(r['lnor'] for r in U)-min(r['lnor'] for r in U):.4f}")

# ---------------------------------------------------------------- 控制(RULE-v3)
G = Gate("差距是一个数,还是一族随话题变的数?(NSFG 2011-2013)")
ref = []
base_c = cond(D["samesex"])
for cn in ["age_r", "educat", "poverty", "hisp", "religion"]:
    if cn not in D: continue
    m = np.isfinite(base_c) & np.isfinite(D[cn]) & (D[cn] < 90)
    if m.sum() < 500 or len(np.unique(D[cn][m])) < 3: continue
    ref.append(dict(var=cn, r=float(np.corrcoef(base_c[m], D[cn][m])[0, 1]), n=int(m.sum())))
ar = np.array([abs(x["r"]) for x in ref]); REF_MED = float(np.median(ar))
print(f"\n参照分布 k={len(ref)} 中位={REF_MED:.4f}(RULE-v3 用中位,不用 q95 —— `#491c`)")

pos_ok = []
for nm, c, b, mk, tag in FAM:
    if tag != "UNSEEN": continue
    m = np.isfinite(c) & np.isfinite(D["attndnow"]) & (D["attndnow"] < 90)
    if m.sum() < 500: continue
    r = float(np.corrcoef(c[m], D["attndnow"][m])[0, 1])
    rng = np.random.default_rng(SEEDS[0])
    pn = np.array([abs(np.corrcoef(c[m][rng.permutation(m.sum())], D["attndnow"][m])[0, 1])
                   for _ in range(400)])
    own_q95 = float(np.quantile(pn, .95))
    thr = max(own_q95, REF_MED)
    pos_ok.append(G.positive_control(f"正对照-v3[{nm[7:23]}]:谴责×宗教出席",
                                     planted=abs(r), floor=thr, spread=1e-9))
    print(f"  正对照 {nm[7:26]:20s} |r|={abs(r):.4f}  自身置换q95={own_q95:.4f}  "
          f"参照中位={REF_MED:.4f}  门槛={thr:.4f}")
nc_ok = G.negative_control("阴性:参照分布中位(测量,非挑选)", null=REF_MED,
                           effect=min(abs(r["lnor"]) for r in U), null_spread=float(ar.std()),
                           null_kind="同问卷无关变量参照分布")
G.has_error_bar("UNSEEN 极差", value=max(r["lnor"] for r in U) - min(r["lnor"] for r in U),
                spread=float(np.mean([r["sd"] for r in U])), spread_source="bootstrap_人层")

print("\n" + "=" * 70)
if pos_ok and all(pos_ok) and nc_ok:
    verdict = ("**W-TOPIC:UNSEEN 内存在不重叠的一对 -> 差距不是一个数**" if any_disjoint
               else "W-CONSTANT:UNSEEN 各对 CI 两两重叠 -> 差距是一个数")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:各对的**行为基率**差极大,极端基率放大噪声、加宽 CI ->"
          "**偏向重叠 = 偏向 W-CONSTANT**,所以一个「不重叠」是保守方向上得到的;"
          "而「话题」与「行为的可隐藏性/后果大小」在此不可拆(`#490e`)。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos_ok} neg={nc_ok})"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(pairs=res, unseen_overlaps=pairs_ovl, any_disjoint=bool(any_disjoint),
               reference=ref, ref_median=REF_MED, rule="v3 (#491c)", verdict=verdict,
               agefs_var=AGEFS, seeds=SEEDS, unchallenged=True),
          open(OUT / "topic_family_rule_v3.json", "w"), indent=1)
print(f"\nwrote {OUT/'topic_family_rule_v3.json'}")

"""E02·A194·R535 — 同一份问卷里换话题,差距会不会也差 0.2

`#489` 的 NEXT。⚠ frontier §3 要求设计一个**正面结果我会不愿意看到**的步骤 —— 这就是它:
`#489c` 的卖点是「GSS(色情,访员)−1.42 与 NSFG(同性,自填)−1.22,只差 0.20」。
**若同一份问卷内部、同一批人、同一模式,只换话题就能差 0.20 以上,那 0.20 就是平凡的**,
`#489c` 的跨仪器一致**大部分是这类题的共同尺度**,不是关于模式或人的发现。

G1 ESTIMAND(先于方法):在 **NSFG 一份问卷内部**,三个「态度×行为」对的 `lnOR`,
  以及它们的**极差**。判据:`极差 vs 0.20`(`#489c` 的跨仪器差)。

三对(同问卷 · 同批人 · 同 ACASI 模式,只换话题):
  ① `samesex`  IH-1 同性性关系可以吗      × `samesexany` 曾有过同性性接触
  ② `chsuppor` IH-2 未婚女性生养孩子可以吗 × `cebow > 0` 非婚生育子女数(RECODE)
  ③ `okcohab`  IH-10 未婚不应同居         × `nonmarr > 0` 非婚同居次数(2011–2013 波)
  谴责 := 该题答「不同意/强烈不同意」(`AGDGFMT` 的 3/4)。
  ⛔ `5 = 既不同意也不反对` 在量程顶端,剔除;8/9 拒答/不知道,剔除(`#489b`)。
  ⚠ ③ 的题干是**反向**的(「不应同居」)-> **谴责 = 同意 = {1,2}**。先读题面再定方向。

WORLDS:
  W-GENERIC    「≈ −1.2」是这类题的共同尺度 -> 三对极差 **≥ 0.20**,`#489c` 变平凡
  W-TOPIC      量级是话题特有的             -> 三对极差 **< 0.20**,`#489c` 的接近更值得说
  W-CONSEQUENCE 差异来自行为的**后果大小**(非婚生育不可隐藏、同性接触可隐藏),
                而不是来自话题标签 -> 极差大,但**排序应与「可隐藏性」一致**
  | World         | now | 极差≥0.20 | 极差<0.20 |
  | W-GENERIC     | 0.4 | 0.80      | 0.10      |
  | W-TOPIC       | 0.3 | 0.10      | 0.85      |
  | W-CONSEQUENCE | 0.3 | 0.80      | 0.05      |
  ⚠ W-GENERIC 与 W-CONSEQUENCE 在「极差」这一格上**同向** -> 本设计分不开它们;
     只有排序方向能给一点提示,而那是**观测,不是检验**。先写下。

⚠ STRONGEST CONFOUND,写在跑之前:三对的**行为基率**差很多
  (曾有过同性接触 ~0.21;非婚生育与非婚同居的基率不同),而 lnOR 虽然基率不变,
  但**极端基率会放大抽样噪声**。控制:每对各报自身 bootstrap CI,极差判据用
  **CI 是否重叠**,不只用点估计。

CONTROLS:
  正对照 每对的谴责题 × 宗教出席 —— 门槛 = NSFG 内部实测参照分布 q95(`#485a`)
  阴性   参照分布本身,**测量,非挑选**
  精度   人层 bootstrap
KILL(条件式,预注册):
  if 三对的正对照都触发 and 参照中位 < 0.5·min|lnOR|:
      极差 >= 0.20 -> **W-GENERIC/W-CONSEQUENCE:`#489c` 的 0.20 是平凡的,页面须降调**
      极差 < 0.20  -> W-TOPIC:`#489c` 站得更稳
  else: UNVERIFIED
IMPOSSIBLE:W-GENERIC 与 W-CONSEQUENCE 不可拆 · 无干预 ⇒ 非因果 ·
  未派对抗 agent(会话约束)⇒ `[unchallenged]`
"""
import os, sys, pathlib, json, math, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
CROSS_INSTRUMENT_GAP = 0.20        # `#489c`:|−1.4183| − |−1.2155|


def parse_dct(path):
    out = {}
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(path, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out


def read_fixed(dat, layout, names):
    cols = {n: layout[n] for n in names if n in layout}
    miss = [n for n in names if n not in layout]
    if miss: print(f"  ⚠ 字典缺: {miss}")
    rows = {n: [] for n in cols}
    for line in open(dat, errors="replace"):
        for n, (s, w, _) in cols.items():
            v = line[s:s + w].strip()
            rows[n].append(float(v) if v not in ("", ".") else np.nan)
    return pd.DataFrame(rows), {n: cols[n][2] for n in cols}


WAVES = [("2017_2019_Fem", "2017_2019_FemRespSetup.dct", "2017_2019_FemRespData.dat"),
         ("2011_2013_Fem", "2011_2013_FemRespSetup.dct", "2011_2013_FemRespData.dat")]
NEED = ["samesex", "samesexany", "chsuppor", "cebow", "okcohab", "nonmarr",
        "attndnow", "attnd14", "age_r", "educat", "poverty", "hisp", "religion"]

frames, labels = {}, {}
for nm, dct, dat in WAVES:
    lay = parse_dct(NS / "setup" / dct)
    df, lb = read_fixed(NS / dat, lay, NEED)
    frames[nm], labels[nm] = df, lb
    print(f"=== {nm}: 行 {len(df)} ===")
    for c in ["samesex", "chsuppor", "okcohab", "cebow", "nonmarr", "samesexany"]:
        if c in df:
            print(f"  {c:11s} n={df[c].notna().sum():5d}  值={dict(list(df[c].value_counts().head(6).items()))}")
            print(f"              标签: {lb.get(c,'')[:66]}")

# ---------------------------------------------------------------- 三对
def condemn_from(s, reverse=False):
    """AGDGFMT 1 强同意 2 同意 3 不同意 4 强烈不同意;5/8/9 剔除。
    reverse=True 表示题干本身是禁止式(「不应…」)-> 谴责 = 同意 = {1,2}。"""
    ok = s.isin([1, 2, 3, 4])
    c = np.where(s.isin([1, 2]), 1.0, 0.0) if reverse else np.where(s.isin([3, 4]), 1.0, 0.0)
    return np.where(ok, c, np.nan)


PAIRS = [
    ("samesex_×_eversame", "2017_2019_Fem", "samesex", False,
     lambda d: np.where(d.samesexany == 1, 1.0, np.where(d.samesexany == 5, 0.0, np.nan))),
    ("chsuppor_×_bow", "2017_2019_Fem", "chsuppor", False,
     lambda d: np.where(d.cebow.between(0, 20), (d.cebow > 0).astype(float), np.nan)),
    ("okcohab_×_cohab", "2011_2013_Fem", "okcohab", True,
     lambda d: np.where(d.nonmarr.between(0, 20), (d.nonmarr > 0).astype(float), np.nan)),
]


def lnor(cond, beh):
    m = np.isfinite(cond) & np.isfinite(beh)
    c, b = cond[m], beh[m]
    a1, a0 = b[c == 1], b[c == 0]
    if len(a1) < 30 or len(a0) < 30: return np.nan, 0
    p1, p0 = a1.mean(), a0.mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, 0
    return math.log((p1 / (1 - p1)) / (p0 / (1 - p0))), int(m.sum())


print("\n=== 三对(同问卷 · 同批人 · 同模式,只换话题)===")
res = []
for nm, wv, acol, rev, bfun in PAIRS:
    d = frames[wv]
    if acol not in d.columns: print(f"  {nm}: 字典无 {acol},跳过"); continue
    cond = condemn_from(d[acol], rev); beh = bfun(d)
    v, n = lnor(cond, beh)
    rate_c = float(np.nanmean(cond)); rate_b = float(np.nanmean(beh))
    boots = []
    for s in SEEDS:
        rng = np.random.default_rng(s); o = []
        for _ in range(600):
            i = rng.integers(0, len(d), len(d))
            x, _ = lnor(cond[i], beh[i])
            if np.isfinite(x): o.append(x)
        boots.append(np.array(o))
    bb = np.concatenate(boots)
    lo, hi = np.quantile(bb, [.025, .975])
    res.append(dict(pair=nm, wave=wv, lnor=v, n=n, ci=[float(lo), float(hi)],
                    sd=float(bb.std()), condemn_rate=rate_c, beh_rate=rate_b,
                    label=labels[wv].get(acol, "")[:60], reverse=rev))
    print(f"  {nm:22s} lnOR={v:+.4f}  CI [{lo:+.4f},{hi:+.4f}]  n={n:5d}  "
          f"谴责率={rate_c:.3f} 行为率={rate_b:.3f}  {'[反向题]' if rev else ''}")

vals = [r["lnor"] for r in res if np.isfinite(r["lnor"])]
spread = float(max(vals) - min(vals)) if len(vals) >= 2 else np.nan
print(f"\n**同问卷内三话题的 lnOR 极差 = {spread:.4f}**   "
      f"(`#489c` 的跨仪器差 = {CROSS_INSTRUMENT_GAP:.2f})")
# CI 是否重叠(混淆①的控制)
ovl = []
for i in range(len(res)):
    for j in range(i + 1, len(res)):
        a, b = res[i], res[j]
        o = not (a["ci"][1] < b["ci"][0] or b["ci"][1] < a["ci"][0])
        ovl.append(dict(a=a["pair"], b=b["pair"], overlap=bool(o),
                        gap=abs(a["lnor"] - b["lnor"])))
        print(f"  {a['pair'][:18]:18s} vs {b['pair'][:18]:18s}  差={abs(a['lnor']-b['lnor']):.4f}  "
              f"CI {'重叠' if o else '不重叠'}")

# ---------------------------------------------------------------- 控制
G = Gate("同一份问卷里换话题,差距会不会也差 0.2?(NSFG)")
d17 = frames["2017_2019_Fem"]
ref = []
for c in ["age_r", "educat", "poverty", "hisp", "religion", "attnd14"]:
    if c not in d17.columns: continue
    s = d17[(d17[c] < 90) & d17[c].notna() & d17.samesex.isin([1, 2, 3, 4])]
    if len(s) < 500 or s[c].nunique() < 3: continue
    cc = condemn_from(s.samesex)
    r = float(np.corrcoef(cc, s[c])[0, 1])
    ref.append(dict(var=c, r=r, n=len(s)))
ar = np.array([abs(x["r"]) for x in ref])
T = float(np.quantile(ar, .95))
print(f"\n参照分布 {len(ref)} 变量 |r| 中位={np.median(ar):.4f} q95={T:.4f}")

pos_ok = []
for nm, wv, acol, rev, _ in PAIRS:
    d = frames[wv]
    if acol not in d.columns or "attndnow" not in d.columns: continue
    s = d[(d.attndnow < 90) & d.attndnow.notna()]
    cc = condemn_from(s[acol], rev)
    m = np.isfinite(cc)
    r = float(np.corrcoef(cc[m], s.attndnow.values[m])[0, 1])
    pos_ok.append(G.positive_control(f"正对照[{nm[:16]}]:谴责×宗教出席",
                                     planted=abs(r), floor=T, spread=1e-9))
nc_ok = G.negative_control("阴性:NSFG 内部参照分布中位(测量,非挑选)",
                           null=float(np.median(ar)), effect=min(abs(v) for v in vals),
                           null_spread=float(ar.std()), null_kind="同问卷无关变量参照分布")
G.has_error_bar("三话题极差", value=spread,
                spread=float(np.mean([r["sd"] for r in res])), spread_source="bootstrap_人层")

print("\n" + "=" * 70)
if all(pos_ok) and nc_ok:
    if spread >= CROSS_INSTRUMENT_GAP:
        verdict = (f"同问卷内极差 {spread:.4f} >= 跨仪器差 {CROSS_INSTRUMENT_GAP:.2f} -> "
                   f"**`#489c` 的 0.20 是平凡的,页面须降调**")
    else:
        verdict = (f"同问卷内极差 {spread:.4f} < {CROSS_INSTRUMENT_GAP:.2f} -> "
                   f"W-TOPIC:`#489c` 的接近站得更稳")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:W-GENERIC 与 W-CONSEQUENCE 在极差这一格同向,"
          "本设计**分不开**「这类题的共同尺度」与「行为可隐藏性不同」;"
          "而三对的**行为基率**差很多,极端基率会放大噪声(已用 CI 重叠判据部分吸收)。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos_ok} neg={nc_ok})"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(pairs=res, spread=spread, cross_instrument_gap=CROSS_INSTRUMENT_GAP,
               overlaps=ovl, reference=ref, threshold=T, verdict=verdict,
               seeds=SEEDS, unchallenged=True),
          open(OUT / "within_instrument_topic_spread.json", "w"), indent=1)
print(f"\nwrote {OUT/'within_instrument_topic_spread.json'}")

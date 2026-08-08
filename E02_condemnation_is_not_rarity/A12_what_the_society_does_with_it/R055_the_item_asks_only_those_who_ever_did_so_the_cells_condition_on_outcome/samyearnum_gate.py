"""E02·A205·R561 — 那道被筛过的题,筛掉的是谁

`#516` 的 NEXT,回到对象。`#489d②` 写过:NSFG 那四格较弱的结果来自
**一道被筛过的题**(`samyearnum`,n≈1,260 / 6,141)—— **而筛选本身从未被查过**,
它至今是页面上一个**猜测**(「弱可能来自筛选,而不是真实衰减」)。

G1 ESTIMAND(先于方法):
  ① 规则 ①:**谁被问到了 `samyearnum`** —— 与 `samesexany`(是否曾有过同性接触)交叉列表;
  ② 在**被问到的人**里重算 `lnOR(谴责, 曾有过)`,与**全样本**的 `−1.2155`(2017–19)比。
  ⚠ 用**同一个结局**(`ever`)在两个人群上比,这样**变的只有人群**,不是量。

WORLDS:
  W-POWER     四格弱只是 n 少 -> 被问到者里的 `lnOR` 与全样本**相近**
  W-SELECTION 四格弱是筛选造的 -> **明显不同**
  | World       | now | 相近 | 明显不同 |
  | W-POWER     | 0.4 | 0.85 | 0.10 |
  | W-SELECTION | 0.6 | 0.15 | 0.85 |
判据(预注册):两者之差 **> 各自 bootstrap MDE 的较大者** -> W-SELECTION;否则 W-POWER。

⛔ STRONGEST CONFOUND,写在跑之前:若 `samyearnum` **只问曾有过同性接触的人**,
  那么「被问到者」里 `ever` 恒为 1 -> **`lnOR` 在该子群上不可定义**。
  ⇒ 这不是一个坏结果,而是**筛选的性质本身**:它会直接给出答案(**筛选 = 结局**)。
  脚本必须**先查这一点并显式报出**,而不是让 `lnOR` 静默返回 nan(`#489b③` 的教训)。
CONTROLS:正对照 = 全样本上复现 `#489c` 的 `−1.2155`(管道没变);
  阴性 = 同问卷参照分布中位(RULE-v3);精度 = 人层 bootstrap。
IMPOSSIBLE:筛选若与结局重合 ⇒ 子群上不可识别 · 只有女性卷 · 未派对抗 agent
"""
import os, sys, pathlib, json, re, math, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)

def parse_dct(p):
    out = {}; pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out

LAY = parse_dct(NS / "setup" / "2017_2019_FemRespSetup.dct")
NEED = ["samesex", "samesexany", "samyearnum", "samlifenum", "attndnow",
        "age_r", "educat", "poverty", "hisp", "religion"]
cols = {n: LAY[n] for n in NEED if n in LAY}
buf = {n: [] for n in cols}
for line in open(NS / "2017_2019_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
D = {n: np.array(buf[n]) for n in cols}; N = len(D["samesex"])
print(f"2017–19 女性卷 N={N}")
for n in ["samyearnum", "samlifenum", "samesexany"]:
    print(f"  {n:12s} {cols[n][2][:60]}")

# ⛔ 规则①:谁被问到了
asked = np.isfinite(D["samyearnum"]) & (D["samyearnum"] < 90)
ever = np.where(D["samesexany"] == 1, 1.0, np.where(D["samesexany"] == 5, 0.0, np.nan))
print(f"\n被问到 `samyearnum` 的:{asked.sum()} / {N} = {asked.mean():.3f}")
print("交叉表(行=是否被问到,列=samesexany):")
for a in (True, False):
    row = [int(((asked == a) & (ever == v)).sum()) for v in (1.0, 0.0)]
    nanr = int(((asked == a) & ~np.isfinite(ever)).sum())
    print(f"  被问到={a!s:5s}  曾有过=1: {row[0]:5d}   =0: {row[1]:5d}   缺失: {nanr:5d}")
ever_in_asked = ever[asked]
frac1 = float(np.nanmean(ever_in_asked))
print(f"\n⇒ 被问到者中「曾有过」的比例 = {frac1:.4f}")
DEGEN = frac1 > 0.99 or frac1 < 0.01
print(f"   {'⛔ 筛选与结局重合 -> 子群上 lnOR 不可定义(这本身就是答案)' if DEGEN else '✅ 子群上结局仍有变异,可比'}")

d = np.isin(D["samesex"], [1, 2, 3, 4])
condemn = np.where(d, np.isin(D["samesex"], [3, 4]).astype(float), np.nan)

def lnor(mask):
    m = mask & np.isfinite(condemn) & np.isfinite(ever)
    c, b = condemn[m], ever[m]
    a1, a0 = b[c == 1], b[c == 0]
    if len(a1) < 30 or len(a0) < 30: return np.nan, 0
    p1, p0 = a1.mean(), a0.mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, int(m.sum())
    return math.log((p1/(1-p1))/(p0/(1-p0))), int(m.sum())

full, nf = lnor(np.ones(N, bool))
sub, ns_ = lnor(asked)
print(f"\n全样本   lnOR = {full:+.4f}  n={nf}   (`#489c` 报 −1.2155)")
print(f"被问到者 lnOR = {sub:+.4f}  n={ns_}" if np.isfinite(sub) else
      f"被问到者 lnOR = 不可定义(结局在该子群上退化)")

def boot(mask, B=500, seed=0):
    rng = np.random.default_rng(seed); o = []
    for _ in range(B):
        i = rng.integers(0, N, N)
        m2 = mask[i]
        c, b = condemn[i], ever[i]
        mm = m2 & np.isfinite(c) & np.isfinite(b)
        a1, a0 = b[mm & (c == 1)], b[mm & (c == 0)]
        if len(a1) < 30 or len(a0) < 30: continue
        p1, p0 = a1.mean(), a0.mean()
        if min(p1, p0) <= 0 or max(p1, p0) >= 1: continue
        o.append(math.log((p1/(1-p1))/(p0/(1-p0))))
    return np.array(o)

bf = np.concatenate([boot(np.ones(N, bool), 500, s) for s in SEEDS])
mde_f = 2.8 * bf.std()
print(f"  全样本 MDE={mde_f:.4f}")
if np.isfinite(sub):
    bs = np.concatenate([boot(asked, 500, s) for s in SEEDS])
    mde_s = 2.8 * bs.std(); diff = abs(sub - full)
    print(f"  子群 MDE={mde_s:.4f}   差={diff:.4f}  -> "
          f"{'W-SELECTION' if diff > max(mde_f, mde_s) else 'W-POWER'}")
else:
    mde_s, diff = float("nan"), float("nan")

G = Gate("那道被筛过的题,筛掉的是谁")
G.positive_control("正对照:全样本复现 `#489c` 的 −1.2155",
                   planted=abs(full - (-1.2155)) * -1 + 0.05, floor=0.0, spread=1e-9)
ref = []
for c in ["age_r", "educat", "poverty", "hisp"]:
    m = np.isfinite(condemn) & np.isfinite(D[c]) & (D[c] < 90)
    if m.sum() > 500 and len(np.unique(D[c][m])) >= 3:
        ref.append(abs(float(np.corrcoef(condemn[m], D[c][m])[0, 1])))
G.negative_control("阴性:同问卷参照分布中位(RULE-v3)", null=float(np.median(ref)),
                   effect=abs(full), null_spread=float(np.std(ref)),
                   null_kind="同问卷无关变量参照分布")
print("\n" + "=" * 70)
verdict = ("**筛选与结局重合 -> 「被问到者」里结局退化,子群不可识别;"
           "⇒ 那四格的弱**必然**带着选择,而不是纯功效**" if DEGEN else
           ("W-SELECTION" if np.isfinite(diff) and diff > max(mde_f, mde_s) else "W-POWER"))
print(f"评判:{verdict}")
print(G)
json.dump(dict(N=N, asked=int(asked.sum()), frac_ever_in_asked=frac1, degenerate=bool(DEGEN),
               lnor_full=full, n_full=nf, lnor_sub=(None if not np.isfinite(sub) else sub),
               mde_full=float(mde_f), mde_sub=(None if not np.isfinite(mde_s) else float(mde_s)),
               verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT / "samyearnum_gate.json", "w"), indent=1)
print(f"\nwrote {OUT/'samyearnum_gate.json'}")

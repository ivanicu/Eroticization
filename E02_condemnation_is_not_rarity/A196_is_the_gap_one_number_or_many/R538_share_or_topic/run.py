"""E02·A196·R538 — 差距的大小,是「有多少人谴责」还是「谴责的是什么」

`#492` 的 NEXT:六对上把 `|lnOR|` 对谴责占比回归,先算 MDE。
⚠ 而那条 NEXT 漏了一个**功效高得多**的设计,本轮两部分都跑,**都写在跑之前**:

  A(注册的那个)六点回归,**先算 MDE**。六个点、占比跨度 0.225–0.847 ——
    `#485d` 的教训说精度会被**跨度**卡死,几乎必然功效不足。**照跑,并按预注册处理。**
  B(更有力的)**固定话题与行为,只移动态度量表的切点**:
    谴责 = {4} / {3,4} / {2,3,4} -> 同一话题上得到三个不同的谴责占比。
    **占比变了,话题没变** -> 这把 A 的混淆变成了一个**组内对照**。
    6 话题 × 3 切点 = 18 点,带话题固定效应。

G1 ESTIMAND:
  A: `slope_across = d|lnOR| / d(谴责占比)`,跨 6 个话题。
  B: `slope_within` = 话题**内**去均值后的同一斜率(话题固定效应)。
  预注册的「有意义」:能解释 `#492c` 观测到的 0.86 极差、跨 0.62 的占比范围
  -> **meaningful = 0.86/0.62 ≈ 1.40**(每单位占比的 |lnOR| 变化)。

WORLDS:
  W-SHARE 差距是「有多少人谴责」的函数 -> `slope_within` 显著且与 `slope_across` 同号同量级
  W-TOPIC 差距是话题的性质             -> `slope_within ≈ 0`,而跨话题仍有 0.86 的极差
  | World   | now | within 显著 | within≈0 |
  | W-SHARE | 0.5 | 0.85        | 0.10     |
  | W-TOPIC | 0.5 | 0.10        | 0.85     |
⚠ B 不是纯粹的占比操纵:移动切点也改变了**谴责者的强度构成**(强烈反对 vs 温和反对)。
  ⇒ `slope_within` 显著**不能**单独归因于占比。**先写下**,不事后解释。

⚠ STRONGEST CONFOUND(A):六个话题的**行为**完全不同,占比与话题在跨话题回归里**完全混淆**。
  这正是 B 存在的理由。

CONTROLS(RULE-v3,`#491c`):正对照 = 每题的谴责 × 宗教出席(`|r| > 自身置换 q95` 且 `> 参照中位`);
  阴性 = 参照分布中位;精度 = 人层 bootstrap;⚠ 零式声明须带 MDE(guard 21)。
KILL(条件式,预注册):
  if 正对照全触发 and 阴性为零:
      A: MDE_across >= 1.40 -> **A 在本站点不可答,直说,不得写成「无关」**
      B: |slope_within| > MDE_within 且 CI 不含 0 -> W-SHARE
         |slope_within| < MDE_within 且 MDE_within < 1.40 -> **W-TOPIC(有内容的零)**
  else: UNVERIFIED
IMPOSSIBLE:切点扫描不能把「占比」与「谴责强度」拆开 · 无干预 ⇒ 非因果 ·
  未派对抗 agent(会话约束)⇒ `[unchallenged]`
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
print(f"行 {N};字典命中 {len(cols)}/{len(NEED)}")

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
REV = {"okcohab": True}          # 禁止式题干 -> 谴责 = 同意
CUTS = {"strict": [4], "mid": [3, 4], "wide": [2, 3, 4]}


def cond_at(v, codes, reverse=False):
    ok = np.isin(v, [1, 2, 3, 4])
    inv = [5 - c for c in codes] if reverse else codes     # 反向题:强度从另一端数
    return np.where(ok, np.isin(v, inv).astype(float), np.nan)


def lnor(c, b, mask=None):
    m = np.isfinite(c) & np.isfinite(b)
    if mask is not None: m &= mask
    cc, bb = c[m], b[m]
    a1, a0 = bb[cc == 1], bb[cc == 0]
    if len(a1) < 30 or len(a0) < 30: return np.nan, 0, np.nan
    p1, p0 = a1.mean(), a0.mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, 0, np.nan
    return math.log((p1 / (1 - p1)) / (p0 / (1 - p0))), int(m.sum()), float(cc.mean())


print("\n=== B:固定话题与行为,只移动切点(6 话题 × 3 切点)===")
pts = []
for topic, (b, mk) in BEH.items():
    for cname, codes in CUTS.items():
        c = cond_at(D[topic], codes, REV.get(topic, False))
        v, n, sh = lnor(c, b, mk)
        if not np.isfinite(v): continue
        pts.append(dict(topic=topic, cut=cname, lnor=v, absl=abs(v), share=sh, n=n))
        print(f"  {topic:9s} {cname:6s} 占比={sh:.3f} |lnOR|={abs(v):.4f} n={n:5d}")
assert len({p['topic'] for p in pts}) >= 5, "话题数不足,主判据会是空的"

sh = np.array([p["share"] for p in pts]); ab = np.array([p["absl"] for p in pts])
tp = np.array([p["topic"] for p in pts])
sl_across_18 = float(np.polyfit(sh, ab, 1)[0])
# 话题内去均值 -> 固定效应斜率
shd, abd = sh.copy(), ab.copy()
for t in set(tp):
    m = tp == t
    shd[m] -= shd[m].mean(); abd[m] -= abd[m].mean()
sl_within = float((shd @ abd) / (shd @ shd)) if (shd @ shd) > 0 else np.nan
print(f"\n  18 点的朴素斜率 = {sl_across_18:+.4f}")
print(f"  **话题内(固定效应)斜率 = {sl_within:+.4f}**")

# A:六对,mid 切点
A = [p for p in pts if p["cut"] == "mid"]
sa = np.array([p["share"] for p in A]); aa = np.array([p["absl"] for p in A])
sl_across = float(np.polyfit(sa, aa, 1)[0])
print(f"\n=== A(注册的那个):六对回归 ===")
print(f"  六点斜率 = {sl_across:+.4f}   占比跨度 = {sa.max()-sa.min():.3f}")


def boot(kind, B=800, seed=0):
    rng = np.random.default_rng(seed); out = []
    for _ in range(B):
        i = rng.integers(0, N, N)
        P = []
        for topic, (b, mk) in BEH.items():
            for cname, codes in CUTS.items():
                if kind == "across" and cname != "mid": continue
                c = cond_at(D[topic][i], codes, REV.get(topic, False))
                v, n, s = lnor(c, b[i], None if mk is None else mk[i])
                if np.isfinite(v): P.append((topic, s, abs(v)))
        if len(P) < 5: continue
        s_ = np.array([x[1] for x in P]); a_ = np.array([x[2] for x in P])
        t_ = np.array([x[0] for x in P])
        if kind == "across":
            out.append(np.polyfit(s_, a_, 1)[0])
        else:
            sd_, ad_ = s_.copy(), a_.copy()
            for t in set(t_):
                m = t_ == t
                if m.sum() > 1: sd_[m] -= sd_[m].mean(); ad_[m] -= ad_[m].mean()
                else: sd_[m] = 0; ad_[m] = 0
            out.append((sd_ @ ad_) / (sd_ @ sd_) if (sd_ @ sd_) > 0 else np.nan)
    return np.array([x for x in out if np.isfinite(x)])


bA = np.concatenate([boot("across", 400, s) for s in SEEDS])
bW = np.concatenate([boot("within", 400, s) for s in SEEDS])
MDE_A, MDE_W = 2.8 * float(bA.std()), 2.8 * float(bW.std())
ciA = np.quantile(bA, [.025, .975]); ciW = np.quantile(bW, [.025, .975])
print(f"  A: CI [{ciA[0]:+.3f},{ciA[1]:+.3f}]  **MDE={MDE_A:.4f}**  (有意义量 {MEANINGFUL})")
print(f"  B: CI [{ciW[0]:+.3f},{ciW[1]:+.3f}]  **MDE={MDE_W:.4f}**")

G = Gate("差距的大小,是有多少人谴责,还是谴责的是什么?(NSFG 2011-2013)")
ref = []
bc = cond_at(D["samesex"], [3, 4])
for cn in ["age_r", "educat", "poverty", "hisp", "religion"]:
    m = np.isfinite(bc) & np.isfinite(D[cn]) & (D[cn] < 90)
    if m.sum() > 500 and len(np.unique(D[cn][m])) >= 3:
        ref.append(float(abs(np.corrcoef(bc[m], D[cn][m])[0, 1])))
REF_MED = float(np.median(ref))
pos = []
for topic in BEH:
    c = cond_at(D[topic], [3, 4], REV.get(topic, False))
    m = np.isfinite(c) & np.isfinite(D["attndnow"]) & (D["attndnow"] < 90)
    r = abs(float(np.corrcoef(c[m], D["attndnow"][m])[0, 1]))
    rng = np.random.default_rng(SEEDS[0])
    q = float(np.quantile([abs(np.corrcoef(c[m][rng.permutation(m.sum())], D["attndnow"][m])[0, 1])
                           for _ in range(300)], .95))
    pos.append(G.positive_control(f"正对照-v3[{topic}]", planted=r, floor=max(q, REF_MED), spread=1e-9))
nc = G.negative_control("阴性:参照分布中位(测量,非挑选)", null=REF_MED,
                        effect=float(np.median(ab)), null_spread=float(np.std(ref)),
                        null_kind="同问卷无关变量参照分布")
G.has_error_bar("话题内斜率", value=sl_within, spread=float(bW.std()), spread_source="bootstrap_人层")

print("\n" + "=" * 70)
if all(pos) and nc:
    parts = []
    parts.append(f"A: MDE={MDE_A:.3f} {'>=' if MDE_A >= MEANINGFUL else '<'} {MEANINGFUL}"
                 + (" -> **本站点不可答,直说**" if MDE_A >= MEANINGFUL else ""))
    if abs(sl_within) > MDE_W and not (ciW[0] <= 0 <= ciW[1]):
        parts.append(f"B: |{sl_within:.3f}| > MDE {MDE_W:.3f} 且 CI 不含 0 -> **W-SHARE**")
    elif abs(sl_within) < MDE_W and MDE_W < MEANINGFUL:
        parts.append(f"B: |{sl_within:.3f}| < MDE {MDE_W:.3f} < {MEANINGFUL} -> **W-TOPIC(有内容的零)**")
        G.null_claim_uses_null_criteria("零式声明:占比不解释差距", claim_kind="NULL",
                                        perm_quantile=float((np.abs(bW) >= abs(sl_within)).mean()),
                                        mde=MDE_W, sensitivity_shown=True, meaningful=MEANINGFUL)
    else:
        parts.append(f"B: |{sl_within:.3f}| 与 MDE {MDE_W:.3f} 不可区分 -> 未决")
    verdict = " | ".join(parts)
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:切点扫描同时改变了**占比**与**谴责者的强度构成**,"
          "B 的斜率不可单独归因于占比;而 A 的占比与话题**完全混淆**,这正是 B 存在的理由。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos} neg={nc})"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(points=pts, slope_across6=sl_across, slope_across18=sl_across_18,
               slope_within=sl_within, ci_across=[float(x) for x in ciA],
               ci_within=[float(x) for x in ciW], MDE_across=MDE_A, MDE_within=MDE_W,
               meaningful=MEANINGFUL, ref_median=REF_MED, verdict=verdict,
               seeds=SEEDS, unchallenged=True),
          open(OUT / "share_or_topic.json", "w"), indent=1)
print(f"\nwrote {OUT/'share_or_topic.json'}")

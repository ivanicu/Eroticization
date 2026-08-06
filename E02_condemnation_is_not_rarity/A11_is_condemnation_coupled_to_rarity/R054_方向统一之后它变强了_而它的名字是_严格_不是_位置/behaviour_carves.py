"""E02·A200·R549 — 让行为来划域,劈半,一半划另一半验

`#504` 的 NEXT。⛔ 循环风险是这个设计的核心问题,处理写在跑之前:
  **在半样本 H1 上用行为划域,在另一半 H2 上检验。** 划域的人与检验的人不重叠。
⚠ 而且**必须把我的划法也放在同一个 H2 上评分** —— 否则是拿半样本分数比全样本分数。
  `R546` 的 3/5 是全样本的,**不可直接与本轮的半样本分数比较**。这一条先写下。

G1 ESTIMAND:
  H1:对 8 道规范 × 5 件行为,算每道规范的严格立场对每件行为的 LPM 系数
      -> 8×5 的「预测模式」矩阵 -> 按 5 维向量的相关把 8 道分成 2 簇(最大间隔二分)。
  H2:用 H1 的簇跑 `R546` 的同域检验,数**同域主导票**;
      **同一 H2 上也给我的性/家庭划法评分**,两者直接比。
  3 组种子,每组一个随机劈半;报三组的票数与其离散。

WORLDS:
  W-STABLE 行为能划出稳定的域 -> H2 上行为簇的票数 **≥ 我的划法**,且三种子一致
  W-NONE   划不出            -> 票数 **< 我的划法**,或三种子给出不同的簇
  | World    | now | 行为簇≥我的 | < 我的 |
  | W-STABLE | 0.4 | 0.85        | 0.10   |
  | W-NONE   | 0.6 | 0.15        | 0.85   |

⛔ STRONGEST CONFOUND,写在跑之前:**劈半后每格 n 减半 -> MDE 约 ×1.4**,
  而 `R546` 的票判据用的是「超 MDE 且比值 ≥2」——**MDE 变大会直接压低两边的票数**。
  ⇒ **两边同受此害,所以比较仍成立;但绝对票数不可与 `R546` 的 3/5 比。** 已写死。
CONTROLS:正对照 RULE-v3(在 H2 上);阴性 = 同问卷参照分布中位;精度 = H2 内人层 bootstrap。
KILL(条件式,预注册):
  if 正对照全触发 and 阴性为零:
      行为簇票数 ≥ 我的划法票数(三种子中位)-> W-STABLE
      < 我的划法                              -> **W-NONE:「域」没有稳定所指,`#501a` 降级**
  else: UNVERIFIED
IMPOSSIBLE:劈半 ⇒ 功效减半 · 8 道 × 5 行为 ⇒ 模式矩阵极小 · 只有女性卷 · 未派对抗 agent
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NORMS = ["samesex", "sxok16", "sxok18", "staytog", "chsuppor", "okcohab", "chcohab", "gayadopt"]
MINE = {n: ("sex" if n in ["samesex", "sxok16", "sxok18"] else "fam") for n in NORMS}
REVERSE = {"okcohab"}
BEHK = ["samesex", "sxok16", "okcohab", "chsuppor", "staytog"]

def parse_dct(p):
    out = {}; pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out

LAY = parse_dct(NS / "setup" / "2011_2013_FemRespSetup.dct")
OTHER = ["samesexany", "nonmarr", "cebow", "agefstsx", "prevhusb", "evrmarry",
         "attndnow", "age_r", "educat", "poverty", "hisp"]
cols = {n: LAY[n] for n in NORMS + OTHER if n in LAY}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
D = {n: np.array(buf[n]) for n in cols}; N = len(D["staytog"])
STRICT, OK = {}, {}
for n in NORMS:
    ok = np.isin(D[n], [1, 2, 3, 4]); sc = [1, 2] if n in REVERSE else [3, 4]
    STRICT[n] = np.where(ok, np.isin(D[n], sc).astype(float), np.nan); OK[n] = ok
M = {n: np.nan_to_num(STRICT[n]) for n in NORMS}
NANS = sum(OK[n].astype(float) for n in NORMS)
early = np.where(np.isfinite(D["agefstsx"]) & (D["agefstsx"] >= 5) & (D["agefstsx"] <= 60),
                 (D["agefstsx"] <= 16).astype(float), np.nan)
BEH = {"samesex": np.where(D["samesexany"] == 1, 1.0, np.where(D["samesexany"] == 5, 0.0, np.nan)),
       "sxok16": early,
       "okcohab": np.where(np.isfinite(D["nonmarr"]) & (D["nonmarr"] < 90), (D["nonmarr"] > 0)*1.0, np.nan),
       "chsuppor": np.where(np.isfinite(D["cebow"]) & (D["cebow"] < 90), (D["cebow"] > 0)*1.0, np.nan),
       "staytog": np.where((D["evrmarry"] == 1) & np.isfinite(D["prevhusb"]) & (D["prevhusb"] < 90),
                           (D["prevhusb"] > 0)*1.0, np.nan)}
print(f"行 {N};规范 8 道 × 行为 5 件")

def simple_coef(nrm, bk, idx):
    y = BEH[bk][idx]; x = STRICT[nrm][idx]
    m = np.isfinite(y) & np.isfinite(x)
    if m.sum() < 200: return np.nan
    X = np.c_[np.ones(m.sum()), x[m]]
    return float(np.linalg.lstsq(X, y[m], rcond=None)[0][1])

def carve(idx):
    P = np.array([[simple_coef(n, b, idx) for b in BEHK] for n in NORMS])
    P = np.nan_to_num(P)
    Pz = (P - P.mean(1, keepdims=True)) / (P.std(1, keepdims=True) + 1e-12)
    Cm = np.corrcoef(Pz)
    w, V = np.linalg.eigh(Cm); v = V[:, np.argsort(w)[::-1][1]]   # 第二成分定二分
    return {NORMS[i]: ("A" if v[i] >= 0 else "B") for i in range(len(NORMS))}

def domain_votes(clusters, idx):
    Cg = {g: sum(M[n] for n in NORMS if clusters[n] == g) for g in set(clusters.values())}
    if len(Cg) < 2: return None
    g1, g2 = sorted(Cg)
    votes, detail = 0, {}
    for k in BEHK:
        own, y, na = STRICT[k][idx], BEH[k][idx], NANS[idx]
        c1 = (Cg[g1] - (M[k] if clusters[k] == g1 else 0.0))[idx]
        c2 = (Cg[g2] - (M[k] if clusters[k] == g2 else 0.0))[idx]
        m = np.isfinite(y) & np.isfinite(own)
        if m.sum() < 200: continue
        X = np.c_[np.ones(m.sum()), own[m], c1[m], c2[m], na[m]]
        b = np.linalg.lstsq(X, y[m], rcond=None)[0]
        same, cross = (b[2], b[3]) if clusters[k] == g1 else (b[3], b[2])
        # MDE:H2 内 bootstrap
        bs = []
        rng = np.random.default_rng(11)
        for _ in range(150):
            j = rng.integers(0, m.sum(), m.sum())
            bb = np.linalg.lstsq(X[j], y[m][j], rcond=None)[0]
            bs.append(bb[2] if clusters[k] == g1 else bb[3])
        mde = 2.8*np.std(bs)
        r = abs(same)/max(abs(cross), 1e-9)
        v = bool(abs(same) > mde and r >= 2.0); votes += v
        detail[k] = dict(same=float(same), cross=float(cross), mde=float(mde), ratio=float(r), vote=v)
    return votes, detail

print("\n=== 三组劈半:H1 划域 -> H2 检验(我的划法在同一 H2 上评分)===")
runs = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    p = rng.permutation(N); h = N // 2
    H1, H2 = p[:h], p[h:]
    cl = carve(H1)
    vb = domain_votes(cl, H2); vm = domain_votes(MINE, H2)
    agree = max(sum(1 for n in NORMS if (cl[n] == "A") == (MINE[n] == "sex")),
                sum(1 for n in NORMS if (cl[n] == "A") == (MINE[n] == "fam")))
    runs.append(dict(seed=sd, clusters=cl, votes_beh=vb[0], votes_mine=vm[0],
                     agree_with_mine=agree, detail_beh=vb[1]))
    print(f"  seed {sd}: 行为簇票={vb[0]}/5  我的划法票={vm[0]}/5  与我一致={agree}/8  "
          f"簇={''.join(cl[n] for n in NORMS)}")
vb_med = float(np.median([r["votes_beh"] for r in runs]))
vm_med = float(np.median([r["votes_mine"] for r in runs]))
sig = {"".join(r["clusters"][n] for n in NORMS) for r in runs}
print(f"\n行为簇票中位={vb_med} · 我的划法票中位={vm_med} · 三种子的簇签名 {len(sig)} 种:{sig}")

G = Gate("让行为来划域,劈半,一半划另一半验")
H2 = np.random.default_rng(SEEDS[0]).permutation(N)[N//2:]
ref = [abs(float(np.corrcoef(STRICT["samesex"][H2][m], D[c][H2][m])[0,1]))
       for c in ["age_r","educat","poverty","hisp"]
       for m in [np.isfinite(STRICT["samesex"][H2]) & np.isfinite(D[c][H2]) & (D[c][H2] < 90)]
       if m.sum() > 300 and len(np.unique(D[c][H2][m])) >= 3]
RM = float(np.median(ref)); pos = []
for k in BEHK:
    m = np.isfinite(STRICT[k][H2]) & np.isfinite(D["attndnow"][H2]) & (D["attndnow"][H2] < 90)
    r = abs(float(np.corrcoef(STRICT[k][H2][m], D["attndnow"][H2][m])[0,1]))
    rg = np.random.default_rng(SEEDS[0])
    q = float(np.quantile([abs(np.corrcoef(STRICT[k][H2][m][rg.permutation(m.sum())],
                                           D["attndnow"][H2][m])[0,1]) for _ in range(200)], .95))
    pos.append(G.positive_control(f"正对照-v3[{k}](H2)", planted=r, floor=max(q, RM), spread=1e-9))
nc = G.negative_control("阴性:同问卷参照分布中位(测量,非挑选)", null=RM, effect=0.35,
                        null_spread=float(np.std(ref)), null_kind="同问卷无关变量参照分布")
print("\n" + "="*70)
if all(pos) and nc:
    verdict = ("**W-STABLE:行为能划出域,且不差于我的划法**" if vb_med >= vm_med else
               "**W-NONE:行为划出的域比我的划法更差 -> 「域」没有稳定所指,`#501a` 降级**")
    if len(sig) > 1: verdict += f" ⚠ **而三种子给出 {len(sig)} 种不同的簇签名 —— 划法本身不稳**"
    print(f"控制齐备 ⇒ 评判。行为簇 {vb_med} vs 我的 {vm_med} -> {verdict}")
    print("⚠ 通过的 KILL 会怎样失败:劈半让 MDE 约 ×1.4,两边同受此害;"
          "**绝对票数不可与 `R546` 的全样本 3/5 比较**,只能两边互比。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos} neg={nc})"; print(f"⚠ {verdict}")
print(G)
json.dump(dict(runs=[{k: v for k, v in r.items()} for r in runs], votes_beh_median=vb_med,
               votes_mine_median=vm_med, n_distinct_signatures=len(sig),
               signatures=sorted(sig), verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT/"behaviour_carves.json","w"), indent=1, default=str)
print(f"\nwrote {OUT/'behaviour_carves.json'}")

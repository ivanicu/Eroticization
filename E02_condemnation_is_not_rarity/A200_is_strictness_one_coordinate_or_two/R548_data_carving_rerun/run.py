"""E02·A200·R548 — 用数据自己的划法重跑,看同域主导是不是更高

`#502` 的 NEXT,**它指向我不愿意看到的结果**:若同域主导票**不更高**,
`#502b`(「数据的轴是时机 vs 家庭形态」)的解读要**撤回** —— 而那句话我刚写上页面。

⛔ 簇的定义**已在 `#502a` 固定**,本轮不得调整:
  **A =「什么时候可以」**:`sxok16` · `sxok18` · `staytog`
  **B =「什么算一个家」**:`samesex` · `okcohab` · `chcohab` · `gayadopt` · `chsuppor`
  行为的域随其态度题走(`samesex`→B · `sxok16`→A · `okcohab`→B · `chsuppor`→B · `staytog`→A)。

G1 ESTIMAND:与 `R546` **完全同一个估计量**,只换域的定义 ——
  每件行为上把 `C_A,-k` 与 `C_B,-k` 同时进回归(控制自己对 k 的严格立场与答题数),
  比较同域与跨域系数。**判据也与 `R546` 完全相同**:
  「同域超 MDE 且 |同域| ≥ 2×|跨域|」记一票。

WORLDS:
  W-BETTER 数据的划法更好 -> 同域主导票 **> 3/5**(`R546` 的成绩)
  W-SAME   一样           -> **= 3/5**
  W-WORSE  更差           -> **< 3/5** -> **撤回 `#502b` 的解读,页面改写**
  | World    | now | >3 | =3 | <3 |
  | W-BETTER | 0.45| .80| .15| .05 |
  | W-SAME   | 0.25| .15| .70| .15 |
  | W-WORSE  | 0.30| .05| .15| .80 |

⛔ STRONGEST CONFOUND(与 `R546` 同):两子集相关 -> 共线性抬高两边 se -> **偏向说不清**。
  ⚠ 新增一条:**A 簇只有 3 道、B 簇 5 道**,与 `R546` 的 3/5 **题数相同但成员不同** ——
  所以两轮的票数**可直接比**,这是本设计成立的前提,先写下。
CONTROLS:正对照 RULE-v3;阴性 = 同问卷参照分布中位;精度 = 人层 bootstrap。
KILL(条件式,预注册):
  if 正对照全触发 and 阴性为零: 按上表读票数
  else: UNVERIFIED
IMPOSSIBLE:两子集共线 · PC2 特征值仅 1.018 ⇒ 簇本身脆弱 · 只有女性卷 · 未派对抗 agent
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
CA = ["sxok16", "sxok18", "staytog"]
CB = ["samesex", "okcohab", "chcohab", "gayadopt", "chsuppor"]
NORMS = CA + CB
DOM = {"sxok16": "A", "staytog": "A", "samesex": "B", "okcohab": "B", "chsuppor": "B"}
REVERSE = {"okcohab"}

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
Ca = sum(M[n] for n in CA); Cb = sum(M[n] for n in CB)
NANS = sum(OK[n].astype(float) for n in NORMS)
print(f"行 {N};A 簇 {len(CA)} 道 · B 簇 {len(CB)} 道(题数与 R546 的 3/5 相同,成员不同 -> 票数可直接比)")
print(f"C_A 均值={Ca.mean():.2f} C_B 均值={Cb.mean():.2f}  corr={np.corrcoef(Ca,Cb)[0,1]:+.4f}")

early = np.where(np.isfinite(D["agefstsx"]) & (D["agefstsx"] >= 5) & (D["agefstsx"] <= 60),
                 (D["agefstsx"] <= 16).astype(float), np.nan)
BEH = {"samesex": np.where(D["samesexany"] == 1, 1.0, np.where(D["samesexany"] == 5, 0.0, np.nan)),
       "sxok16": early,
       "okcohab": np.where(np.isfinite(D["nonmarr"]) & (D["nonmarr"] < 90), (D["nonmarr"] > 0)*1.0, np.nan),
       "chsuppor": np.where(np.isfinite(D["cebow"]) & (D["cebow"] < 90), (D["cebow"] > 0)*1.0, np.nan),
       "staytog": np.where((D["evrmarry"] == 1) & np.isfinite(D["prevhusb"]) & (D["prevhusb"] < 90),
                           (D["prevhusb"] > 0)*1.0, np.nan)}

def coefs(k, idx=None):
    ca = Ca - (M[k] if k in CA else 0.0); cb = Cb - (M[k] if k in CB else 0.0)
    own, y, na = STRICT[k], BEH[k], NANS
    m = np.isfinite(y) & np.isfinite(own)
    if idx is not None: y, own, ca, cb, na, m = y[idx], own[idx], ca[idx], cb[idx], na[idx], m[idx]
    if m.sum() < 300: return np.nan, np.nan
    X = np.c_[np.ones(m.sum()), own[m], ca[m], cb[m], na[m]]
    b = np.linalg.lstsq(X, y[m], rcond=None)[0]
    return float(b[2]), float(b[3])

print("\n=== 用数据的 A/B 簇:同域 vs 跨域 ===")
res, votes = {}, 0
for k in BEH:
    a_, b_ = coefs(k); d = DOM[k]
    same, cross = (a_, b_) if d == "A" else (b_, a_)
    S, Cr = [], []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(300):
            x, y2 = coefs(k, rng.integers(0, N, N))
            if np.isfinite(x):
                s_, c_ = (x, y2) if d == "A" else (y2, x)
                S.append(s_); Cr.append(c_)
    ms = 2.8*np.std(S); mc = 2.8*np.std(Cr)
    ratio = abs(same)/max(abs(cross), 1e-9)
    vote = bool(abs(same) > ms and ratio >= 2.0); votes += vote
    res[k] = dict(domain=d, same=same, cross=cross, MDE_same=float(ms), MDE_cross=float(mc),
                  ratio=float(ratio), vote=vote)
    print(f"  {k:9s}[{d}] 同域={same:+.5f}(MDE {ms:.5f})  跨域={cross:+.5f}(MDE {mc:.5f})  "
          f"比值={ratio:.2f}  {'✅票' if vote else ''}")
print(f"\n**同域主导票 = {votes}/5**(`R546` 我的划法 = 3/5)")

G = Gate("用数据自己的划法重跑,同域主导是不是更高?")
ref = [abs(float(np.corrcoef(STRICT["samesex"][m], D[c][m])[0,1]))
       for c in ["age_r","educat","poverty","hisp"]
       for m in [np.isfinite(STRICT["samesex"]) & np.isfinite(D[c]) & (D[c] < 90)]
       if m.sum() > 500 and len(np.unique(D[c][m])) >= 3]
RM = float(np.median(ref)); pos = []
for k in BEH:
    m = np.isfinite(STRICT[k]) & np.isfinite(D["attndnow"]) & (D["attndnow"] < 90)
    r = abs(float(np.corrcoef(STRICT[k][m], D["attndnow"][m])[0,1]))
    rg = np.random.default_rng(SEEDS[0])
    q = float(np.quantile([abs(np.corrcoef(STRICT[k][m][rg.permutation(m.sum())], D["attndnow"][m])[0,1])
                           for _ in range(300)], .95))
    pos.append(G.positive_control(f"正对照-v3[{k}]", planted=r, floor=max(q, RM), spread=1e-9))
nc = G.negative_control("阴性:同问卷参照分布中位(测量,非挑选)", null=RM,
                        effect=float(np.mean([abs(res[k]["same"]) for k in res]))*20,
                        null_spread=float(np.std(ref)), null_kind="同问卷无关变量参照分布")
print("\n" + "="*70)
if all(pos) and nc:
    verdict = ("**W-BETTER:数据的划法更好**" if votes > 3 else
               "W-SAME:与我的划法一样(3/5)" if votes == 3 else
               "**W-WORSE:数据的划法更差 -> 撤回 `#502b` 的解读,页面改写**")
    print(f"控制齐备 ⇒ 评判。{votes}/5 vs 3/5 -> {verdict}")
    print("⚠ 通过的 KILL 会怎样失败:两簇仍相关,共线性偏向说不清;"
          "而 PC2 的特征值只有 1.018 —— 簇本身就脆弱,票数差 1 不该被当成结构差别。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos} neg={nc})"; print(f"⚠ {verdict}")
print(G)
json.dump(dict(cluster_A=CA, cluster_B=CB, corr_AB=float(np.corrcoef(Ca,Cb)[0,1]),
               results=res, votes=votes, votes_R546=3, verdict=verdict,
               seeds=SEEDS, unchallenged=True),
          open(OUT/"data_carving_rerun.json","w"), indent=1)
print(f"\nwrote {OUT/'data_carving_rerun.json'}")

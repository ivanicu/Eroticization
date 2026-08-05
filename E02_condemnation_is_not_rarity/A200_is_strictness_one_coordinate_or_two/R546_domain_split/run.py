"""E02·A200·R546 — 「严格」是一个坐标,还是按域分块的两个

`#500` 的 NEXT。划分与理由**写在跑之前**:
  **性规范**(直接关于性行为本身,3 道):`samesex` 同性性关系 · `sxok16` · `sxok18` 未婚者性行为
  **家庭规范**(关于婚姻与生养的安排,5 道):`staytog` 离婚 · `chsuppor` 未婚生育 ·
     `okcohab` 未婚同居 · `chcohab` 同居生养 · `gayadopt` 同性收养
  行为同样分域:**性**(曾有同性接触 · 初次性交 ≤16)· **家庭**(非婚同居 · 非婚生育 · 曾离婚)

G1 ESTIMAND:每件行为 k 上,**把 `C_sex,-k` 与 `C_fam,-k` 同时放进回归**
  (并控制自己对 k 的严格立场与答题数),看**同域**与**跨域**两个系数。

WORLDS:
  W-DOMAIN 「严格」按域分块 -> 同域系数为负且超 MDE,跨域 ≈0
  W-SINGLE 它是一个坐标     -> 两个系数量级相当
  | World    | now | 同域>跨域 | 相当 |
  | W-DOMAIN | 0.5 | 0.85      | 0.10 |
  | W-SINGLE | 0.5 | 0.10      | 0.85 |

⛔ STRONGEST CONFOUND,写在跑之前:**两个子集高度相关**(严格的人到处严格)->
  共线性抬高两个系数的 se -> **偏向「两个都看不见」**,即**偏向说不清**,不是偏向任一世界。
  已在输出里报 `corr(C_sex, C_fam)`,并对每个系数各报 MDE。
⚠ 题数不同(3 vs 5)-> **每个系数各算自己的 MDE**,不共用。
CONTROLS:正对照 RULE-v3;阴性 = 同问卷参照分布中位;SHAM = 极端应答;精度 = 人层 bootstrap。
KILL(条件式,预注册):
  if 正对照全触发 and 阴性为零:
      ≥3/5 行为上「同域系数超 MDE 且 |同域| ≥ 2×|跨域|」 -> W-DOMAIN
      ≥3/5 行为上两者量级比在 [0.5, 2] 之间               -> W-SINGLE
      否则                                               -> 未决(多半是共线性)
  else: UNVERIFIED
IMPOSSIBLE:两子集共线 ⇒ 分辨力有限 · NSFG 无羞耻变量 · 只有女性卷 · 未派对抗 agent ⇒ [unchallenged]
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
SEXN = ["samesex", "sxok16", "sxok18"]
FAMN = ["staytog", "chsuppor", "okcohab", "chcohab", "gayadopt"]
NORMS = SEXN + FAMN
REVERSE = {"okcohab"}
DOMAIN = {"samesex": "sex", "sxok16": "sex", "okcohab": "fam", "chsuppor": "fam", "staytog": "fam"}

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
Csex = sum(M[n] for n in SEXN); Cfam = sum(M[n] for n in FAMN)
NANS = sum(OK[n].astype(float) for n in NORMS)
EXTR = sum(np.where(OK[n], np.isin(D[n], [1, 4]).astype(float), 0.0) for n in NORMS)
print(f"行 {N};性规范 {len(SEXN)} 道 家庭规范 {len(FAMN)} 道")
print(f"C_sex 均值={Csex.mean():.2f}  C_fam 均值={Cfam.mean():.2f}  "
      f"⛔ corr(C_sex,C_fam) = {np.corrcoef(Csex,Cfam)[0,1]:+.4f}(共线性偏向「说不清」)")

early = np.where(np.isfinite(D["agefstsx"]) & (D["agefstsx"] >= 5) & (D["agefstsx"] <= 60),
                 (D["agefstsx"] <= 16).astype(float), np.nan)
BEH = {"samesex": np.where(D["samesexany"] == 1, 1.0, np.where(D["samesexany"] == 5, 0.0, np.nan)),
       "sxok16": early,
       "okcohab": np.where(np.isfinite(D["nonmarr"]) & (D["nonmarr"] < 90), (D["nonmarr"] > 0)*1.0, np.nan),
       "chsuppor": np.where(np.isfinite(D["cebow"]) & (D["cebow"] < 90), (D["cebow"] > 0)*1.0, np.nan),
       "staytog": np.where((D["evrmarry"] == 1) & np.isfinite(D["prevhusb"]) & (D["prevhusb"] < 90),
                           (D["prevhusb"] > 0)*1.0, np.nan)}

def coefs(k, idx=None, sham=False):
    cs = Csex - (M[k] if k in SEXN else 0.0); cf = Cfam - (M[k] if k in FAMN else 0.0)
    if sham: cs = cf = EXTR - M[k]
    own, y, na = STRICT[k], BEH[k], NANS
    m = np.isfinite(y) & np.isfinite(own)
    if idx is not None: y, own, cs, cf, na, m = y[idx], own[idx], cs[idx], cf[idx], na[idx], m[idx]
    if m.sum() < 300: return np.nan, np.nan
    X = np.c_[np.ones(m.sum()), own[m], cs[m], cf[m], na[m]] if not sham else \
        np.c_[np.ones(m.sum()), own[m], cs[m], na[m]]
    b = np.linalg.lstsq(X, y[m], rcond=None)[0]
    return (float(b[2]), float(b[3])) if not sham else (float(b[2]), np.nan)

print("\n=== 主:同域 vs 跨域(同时进回归)===")
res, votes = {}, {"dom": 0, "sing": 0}
for k in BEH:
    bs_, bf_ = coefs(k); dom = DOMAIN[k]
    same, cross = (bs_, bf_) if dom == "sex" else (bf_, bs_)
    B = {"s": [], "c": []}
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(300):
            a, b2 = coefs(k, rng.integers(0, N, N))
            if np.isfinite(a):
                s_, c_ = (a, b2) if dom == "sex" else (b2, a)
                B["s"].append(s_); B["c"].append(c_)
    ms, mc = 2.8*np.std(B["s"]), 2.8*np.std(B["c"])
    ratio = abs(same)/max(abs(cross), 1e-9)
    res[k] = dict(domain=dom, same=same, cross=cross, MDE_same=float(ms), MDE_cross=float(mc),
                  ratio=float(ratio))
    if abs(same) > ms and ratio >= 2.0: votes["dom"] += 1
    if 0.5 <= ratio <= 2.0: votes["sing"] += 1
    print(f"  {k:9s}[{dom}] 同域={same:+.5f}(MDE {ms:.5f})  跨域={cross:+.5f}(MDE {mc:.5f})  "
          f"比值={ratio:.2f}")
print(f"\n同域主导票 {votes['dom']}/5;量级相当票 {votes['sing']}/5")

G = Gate("「严格」是一个坐标,还是按域分块的两个?")
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
    if votes["dom"] >= 3: verdict = f"同域主导 {votes['dom']}/5 -> **W-DOMAIN:「严格」按域分块**"
    elif votes["sing"] >= 3: verdict = f"量级相当 {votes['sing']}/5 -> **W-SINGLE:它是一个坐标**"
    else: verdict = f"同域 {votes['dom']}/5、相当 {votes['sing']}/5 -> **未决(多半是共线性)**"
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:两子集高度相关,共线性抬高两个 se ->"
          "**偏向说不清**;而域的划分是我写的,换一种划法可能换一个答案。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos} neg={nc})"; print(f"⚠ {verdict}")
print(G)
json.dump(dict(sex_norms=SEXN, fam_norms=FAMN, corr_sex_fam=float(np.corrcoef(Csex,Cfam)[0,1]),
               results=res, votes=votes, verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT/"domain_split.json","w"), indent=1)
print(f"\nwrote {OUT/'domain_split.json'}")

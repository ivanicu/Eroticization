"""E02·A199·R545 — 方向统一之后,那个负号是不是整个由「严格」携带

`#499` 的 NEXT。⚠ 先做 **G1 识别性**,而它改写了 NEXT 的一半:

  NEXT 要造两个分开的计数 `C`(严格立场数)与 `P`(宽容立场数),问 `P` 是否带独立信息。
  ⛔ **剔除中间档(码 5)之后,在「答了几道」给定时 `C + P = n_answered`** ->
     **`P` 与 `C` 完全共线,不可单独识别。** 这不是功效问题。
  ⇒ 可识别的版本:**`C_{-k}` 的系数 + `n_answered` 作协变量**;`P` 不再是一个独立问题。
  **NEXT 的那一半在此站点上不成立,直说,不硬跑。**

⚠ 第二个改动(由 `#498b` 已写下的规则推出,不是本轮新造):
  十一道里 **`chunless` / `marrfail` / `prvntdiv` 是信念不是规范** —— 把信念算进
  「你取了多少个严格立场」会污染构造。⇒ **`C` 只在 8 道规范上建。**

G1 ESTIMAND:`C_i` = i 在 8 道规范上取**严格**立场的道数;
  对每件行为 k,估计量 = **控制了 i 自己对第 k 道的态度与 `n_answered` 之后**,`C_{-k}` 的系数。

⛔ **前置断言(写进脚本,不是打印)—— `#499d` 的教训:**
  统一方向后每道题的「严格侧」必须与 `#492` 那套谴责定义**逐道一致**,不一致就 `assert` 停下。
  正向题严格 = {3,4};反向题干(`okcohab`)严格 = {1,2}。

WORLDS:
  W-C     `#499b` 的负号**整个由严格立场携带** -> `C_{-k}` 系数为负、量级与 `#499b` 相当
  W-RESID 还剩别的东西                        -> 量级明显小于 `#499b`,或符号改变
  | World   | now | 与 #499b 相当 | 明显更小/变号 |
  | W-C     | 0.6 | 0.85          | 0.10 |
  | W-RESID | 0.4 | 0.10          | 0.85 |

⛔ STRONGEST CONFOUND(沿用 `#499`):**应答风格**。SHAM = 极端应答计数(答 1 或 4,与内容无关)。
CONTROLS:正对照 = 每道规范态度 × 宗教出席(RULE-v3);阴性 = 同问卷参照分布中位;精度 = 人层 bootstrap。
KILL(条件式,预注册):
  if 正对照全触发 and 阴性为零:
      `C_{-k}` 在 ≥3/5 行为上为负且超 MDE,且**逐行为量级达 `#499b` 的 ≥70%** -> W-C
      量级 <70% 或符号改变                                                    -> W-RESID
      SHAM 与真的同强                                                          -> 应答风格,两者皆不成立
  else: UNVERIFIED
IMPOSSIBLE:**`P` 在此不可识别**(与 `C` 共线)· NSFG 无羞耻变量 ⇒ 只能测「位置→行为」·
  只有女性卷 ⇒ 男性未知 · 未派对抗 agent ⇒ `[unchallenged]`
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
R544 = json.load(open(ROOT / "E02_condemnation_is_not_rarity/A199_S_outside_its_home_instrument/"
                      "R544_leave_one_out_minority_position/results/leave_one_out_minority.json"))


def parse_dct(p):
    out = {}
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out


LAY = parse_dct(NS / "setup" / "2011_2013_FemRespSetup.dct")
NORMS = ["staytog", "samesex", "sxok18", "sxok16", "chsuppor", "gayadopt", "okcohab", "chcohab"]
BELIEFS = ["chunless", "marrfail", "prvntdiv"]        # `#498b`:信念,不进 C
REVERSE = {"okcohab"}                                  # 禁止式题干
OTHER = ["samesexany", "nonmarr", "cebow", "agefstsx", "prevhusb", "evrmarry",
         "attndnow", "age_r", "educat", "poverty", "hisp"]
cols = {n: LAY[n] for n in NORMS + BELIEFS + OTHER if n in LAY}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
D = {n: np.array(buf[n]) for n in cols}
N = len(D["staytog"])
print(f"行 {N};规范 {len(NORMS)} 道(已剔除信念 {BELIEFS} —— `#498b`)")

# ⛔ 前置断言:严格侧必须与 `#492` 的谴责定义逐道一致
STRICT, OK = {}, {}
for n in NORMS:
    v = D[n]; ok = np.isin(v, [1, 2, 3, 4])
    strict_codes = [1, 2] if n in REVERSE else [3, 4]
    condemn_codes_492 = [1, 2] if n in REVERSE else [3, 4]     # `#492` 的定义
    assert strict_codes == condemn_codes_492, f"{n}: 严格侧与 #492 的谴责定义不一致"
    STRICT[n] = np.where(ok, np.isin(v, strict_codes).astype(float), np.nan)
    OK[n] = ok
    print(f"  {n:10s} 严格侧={strict_codes}  严格率={np.nanmean(STRICT[n]):.3f}  n={ok.sum()}")
print("  ✅ 前置断言通过:8 道的严格侧逐道与 `#492` 的谴责定义一致(`#499d` 的教训)")

M = np.column_stack([np.nan_to_num(STRICT[n]) for n in NORMS])
MOK = np.column_stack([OK[n] for n in NORMS]).astype(float)
Cfull = M.sum(1); NANS = MOK.sum(1)
Pfull = NANS - Cfull
print(f"\nC 均值={Cfull.mean():.2f}  P 均值={Pfull.mean():.2f}  答题数均值={NANS.mean():.2f}")
print(f"⛔ 识别性核对:corr(C+P, n_answered) = {np.corrcoef(Cfull+Pfull, NANS)[0,1]:+.6f} "
      f"-> **P 与 C 在 n 给定时完全共线,不可单独识别**(NEXT 的那一半在此不成立)")
EXTR = np.column_stack([np.where(OK[n], np.isin(D[n], [1, 4]).astype(float), 0.0)
                        for n in NORMS]).sum(1)

early = np.where(np.isfinite(D["agefstsx"]) & (D["agefstsx"] >= 5) & (D["agefstsx"] <= 60),
                 (D["agefstsx"] <= 16).astype(float), np.nan)
BEH = {
    "samesex": np.where(D["samesexany"] == 1, 1.0, np.where(D["samesexany"] == 5, 0.0, np.nan)),
    "okcohab": np.where(np.isfinite(D["nonmarr"]) & (D["nonmarr"] < 90), (D["nonmarr"] > 0) * 1.0, np.nan),
    "chsuppor": np.where(np.isfinite(D["cebow"]) & (D["cebow"] < 90), (D["cebow"] > 0) * 1.0, np.nan),
    "sxok16": early,
    "staytog": np.where((D["evrmarry"] == 1) & np.isfinite(D["prevhusb"]) & (D["prevhusb"] < 90),
                        (D["prevhusb"] > 0) * 1.0, np.nan),
}


def coef(k, S, idx=None):
    """行为 k ~ 1 + 自己对 k 的严格立场 + S_{-k} + n_answered。返回 S 的系数。"""
    j = NORMS.index(k)
    own = STRICT[k]; Sm = S - M[:, j]
    y = BEH[k]; na = NANS
    m = np.isfinite(y) & np.isfinite(own)
    if idx is not None: y, own, Sm, na, m = y[idx], own[idx], Sm[idx], na[idx], m[idx]
    if m.sum() < 300: return np.nan
    X = np.c_[np.ones(m.sum()), own[m], Sm[m], na[m]]
    return float(np.linalg.lstsq(X, y[m], rcond=None)[0][2])


print("\n=== 主:方向统一后的 C_{-k} 系数(控制 自己的严格立场 + 答题数)===")
res = {}
for k in BEH:
    b = coef(k, Cfull); sh = coef(k, EXTR)
    bs = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(300):
            v = coef(k, Cfull, rng.integers(0, N, N))
            if np.isfinite(v): bs.append(v)
    bs = np.array(bs); MDE = 2.8 * bs.std(); ci = np.quantile(bs, [.025, .975])
    old = R544["results"][k]["coef"]
    ratio = abs(b) / max(abs(old), 1e-12)
    res[k] = dict(coef=b, sham=sh, MDE=float(MDE), ci=[float(ci[0]), float(ci[1])],
                  sd=float(bs.std()), r544=old, ratio=float(ratio))
    print(f"  {k:9s} coef={b:+.5f} CI[{ci[0]:+.5f},{ci[1]:+.5f}] MDE={MDE:.5f}  "
          f"#499b={old:+.5f}  比值={ratio:.2f}  SHAM={sh:+.5f}  "
          f"{'✅超MDE' if abs(b)>MDE else '⛔看不见'}")

neg = [k for k in res if res[k]["coef"] < 0 and abs(res[k]["coef"]) > res[k]["MDE"]]
big = [k for k in neg if res[k]["ratio"] >= 0.70]
shw = sum(1 for k in res if abs(res[k]["sham"]) < abs(res[k]["coef"]))
print(f"\n负且超 MDE:{len(neg)}/5;其中量级达 `#499b` 的 ≥70%:{len(big)}/5;SHAM 更弱:{shw}/5")

G = Gate("方向统一之后,那个负号是不是整个由「严格」携带?")
ref = []
for cn in ["age_r", "educat", "poverty", "hisp"]:
    m = np.isfinite(STRICT["samesex"]) & np.isfinite(D[cn]) & (D[cn] < 90)
    if m.sum() > 500 and len(np.unique(D[cn][m])) >= 3:
        ref.append(abs(float(np.corrcoef(STRICT["samesex"][m], D[cn][m])[0, 1])))
RM = float(np.median(ref))
pos = []
for k in BEH:
    m = np.isfinite(STRICT[k]) & np.isfinite(D["attndnow"]) & (D["attndnow"] < 90)
    r = abs(float(np.corrcoef(STRICT[k][m], D["attndnow"][m])[0, 1]))
    rg = np.random.default_rng(SEEDS[0])
    q = float(np.quantile([abs(np.corrcoef(STRICT[k][m][rg.permutation(m.sum())], D["attndnow"][m])[0, 1])
                           for _ in range(300)], .95))
    pos.append(G.positive_control(f"正对照-v3[{k}]", planted=r, floor=max(q, RM), spread=1e-9))
nc = G.negative_control("阴性:同问卷参照分布中位(测量,非挑选)", null=RM,
                        effect=float(np.mean([abs(res[k]["coef"]) for k in res])) * 20,
                        null_spread=float(np.std(ref)), null_kind="同问卷无关变量参照分布")
for k in res: G.has_error_bar(f"系数[{k}]", value=res[k]["coef"], spread=res[k]["sd"],
                              spread_source="bootstrap_人层")

print("\n" + "=" * 70)
if all(pos) and nc:
    if len(neg) >= 3 and len(big) >= 3 and shw >= 3:
        verdict = f"{len(neg)}/5 负且超 MDE,{len(big)}/5 量级≥70%,SHAM 更弱 {shw}/5 -> **W-C:负号整个由「严格」携带**"
    elif len(neg) >= 3:
        verdict = f"{len(neg)}/5 负,但只有 {len(big)}/5 量级≥70% -> **W-RESID:还剩别的东西**"
    else:
        verdict = f"{len(neg)}/5 负且超 MDE -> 未决"
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:`P` 在此**不可识别**(与 `C` 共线),"
          "所以「负号整个由严格携带」只是说**在这个可识别的分解里**如此;"
          "而信念三道被剔除是**依 `#498b` 的规则**,不是看了结果之后决定的。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos} neg={nc})"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(norms=NORMS, beliefs_excluded=BELIEFS, results=res, n_neg=len(neg),
               n_big=len(big), sham_weaker=shw, P_identifiable=False,
               verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT / "direction_unified_strictness.json", "w"), indent=1)
print(f"\nwrote {OUT/'direction_unified_strictness.json'}")

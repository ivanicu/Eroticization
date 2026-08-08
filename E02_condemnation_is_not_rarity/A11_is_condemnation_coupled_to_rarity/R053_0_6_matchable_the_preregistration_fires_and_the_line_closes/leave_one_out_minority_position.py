"""E02·A199·R544 — 「在别的规范上站在少数一边」能不能预测这一件事你做没做

`#498` 的 NEXT。⛔ 先兑现它的第一条:**NSFG 里没有羞耻/自我评价变量。**
字典 3091 个变量,搜 `shame|embarrass|guilt` 只命中 **`embarras`**
(II-4「谈避孕套会不会让 R 和伴侣尴尬」—— 工具性的,不是关于自己位置的羞耻);
`regret|stigma|tell anyone` **全 0**。⇒ **NSFG 结构上不能测羞耻。**

⭐ 但 `S` 的同构物**不需要羞耻**。E01 的 `S` 本来就是**留一构造**,
   而 NSFG 的 IH 段十一道规范(`#497a` 已验:跨度全部合格)给了它一个外部对应物:

G1 ESTIMAND(先于方法):
  `Smin_i,-k` := person i 在**除第 k 道之外的其余十道**规范上,站在**少数一边**的次数(0–10)。
  对每一件行为 k,估计量 = **在控制了 i 自己对第 k 道规范的态度之后,`Smin_-k` 的系数**。
  ⇒ 这**不是**同义反复:同一道规范上的态度已被控制,问的是**别处的少数位置**是否还带信息。
  「少数一边」由**该题的实际分布**定义(占比 < 0.5 的那一侧)—— 是定义,不是结果。

五对(态度 -> 对应行为,均在 NSFG 2011–13 女性卷内):
  `samesex`→曾有同性接触 · `okcohab`→曾非婚同居 · `chsuppor`→非婚生育 ·
  `sxok16`→初次性交 ≤16 · `staytog`→曾离婚(限已婚过者)

WORLDS:
  W-GENERAL 「整体上站在少数一边」是一个人层坐标,带跨话题的信息 -> 系数**正**且超 MDE
  W-SPECIFIC 只有对应那道态度重要                                -> 系数 **≈ 0**
  | World      | now | 系数正 | ≈0 |
  | W-GENERAL  | 0.5 | 0.85   | 0.10 |
  | W-SPECIFIC | 0.5 | 0.10   | 0.85 |

⛔ STRONGEST CONFOUND,写在跑之前:**应答风格**。
  一个对什么都选「强烈」的人,`Smin` 会高,而她回答行为题的方式也可能不同。
  ⇒ **SHAM**:用**极端应答计数**(答 1 或 4 的道数,与内容无关)造一个假 `Smin`,
     同样跑一遍。**若假的也一样正,真的就不是「少数位置」在起作用。**

CONTROLS:正对照 = 每道规范的态度 × 宗教出席(RULE-v3);阴性 = 同问卷参照分布中位;
  精度 = 人层 bootstrap;多重性 = 五个行为的族内最大值零。
KILL(条件式,预注册):
  if 正对照全触发 and 阴性为零:
      `Smin_-k` 系数在 ≥3/5 行为上为正且超 MDE **且 SHAM 明显更弱** -> W-GENERAL
      系数 ≈0(全部 < MDE)且 MDE 够小                              -> W-SPECIFIC(有内容的零)
      真假同强                                                      -> **应答风格,两个世界都不成立**
  else: UNVERIFIED
IMPOSSIBLE:NSFG 无羞耻变量 ⇒ **E01 的「位置→羞耻」在此结构上不可复现**,
  本轮只能测「位置→行为」· 无干预 ⇒ 非因果 · 未派对抗 agent ⇒ `[unchallenged]`
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


def parse_dct(p):
    out = {}
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out


LAY = parse_dct(NS / "setup" / "2011_2013_FemRespSetup.dct")
NORMS = ["staytog", "samesex", "sxok18", "sxok16", "chunless", "chsuppor",
         "gayadopt", "okcohab", "marrfail", "chcohab", "prvntdiv"]
OTHER = ["samesexany", "nonmarr", "cebow", "agefstsx", "prevhusb", "evrmarry",
         "attndnow", "age_r", "educat", "poverty", "hisp"]
cols = {n: LAY[n] for n in NORMS + OTHER if n in LAY}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
D = {n: np.array(buf[n]) for n in cols}
N = len(D["staytog"])
print(f"行 {N};规范 {len(NORMS)} 道")

# 「少数一边」:把每道题二分为 {1,2} vs {3,4},取占比 < 0.5 的那一侧为少数
MIN_SIDE, valid = {}, {}
for n in NORMS:
    v = D[n]; ok = np.isin(v, [1, 2, 3, 4])
    lowhalf = np.where(ok, np.isin(v, [1, 2]).astype(float), np.nan)
    share = np.nanmean(lowhalf)
    MIN_SIDE[n] = (share < 0.5)          # True: {1,2} 是少数
    valid[n] = ok
    print(f"  {n:10s} 答{{1,2}}占比={share:.3f} -> 少数一边 = {'{1,2}' if MIN_SIDE[n] else '{3,4}'}")

MINO = np.zeros((N, len(NORMS))); MINOK = np.zeros((N, len(NORMS)), bool)
EXTR = np.zeros((N, len(NORMS)))     # SHAM:极端应答(答 1 或 4)
for j, n in enumerate(NORMS):
    v = D[n]; ok = valid[n]
    side = np.isin(v, [1, 2]) if MIN_SIDE[n] else np.isin(v, [3, 4])
    MINO[:, j] = np.where(ok, side.astype(float), 0.0); MINOK[:, j] = ok
    EXTR[:, j] = np.where(ok, np.isin(v, [1, 4]).astype(float), 0.0)
print(f"\nSmin(全 11 道)均值={MINO.sum(1).mean():.2f}  SHAM 极端应答均值={EXTR.sum(1).mean():.2f}  "
      f"corr={np.corrcoef(MINO.sum(1), EXTR.sum(1))[0,1]:+.3f}")

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


def coef_lo(k, S, idx=None):
    """行为 k ~ 1 + 自己对第 k 道的态度(二分) + S_{-k};返回 S 的系数。"""
    j = NORMS.index(k)
    own = np.where(valid[k], MINO[:, j], np.nan)
    Sm = (S.sum(1) - S[:, j])
    y = BEH[k]
    m = np.isfinite(y) & np.isfinite(own)
    if idx is not None:
        y, own, Sm, m = y[idx], own[idx], Sm[idx], m[idx]
    if m.sum() < 300: return np.nan
    X = np.c_[np.ones(m.sum()), own[m], Sm[m]]
    return float(np.linalg.lstsq(X, y[m], rcond=None)[0][2])


print("\n=== 主:留一少数位置计数的系数(已控制自己对该道的态度)===")
res = {}
for k in BEH:
    b = coef_lo(k, MINO); sh = coef_lo(k, EXTR)
    bs = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(300):
            v = coef_lo(k, MINO, rng.integers(0, N, N))
            if np.isfinite(v): bs.append(v)
    bs = np.array(bs); MDE = 2.8 * bs.std(); ci = np.quantile(bs, [.025, .975])
    res[k] = dict(coef=b, sham=sh, MDE=float(MDE), ci=[float(ci[0]), float(ci[1])], sd=float(bs.std()))
    print(f"  {k:9s} coef={b:+.5f}  CI[{ci[0]:+.5f},{ci[1]:+.5f}]  MDE={MDE:.5f}  "
          f"SHAM(极端应答)={sh:+.5f}  {'✅超MDE' if abs(b)>MDE else '⛔看不见'}")

npos = sum(1 for k in res if res[k]["coef"] > 0 and abs(res[k]["coef"]) > res[k]["MDE"])
sham_weaker = sum(1 for k in res if abs(res[k]["sham"]) < abs(res[k]["coef"]))
print(f"\n正且超 MDE 的行为数 = {npos}/5;SHAM 更弱的行为数 = {sham_weaker}/5")

G = Gate("在别的规范上站在少数一边,能不能预测这一件事你做没做?")
ref = []
base = np.where(valid["samesex"], MINO[:, NORMS.index("samesex")], np.nan)
for cn in ["age_r", "educat", "poverty", "hisp"]:
    m = np.isfinite(base) & np.isfinite(D[cn]) & (D[cn] < 90)
    if m.sum() > 500 and len(np.unique(D[cn][m])) >= 3:
        ref.append(abs(float(np.corrcoef(base[m], D[cn][m])[0, 1])))
RM = float(np.median(ref))
pos = []
for k in BEH:
    j = NORMS.index(k); c = np.where(valid[k], MINO[:, j], np.nan)
    m = np.isfinite(c) & np.isfinite(D["attndnow"]) & (D["attndnow"] < 90)
    r = abs(float(np.corrcoef(c[m], D["attndnow"][m])[0, 1]))
    rg = np.random.default_rng(SEEDS[0])
    q = float(np.quantile([abs(np.corrcoef(c[m][rg.permutation(m.sum())], D["attndnow"][m])[0, 1])
                           for _ in range(300)], .95))
    pos.append(G.positive_control(f"正对照-v3[{k}]", planted=r, floor=max(q, RM), spread=1e-9))
nc = G.negative_control("阴性:同问卷参照分布中位(测量,非挑选)", null=RM,
                        effect=float(np.mean([abs(res[k]["coef"]) for k in res])) * 20,
                        null_spread=float(np.std(ref)), null_kind="同问卷无关变量参照分布")
for k in res: G.has_error_bar(f"系数[{k}]", value=res[k]["coef"], spread=res[k]["sd"],
                              spread_source="bootstrap_人层")

print("\n" + "=" * 70)
if all(pos) and nc:
    if npos >= 3 and sham_weaker >= 3:
        verdict = f"{npos}/5 为正且超 MDE,且 SHAM 在 {sham_weaker}/5 上更弱 -> **W-GENERAL**"
    elif npos >= 3:
        verdict = f"{npos}/5 为正,但 SHAM 只在 {sham_weaker}/5 上更弱 -> **应答风格未被排除**"
    elif all(abs(res[k]["coef"]) < res[k]["MDE"] for k in res):
        verdict = "全部 < MDE -> W-SPECIFIC(只有对应那道态度重要);须带 MDE 陈述"
    else:
        verdict = f"{npos}/5 -> 未决"
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:「少数一边」由样本自身的分布定义,"
          "换一个人群会换一条边;而 NSFG 只有女性卷,**男性上是否成立本轮不可知**。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pos} neg={nc})"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(minority_side={k: bool(v) for k, v in MIN_SIDE.items()}, results=res,
               n_pos=npos, sham_weaker=sham_weaker, verdict=verdict,
               nsfg_has_shame=False, seeds=SEEDS, unchallenged=True),
          open(OUT / "leave_one_out_minority.json", "w"), indent=1)
print(f"\nwrote {OUT/'leave_one_out_minority.json'}")

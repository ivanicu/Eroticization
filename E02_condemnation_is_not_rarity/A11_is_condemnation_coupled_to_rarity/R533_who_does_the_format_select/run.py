"""E02·A193·R533 — 换一道题,基率移动十分之一;而它挑走的是不是同一批人

`#487` 的 NEXT:两臂内分别把行为对人口学回归,**比较系数向量而不是基率**。

⚠ 先补一个 `#487` 的 NEXT 漏掉的世界(meta-separator):
  两臂不只**格式**不同,**所指对象**也不同(「X 级电影」vs「色情内容」)。
  所以 `{A 门槛平移, B 格式改变了在测谁}` 这个二分**不完整**,还有
  **W-OBJECT:两道题问的本来就是不同的东西**,而它与 B 在数据上无法分开。
  ⇒ 「不一致」**仍然不可归因**。
  ⇒ 因此本轮的可决策产出改写为:**跨题面比较安不安全** ——
     只有「不一致」能回答它,而「一致」只给一个分辨率界。

G1 ESTIMAND(先于方法):
  两臂各自把「行为阳性」对**标准化人口学**(age·sex·degree·attend·polviews)做线性概率回归,
  取**去截距的系数向量**,估计量 = 两臂系数向量的**余弦相似度**。
  **方向 = 「挑走了谁」,与基率无关;模长 = 「挑走了多少」,与基率有关。** 只判方向。

⚠ 零不是发明的:**臂内随机劈半**,算臂内余弦分布 -> 这是**实测的分辨率底**
  (`#482c`/`#485a` 的教训:测参照,别挑门槛)。

⚠ STRONGEST CONFOUND,写在跑之前:**基率差本身会拉伸 LPM 系数**。
  余弦对尺度不变,方向受保护;但若关系非线性,基率可以**旋转**向量。
  控制(SHAM,同操作减去被研究的成分):把频率臂在**与二值臂基率匹配**的切点上二分
  —— 匹配基率后,机械成分被移除,方向差若仍在,就不是基率造的。

WORLDS:
  W-SAME      两道题挑走同一批人 -> 臂间余弦落在**臂内劈半分布**之内
  W-DIFFERENT 两道题挑走不同的人 -> 臂间余弦**低于**臂内底(< q05)
  | World       | now | 低于底 | 落在底内 |
  | W-SAME      | 0.5 | 0.10   | 0.85     |
  | W-DIFFERENT | 0.5 | 0.85   | 0.10     |
  ⚠ 按上文,只有「低于底」这一格是可归因于「跨题面比较不安全」这个**决策**的;
     「落在底内」只给一个界,**不授权 W-SAME 作为机制**。

CONTROLS:
  正对照 谴责×礼拜出席,两臂各一次,门槛 = 同问卷参照分布 q95(`#485a`)
  安慰剂 谴责×星座,两臂各一次
  底     臂内随机劈半余弦(两臂各算,取更宽的那个作底)
  精度   人层 bootstrap

KILL(条件式,预注册):
  if 正对照两臂都触发 and 安慰剂两臂都为零:
      臂间余弦 < 臂内底 q05 -> **W-DIFFERENT:跨题面比较不安全,`#486b` 必须带「题面固定」**
      臂间余弦 >= 底 q05    -> 只报界;**不授权 W-SAME**
  else: UNVERIFIED

IMPOSSIBLE:格式与对象不可拆(GSS 设计)· 无干预 ⇒ 非因果 ·
  未派对抗 agent(会话约束)⇒ `[unchallenged]` · 两臂共享自报方法 ⇒ 一致方向不可授权
"""
import os, sys, pathlib, json, math, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
from lib.gates import Gate, check_columns

SEEDS = [20260805, 7, 991]
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
REF529 = json.load(open(ROOT / "E02_condemnation_is_not_rarity/"
                        "A192_does_condemnation_suppress_the_report_or_the_act/"
                        "R529_gss_gap_over_51_years/results/gss_gap_over_51_years.json"))["reference"]
T_POS = float(np.quantile([abs(x["r"]) for x in REF529], .95))
PRED = ["age", "sex", "degree", "attend", "polviews"]

cols = ["year", "zodiac", "pornlaw", "xmovie", "xmoviey"] + PRED
it = pd.read_stata(DTA, iterator=True); have = set(it.variable_labels())
d = pd.read_stata(DTA, columns=[c for c in cols if c in have], convert_categoricals=False)
check_columns(d, where="R533")
d = d[d.year.isin([2022, 2024])].copy()
d["arm"] = np.where(d.xmovie.notna(), "xmovie", np.where(d.xmoviey.notna(), "xmoviey", None))
d = d[d.arm.notna()].copy()
d["condemn"] = np.where(d.pornlaw.notna(), (d.pornlaw == 1).astype(float), np.nan)
d = d.dropna(subset=PRED).copy()
for c in PRED:                                     # 标准化,两臂**合并**标准化以保尺度一致
    d[c] = (d[c] - d[c].mean()) / d[c].std()
print(f"n={len(d)}  两臂 {d.arm.value_counts().to_dict()}")

MINC = 1.0                                          # xmoviey: 1 = never(`#487` 已核)
d["posA"] = np.where(d.arm == "xmovie", (d.xmovie == 1).astype(float), np.nan)
d["posB"] = np.where(d.arm == "xmoviey", (d.xmoviey > MINC).astype(float), np.nan)
A = d[d.arm == "xmovie"]; B = d[d.arm == "xmoviey"]
baseA = float(A.posA.mean()); baseB = float(B.posB.mean())
print(f"基率  xmovie={baseA:.4f}   xmoviey(非从不)={baseB:.4f}")

# 基率匹配的切点(SHAM):在频率臂上取分位,使阳性率 ≈ baseA
q = B.xmoviey.quantile(1 - baseA)
d["posB_matched"] = np.where(d.arm == "xmoviey", (d.xmoviey > q).astype(float), np.nan)
A = d[d.arm == "xmovie"]; B = d[d.arm == "xmoviey"]      # ⚠ 重新切片:posB_matched 是刚加的列
baseBm = float(B.posB_matched.mean())
print(f"基率匹配切点 = >{q:.1f}  -> 匹配后 xmoviey 阳性={baseBm:.4f}(目标 {baseA:.4f})")


def coefs(sub, ycol):
    y = sub[ycol].values
    X = np.c_[np.ones(len(sub)), sub[PRED].values]
    b = np.linalg.lstsq(X, y, rcond=None)[0][1:]     # 去截距
    return b


def cos(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


bA = coefs(A, "posA"); bB = coefs(B, "posB"); bBm = coefs(B, "posB_matched")
print("\n=== 标准化系数向量(线性概率)===")
print(f"{'pred':10s} {'xmovie':>9s} {'xmoviey':>9s} {'xmoviey匹配':>12s}")
for i, p in enumerate(PRED):
    print(f"{p:10s} {bA[i]:+9.4f} {bB[i]:+9.4f} {bBm[i]:+12.4f}")
print(f"{'‖b‖':10s} {np.linalg.norm(bA):9.4f} {np.linalg.norm(bB):9.4f} {np.linalg.norm(bBm):12.4f}")

c_main = cos(bA, bB); c_matched = cos(bA, bBm)
print(f"\n臂间余弦(原切点) = {c_main:+.4f}")
print(f"臂间余弦(基率匹配 SHAM) = {c_matched:+.4f}")

# ---------------------------------------------------------------- 实测底:臂内随机劈半
def within_arm_cos(sub, ycol, n=600, seed=0):
    rng = np.random.default_rng(seed); out = []
    idx = np.arange(len(sub))
    for _ in range(n):
        p = rng.permutation(idx); h = len(idx) // 2
        u = coefs(sub.iloc[p[:h]], ycol); v = coefs(sub.iloc[p[h:]], ycol)
        c = cos(u, v)
        if np.isfinite(c): out.append(c)
    return np.array(out)


floors = {}
for nm, sub, yc in [("xmovie", A, "posA"), ("xmoviey", B, "posB")]:
    cc = np.concatenate([within_arm_cos(sub, yc, 600, s) for s in SEEDS])
    floors[nm] = dict(median=float(np.median(cc)), q05=float(np.quantile(cc, .05)),
                      q25=float(np.quantile(cc, .25)))
    print(f"臂内劈半底[{nm}]: 中位={np.median(cc):+.4f}  q25={np.quantile(cc,.25):+.4f}  "
          f"q05={np.quantile(cc,.05):+.4f}")
FLOOR = min(floors[k]["q05"] for k in floors)       # 取更宽的底(更保守)
print(f"⇒ 采用的底 q05 = {FLOOR:+.4f}(两臂中更低者,更保守)")

# 臂间余弦的精度
def boot_cos(n=500, seed=0, ycolB="posB"):
    rng = np.random.default_rng(seed); out = []
    for _ in range(n):
        aa = A.iloc[rng.integers(0, len(A), len(A))]
        bb = B.iloc[rng.integers(0, len(B), len(B))]
        c = cos(coefs(aa, "posA"), coefs(bb, ycolB))
        if np.isfinite(c): out.append(c)
    return np.array(out)


bc = np.concatenate([boot_cos(500, s) for s in SEEDS])
lo_ci, hi_ci = np.quantile(bc, [.025, .975])
print(f"臂间余弦 95% CI [{lo_ci:+.4f}, {hi_ci:+.4f}]  sd={bc.std():.4f}  "
      f"seed_spread={np.std([boot_cos(500,s).mean() for s in SEEDS]):.5f}")

# ---------------------------------------------------------------- G4 规格曲线
print("\n=== G4 规格曲线(全格公布)===")
spec = []
for cut in [1, 2, 3, 4, 5, 6]:
    d["_p"] = np.where(d.arm == "xmoviey", (d.xmoviey > cut).astype(float), np.nan)
    bb2 = coefs(d[d.arm == "xmoviey"], "_p")
    spec.append(dict(kind="cut", cut=cut, cos=cos(bA, bb2),
                     base=float(d[d.arm == "xmoviey"]._p.mean())))
    print(f"  切点 >{cut}  基率={spec[-1]['base']:.3f}  余弦={spec[-1]['cos']:+.4f}")
for nm, ps in [("no_polviews", [p for p in PRED if p != "polviews"]),
               ("no_attend", [p for p in PRED if p != "attend"]),
               ("age_sex_only", ["age", "sex"])]:
    g = PRED[:]; PRED[:] = ps
    spec.append(dict(kind="predset", name=nm, cos=cos(coefs(A, "posA"), coefs(B, "posB"))))
    print(f"  预测集 {nm:14s} 余弦={spec[-1]['cos']:+.4f}")
    PRED[:] = g
spec.append(dict(kind="sham", name="base_rate_matched", cos=c_matched))
print(f"  SHAM 基率匹配        余弦={c_matched:+.4f}")
cs = [s["cos"] for s in spec if np.isfinite(s["cos"])]
below = sum(1 for c in cs if c < FLOOR)
print(f"\nspec: {len(cs)} 格,**低于底的格数 = {below}/{len(cs)}**;"
      f" 余弦范围 [{min(cs):+.4f}, {max(cs):+.4f}]")

# ---------------------------------------------------------------- 控制
G = Gate("换一道题,挑走的是不是同一批人?(GSS 2022/2024 随机题面臂)")
oks = {}
for nm, arm in [("xmovie", A), ("xmoviey", B)]:
    s = arm.dropna(subset=["attend", "condemn"])
    r = float(np.corrcoef(s.condemn, s.attend)[0, 1])
    oks[f"pos_{nm}"] = G.positive_control(f"正对照[{nm}]:谴责×礼拜出席(门槛=参照 q95)",
                                          planted=abs(r), floor=T_POS, spread=1e-9)
    z = arm.dropna(subset=["zodiac", "condemn"])
    zr = float(np.corrcoef(z.condemn, z.zodiac)[0, 1])
    rngz = np.random.default_rng(SEEDS[0])
    zn = np.array([np.corrcoef(z.condemn.values[rngz.permutation(len(z))], z.zodiac)[0, 1]
                   for _ in range(300)])
    oks[f"plac_{nm}"] = G.negative_control(f"安慰剂[{nm}]:谴责×星座", null=zr, effect=r,
                                           null_spread=float(zn.std()), null_kind="个体层标签置换")
G.has_error_bar("臂间余弦", value=c_main, spread=float(bc.std()), spread_source="bootstrap_人层")
G.asserted("底不是发明的", True, "臂内随机劈半 3600 次,取两臂更保守的 q05", kind="control")

all_pos = oks["pos_xmovie"] and oks["pos_xmoviey"]
all_plac = oks["plac_xmovie"] and oks["plac_xmoviey"]
print("\n" + "=" * 70)
if all_pos and all_plac:
    if c_main < FLOOR:
        verdict = (f"臂间余弦 {c_main:+.4f} < 臂内底 q05 {FLOOR:+.4f} -> "
                   f"**W-DIFFERENT:两道题挑走的不是同一批人,跨题面比较不安全**")
    else:
        verdict = (f"臂间余弦 {c_main:+.4f} >= 臂内底 q05 {FLOOR:+.4f} -> "
                   f"在本设计分辨率内**分不开**;只报界,**不授权 W-SAME**")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:格式与所指对象在两臂间同时变了,"
          "任何方向差**不可归因于任一个**;而余弦只看方向,两臂在**模长**上的差"
          "(挑走多少)本设计不判。")
else:
    verdict = f"UNVERIFIED —— 控制未齐 {oks}"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(coef_xmovie=dict(zip(PRED, bA.tolist())),
               coef_xmoviey=dict(zip(PRED, bB.tolist())),
               coef_xmoviey_matched=dict(zip(PRED, bBm.tolist())),
               norm=dict(xmovie=float(np.linalg.norm(bA)), xmoviey=float(np.linalg.norm(bB))),
               cos_main=c_main, cos_matched=c_matched, ci=[float(lo_ci), float(hi_ci)],
               floors=floors, floor_used=FLOOR, base=dict(xmovie=baseA, xmoviey=baseB,
                                                          xmoviey_matched=baseBm),
               spec=[{k: (None if isinstance(v, float) and not np.isfinite(v) else v)
                      for k, v in s.items()} for s in spec],
               cells_below_floor=below, controls=oks, verdict=verdict,
               seeds=SEEDS, unchallenged=True),
          open(OUT / "who_does_the_format_select.json", "w"), indent=1)
print(f"\nwrote {OUT/'who_does_the_format_select.json'}")

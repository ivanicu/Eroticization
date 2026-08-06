"""E02·A13·R643 —— 那把「唯一的尺子」,是性道德,还是宗教?(并顺手量了那句「同一批人」)

BASIN RULE:`#600`–`#606` 连续七轮全是工具与归档。本轮**故意攻击本项目唯一还站着的旗舰主张**,
且它的**阳性结果是我不希望的**:若成立,「人层只有一把性道德的尺子」降级为「宗教承诺的再描述」。

⚠ 第一版设计死了,而它死得有用:我要求性题与参照题**同时**非缺失,联合 n **34 年全部为 0** ——
  GSS 是**分票问卷**。`R579` 用的是**逐对**可用子样本,所以本轮沿用逐对设计。
  ⚠ 而这暴露了页面上一句话:「against 0.093 for pairs ... **answered by the same people**」——
     `homosex×premarsx` 与 `cappun×grass` 的受访者**零重叠,20 年全部**。本轮把它量出来。

G1 ESTIMAND(先于方法):**超出量** `E = median|ρ|(性题对) − median|ρ|(参照对)`,同年内、逐对子样本。
  关心的量是 **`retain = E_partial / E_raw`** —— 偏掉宗教之后超出量保住了多少。
  ⚠ 不是「性题相关掉了多少」:偏相关会**同时**抽掉参照对的方差,只看一端 = 「两个有界分数之差」。

SCOPE:population = GSS 受访者(逐年);instrument = **GSS 问卷,同一批访员**(硬规则②:全程 route через GSS);
  baseline = 同问卷非性态度题对;regime = 每对子样本 n ≥ 200 的年份。

WORLDS:A 它是道德坐标(retain 高)· B 它**就是**宗教(retain ≈ 0)· C 宗教是若干输入之一(中间)
CONTROLS:正对照 = 偏掉一道性题本身,其余三题的超出量必须塌(且 g=0 时不通过)·
  安慰剂 = 偏掉 `zodiac`,必须不动 ·
  **offset_control**(「这个零该是零吗?」**不该**)= 边际相关匹配的假变量,null 种类 = **机械衰减**
KILL(条件式,预注册):控制齐备才评判 —— 校正后 retain < 0.25 -> W-B;> 0.70 -> W-A;之间 -> W-C
IMPOSSIBLE:无干预 ⇒ 非因果 · 宗教与性道德一生里互为因果,横截面**结构上无法定序** ·
  单一仪器 ⇒ cross-dataset 做不到(需要另一份同时问四道性题与宗教强度的问卷)· `[unchallenged]`
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate

SEEDS = [20260806, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
MINN = 200

SEX = {"homosex": [1,2,3,4], "premarsx": [1,2,3,4], "xmarsex": [1,2,3,4], "teensex": [1,2,3,4]}
REF = {"cappun": [1,2], "grass": [1,2], "letdie1": [1,2], "suicide1": [1,2], "fepol": [1,2]}
RELB = {"仅出席": ["attend"], "出席+强度+原教旨": ["attend","reliten","fund"],
        "全五项": ["attend","reliten","fund","pray","bible"]}
ALLREL = ["attend","reliten","fund","pray","bible"]

cols = ["year","zodiac"] + list(SEX) + list(REF) + ALLREL
g = pd.read_stata(DTA, columns=cols, convert_categoricals=False)
print(f"载入 n={len(g)}  年份 {int(g.year.min())}-{int(g.year.max())}")


def resid(v, X):
    r = rankdata(v).astype(float)
    if X is None or X.shape[1] == 0: return r
    Z = np.column_stack([rankdata(X[:, j]) for j in range(X.shape[1])]).astype(float)
    Z = np.column_stack([np.ones(len(r)), (Z - Z.mean(0)) / (Z.std(0) + 1e-12)])
    b, *_ = np.linalg.lstsq(Z, r, rcond=None)
    return r - Z @ b


# ⚠ 第一版死在这里,而安慰剂抓住了它:偏相关时子样本随 part 收缩,
#   进入中位数的**题对集合本身在变** —— 分母在动,于是偏掉星座也抹掉 82% 的超出量。
#   修法:**行钉死**。原始、偏宗教、安慰剂、offset 全部跑在**同一批行**上 —— FIXED 给出那批行。
FIXED = ["attend","reliten","fund","pray","bible","zodiac"]

def rho_pair(d, a, b, va, vb, part, lock=True):
    """一对题的 |ρ|。lock=True 时子样本 = 该对可用 ∩ FIXED 全部非缺失(与 part 无关)。"""
    m = d[a].isin(va) & d[b].isin(vb)
    if lock:
        for p in FIXED:
            if p in d.columns: m &= d[p].notna()
    else:
        for p in part: m &= d[p].notna()
    # ⚠ 被偏掉的变量若在锁定行里有缺失,lstsq 会炸(实测 SVD did not converge)——
    #   所以 part 也必须非缺失,而这一步**对所有规格一视同仁**,不是给某一格开的后门。
    for pcol in part:
        if pcol in d.columns: m &= d[pcol].notna()
    if m.sum() < MINN: return np.nan, int(m.sum())
    X = d.loc[m, part].values if part else None
    x, y = resid(d[a][m].values, X), resid(d[b][m].values, X)
    if x.std() < 1e-12 or y.std() < 1e-12: return np.nan, int(m.sum())
    return abs(float(np.corrcoef(rankdata(x), rankdata(y))[0,1])), int(m.sum())


def med(d, spec, part, lock=True):
    out, ns = [], []
    for a, b in itertools.combinations(spec, 2):
        r, n = rho_pair(d, a, b, spec[a], spec[b], part, lock)
        if np.isfinite(r): out.append(r); ns.append(n)
    return (float(np.median(out)), len(out), int(np.median(ns))) if out else (np.nan, 0, 0)


# ── ① 先把那句「同一批人」量出来 ────────────────────────────────────
print("\n=== ① 页面说参照对是「同一批人」回答的 —— 逐对之对量重叠 ===")
ov = []
for y in sorted(g.year.unique()):
    d = g[g.year == y]
    for (a,b) in itertools.combinations(SEX, 2):
        ms = d[a].isin(SEX[a]) & d[b].isin(SEX[b])
        if ms.sum() < MINN: continue
        for (c,e) in itertools.combinations(REF, 2):
            mr = d[c].isin(REF[c]) & d[e].isin(REF[e])
            if mr.sum() < MINN: continue
            ov.append(dict(year=int(y), sexpair=f"{a}×{b}", refpair=f"{c}×{e}",
                           n_sex=int(ms.sum()), n_ref=int(mr.sum()), overlap=int((ms & mr).sum())))
OV = pd.DataFrame(ov)
share0 = float((OV.overlap == 0).mean())
print(f"  {len(OV)} 个「性题对 × 参照对」组合 · **重叠为 0 的占 {share0*100:.1f}%** · "
      f"重叠中位 {int(OV.overlap.median())} · 最大 {int(OV.overlap.max())}")
print(f"  ⇒ 那句「same people」对 {(1-share0)*100:.1f}% 的组合成立,对 {share0*100:.1f}% **不成立**")

# ── ② 主网格 ─────────────────────────────────────────────────────
rows = []
for y in sorted(g.year.unique()):
    d = g[g.year == y]
    s0, ks, ns0 = med(d, SEX, []); r0, kr, nr0 = med(d, REF, [])
    if not (np.isfinite(s0) and np.isfinite(r0)): continue
    E0 = s0 - r0
    if E0 <= 1e-9: 
        rows.append(dict(year=int(y), block="—", E_raw=E0, E_par=np.nan, retain=np.nan,
                         sex_raw=s0, ref_raw=r0, k_sex=ks, k_ref=kr, n_sex=ns0, n_ref=nr0)); continue
    for bname, blk in RELB.items():
        blk = [c for c in blk if d[c].notna().sum() > MINN]
        if not blk: continue
        s1, _, ns1 = med(d, SEX, blk); r1, _, nr1 = med(d, REF, blk)
        if not (np.isfinite(s1) and np.isfinite(r1)): continue
        rows.append(dict(year=int(y), block=bname, blk=",".join(blk), E_raw=E0, E_par=s1-r1,
                         retain=(s1-r1)/E0, sex_raw=s0, sex_par=s1, ref_raw=r0, ref_par=r1,
                         k_sex=ks, k_ref=kr, n_sex=ns1, n_ref=nr1))
G_ = pd.DataFrame(rows)
GG = G_[G_.retain.notna()]
print(f"\n=== ② G3 整张网格:{len(GG)} 格(年份 {GG.year.nunique()} × 宗教块 {GG.block.nunique()})===")
print(GG.groupby("block").agg(格=("retain","size"), 年=("year","nunique"),
      E_raw中位=("E_raw","median"), E_par中位=("E_par","median"),
      retain中位=("retain","median"), retain最小=("retain","min"), retain最大=("retain","max")).round(4).to_string())
RET = float(GG.retain.median())
print(f"\n**整张网格 retain 中位 = {RET:.4f}** · 四分位 [{GG.retain.quantile(.25):.4f}, {GG.retain.quantile(.75):.4f}]")
print(f"  **retain < 0.25 的格 {int((GG.retain<0.25).sum())}/{len(GG)} · > 0.70 的格 {int((GG.retain>0.70).sum())}/{len(GG)}**")
print(f"  性题中位 {GG.sex_raw.median():.4f} -> 偏掉宗教后 {GG.sex_par.median():.4f} · "
      f"参照中位 {GG.ref_raw.median():.4f} -> {GG.ref_par.median():.4f}")
print("  ⚠ 参照也在下降 —— 这正是为什么估计量是**超出量**而不是相关本身")

# ── ③ 正对照 ─────────────────────────────────────────────────────
print("\n=== ③ 正对照:偏掉四道性题里的一道,其余三题的超出量必须塌 ===")
pos = []
for y in sorted(g.year.unique()):
    d = g[g.year == y]
    for drop in SEX:
        rest = {k: v for k, v in SEX.items() if k != drop}
        s0, k0, _ = med(d, rest, []); r0, _, _ = med(d, REF, [])
        if not (np.isfinite(s0) and np.isfinite(r0)) or k0 < 3 or s0 - r0 <= 1e-9: continue
        # 正对照偏的是性题本身(不在 FIXED 里)-> 原始也必须限制在同一批行,否则又是分母在动
        s0b, _, _ = med(d, rest, [], lock=True); r0b, _, _ = med(d, REF, [], lock=True)
        s1, _, _ = med(d, rest, [drop]); r1, _, _ = med(d, REF, [drop])
        if np.isfinite(s0b) and np.isfinite(r0b) and s0b-r0b > 1e-9: s0, r0 = s0b, r0b
        if np.isfinite(s1) and np.isfinite(r1): pos.append((s1-r1)/(s0-r0))
POS = float(np.median(pos)) if pos else np.nan
print(f"  偏掉一道性题 -> retain 中位 = **{POS:.4f}**({len(pos)} 格)")
print(f"  g=0(什么都不偏)-> retain ≡ 1.0000 ⇒ 判据「retain 明显 < 1」在 g=0 **不通过** ✅ 正对照可以失败")

# ── ④ 安慰剂 ─────────────────────────────────────────────────────
pla = []
for y in sorted(g.year.unique()):
    d = g[g.year == y]
    if d.zodiac.notna().sum() < MINN: continue
    s0,_,_ = med(d, SEX, []); r0,_,_ = med(d, REF, [])
    if not (np.isfinite(s0) and np.isfinite(r0)) or s0-r0 <= 1e-9: continue
    s1,_,_ = med(d, SEX, ["zodiac"]); r1,_,_ = med(d, REF, ["zodiac"])
    if np.isfinite(s1) and np.isfinite(r1): pla.append((s1-r1)/(s0-r0))
PLA = float(np.median(pla)) if pla else np.nan
print(f"\n=== ④ 安慰剂:偏掉星座 -> retain 中位 = **{PLA:.4f}**({len(pla)} 格),必须 ≈ 1 ===")

# ── ⑤ offset:边际相关匹配的假变量 ──────────────────────────────────
print("\n=== ⑤ offset(不是 negative):偏掉任何与两端都相关的变量都有**机械衰减** ===")
off = []
for y in sorted(g.year.unique()):
    d = g[g.year == y].copy()
    blk = [c for c in RELB["出席+强度+原教旨"] if d[c].notna().sum() > MINN]
    if not blk: continue
    s0,_,_ = med(d, SEX, []); r0,_,_ = med(d, REF, [])
    if not (np.isfinite(s0) and np.isfinite(r0)) or s0-r0 <= 1e-9: continue
    # 目标边际:宗教块 PC1 与每道题的 |ρ|(在该题自己的子样本上)
    mb = d[blk].notna().all(axis=1)
    Z = np.column_stack([rankdata(d.loc[mb, c].values) for c in blk]).astype(float)
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-12)
    pc = np.zeros(len(d)); pc[:] = np.nan
    pc[mb.values] = np.linalg.svd(Z, full_matrices=False)[0][:, 0]
    tgt = {}
    for c in list(SEX) + list(REF):
        m = d[c].isin((SEX | REF)[c]) & np.isfinite(pc)
        tgt[c] = abs(float(np.corrcoef(rankdata(d[c][m]), rankdata(pc[m.values]))[0,1])) if m.sum() > MINN else 0.0
    sub = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd + int(y))
        # 每道题各造一份「只与它相关 tgt[c]」的噪声,取均值 -> 边际匹配,**共同结构被打散**
        parts = []
        for c in list(SEX) + list(REF):
            m = d[c].isin((SEX | REF)[c]).values
            v = np.zeros(len(d))
            if m.sum() > 2:
                yy = rankdata(d[c].values[m]).astype(float); yy = (yy - yy.mean())/(yy.std()+1e-12)
                v[m] = tgt[c]*yy + np.sqrt(max(0., 1-tgt[c]**2))*rng.standard_normal(m.sum())
            v[~m] = rng.standard_normal((~m).sum())
            parts.append(v)
        d["_z"] = np.mean(parts, axis=0)
        s1,_,_ = med(d, SEX, ["_z"]); r1,_,_ = med(d, REF, ["_z"])
        if np.isfinite(s1) and np.isfinite(r1): sub.append((s1-r1)/(s0-r0))
    if sub: off.append(float(np.mean(sub)))
OFF = float(np.median(off)) if off else np.nan
print(f"  边际匹配假变量 -> retain 中位 = **{OFF:.4f}**({len(off)} 年 × {len(SEEDS)} 种子)")
print(f"  ⇒ **机械衰减本身带走 {(1-OFF)*100:.1f}% 的超出量,而这与宗教无关。**")

# ── ⑥ 条件式 KILL ────────────────────────────────────────────────
G = Gate("那把唯一的尺子,是性道德还是宗教?(GSS,人层,同年内,逐对子样本)")
spread = float(GG.retain.std())
pos_ok = G.positive_control("正对照:偏掉一道性题本身", planted=1.0-POS, floor=0.10, spread=spread)
pla_ok = G.negative_control("安慰剂:偏掉星座", null=abs(1.0-PLA), effect=1.0-RET,
                            null_spread=spread, null_kind="与两端都无关的个体层标签")
G.offset_control("offset:边际相关匹配的假变量", effect=1.0-RET, offset=1.0-OFF, spread=spread,
                 null_kind="偏掉任何与两端都相关的变量都有的机械衰减,与宗教无关")
bs = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(400):
        bs.append(float(np.median(GG.retain.values[rng.integers(0, len(GG), len(GG))])))
G.has_error_bar("retain(整张网格)", value=RET, spread=float(np.std(bs)), spread_source="bootstrap_人层")
CORR = RET/OFF if OFF and OFF > 1e-9 else np.nan
print(f"\n**校正后 retain = {RET:.4f} / {OFF:.4f} = {CORR:.4f}**")
if pos_ok and pla_ok:
    verdict = ("W-B:那把尺子**就是宗教**,旗舰主张降级" if CORR < 0.25 else
               "W-A:**宗教不是它** —— 偏掉宗教后超出量基本保住" if CORR > 0.70 else
               "W-C:**宗教是若干输入之一** —— 超出量明显下降但远离零")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
    print("⚠ 通过的 KILL 会怎样失败:偏相关**不定序** —— 一个人先虔诚后严格,还是先严格后虔诚,"
          "这份横截面数据结构上答不了。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · 安慰剂 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(overlap_cells=len(OV), overlap_zero_share=share0,
               overlap_median=float(OV.overlap.median()), overlap_max=int(OV.overlap.max()),
               grid=GG.to_dict("records"), retain_median=RET,
               retain_iqr=[float(GG.retain.quantile(.25)), float(GG.retain.quantile(.75))],
               cells=len(GG), below_025=int((GG.retain<0.25).sum()), above_070=int((GG.retain>0.70).sum()),
               positive_retain=POS, placebo_retain=PLA, offset_retain=OFF, corrected_retain=CORR,
               verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT/"is_the_ruler_religion.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'is_the_ruler_religion.json'}")

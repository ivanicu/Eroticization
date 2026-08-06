"""E02·A12·R653 —— 一个行为、一条尺子、三圈「对谁」:谴责是关于行为,还是关于对谁做?

`#616` 的 NEXT。**这是本项目第一次拿到「同一行为 × 同一量表 × 三个社会距离」的设计。**

⚠ **硬规则①当场抓到最强混淆,而它必须在同一轮里被控制:**
   `SCCS781`(本地社群)**没有 `1=Valued` 这一档**,取值只有 2–4;而 `782/783` 有 1–4。
   ⇒ **本地那一项在结构上就够不到最低分,这会机械地造出我预期的那个梯度。**
   ⇒ 控制:把 `782/783` 的 1 与 2 合并,让三项都变成三档,**两个口径都报**。

⚠ **而这一轮自带一个跨编码团队的复制(硬规则④),在同一份数据内部:**
   **Ross 1983**(`781/782/783`,1=Valued…4=Disapproved)与
   **Lang 1998**(`1768/1769/1770`,1=rejected…3=appreciated,**极性相反**)——
   **两个团队、两个年代、同一个三距离结构。**
   `SCCS1770` 用两位码(十位=程度,个位=限定),`88=无接触` 必须剔除。

INSTRUMENT(硬规则②):SCCS/D-PLACE,186 个社会,**民族志编码**。
  主威胁 = **同一批编码者的自身耦合**,`#530` 已量为 **+0.187 ± 0.152**。

G1 ESTIMAND(先于方法):`ρ(社会距离, 可接受度)`,**观测单位 = (社会 × 距离)**,
  聚类单位 = 社会,并按 `#529` 的 **10°×10° 经纬网格块** 做 block bootstrap(Galton 问题的部分修正)。
  **同时报三点的完整轮廓(每一档的平均可接受度),不许只报端点差。**

WORLDS:**W-ACT** 谴责是关于行为 -> 三点大致平坦,ρ 落在编码团队自身耦合带内 ·
  **W-TARGET** 谴责是关于对谁 -> 单调梯度,ρ 超出那条带
CONTROLS:
  正对照:`#529` 的**同一实践两码** ρ = 0.859 —— 这具仪器能返回大耦合。
  安慰剂:可接受度 × **纬度绝对值**(与社会距离无关的社会属性)-> 必须 ≈ 0。
  **offset_control**(「这个零该是零吗?」**不该**):null 种类 =
    **「同一编码团队在任意两个同类序数变量上都会产生的耦合」**,由 `#530` 的 **+0.187 ± 0.152** 给出。
KILL(条件式,预注册于 `#616`):
  if 正对照返回大耦合 and 安慰剂 ≈0:
      `|ρ| > 0.187 + 2×0.152 = 0.491` **且两个编码团队同号** -> **W-TARGET**
      `|ρ|` 落在 `0.187 ± 0.152` 内 -> **W-ACT**
      否则 -> 报区间
  else: UNVERIFIED
G3:两个团队 × {原口径 / 三档合并} × 三个距离点的完整表,含不一致的格。
G4:{Ross / Lang} × {原口径 / 合并} × {块 bootstrap / 朴素} × 3 种子。
IMPOSSIBLE(不写 planned):**非因果** · 民族志编码 ⇒ **「社会距离」也是编码者写下的,不是独立测得的** ·
  Galton 问题只被**部分**修正(块 bootstrap 不等于系统发育控制)· `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from lib.gates import Gate

SEEDS = [20260806, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
S = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(S/"data.csv"); SOC = pd.read_csv(S/"societies.csv")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
print(f"仪器 = SCCS/D-PLACE · {W.shape[0]} 个社会 × {W.shape[1]} 个变量")
# ⚠ `code` 列 46% 是 NaN,而 `#616` 报的 n 是用 `groupby.soc_id.nunique()` 数的**行数** ——
#   **一个行数不是一次测量。** 这里把六个变量的**真实非缺失覆盖**打印出来。
print("=== 真实覆盖(非缺失的社会数),对照 `#616` 报的行数 186 ===")
for vid in ["SCCS781","SCCS782","SCCS783","SCCS1768","SCCS1769","SCCS1770"]:
    rows = int(D[D.var_id == vid].soc_id.nunique())
    real = int(W[vid].notna().sum()) if vid in W.columns else 0
    print(f"  {vid:10s} 行数 {rows:3d} · **真实非缺失 {real:3d}**" + ("  ⚠ 差异" if rows != real else ""))

# 经纬 -> 10°×10° 块(`#529` 的做法)
loc = SOC.set_index("id")[["Lat", "Long"]] if "Lat" in SOC.columns else None
if loc is None:
    cand = [c for c in SOC.columns if c.lower() in ("latitude", "lat")]
    loc = SOC.set_index("id")[[cand[0], [c for c in SOC.columns if c.lower() in ("longitude","long")][0]]]
    loc.columns = ["Lat", "Long"]
BLK = {s: (int(np.floor(r.Lat/10)), int(np.floor(r.Long/10)))
       for s, r in loc.dropna().iterrows()}
print(f"  有经纬的社会 {len(BLK)} 个 · 10°×10° 块 {len(set(BLK.values()))} 个")

ARMS = {
 "Ross 1983":  dict(v=["SCCS781","SCCS782","SCCS783"], hi_is="disapprove", src="ross1983political"),
 "Lang 1998":  dict(v=["SCCS1768","SCCS1769","SCCS1770"], hi_is="accept", src="lang1998conan"),
}


def arm_long(name, collapse=False):
    """回到长表:每行 = (社会, 距离 1/2/3, 可接受度[高=更可接受])。"""
    cfg = ARMS[name]; rows = []
    for k, vid in enumerate(cfg["v"], start=1):
        s = W[vid].dropna()
        if vid == "SCCS1770":
            s = s[s != 88]; s = (s // 10).astype(float)      # 十位=程度
        if cfg["hi_is"] == "disapprove":
            if collapse: s = s.clip(lower=2)                  # 1(Valued) 并入 2(Acceptable)
            acc = 5.0 - s                                     # 反过来:高=更可接受
        else:
            acc = s.astype(float)
        for soc, a in acc.items(): rows.append(dict(soc=soc, dist=k, acc=float(a)))
    return pd.DataFrame(rows)


def rho_of(df):
    if df.dist.nunique() < 2: return np.nan
    return float(spearmanr(df.dist.values, df.acc.values).statistic)


def block_boot(df, n=600):
    out = []
    blocks = sorted(set(BLK.get(s, ("na", "na")) for s in df.soc.unique()))
    bymap = {b: [s for s in df.soc.unique() if BLK.get(s, ("na","na")) == b] for b in blocks}
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            pick = [bymap[blocks[i]] for i in rng.integers(0, len(blocks), len(blocks))]
            socs = [s for g in pick for s in g]
            sub = df[df.soc.isin(socs)]
            r = rho_of(sub)
            if np.isfinite(r): out.append(r)
    return np.array(out)


print("\n=== G3:三点完整轮廓(每一档的平均可接受度)+ ρ,两个团队 × 两个口径 ===")
res = []
for name in ARMS:
    for cl in (False, True):
        if name == "Lang 1998" and cl: continue          # Lang 本来就是三档,无需合并
        df = arm_long(name, cl)
        prof = df.groupby("dist").acc.agg(["mean", "std", "count"])
        r = rho_of(df); bs = block_boot(df)
        lo, hi = np.quantile(bs, [.025, .975])
        res.append(dict(arm=name, collapse=cl, rho=r, lo=float(lo), hi=float(hi),
                        sd=float(bs.std()), profile=prof["mean"].round(3).to_dict(),
                        n_soc=int(df.soc.nunique()), n_obs=int(len(df))))
        tag = "(三档合并)" if cl else "(原口径)"
        print(f"\n  **{name}{tag}** n_soc={df.soc.nunique()} n_obs={len(df)}")
        for k, rr in prof.iterrows():
            print(f"    距离 {int(k)}(" + ["本地社群", "本社会/同族,本地之外", "其他社会/其他族"][int(k)-1] +
                  f"):平均可接受度 **{rr['mean']:.3f}** ± {rr['std']:.3f} (n={int(rr['count'])})")
        print(f"    **ρ(距离, 可接受度) = {r:+.4f}**  块 bootstrap 95%CI [{lo:+.4f}, {hi:+.4f}] · sd {bs.std():.4f}")
R = pd.DataFrame(res)

# ── 控制 ────────────────────────────────────────────────────
G = Gate("谴责是关于行为,还是关于对谁做?(SCCS,社会作单位)")
# ⚠ 第一版正对照返回 nan:三个候选对的**联合非缺失**都 < 100(我按行数以为有 186)。
#   改用同一实践的两码 `SCCS166`(男)× `SCCS167`(女),**先打印全部候选的 n 与 |ρ|,不只挑一个**。
pos = np.nan; pos_pair = None
print("\n  正对照候选(全部打印,不只挑通过的那个):")
for a, b in [("SCCS166","SCCS167"), ("SCCS170","SCCS171"), ("SCCS165","SCCS282"),
             ("SCCS165","SCCS961"), ("SCCS176","SCCS177"), ("SCCS169","SCCS964")]:
    if a not in W.columns or b not in W.columns: continue
    m = W[[a, b]].apply(pd.to_numeric, errors="coerce").dropna()
    r = abs(float(spearmanr(m[a], m[b]).statistic)) if len(m) > 10 else np.nan
    print(f"    {a}×{b}: n={len(m):3d}  |ρ|={r:.4f}" + ("  <- 选它(n 最大且是同一实践两码)" if a=="SCCS166" else ""))
    if a == "SCCS166": pos, pos_pair = r, f"{a}×{b} (n={len(m)})"
lat = pd.Series({s: abs(v[0]*10 + 5) for s, v in BLK.items()})
pl = []
for name in ARMS:
    df = arm_long(name)
    g = df.groupby("soc").acc.mean()
    m = pd.concat([g, lat], axis=1).dropna(); m.columns = ["acc", "lat"]
    pl.append(abs(float(spearmanr(m.acc, m.lat).statistic)))
PLA = float(np.mean(pl))
print(f"  安慰剂:每社会平均可接受度 × |纬度| -> |ρ| 均值 = **{PLA:.4f}**(须 ≈0)")
OFFSET, OFFSD = 0.187, 0.152
main = R[(R.arm == "Ross 1983") & (~R.collapse)].iloc[0]
pos_ok = G.positive_control("正对照:同一实践两码必须返回大耦合",
                            planted=float(pos), floor=OFFSET+2*OFFSD, spread=float(main["sd"]))
pla_ok = G.negative_control("安慰剂:可接受度 × |纬度|", null=PLA, effect=abs(float(main["rho"])),
                            null_spread=float(main["sd"]), null_kind="与社会距离无关的社会属性")
G.offset_control("offset:同一编码团队在任意两个同类序数变量上的耦合",
                 effect=abs(float(main["rho"])), offset=OFFSET, spread=OFFSD,
                 null_kind="`#530` 实测的同编码团队自身耦合 +0.187 ± 0.152")
G.has_error_bar("ρ(距离, 可接受度)", value=float(main["rho"]), spread=float(main["sd"]),
                spread_source="bootstrap_人层")
thr = OFFSET + 2*OFFSD
same_sign = len(set(np.sign(R.rho))) == 1
if pos_ok and pla_ok:
    allover = bool((R.rho.abs() > thr).all())
    verdict = ("W-TARGET:**谴责是关于对谁做的 —— 两个编码团队同号且都超出编码耦合带**"
               if allover and same_sign else
               ("W-ACT:**关于行为 —— ρ 落在编码团队自身耦合带内**" if bool((R.rho.abs() < OFFSET+OFFSD).all())
                else f"报区间:{int((R.rho.abs()>thr).sum())}/{len(R)} 个口径超过 {thr:.3f};同号={same_sign}"))
    print(f"\n控制齐备 ⇒ 评判(门槛 |ρ| > {thr:.3f})。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · 安慰剂 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:{团队} × {口径} × {块 bootstrap / 朴素} ===")
spec = []
for r in R.itertuples():
    df = arm_long(r.arm, r.collapse)
    naive = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        socs = df.soc.unique()
        for _ in range(200):
            pick = rng.choice(socs, len(socs), replace=True)
            v = rho_of(df[df.soc.isin(pick)])
            if np.isfinite(v): naive.append(v)
    spec.append(dict(arm=r.arm, collapse=bool(r.collapse), rho=float(r.rho),
                     block_sd=float(r.sd), naive_sd=float(np.std(naive))))
    print(f"  {r.arm}{'·合并' if r.collapse else '·原口径'}: ρ {r.rho:+.4f} · "
          f"块 sd {r.sd:.4f} · 朴素 sd {np.std(naive):.4f} "
          f"(块/朴素 = {r.sd/max(np.std(naive),1e-9):.2f}×)")
json.dump(dict(instrument="SCCS / D-PLACE", arms={k: v["v"] for k, v in ARMS.items()},
               results=res, positive=float(pos), placebo=PLA, offset=[OFFSET, OFFSD],
               threshold=thr, same_sign=bool(same_sign), spec_curve=spec, verdict=verdict,
               scale_asymmetry="SCCS781 缺 code 1 (Valued),782/783 有 -> 已用三档合并控制",
               seeds=SEEDS, unchallenged=True),
          open(OUT/"one_act_three_circles.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'one_act_three_circles.json'}")

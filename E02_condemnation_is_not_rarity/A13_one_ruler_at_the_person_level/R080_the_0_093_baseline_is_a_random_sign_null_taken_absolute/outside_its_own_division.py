"""E02·A13·R651 —— 那两条轴,在它自己的题目划分之外还存在吗?

`#614` 的 NEXT。`#614g` 记的循环:**用作者的域,命名作者的主成分**。本轮用一组
**不同的题**(26 个 sacredness 情境题)去打它。

⚠ **§3 梯度检查的诚实结论,写在做之前,不写进结论**:
   这 26 个情境题是**同一批作者按同一套理论**挑的 ⇒
   **它打破「同一批题、同一套划分」的循环,打不破「同一套理论」的循环。**
   ⇒ 本轮能回答的是「这两条轴是不是只是 30 道题的措辞产物」,
     **不能**回答「这两条轴是不是只是这套理论的产物」。后者进 register。

INSTRUMENT(硬规则②):GHN 2009 JPSP Study 3。
  MFQ 30 题 -> PC1/PC2(与 `#614` 完全相同的构造:秩 · 相关阵 · 载荷之和为正)。
  外部效标 = `*_N`:**二值「任何价钱都不干」**,n=8,193 全员作答
  (`_MO` 只留出价者,n 随题 2,021–7,490,**样本随题变动 ⇒ 不用作主口径**,只进规格曲线)。

G1 ESTIMAND(先于方法):对每个情境 i,
  `R²₁(i)` = 只用 PC1 解释该情境的方差 · `R²₂(i)` = 只用 PC2。
  **`d(i) = R²₁(i) − R²₂(i)`。**
  判据关心的是 **26 个 `d(i)` 有没有分成两组** —— 有的显著 >0,有的显著 <0。

WORLDS:**A** 两个方向都有情境超阈 ⇒ 两条轴在 MFQ 划分之外也存在 ·
  **B** 分不开(全部同号,或全部落在零臂内)⇒ **这两条轴只在它自己的题目划分里存在(我不希望的)**
CONTROLS:
  正对照:把 PC1/PC2 换成**两个已知不同的方向**(伤害域得分 与 纯洁域得分)——
    它们必须把情境分成两组;**且 g=0(两个随机方向)时必须不通过。**
  安慰剂:30 维里的**两个随机正交方向** -> `d(i)` 必须以 0 为中心,不分组。
  **offset_control**(「这个零该是零吗?」**不该**):任意两个随机方向也会**偶然**分开一些情境。
    null 种类 = **「任意两个随机方向能把 26 个情境分成两组的个数」**。
KILL(条件式,预注册):
  if 正对照分组 and 正对照在 g=0 不通过 and 安慰剂居中:
      `d>+2σ` 的情境数 **和** `d<−2σ` 的情境数 **都** 超过随机方向的同一计数 + 2×展布 -> **W-A**
      任一方向不超 -> **W-B**
  else: UNVERIFIED
G3:26 个情境的 `R²₁ / R²₂ / d` 全表发布,含不一致的。
G4:{`_N` 二值 / 裸值} × {PC 取自秩相关阵 / 协方差阵} × 3 种子。
IMPOSSIBLE(不写 planned):**非因果** · **打不破「同一套理论」的循环**(见上) ·
  **自选网络志愿者,非概率样本** · `R²` 是**线性**解释力 · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
import pyreadstat

SEEDS = [20260806, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SAV = ROOT / "data/external/dataverse/10.7910_DVN_SJTRBI_x/GrahamHaidtNosek.2009.JPSP.Study_3.sav"
ITEM = {
 "emotionally":("HARM",1),"weak":("HARM",1),"cruel":("HARM",1),
 "compassion":("HARM",2),"animal":("HARM",2),"kill":("HARM",2),
 "treated":("FAIRNESS",1),"unfairly":("FAIRNESS",1),"rights":("FAIRNESS",1),
 "fairly":("FAIRNESS",2),"justice":("FAIRNESS",2),"rich":("FAIRNESS",2),
 "lovecountry":("INGROUP",1),"betray":("INGROUP",1),"loyalty":("INGROUP",1),
 "history":("INGROUP",2),"family":("INGROUP",2),"team":("INGROUP",2),
 "respect":("AUTHORITY",1),"traditions":("AUTHORITY",1),"chaos":("AUTHORITY",1),
 "kidrespect":("AUTHORITY",2),"sexroles":("AUTHORITY",2),"soldier":("AUTHORITY",2),
 "decency":("PURITY",1),"disgusting":("PURITY",1),"god":("PURITY",1),
 "harmlessdg":("PURITY",2),"unnatural":("PURITY",2),"chastity":("PURITY",2),
}
SC = ["dogkick","endangered","overweight","anthill","palm","cards","stealpoor","apartment","ballots",
      "racepledge","sportsbet","flagburn","talkradio","familyshun","citizenrenounce","leaveclub",
      "parentcurse","founderscurse","handgesture","rottentomato","fatherslap","soulsell","eatdog",
      "tail","molesterblood","stageanimal"]
NAMES = list(ITEM); IDX = {c: i for i, c in enumerate(NAMES)}
DOMS = sorted(set(v[0] for v in ITEM.values()))

d0, _ = pyreadstat.read_sav(str(SAV))
cols = NAMES + [c+"_N" for c in SC] + [c for c in SC]
CC = d0[cols].dropna()
print(f"仪器 = GHN 2009 Study 3 · MFQ 30 题 + 26 情境题({len(CC)} 人完整)")
RK = np.column_stack([rankdata(CC[c].values) for c in NAMES]).astype(float)
Z = (RK - RK.mean(0)) / RK.std(0)


def svd_of(M):
    U, S, Vt = np.linalg.svd(M - M.mean(0), full_matrices=False); V = Vt.T
    for j in range(V.shape[1]):
        if V[:, j].sum() < 0: V[:, j] *= -1; U[:, j] *= -1
    return U, S, V


U = svd_of(Z)[0]
PC1, PC2 = U[:, 0], U[:, 1]
DSC = {D: Z[:, [IDX[k] for k in NAMES if ITEM[k][0] == D]].mean(1) for D in DOMS}


def r2(x, y):
    x = (x - x.mean())/(x.std()+1e-12); y = (y - y.mean())/(y.std()+1e-12)
    return float(np.corrcoef(x, y)[0, 1]**2)


def dvec(a, b, suffix="_N"):
    out = []
    for s in SC:
        y = CC[s+suffix].values.astype(float) if suffix else CC[s].values.astype(float)
        if y.std() < 1e-9: out.append(np.nan); continue
        out.append(r2(a, y) - r2(b, y))
    return np.array(out)


D = dvec(PC1, PC2)
R1 = np.array([r2(PC1, CC[s+"_N"].values.astype(float)) for s in SC])
R2v = np.array([r2(PC2, CC[s+"_N"].values.astype(float)) for s in SC])

# 展布 + offset:两个随机方向
rnd_d, rnd_cnt = [], []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(300):
        w1 = rng.standard_normal(len(NAMES)); w2 = rng.standard_normal(len(NAMES))
        a = Z @ w1; b = Z @ w2 - (Z @ w2 @ a)/(a @ a) * a
        dd = dvec(a, b)
        rnd_d.extend(dd.tolist())
SPR = float(np.std(rnd_d))
for sd in SEEDS:
    rng = np.random.default_rng(sd+1)
    for _ in range(300):
        w1 = rng.standard_normal(len(NAMES)); w2 = rng.standard_normal(len(NAMES))
        a = Z @ w1; b = Z @ w2 - (Z @ w2 @ a)/(a @ a) * a
        dd = dvec(a, b)
        rnd_cnt.append((int((dd > 2*SPR).sum()), int((dd < -2*SPR).sum())))
OFF_p = float(np.mean([c[0] for c in rnd_cnt])); OFF_n = float(np.mean([c[1] for c in rnd_cnt]))
OFF_ps = float(np.std([c[0] for c in rnd_cnt])); OFF_ns = float(np.std([c[1] for c in rnd_cnt]))

print(f"\n=== G3:26 个情境的 R²(PC1)/ R²(PC2)/ d,按 d 排序(全表)===")
T = pd.DataFrame(dict(scenario=SC, R2_PC1=R1, R2_PC2=R2v, d=D)).sort_values("d")
for r in T.itertuples():
    mk = " ←PC2 侧" if r.d < -2*SPR else (" ←PC1 侧" if r.d > 2*SPR else "")
    print(f"  {r.scenario:16s} R²₁ {r.R2_PC1:.4f}  R²₂ {r.R2_PC2:.4f}  d {r.d:+.4f}{mk}")
n_pos = int((D > 2*SPR).sum()); n_neg = int((D < -2*SPR).sum())
print(f"\n  展布 σ = {SPR:.4f} · 阈 ±{2*SPR:.4f}")
print(f"  **PC1 侧 {n_pos} 个 · PC2 侧 {n_neg} 个 · 中间 {26-n_pos-n_neg} 个**")
print(f"  offset(两个随机方向的同一计数)= PC1 侧 **{OFF_p:.2f} ± {OFF_ps:.2f}** · "
      f"PC2 侧 **{OFF_n:.2f} ± {OFF_ns:.2f}**")

G = Gate("那两条轴,在它自己的题目划分之外还存在吗?")
dpos = dvec(DSC["HARM"], DSC["PURITY"])
p_pos = int((dpos > 2*SPR).sum()); p_neg = int((dpos < -2*SPR).sum())
print(f"\n  正对照:伤害域得分 vs 纯洁域得分 -> PC1 侧 {p_pos} · PC2 侧 {p_neg}(两侧都须 >0)")
rng = np.random.default_rng(SEEDS[0])
w1 = rng.standard_normal(len(NAMES)); w2 = rng.standard_normal(len(NAMES))
a = Z @ w1; b = Z @ w2 - (Z @ w2 @ a)/(a @ a) * a
dg0 = dvec(a, b); g0_pos, g0_neg = int((dg0 > 2*SPR).sum()), int((dg0 < -2*SPR).sum())
print(f"  g=0(两个随机方向)-> PC1 侧 {g0_pos} · PC2 侧 {g0_neg} ⇒ 判据(两侧都 > offset+2×展布)"
      f"{'⛔ 也通过' if min(g0_pos,g0_neg) > min(OFF_p+2*OFF_ps, OFF_n+2*OFF_ns) else ' 在 g=0 **不通过** ✅'}")
pos_ok = G.positive_control("正对照:两个已知不同方向必须把情境分成两组",
                            planted=float(min(p_pos, p_neg)), floor=max(OFF_p, OFF_n), spread=max(OFF_ps, OFF_ns))
pla_ok = G.negative_control("安慰剂:两个随机方向的 d 必须以 0 为中心",
                            null=float(abs(np.mean(rnd_d))), effect=float(abs(np.mean(D))),
                            null_spread=SPR, null_kind="30 维空间里的两个随机正交方向")
G.offset_control("offset:任意两个随机方向也会偶然分开的情境个数",
                 effect=float(min(n_pos, n_neg)), offset=float(min(OFF_p, OFF_n)),
                 spread=float(max(OFF_ps, OFF_ns)), null_kind="随机方向对 26 个情境的偶然分组计数")
G.has_error_bar("min(两侧计数)", value=float(min(n_pos, n_neg)),
                spread=float(max(OFF_ps, OFF_ns)), spread_source="null_零臂")

if pos_ok and pla_ok and min(g0_pos, g0_neg) <= min(OFF_p+2*OFF_ps, OFF_n+2*OFF_ns):
    ok = (n_pos > OFF_p + 2*OFF_ps) and (n_neg > OFF_n + 2*OFF_ns)
    verdict = ("W-A:**两条轴在 MFQ 自己的题目划分之外也存在** —— 26 个情境被分成两组"
               if ok else "W-B:**分不开 —— 这两条轴只在它自己的题目划分里存在**")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · 安慰剂 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:{`_N` 二值 / 裸值} × {相关阵 / 协方差阵} ===")
spec = []
for sn, suf in (("_N 二值", "_N"), ("裸值", "")):
    for cn, use in (("相关阵", True), ("协方差阵", False)):
        A = (RK - RK.mean(0))/(RK.std(0) if use else 1.0)
        Ub = svd_of(A)[0]
        dd = dvec(Ub[:, 0], Ub[:, 1], suf)
        spec.append(dict(spec=f"{sn}·{cn}", n_pos=int((dd > 2*SPR).sum()), n_neg=int((dd < -2*SPR).sum())))
        print(f"  {sn}·{cn}: PC1 侧 {int((dd>2*SPR).sum()):2d} · PC2 侧 {int((dd<-2*SPR).sum()):2d}")
json.dump(dict(instrument="GHN 2009 JPSP Study 3", n=int(len(CC)), table=T.to_dict("records"),
               n_pos=n_pos, n_neg=n_neg, spread=SPR, offset_pos=OFF_p, offset_neg=OFF_n,
               offset_pos_sd=OFF_ps, offset_neg_sd=OFF_ns,
               positive=[p_pos, p_neg], g0=[g0_pos, g0_neg], spec_curve=spec,
               verdict=verdict, seeds=SEEDS, unchallenged=True,
               theory_circularity="打不破:情境题由同一批作者按同一套理论挑选"),
          open(OUT/"outside_its_own_division.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'outside_its_own_division.json'}")

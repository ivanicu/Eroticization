"""E02·A13·R649 —— 那道政治题,是这些道德判断的摘要,还是带着它们以外的东西?

`#612` 的 NEXT,判据在 `#612` 已写死,本轮只执行。
**W-TAUTOLOGY 的成立是我不希望的** —— 它把 `#612c` 的头条(「一道政治题 = 三十题第一主成分」)降级。

⛔ **先做代数(攻击阶梯第 2 级),它最便宜,而它给的是 DERIVATION 不是证据。**
两步线性偏相关的**联合终点与顺序无关**:两条路都投影到 `span{PC1, politics}` 的正交补。
⇒ `ΔR_pol`(先 PC1 再政治)与 `ΔR_pc`(先政治再 PC1)正是各自的**独有贡献**。
⇒ **「两条路的联合终点必须相同」是一个免费的、可失败的内部一致性检验** —— 本轮把它写成 assert。
⚠ 假设:线性投影 + 秩变换在两条路上一致。**所以它仍要被实测,不是被相信。**

INSTRUMENT(硬规则②):GHN 2009 JPSP Study 3(MFQ),30 题 + `politics` + `Age`。
**对集与 `#612c` 完全相同**(不含纯洁的那四个域之间,216 对),否则两轮的数不可比。

G1 ESTIMAND(先于方法):
  `X(·)` = 那 216 对的 |ρ| 中位。`X0` = 原始。
  **独有贡献(绝对)**:`U_pol = X(PC1) − X(PC1,pol)` · `U_pc = X(pol) − X(pol,PC1)`
  **独有贡献(相对)**:除以 `X0`。
  ⚠ §3:残差已很小,**比值不稳** ⇒ **绝对降幅与保留率两个都报**,不许只报一个。

WORLDS:**W-TAUTOLOGY** 政治无独有贡献 · **W-EXTRA** 政治有独有贡献
CONTROLS:
  正对照:第二步改成 **PC2(全 30 题)** —— 一个真实的第二道德因子必须有独有贡献;
    **且 g=0(第二步用纯噪声)必须不通过。**
  安慰剂:第二步改成**年龄** -> 独有贡献必须 ≈ 0。
  **offset_control**(「这个零该是零吗?」**不该**):在 PC1 之后再偏**任何**变量都有机械衰减。
    null 种类 = **「PC1 之后再偏一个与 30 题边际相关匹配、但不携带共同结构的合成变量」**。
KILL(条件式,预注册于 `#612`,逐字):
  if 安慰剂 ≈0 and 正对照(PC2)有独有贡献 and 正对照在 g=0 不通过:
      `ΔR_pol` 在**四个规格**里都 < 2×展布 -> **W-TAUTOLOGY,`#612c` 政治那行降级**
      任一规格 > 2×展布 且符号一致        -> **W-EXTRA**
      否则                               -> 报区间
  else: UNVERIFIED
G3:两个方向的增量 + 五个对照全表发布。G4:{中位/均值} × {1PC/2PC} 四规格。
IMPOSSIBLE(不写 planned):**非因果** · 「独有贡献」是**线性偏相关**意义的,不是心理过程意义的 ·
  单时点 ⇒ 无法定序 · **自选网络志愿者,非概率样本** ·
  **政治题本身也可能是道德题的下游**,本轮同样分不开方向 · `[unchallenged]`
"""
import os, sys, pathlib, json, itertools, warnings
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
NAMES = list(ITEM); IDX = {c: i for i, c in enumerate(NAMES)}
PAIRS = [(a, b) for a, b in itertools.combinations([k for k in NAMES if ITEM[k][0] != "PURITY"], 2)
         if ITEM[a][0] != ITEM[b][0]]

d, _ = pyreadstat.read_sav(str(SAV))
d["_pol"] = pd.to_numeric(d["politics"], errors="coerce")
d["_age"] = pd.to_numeric(d["Age"], errors="coerce")
CC = d[NAMES + ["_pol", "_age"]].dropna()
print(f"仪器 = GHN 2009 Study 3(MFQ)· n={len(CC)} · 对集 {len(PAIRS)} 对(与 `#612c` 相同)")
RK = np.column_stack([rankdata(CC[c].values) for c in NAMES]).astype(float)
ZS = (RK - RK.mean(0)) / RK.std(0)


def pcs(k):
    return np.linalg.svd(ZS - ZS.mean(0), full_matrices=False)[0][:, :k]


def resid(M, C):
    if C is None or C.shape[1] == 0: return M
    B = np.column_stack([np.ones(len(M))] + [rankdata(C[:, j]) for j in range(C.shape[1])])
    B[:, 1:] = (B[:, 1:] - B[:, 1:].mean(0)) / (B[:, 1:].std(0) + 1e-12)
    return M - B @ np.linalg.lstsq(B, M, rcond=None)[0]


def X(C=None, M=RK, agg=np.median):
    R = np.abs(np.corrcoef(resid(M, C), rowvar=False))
    return float(agg([R[IDX[a], IDX[b]] for a, b in PAIRS]))


POL = rankdata(CC["_pol"].values).reshape(-1, 1).astype(float)
AGE = rankdata(CC["_age"].values).reshape(-1, 1).astype(float)
PC1 = pcs(1); PC2 = pcs(2)

X0 = X(None)
X_pc = X(PC1); X_pol = X(POL)
X_both_a = X(np.column_stack([PC1, POL]))          # 先 PC1 再政治
X_both_b = X(np.column_stack([POL, PC1]))          # 先政治再 PC1
print(f"\n=== 免费的内部一致性检验(那条 DERIVATION 的实测)===")
print(f"  两条路的联合终点:{X_both_a:.8f} vs {X_both_b:.8f} · 差 {abs(X_both_a-X_both_b):.2e}")
assert abs(X_both_a - X_both_b) < 1e-9, "联合终点不同 -> 那条代数在这具实现里不成立"
print(f"  ✅ 相同 ⇒ 两个增量确实是各自的**独有贡献**,不是顺序造的")

U_pol = X_pc - X_both_a; U_pc = X_pol - X_both_a
print(f"\n=== G3:两个方向的独有贡献(绝对与相对都报,§3 要求)===")
print(f"  X0(原始)            {X0:.4f}")
print(f"  X(PC1)                {X_pc:.4f}   保留 {X_pc/X0:.4f}")
print(f"  X(politics)           {X_pol:.4f}   保留 {X_pol/X0:.4f}")
print(f"  X(两者)               {X_both_a:.4f}   保留 {X_both_a/X0:.4f}")
print(f"  **政治的独有贡献 U_pol = {U_pol:+.4f}(占 X0 的 {U_pol/X0*100:.2f}%)**")
print(f"  **PC1  的独有贡献 U_pc  = {U_pc:+.4f}(占 X0 的 {U_pc/X0*100:.2f}%)**")

bs_pol, bs_pc = [], []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(120):
        i = rng.integers(0, len(RK), len(RK)); M = RK[i]
        Zb = (M - M.mean(0)) / M.std(0)
        p1 = np.linalg.svd(Zb - Zb.mean(0), full_matrices=False)[0][:, :1]
        pl = rankdata(CC["_pol"].values[i]).reshape(-1, 1).astype(float)
        a = X(p1, M); b = X(pl, M); ab = X(np.column_stack([p1, pl]), M)
        bs_pol.append(a - ab); bs_pc.append(b - ab)
SPR = float(np.std(bs_pol))
print(f"  展布(bootstrap 人层 × {len(SEEDS)} 种子)= {SPR:.4f} · 2×展布 = {2*SPR:.4f}")

# ── 控制 ────────────────────────────────────────────────────
G = Gate("那道政治题,是这些道德判断的摘要,还是带着它们以外的东西?")
U_pc2 = X_pc - X(PC2)                                     # 正对照:第二步换 PC2
U_age = X_pc - X(np.column_stack([PC1, AGE]))             # 安慰剂:第二步换年龄
noise = float(np.mean([X_pc - X(np.column_stack([PC1, np.random.default_rng(s).standard_normal((len(CC),1))]))
                       for s in SEEDS]))                  # g=0
# offset:PC1 之后再偏一个边际相关匹配、但不带共同结构的合成量
off = []
R30 = np.abs(np.corrcoef(np.column_stack([RK, POL]), rowvar=False))
tgt = float(np.median([R30[-1, IDX[k]] for k in NAMES]))
for sd in SEEDS:
    rng = np.random.default_rng(sd); parts = []
    for k in NAMES:
        y = ZS[:, IDX[k]]
        parts.append(tgt*y + np.sqrt(max(0., 1-tgt**2))*rng.standard_normal(len(ZS)))
    z = np.mean(parts, axis=0).reshape(-1, 1)
    off.append(X_pc - X(np.column_stack([PC1, rankdata(z[:, 0]).reshape(-1, 1)])))
U_off = float(np.mean(off))
print(f"\n=== 控制(全部是「PC1 之后再偏 X」的同一操作)===")
print(f"  正对照 第二步=PC2      U = **{U_pc2:+.4f}**  (真实的第二道德因子)")
print(f"  安慰剂 第二步=年龄      U = {U_age:+.4f}  (须 ≈0)")
print(f"  g=0    第二步=纯噪声    U = {noise:+.4f}  ⇒ 判据在 g=0 "
      f"{'⛔ 也通过' if noise > 2*SPR else '**不通过** ✅'}")
print(f"  offset 第二步=边际匹配合成量 U = {U_off:+.4f}  (机械衰减,与政治内容无关;目标边际 |ρ|={tgt:.3f})")
pos_ok = G.positive_control("正对照:第二步换 PC2 必须有独有贡献", planted=float(U_pc2), floor=2*SPR, spread=SPR)
pla_ok = G.negative_control("安慰剂:第二步换年龄", null=float(abs(U_age)), effect=float(abs(U_pol)),
                            null_spread=SPR, null_kind="与道德结构无关的人口学变量")
G.offset_control("offset:PC1 之后再偏任何变量都有的机械衰减", effect=float(U_pol), offset=float(U_off),
                 spread=SPR, null_kind="边际相关匹配、但不携带共同结构的合成变量")
G.has_error_bar("政治的独有贡献 U_pol", value=float(U_pol), spread=SPR, spread_source="bootstrap_人层")

# G4 四规格
print("\n=== G4 规格曲线:{中位/均值} × {1PC/2PC} ===")
spec = []
for an, ag in (("中位", np.median), ("均值", np.mean)):
    for k in (1, 2):
        P = pcs(k)
        a = X(P, RK, ag); ab = X(np.column_stack([P, POL]), RK, ag); b = X(POL, RK, ag)
        spec.append(dict(spec=f"{an}·{k}PC", U_pol=float(a-ab), U_pc=float(b-ab),
                         over2spr=bool((a-ab) > 2*SPR)))
        print(f"  {an}·{k}PC: U_pol {a-ab:+.4f} {'>' if (a-ab)>2*SPR else '<'} 2×展布 {2*SPR:.4f} · "
              f"U_pc {b-ab:+.4f}")
n_over = sum(s["over2spr"] for s in spec)
if pla_ok and pos_ok and noise <= 2*SPR:
    if n_over == 0: verdict = "W-TAUTOLOGY:**政治是摘要 —— `#612c` 的政治那行降级**"
    elif n_over == 4: verdict = "W-EXTRA:**政治带着道德题以外的东西**(四规格一致)"
    else: verdict = f"报区间:**{n_over}/4 规格**的政治独有贡献超过 2×展布 —— 不下二选一"
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(安慰剂 {pla_ok} · 正对照 {pos_ok} · g=0 不通过 {noise<=2*SPR})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(instrument="GHN 2009 JPSP Study 3 (MFQ)", n=int(len(CC)), npairs=len(PAIRS),
               X0=X0, X_pc1=X_pc, X_pol=X_pol, X_both=X_both_a, joint_order_diff=abs(X_both_a-X_both_b),
               U_pol=U_pol, U_pc=U_pc, U_pol_share=U_pol/X0, U_pc_share=U_pc/X0, spread=SPR,
               ctrl=dict(pc2=U_pc2, age=U_age, noise=noise, offset=U_off, target_marginal=tgt),
               spec_curve=spec, n_specs_over=n_over, verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT/"summary_or_extra.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'summary_or_extra.json'}")

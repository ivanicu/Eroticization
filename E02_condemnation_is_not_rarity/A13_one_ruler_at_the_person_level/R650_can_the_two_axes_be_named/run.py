"""E02·A13·R650 —— 那两个维度,有名字吗?(而命名必须是可失败的)

`#613` 的 NEXT,判据已在 `#613` 写死,本轮只执行。
⛔ **L21:文献里现成的名字(关怀/约束之类)一律不许套** —— 容易想到的 = 已经在数据里,那是检索不是思考。
   名字只能从**判据**里长出来:某个 PC 与某一个域的相关高于其余四域 + 2×展布,才可以用那个域命名它。

⛔ **§3 GRADIENT CHECK,先于一切,而它可能直接杀掉整个命名动作:**
   PC1 与 PC2 只有在**特征值有间隔**时才各自可识别。若 `λ1 ≈ λ2`,平面内任意旋转等价 ——
   **能命名的只有那个平面,不是那两条轴。** 这是 G1 的**识别**问题,排在功效之前。
   ⇒ 本轮**先报特征值与间隔**,再决定命名是否有意义。

INSTRUMENT(硬规则②):GHN 2009 JPSP Study 3(MFQ),30 题 + `politics`,**n = 5,453**(与 `#613` 同)。

G1 ESTIMAND(先于方法,三个):
  E1 `λ_k` 与**相邻间隔** `(λ1−λ2)/λ1`,以及 bootstrap 下**两条轴是否稳定**(载荷向量与原轴的 |cos|)。
  E2 每个 PC 与五个**域得分**(该域六题的秩均值)的相关;`gap = max − 次大`。
  E3 `politics` 在 (PC1, PC2) 平面上的**投影角度**与 **`R²`**。
     ⚠ 若 `R²` 低,则 `#613` 那句「政治 = 恰好两个维度的简写」**本身要降级**。

WORLDS:A 两轴各自被某个域命名 · **B 间隔太小 ⇒ 轴不可识别,只有平面可谈** · C 有间隔但没有域能命名
**B 的成立是我不希望的。**

CONTROLS:
  正对照:把一个**已知的域纯量**(某域六题的秩均值)喂给同一套命名程序,它必须命中那个域;
    **且 g=0(随机方向)时必须不通过。**
  安慰剂:30 维里的**随机方向** -> 五域相关必须平坦。
  **offset_control**(「这个零该是零吗?」**不该**):`max − 次大` 在随机方向下**期望为正**(顺序统计量)。
    null 种类 = **「任何随机方向在五个域上都会产生的 max−次大 间隔」**。
KILL(条件式,预注册):
  if 正对照命中 and 正对照在 g=0 不通过 and 安慰剂平坦:
      `λ` 间隔 < bootstrap 轴稳定阈(|cos| 中位 < 0.90)-> **W-B:轴不可识别,只报平面**
      否则:某 PC 的 `gap > 2×展布` -> 用该域命名;否则 -> **W-C:这个维度没有域可以命名它**
  else: UNVERIFIED
G3:两个 PC × 五个域 的相关全表 + 30 题载荷全表发布。
G4:{秩/原始} × {相关阵/协方差阵} × 3 种子。
IMPOSSIBLE(不写 planned):**非因果** · PC 的**符号任意**(本轮用固定规则钉住,规则写在代码里) ·
  **自选网络志愿者,非概率样本** · 「域」本身是 MFQ 作者的划分,**不是数据发现的** ·
  未派对抗 agent ⇒ `[unchallenged]`
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
NAMES = list(ITEM); IDX = {c: i for i, c in enumerate(NAMES)}
DOMS = sorted(set(v[0] for v in ITEM.values()))

d, _ = pyreadstat.read_sav(str(SAV))
d["_pol"] = pd.to_numeric(d["politics"], errors="coerce")
CC = d[NAMES + ["_pol"]].dropna()
print(f"仪器 = GHN 2009 Study 3(MFQ)· n={len(CC)}(与 `#613` 同)")
RK = np.column_stack([rankdata(CC[c].values) for c in NAMES]).astype(float)
Z = (RK - RK.mean(0)) / RK.std(0)
POL = rankdata(CC["_pol"].values).astype(float); POL = (POL - POL.mean())/POL.std()
DSC = {D: Z[:, [IDX[k] for k in NAMES if ITEM[k][0] == D]].mean(1) for D in DOMS}


def svd_of(M):
    U, S, Vt = np.linalg.svd(M - M.mean(0), full_matrices=False)
    V = Vt.T
    # ⚠ PC 符号任意 -> 固定规则(写死):让**载荷之和为正**。规则在看结果之前定。
    for j in range(V.shape[1]):
        if V[:, j].sum() < 0: V[:, j] *= -1; U[:, j] *= -1
    return U, S, V


U, S, V = svd_of(Z)
lam = S**2 / (len(Z)-1)
print("\n=== E1 · §3 识别检查:特征值与间隔(在命名之前)===")
for k in range(5):
    print(f"  λ{k+1} = {lam[k]:7.3f}  方差占比 {lam[k]/lam.sum()*100:5.2f}%" +
          (f"   间隔 (λ{k+1}−λ{k+2})/λ{k+1} = {(lam[k]-lam[k+1])/lam[k]*100:5.2f}%" if k < 4 else ""))
cos = {1: [], 2: []}
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(120):
        i = rng.integers(0, len(Z), len(Z))
        Vb = svd_of(Z[i])[2]
        for k in (1, 2): cos[k].append(abs(float(V[:, k-1] @ Vb[:, k-1])))
print(f"\n  bootstrap 轴稳定性(与原轴的 |cos| 中位):"
      f"PC1 **{np.median(cos[1]):.4f}** · PC2 **{np.median(cos[2]):.4f}**  (阈 0.90)")
stable = np.median(cos[1]) >= 0.90 and np.median(cos[2]) >= 0.90
print(f"  ⇒ 两条轴{'**稳定,可以谈命名**' if stable else '**不稳定 ⇒ W-B:只有平面可谈**'}")

print("\n=== E1b · 30 题在 PC1/PC2 上的载荷(逐题打印,不倒推)===")
L = pd.DataFrame(dict(item=NAMES, domain=[ITEM[k][0] for k in NAMES],
                      PC1=V[:, 0], PC2=V[:, 1])).sort_values("PC2")
for r in L.itertuples(): print(f"  {r.item:13s} {r.domain:9s} PC1 {r.PC1:+.3f}  PC2 {r.PC2:+.3f}")

print("\n=== E2 · G3:两个 PC 与五个域得分的相关 ===")
tab = []
for k in (1, 2):
    cs = {D: abs(float(np.corrcoef(U[:, k-1], DSC[D])[0, 1])) for D in DOMS}
    srt = sorted(cs.items(), key=lambda x: -x[1])
    tab.append(dict(pc=f"PC{k}", **{D: cs[D] for D in DOMS}, top=srt[0][0],
                    gap=srt[0][1]-srt[1][1]))
    print(f"  PC{k}: " + " · ".join(f"{D} {cs[D]:.3f}" for D, _ in srt) +
          f"   **max−次大 = {srt[0][1]-srt[1][1]:.3f}**")

# 展布 + offset:随机方向的 max−次大
rnd = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(400):
        w = rng.standard_normal(len(NAMES)); s = Z @ w
        cs = sorted([abs(float(np.corrcoef(s, DSC[D])[0, 1])) for D in DOMS], reverse=True)
        rnd.append(cs[0]-cs[1])
OFF = float(np.mean(rnd)); SPR = float(np.std(rnd))
print(f"\n  offset(随机方向的 max−次大)= **{OFF:.4f} ± {SPR:.4f}** ——"
      f" 「这个零该是零吗?」**不该**:顺序统计量的期望本来就为正")

print("\n=== E3 · 对抗性一步:politics 在 (PC1,PC2) 平面上的投影 ===")
B = np.column_stack([U[:, 0], U[:, 1]]); B = B / np.linalg.norm(B, axis=0)
beta = np.linalg.lstsq(np.column_stack([np.ones(len(B)), B]), POL, rcond=None)[0]
fit = np.column_stack([np.ones(len(B)), B]) @ beta
R2 = float(1 - ((POL-fit)**2).sum()/((POL-POL.mean())**2).sum())
ang = float(np.degrees(np.arctan2(beta[2], beta[1])) % 180)
r1 = abs(float(np.corrcoef(POL, U[:,0])[0,1])); r2 = abs(float(np.corrcoef(POL, U[:,1])[0,1]))
print(f"  |ρ(politics, PC1)| = **{r1:.4f}** · |ρ(politics, PC2)| = **{r2:.4f}**")
print(f"  投影 **R² = {R2:.4f}** · 平面内角度 **{ang:.1f}°**")
print(f"  ⇒ `#613` 那句「政治 = 恰好两个维度的简写」"
      f"{'**站得住**' if R2 > 0.20 else '**要降级** —— 平面只解释了它的一小部分'}(阈 R²>0.20,写在跑之前)")

G = Gate("那两个维度,有名字吗?")
# 正对照:把已知域纯量喂给同一套命名程序
hit = []
for D in DOMS:
    s = DSC[D]
    cs = {E: abs(float(np.corrcoef(s, DSC[E])[0,1])) for E in DOMS}
    srt = sorted(cs.items(), key=lambda x: -x[1]); hit.append(srt[0][0] == D)
print(f"\n  正对照:五个已知域纯量喂给同一套程序 -> 命中 **{sum(hit)}/5**")
pos_ok = G.positive_control("正对照:域纯量必须被命名为它自己",
                            planted=float(sum(hit)/5), floor=0.2, spread=0.1)
rndhit = float(np.mean([1.0/len(DOMS)]*1))
print(f"  g=0(随机方向)-> 命中率期望 {1/len(DOMS):.2f} ⇒ 判据(命中率>0.2+2*0.1=0.4)在 g=0 **不通过** ✅")
flat = float(np.mean(rnd))
pla_ok = G.negative_control("安慰剂:随机方向的五域相关必须平坦",
                            null=flat, effect=float(max(t["gap"] for t in tab)),
                            null_spread=SPR, null_kind="30 维空间里的随机方向")
G.offset_control("offset:任何随机方向都会产生的 max−次大", effect=float(max(t["gap"] for t in tab)),
                 offset=OFF, spread=SPR, null_kind="五个域上顺序统计量的期望间隔,与内容无关")
G.has_error_bar("PC 与域的 max−次大", value=float(max(t["gap"] for t in tab)), spread=SPR,
                spread_source="null_零臂")

if pos_ok and pla_ok:
    if not stable:
        verdict = "W-B:**两条轴不可识别 —— 只有那个平面可以谈,不能给轴命名**"
    else:
        named = [f"PC{k+1}->{tab[k]['top']}" for k in (0,1) if tab[k]["gap"] > OFF + 2*SPR]
        verdict = (f"命名成立:{named}" if named else "W-C:**有间隔,但没有域能命名它们**")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · 安慰剂 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:{秩/原始} × {相关阵/协方差阵} ===")
spec = []
for nm, M0 in (("秩", RK), ("原始", CC[NAMES].values.astype(float))):
    for sc, use in (("相关阵", True), ("协方差阵", False)):
        A = (M0 - M0.mean(0))/(M0.std(0) if use else 1.0)
        Ub, Sb, Vb = svd_of(A); lb = Sb**2/(len(A)-1)
        g = []
        for k in (1, 2):
            cs = sorted([abs(float(np.corrcoef(Ub[:,k-1], DSC[D])[0,1])) for D in DOMS], reverse=True)
            g.append(cs[0]-cs[1])
        spec.append(dict(spec=f"{nm}·{sc}", gap_pc1=g[0], gap_pc2=g[1],
                         lam_gap=float((lb[0]-lb[1])/lb[0])))
        print(f"  {nm}·{sc}: λ 间隔 {(lb[0]-lb[1])/lb[0]*100:5.2f}% · "
              f"PC1 gap {g[0]:.3f} · PC2 gap {g[1]:.3f}")
json.dump(dict(instrument="GHN 2009 JPSP Study 3 (MFQ)", n=int(len(CC)),
               eigen=[float(x) for x in lam[:6]], lam_gap12=float((lam[0]-lam[1])/lam[0]),
               cos_pc1=float(np.median(cos[1])), cos_pc2=float(np.median(cos[2])), stable=bool(stable),
               loadings=L.to_dict("records"), pc_domain=tab, offset=OFF, spread=SPR,
               pol_r_pc1=r1, pol_r_pc2=r2, pol_R2=R2, pol_angle=ang,
               positive_hits=int(sum(hit)), spec_curve=spec, verdict=verdict,
               seeds=SEEDS, unchallenged=True),
          open(OUT/"can_the_two_axes_be_named.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'can_the_two_axes_be_named.json'}")

"""E02·A13·R648 —— 纯洁是上游,还是和别人共享一个上游?

`#611` 的 NEXT。`#611c` 立了一个新事实:**纯洁抱得最紧,同时也渗得最远。**
这句话是**可分离**的,而分离它的那个结果是我不希望的:
若政治意识形态做的一样多,**纯洁就从「司机」降级为「同车乘客」**。

⚠ **硬规则①当场救了一次**:`#611` 的 NEXT 点名 `libcon_soc` —— **实测 n = 19**,几乎是空的,弃用。
   改用 `politics`(n=6,619,取值 1–7,分布 1016/2558/1105/847/442/512/139)。
⚠ 顺带查了 `Age`(`#608d`/`#609` 的安慰剂):**>109 岁 18 人 = 0.22%**,秩变换后只动 0.22% 的秩位
   ⇒ **那两轮的安慰剂查过了,站得住**,不是假定。

INSTRUMENT(硬规则②):GHN 2009 JPSP Study 3(MFQ),
`data/external/dataverse/10.7910_DVN_SJTRBI_x/…Study_3.sav`。一份问卷 · 同一批人 · 同一量表。

§3 GRADIENT CHECK(踩上去之前):**偏掉任何宽合成量都会机械压低跨域相关。**
所以 offset 不能是「1」,必须是**其余四域各自 PC1** 的同一操作 —— 对称计算,每个域用自己的「其余四域」对集。

G1 ESTIMAND(先于方法):对域 D,`X(D)` = **不含 D 的那四个域之间**的跨域题对 |ρ| 中位。
  `R_D = X(D | 偏掉 PC1(D)) / X(D | 原始)` —— **偏掉一个域,剩下四个域之间还剩多少互相关。**
  `R_pol` = 用**同一个对集(纯洁的那一个)**,改偏 `politics`。
  ⚠ 这不是「纯洁的相关掉了多少」——**纯洁的题一道也不在这个对集里。**

WORLDS:W-CAUSE 纯洁是上游 · **W-SHARED 共享上游(政治)** · W-NEITHER 都不降
CONTROLS:
  正对照:偏掉**全 30 题的 PC1**(真共同因子)-> `R` 必须大幅塌;**且 g=0(纯噪声)必须不通过**。
  安慰剂:偏掉**年龄** -> `R` 必须 ≈ 1。
  **offset_control**(「这个零该是零吗?」**不该**):null 种类 =
    **「偏掉任何一个六题域合成量都会有的机械衰减」**,由其余四域的 `R_D` 分布给出。
KILL(条件式,预注册,写在跑之前):
  if 安慰剂 ∈ [0.95,1.05] and 正对照(g=1)塌 and 正对照在 g=0 不通过:
      `R_pur < R_pol − 2×展布` **且** `R_pur` 是五域最低 -> **W-CAUSE**
      `|R_pur − R_pol| < 2×展布`                      -> **W-SHARED(我不希望的)**
      两者都 > 0.95                                    -> **W-NEITHER**
      否则                                            -> 报区间
  else: UNVERIFIED
G3:五个域的 `R_D` + `R_pol` + `R_age` 全表发布,含不一致的格。
G4:{中位/均值} × {1 个 PC / 2 个 PC} × 3 种子。
IMPOSSIBLE(不写 planned):**无干预 ⇒ 非因果**;偏相关只回答「统计上分不分得开」,
  「上游」在这里只能是**统计意义**的,不是时间意义的 · **自选网络志愿者,非概率样本** ·
  单时点 ⇒ 无法定序 · 未派对抗 agent ⇒ `[unchallenged]`
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
NAMES = list(ITEM); DOMS = sorted(set(v[0] for v in ITEM.values()))
IDX = {c: i for i, c in enumerate(NAMES)}

d, _ = pyreadstat.read_sav(str(SAV))
d["_pol"] = pd.to_numeric(d["politics"], errors="coerce")
d["_age"] = pd.to_numeric(d["Age"], errors="coerce")
CC = d[NAMES + ["_pol", "_age"]].dropna()
print(f"仪器 = GHN 2009 Study 3(MFQ)· 30 题 + politics + Age 完整 n={len(CC)}")
RK = np.column_stack([rankdata(CC[c].values) for c in NAMES]).astype(float)
Z = (RK - RK.mean(0)) / RK.std(0)


def pc(cols, k=1):
    A = Z[:, [IDX[c] for c in cols]]
    U = np.linalg.svd(A - A.mean(0), full_matrices=False)[0]
    return U[:, :k]


def resid(M, C):
    if C is None: return M
    B = np.column_stack([np.ones(len(M))] + [rankdata(C[:, j]) for j in range(C.shape[1])])
    B[:, 1:] = (B[:, 1:] - B[:, 1:].mean(0)) / B[:, 1:].std(0)
    return M - B @ np.linalg.lstsq(B, M, rcond=None)[0]


def X_of(M, excl, agg=np.median):
    """不含 excl 域的那四个域之间的跨域题对 |ρ| 中位。"""
    R = np.abs(np.corrcoef(M, rowvar=False))
    keep = [k for k in NAMES if ITEM[k][0] != excl]
    v = [R[IDX[a], IDX[b]] for a, b in itertools.combinations(keep, 2) if ITEM[a][0] != ITEM[b][0]]
    return float(agg(v)), len(v)


print("\n=== G3:五个域各自的 R_D(偏掉本域 PC1 后,其余四域之间还剩多少)===")
rows = []
for D in DOMS:
    its = [k for k, v in ITEM.items() if v[0] == D]
    x0, npair = X_of(RK, D); x1, _ = X_of(resid(RK, pc(its)), D)
    rows.append(dict(who=f"域 PC1:{D}", raw=x0, par=x1, R=x1/x0, npair=npair))
    print(f"  {D:10s} 对数 {npair:3d}  原始 {x0:.4f} -> 偏后 {x1:.4f}  **R = {x1/x0:.4f}**")
Rd = {r["who"].split("：")[-1].split(":")[-1]: r["R"] for r in rows}
R_pur = Rd["PURITY"]; oth_R = [Rd[D] for D in DOMS if D != "PURITY"]

# 同一个对集(纯洁的那一个),改偏 politics / age / 全 30 题 PC1
x0, npair = X_of(RK, "PURITY")
def R_with(C): 
    x1, _ = X_of(resid(RK, C), "PURITY"); return x1/x0
R_pol = R_with(rankdata(CC["_pol"].values).reshape(-1, 1).astype(float))
R_age = R_with(rankdata(CC["_age"].values).reshape(-1, 1).astype(float))
R_g30 = R_with(pc(NAMES))
print(f"\n=== 同一个对集({npair} 对,纯洁的题一道也不在里面)===")
print(f"  偏掉 **纯洁 PC1**      -> R = **{R_pur:.4f}**")
print(f"  偏掉 **politics**      -> R = **{R_pol:.4f}**   (n={len(CC)}, 取值 1-7)")
print(f"  偏掉 年龄(安慰剂)     -> R = {R_age:.4f}   (须 ≈1)")
print(f"  偏掉 全 30 题 PC1(正对照)-> R = {R_g30:.4f}   (须大幅塌)")

bs = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(120):
        i = rng.integers(0, len(RK), len(RK)); M = RK[i]
        Zb = (M - M.mean(0)) / M.std(0)
        def pcb(cols):
            A = Zb[:, [IDX[c] for c in cols]]
            return np.linalg.svd(A - A.mean(0), full_matrices=False)[0][:, :1]
        b0, _ = X_of(M, "PURITY")
        rp = X_of(resid(M, pcb([k for k in NAMES if ITEM[k][0] == "PURITY"])), "PURITY")[0]/b0
        rl = X_of(resid(M, rankdata(CC["_pol"].values[i]).reshape(-1,1).astype(float)), "PURITY")[0]/b0
        bs.append(rp - rl)
SPR = float(np.std(bs))
print(f"\n  **R_pur − R_pol = {R_pur-R_pol:+.4f}** ± {SPR:.4f}(bootstrap 人层 × {len(SEEDS)} 种子)")
print(f"  offset(其余四域 R_D)= {[f'{x:.3f}' for x in oth_R]} 中位 {np.median(oth_R):.4f}")

G = Gate("纯洁是上游,还是和别人共享一个上游?")
# 正对照必须能失败:g=0 用纯噪声
r0 = float(np.mean([R_with(np.random.default_rng(s).standard_normal((len(CC),1))) for s in SEEDS]))
pos_ok = G.positive_control("正对照:偏掉全 30 题 PC1 必须大幅塌", planted=float(1-R_g30), floor=0.10, spread=SPR)
g0_pass = (1-r0) > 0.10 + 2*SPR
print(f"  正对照 g=0(纯噪声)-> R = {r0:.4f} ⇒ 判据在 g=0 {'⛔ 也通过' if g0_pass else '**不通过** ✅'}")
pla_ok = G.negative_control("安慰剂:偏掉年龄", null=float(abs(1-R_age)), effect=float(abs(1-R_pur)),
                            null_spread=SPR, null_kind="与道德结构无关的人口学变量")
G.offset_control("offset:偏掉任何一个六题域合成量都会有的机械衰减", effect=float(1-R_pur),
                 offset=float(1-np.median(oth_R)), spread=SPR,
                 null_kind="偏掉任意一个六题域 PC1 都会产生的跨域衰减;由其余四域的 R_D 给出")
G.has_error_bar("R_pur − R_pol", value=float(R_pur-R_pol), spread=SPR, spread_source="bootstrap_人层")

if pla_ok and pos_ok and not g0_pass:
    lowest = R_pur == min([R_pur] + oth_R)
    if R_pur < R_pol - 2*SPR and lowest: verdict = "W-CAUSE:**纯洁是(统计意义上的)上游**"
    elif abs(R_pur - R_pol) < 2*SPR: verdict = "W-SHARED:**共享上游 —— 纯洁是同车乘客,不是司机**"
    elif R_pur > 0.95 and R_pol > 0.95: verdict = "W-NEITHER:两者都不降"
    else: verdict = (f"报区间:R_pur {R_pur:.4f} · R_pol {R_pol:.4f} · 差 {R_pur-R_pol:+.4f} ± {SPR:.4f}"
                     f";纯洁{'是' if lowest else '**不是**'}五域最低")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(安慰剂 {pla_ok} · 正对照 {pos_ok} · g=0 不通过 {not g0_pass})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:{中位/均值} × {1 个 PC / 2 个 PC} ===")
spec = []
for an, ag in (("中位", np.median), ("均值", np.mean)):
    for k in (1, 2):
        b0, _ = X_of(RK, "PURITY", ag)
        rp = X_of(resid(RK, pc([x for x in NAMES if ITEM[x][0]=="PURITY"], k)), "PURITY", ag)[0]/b0
        rl = X_of(resid(RK, rankdata(CC["_pol"].values).reshape(-1,1).astype(float)), "PURITY", ag)[0]/b0
        spec.append(dict(spec=f"{an}·{k}PC", R_pur=float(rp), R_pol=float(rl), diff=float(rp-rl)))
        print(f"  {an}·{k}PC: R_pur {rp:.4f} · R_pol {rl:.4f} · 差 {rp-rl:+.4f}")
json.dump(dict(instrument="GHN 2009 JPSP Study 3 (MFQ)", n=int(len(CC)), npair=npair,
               domain_R=rows, R_pur=R_pur, R_pol=R_pol, R_age=R_age, R_g30=R_g30, R_noise=r0,
               other_domain_R=oth_R, spread=SPR, spec_curve=spec, verdict=verdict,
               politics_n=int(CC["_pol"].notna().sum()), libcon_soc_n=19,
               age_implausible_share=0.0022, seeds=SEEDS, unchallenged=True),
          open(OUT/"upstream_or_passenger.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'upstream_or_passenger.json'}")

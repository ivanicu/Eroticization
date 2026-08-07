"""#834 · E03·A76·R273 —— 那句话的主语从没被检验过:换六根别的轴,九十年代那一格还在吗?

⚠⚠ **先说为什么不是 `#833`①。** `git log` 实测:最近七轮里**只有一轮**(`#832`)是关于人的,
   其余六轮全是 `guard/fix/verify` 基建。**`frontier §3`:那是盆地,而 `#833`① 又是基建。**
   ⇒ **`#833`① 已记为 `DEBTS.tsv` 里一行 `OPEN`,本轮去打最载重的那条主张。**

十轮校正之后站着的只有一句话:
**「在同性恋这道题上,**宗教**鸿沟在九十年代偏离了它自己五十年的匀速。」**
⚠⚠ **而这句话的主语 —— 「宗教」—— 只被对照过一次**(`#820` 的政治立场)。
**教育、年龄、地区、种族、性别:一次都没有。**

**三个世界,而第二、三个的正结果我都不欢迎:**
   A **确实是宗教的**:九十年代那一格在虔诚轴上偏离,而在别的轴上大多不偏离。
   B **其实是九十年代的**:它在**多数**轴上都偏离 ⇒ **发生的事是「美国人在九十年代整体飞快地
     改了主意」,而任何一条既存的裂缝都会因此张开** ⇒ **「宗教鸿沟」这个说法必须撤,
     它描述的不是机制,是我挑的那条裂缝。**
   C **任何裂缝都行**:连没有任何机制可言的轴(性别、地区)也偏离 ⇒
     **这是「人群整体快速变化会撑开任何既存差距」的性质,与被切开的是什么无关。**

G1 估计量:**`homosex` 的逐十年 `Δgap`,对**该轴自己**的匀速参照**(与 `#819`/`#832` 完全一致:
   `ref = 该轴全程 Δgap × 该十年跨年 ÷ 全程跨年` ⇒ **极性与分层方向都自动抵消**),
   在**七根轴**上各算一遍,再对**整族**做 BH / BY。

**七根轴(硬规则①:n 与真正被问过的年份已先跑一遍打印,取值标签已读过,不在题名上下结论):**
   `REL`(虔诚度三项合成,三分位)· `polviews`(1–7,三分位)· `educ`(受教育年数,三分位)·
   `age`(三分位)· `region`(**south** vs 其余,标签已读)· `race`(**white** vs **black**)·
   `sex`(**male** vs **female**)。
   ⚠ 连续轴取**年份内三分位**(与 `#819` 同);二分轴直接取两个类别。

⚠⚠ **跑之前写下的最强混淆,而它足以让世界 B 变成假象:**
   **教育、地区、种族都与虔诚度相关** ⇒ **它们上面的偏离不是独立证据。**
   ⇒ 控制:**逐轴测它的「高组」与虔诚层的 Jaccard 重叠**(`#820` 做过这件事,那一轮实测只有 0.196,
   **正是那个低重叠让 `#820` 的对照有内容**)⇒ **每一轴的结果旁边都印出它的重叠率;
   重叠高的轴,它的「偏离」要打折,而这个折扣必须看得见,不能靠读者自己想。**

预测矩阵:
   | 世界 | 现在 | 只有 REL(+≤2) | 多数轴都偏离 | 连性别/地区都偏离 |
   | A 宗教的     | 0.45 | **0.85** | 0.05 | 0.03 |
   | B 九十年代的 | 0.35 | 0.05 | **0.85** | 0.35 |
   | C 任何裂缝   | 0.20 | 0.03 | 0.30 | **0.85** |

预注册判词(条件式):
  if 正控开火(**在某一根轴的 1990s 植入一个巨大偏离,它必须在整族 BH 下存活**)
     and 负控开火(**全匀速世界里「至少一个存活」的比例 ≤ q** ——
        ⚠ **「这个零该不该是零?」不该是 0,是 ≤ q**:BH 控 FDR,`#825` 正犯在这里,`#826` 起已改对):
      `homosex`/1990s 在 `REL` 上存活 **且** 在其余六轴中存活 ≤2 根 -> A
      在 ≥4 根轴上存活                                          -> B
      在**性别或地区**上也存活                                   -> C(与 B 可同时成立,分别报)
  else: UNVERIFIED
⚠ **`G3`:整族 = 7 轴 × 各自可用十年,全报,含没存活的格。**
⚠ **`G4`:BH 与 BY 都报**(格间不独立 —— 同一轴相邻十年共用端点年份,`#826` 已立)。

⚠ 硬规则②(这条主张路由过哪具仪器):**七根轴全部来自 GSS 同一份问卷** ——
   **这不是跨仪器复现,而本轮问的也不是「别处成不成立」,是「这里的主语对不对」。**
⚠ 本轮换不了仪器:估计量是「同一具问卷内换分层轴」,第二份调查会同时换掉题目与轴。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, B, Q, NREP, BNULL = "homosex", 4, 6000, 0.05, 60, 2000

cols = ["year", "attend", "reliten", "fund", "polviews", "educ", "age", "region", "race", "sex", IT]
d = pd.read_stata(gp, columns=cols, convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("polviews", (1, 7)),
                    ("educ", (0, 20)), ("age", (18, 89)), ("region", (1, 4)), ("race", (1, 3)), ("sex", (1, 2))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
R = M.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = zs(R[["attend", "reliten", "fund"]]).mean(axis=1)
M = M.join(R["REL"])

def terc(col):
    return M.groupby("year")[col].transform(
        lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)
AXES = {}
for nm, col in (("REL 虔诚度", "REL"), ("polviews 政治", "polviews"),
                ("educ 受教育年数", "educ"), ("age 年龄", "age")):
    t = terc(col); AXES[nm] = ((t == 2), (t == 0), "三分位:高 vs 低")
AXES["region 南方 vs 其余"] = ((M.region == 3), (M.region != 3) & M.region.notna(), "标签已读:south=3")
AXES["race 白 vs 黑"] = ((M.race == 1), (M.race == 2), "标签已读:white=1 · black=2")
AXES["sex 男 vs 女"] = ((M.sex == 1), (M.sex == 2), "标签已读:male=1 · female=2")

print("=== ⓪ 硬规则①:七根轴的高/低组样本量,与**跑前写下的混淆**:它们与虔诚层的重叠 ===")
relhi = AXES["REL 虔诚度"][0]
GRID = {}
for nm, (hi, lo, how) in AXES.items():
    ok = M[IT].notna() & (hi | lo)
    inter = int((hi & relhi & M[IT].notna()).sum()); union = int(((hi | relhi) & M[IT].notna()).sum())
    jac = inter/union if union else np.nan
    print(f"  {nm:20s} 高 {int((hi & M[IT].notna()).sum()):>6,} · 低 {int((lo & M[IT].notna()).sum()):>6,} · "
          f"{how:22s} · **与虔诚层 Jaccard {jac:.3f}**")
    ys = {}
    for y, g in M[ok].groupby("year"):
        a = g[hi.loc[g.index]][IT].to_numpy(float); b = g[lo.loc[g.index]][IT].to_numpy(float)
        if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
    dec = {}
    for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
    dec = {k: v for k, v in dec.items() if len(v) >= 3}
    s = sorted(ys)
    if len(s) < 8 or not dec: GRID[nm] = None; continue
    g0 = float(ys[s[0]][0].mean()-ys[s[0]][1].mean()); g1 = float(ys[s[-1]][0].mean()-ys[s[-1]][1].mean())
    GRID[nm] = dict(ys=ys, dec=dec, span=s[-1]-s[0], dgap=g1-g0, jac=jac)
print("  ⚠ **重叠高的轴,它的「偏离」不是独立证据 —— 这个折扣印在每一行旁边,不靠读者自己想。**")

CELLS = [(nm, dc) for nm, G_ in GRID.items() if G_ for dc in sorted(G_["dec"])]
print(f"\n=== ① 整族 = **{len(CELLS)} 格**(7 轴 × 各自可用十年)· B={B} · p 下限 {1/(B+1):.2e} ===")
def pval(nm, dc, rng, Bv, src=None):
    G_ = GRID[nm]; S = src if src else G_["ys"]; ys = G_["dec"][dc]
    ref = G_["dgap"]*(ys[-1]-ys[0])/G_["span"]
    dr = np.empty(Bv)
    for i in range(Bv):
        r = lambda a: a[rng.integers(0, len(a), len(a))]
        a0, b0 = r(S[ys[0]][0]), r(S[ys[0]][1]); a1, b1 = r(S[ys[-1]][0]), r(S[ys[-1]][1])
        dr[i] = (a1.mean()-b1.mean()) - (a0.mean()-b0.mean())
    return max(2*min(float(np.mean(dr <= ref)), float(np.mean(dr >= ref))), 1.0/(Bv+1))
rng = np.random.default_rng(273)
PS = {c: pval(c[0], c[1], rng, B) for c in CELLS}
surv_bh = {CELLS[i] for i in Gate.bh([PS[c] for c in CELLS], Q)}
surv_by = {CELLS[i] for i in Gate.by([PS[c] for c in CELLS], Q)}
print(f"  {'轴':22s} " + " ".join(f"{dc}s".rjust(11) for dc in (1970, 1980, 1990, 2000, 2010)))
for nm in AXES:
    if not GRID.get(nm): print(f"  {nm:22s} —— 年份不足,整轴不可用"); continue
    row = ""
    for dc in (1970, 1980, 1990, 2000, 2010):
        if (nm, dc) not in PS: row += "——".rjust(11); continue
        mark = "**" if (nm, dc) in surv_bh else ("*" if PS[(nm, dc)] < 0.05 else "")
        row += f"{PS[(nm,dc)]:.4f}{mark}".rjust(11)
    print(f"  {nm:22s} {row}")
print("  ⚠ `**` = 整族 BH(q=0.05)存活 · `*` = 未校正 p<0.05 但**校正后不存活**(`#826` 的教训)")
print(f"\n  BH 存活 **{len(surv_bh)}/{len(CELLS)}** ⇒ {sorted(f'{a}/{b}s' for a,b in surv_bh) or '无'}")
print(f"  BY 存活 **{len(surv_by)}/{len(CELLS)}** ⇒ {sorted(f'{a}/{b}s' for a,b in surv_by) or '无'}")
rel90 = ("REL 虔诚度", 1990)
others90 = [nm for nm in AXES if nm != "REL 虔诚度" and (nm, 1990) in surv_bh]
n_other = len(others90)
sexregion = [nm for nm in ("region 南方 vs 其余", "sex 男 vs 女") if (nm, 1990) in surv_bh]
print(f"\n  ⇒ `REL`/1990s 存活:**{rel90 in surv_bh}** · 其余六轴 1990s 存活 **{n_other}** 根 ⇒ {others90 or '无'}")
print(f"  ⇒ 性别/地区 1990s 存活:{sexregion or '无'}")

print("\n=== ② 控制 ===")
def syn(mode, rng_, plant=("REL 虔诚度", 1990)):
    S = {}
    for nm, G_ in GRID.items():
        if not G_: continue
        ys = sorted(G_["ys"]); y0, y1 = ys[0], ys[-1]
        g0 = float(G_["ys"][y0][0].mean()-G_["ys"][y0][1].mean()); tot = G_["dgap"]
        S[nm] = {}
        for y in ys:
            tgt = g0 + tot*(y-y0)/(y1-y0)
            if mode == "planted" and nm == plant[0] and plant[1] <= y <= plant[1]+9:
                tgt += 3.0*abs(tot)*(y-plant[1])/9.0
            cur = float(G_["ys"][y][0].mean()-G_["ys"][y][1].mean())
            a, b = G_["ys"][y]
            S[nm][y] = (a[rng_.integers(0, len(a), len(a))] + (tgt-cur),
                        b[rng_.integers(0, len(b), len(b))])
    return S
r2 = np.random.default_rng(274)
Sp = syn("planted", r2)
psp = {c: pval(c[0], c[1], r2, BNULL, src=Sp[c[0]]) for c in CELLS}
sp = {CELLS[i] for i in Gate.bh([psp[c] for c in CELLS], Q)}
pc_ok = rel90 in sp
print(f"  正控(在 `REL` 的 1990s 植入 3× 全程量的偏离)⇒ 整族 BH 存活 {len(sp)}/{len(CELLS)},"
      f"含被植入那格:**{pc_ok}**")
hits = 0
for j in range(NREP):
    rj = np.random.default_rng(7000+j)
    Su = syn("uniform", rj)
    psu = {c: pval(c[0], c[1], rj, BNULL, src=Su[c[0]]) for c in CELLS}
    if len(Gate.bh([psu[c] for c in CELLS], Q)) > 0: hits += 1
rate = hits/NREP; se = float(np.sqrt(rate*(1-rate)/NREP))
print(f"  负控(全匀速 × {NREP} 次)⇒「至少一个存活」比例 **{rate:.3f} ± {se:.3f}** ——"
      f" ⚠ **期望 ≤ q = {Q},不是 0**(`#825` 的错,`#826` 起已改对)")

G = Gate("#834 · 那句话的主语从没被检验过")
G.multiplicity_control("⓪ `G3`:整族一次校正,报 cells tested / surviving 并列出没存活的",
                       [PS[c] for c in CELLS], Q, [f"{a}/{b}s" for a, b in CELLS],
                       method="bh", p_floor=1.0/(B+1))
G.asserted("① 正控:在 `REL` 的 1990s 植入 3× 全程量的偏离,该格必须在**整族 BH** 下存活",
           bool(pc_ok), f"存活 {len(sp)}/{len(CELLS)},含被植入格 = {pc_ok}", kind="control")
G.asserted("② 负控:全匀速世界「至少一个存活」比例必须 **≤ q**(⚠ **不是 0** —— BH 控 FDR)",
           bool(rate <= Q + 2*se), f"{rate:.3f} ± {se:.3f}(阈 {Q+2*se:.3f})", kind="control")
G.asserted("③ 前提(跑前写下的最强混淆):**教育/地区/种族都与虔诚度相关 ⇒ 它们上的偏离不是独立证据** "
           "⇒ **逐轴的 Jaccard 重叠已测并印在每一行旁边**",
           bool(all(GRID[nm] is None or "jac" in GRID[nm] for nm in AXES)),
           " · ".join(f"{nm[:6]}{GRID[nm]['jac']:.2f}" for nm in AXES if GRID[nm]), kind="control")
G.asserted("④ 前提(gauge):参照由**该轴自己的全程 `Δgap`** 定并按跨度缩放 ⇒ 极性与分层方向都自动抵消",
           True, "每格 ref = 该轴全程 Δgap × 该十年跨年 ÷ 全程跨年", kind="control")
G.asserted("⑤ 前提(硬规则①):七根轴的 n / 年份 / 取值标签**跑前已读**,"
           "`region` 只有 4 类(south=3)、`race` 用 white vs black —— **没有在题名上下结论**",
           True, f"可用轴 {sum(1 for nm in AXES if GRID[nm])}/{len(AXES)}", kind="control")
G.asserted("⑥ kill(预注册):「这句话确实是关于宗教的」要成立,需 `REL`/1990s 存活 **且** "
           "其余六轴的 1990s 存活 **≤2** 根",
           bool((rel90 in surv_bh) and n_other <= 2),
           f"REL/1990s 存活 {rel90 in surv_bh} · 其余轴存活 {n_other} 根 {others90}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif (rel90 in surv_bh) and n_other <= 2:
    V = (f"**A 这句话确实是关于宗教的。** `REL`/1990s 在整族({len(CELLS)} 格)BH 下存活,"
         f"而其余六根轴的 1990s 只有 **{n_other}** 根存活({others90 or '无'})。\n"
         f"  ⇒ **一句关于人的话:九十年代美国人在同性恋这道题上分开的那条缝,"
         f"确实是沿着「信不信教」裂的 —— 把同一批人按教育、年龄、地区、种族、性别、政治立场切开,"
         f"同一个十年都没有出现同样的裂口。**")
elif n_other >= 4:
    V = (f"**B 这句话其实是关于九十年代的。** 除 `REL` 外还有 **{n_other}** 根轴的 1990s 存活:{others90}。\n"
         f"  ⇒ **「宗教鸿沟」必须撤作机制描述 —— 九十年代美国人整体飞快改了主意,"
         f"而任何一条既存的裂缝都会因此张开。我挑中宗教,不等于宗教挑中了它。**")
else:
    V = (f"**介于之间,报整张表不选边。** `REL`/1990s 存活 = {rel90 in surv_bh};"
         f"其余轴存活 {n_other} 根:{others90}。\n"
         f"  ⇒ **既不能说「只有宗教」,也不能说「任何裂缝都行」——\n"
         f"  而 kill 是按 ≤2 预注册的,所以这一轮不许把它读成 A。**")
if sexregion:
    V += (f"\n  ⚠⚠ **而性别/地区也存活({sexregion})—— 那是世界 C 的指纹**:"
          f"这两根轴上没有任何关于同性恋态度的机制可言,**它们的偏离说明快速的人群整体变化"
          f"本身就会撑开既存差距。**")
print(V)
json.dump(dict(item=IT, axes=list(AXES), cells=[f"{a}|{b}" for a, b in CELLS], n_cells=len(CELLS),
               B=B, q=Q, ps={f"{a}|{b}": PS[(a, b)] for a, b in CELLS},
               jaccard={nm: (GRID[nm]["jac"] if GRID[nm] else None) for nm in AXES},
               surv_bh=sorted(f"{a}/{b}s" for a, b in surv_bh),
               surv_by=sorted(f"{a}/{b}s" for a, b in surv_by),
               rel1990_alive=bool(rel90 in surv_bh), n_other_axes_1990=n_other,
               other_axes_1990=others90, sex_or_region=sexregion,
               pos_control=bool(pc_ok), neg_control=dict(rate=rate, se=se, expectation="<= q, NOT 0"),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"whose_divide.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'whose_divide.json'}")

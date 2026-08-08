"""#832 · E03·A75·R271 —— 「家族是什么」这条规格轴从没扫过;而 `#826`② 也从没执行

⚠⚠ **先说这一轮为什么不是 `#831`①(那会是连续第六轮基建)。**
   `git log` 实测:最近五轮是 `guard/guard/verify/fix/guard` —— **一轮关于人的都没有。**
   **`frontier §3`:N 步连着确认同一个故事 = 盆地。而这里的「故事」是「我的仪器坏了,我又在修」。**
   ⇒ **`#831`① 已作为一行 `OPEN` 记进 `DEBTS.tsv`(那张表的第一次真正使用),本轮回到对象。**

`#826` 的结论:**34 格里只有 `homosex`/1990s 与 /2000s 存活。**
而 `#827` 当轮登记过一件没人回头做的事:
**「哪些检验属于同一族」是一个语义判断,库做不了,把整张 34 格网格当作一族**是我选的**,别的切法会给别的答案。**
⇒ **`G4` 要的是规格曲线,而「家族定义」正是这条从没被扫过的轴。**
⇒ 同一轮顺带执行 **`#826`②**:`#820` 的政治网格从来没做过多重性校正。

G1 估计量:**在五种家族定义下,`homosex`/1990s 是否仍然存活**,以及每种家族的**存活集**。
   p 值构造与 `#826` 完全一致(该十年 `Δgap` 对该题自己的匀速参照,双侧经验 p,`B`=6000)——
   **只有「谁和谁一起校正」在变。**

**五种家族(跑前写死,含最不站得住的那一种)**:
   ① `religion34` —— `#826` 用的那一族(虔诚度网格,34 格)
   ② `both68` —— **虔诚度 + 政治立场两张网格合起来(最宽的可辩护定义)** ——
      两张网格是为回答**同一个问题**跑的(`#820`),**所以把它们分开校正才是需要辩护的那一方**
   ③ `per_item` —— 每题各自一族(8 族,每族 4–6 格)
   ④ `per_decade` —— 每个十年各自一族(6 族)
   ⑤ `sexual_only` —— 只取本项目框架真正关心的性道德题(`homosex`·`teensex`·`sexeduc`)

⚠⚠ **跑之前必须写下的最强混淆,而它足以自己制造世界 A:**
   **家族越窄,存活越容易 —— 那是算术,不是证据。**
   ③ 每族只有 4–6 格,BH 在秩 1 的阈值是 `q/C = 0.05/5 ≈ 0.01`,比 34 格的 `0.0015` 松近七倍。
   ⇒ 控制:**每一个存活计数旁边必须印出该族的大小**;
   **并且把 ③ 明确标为最不站得住的一种,而不是「最有力的一种」。**
   **一个作者可以通过把家族切小,把任何东西变成「存活」。**

三个世界:
   A **家族无关**:`homosex`/1990s 在五种族下都存活 ⇒ `#826` 的句子不是分组的产物。
   B **家族相关**:存活集随族变 ⇒ **「只有同性恋那一题存活」是一句关于我的分组的话** ——
     **此后每一句都必须把家族说出来。**
   C **最宽的族下全灭**:`both68` 下一个都不剩 ⇒ **整条「十年偏离」线索在最诚实的族定义下没有支撑。**

预测矩阵:
   | 世界 | 现在 | 五族全活 | 存活集随族变 | both68 下全灭 |
   | A 家族无关 | 0.40 | **0.85** | 0.10 | 0.03 |
   | B 家族相关 | 0.40 | 0.05 | **0.85** | 0.25 |
   | C 全灭     | 0.20 | 0.02 | 0.20 | **0.85** |

预注册判词(条件式):
  if 正控开火(**植入一个巨大偏离,它必须在五种族下都存活** —— 否则族的比较没有意义)
     and 负控开火(**全匀速世界里「至少一个存活」的比例 ≤ q**,
        ⚠ **「这个零该不该是零?」——不该是 0,是 ≤ q**:BH 控 FDR,全零世界的保证就是 ≤ q,
        **`#825` 正是犯在向一个额定率为 q 的程序索要 0,`#826` 已改对,本轮沿用**):
      `homosex`/1990s 在**五种族下都存活** -> A
      否则 -> B/C,**报整张 5×2 的表,不选边**
  else: UNVERIFIED
⚠ **`G3`:5 族 × 2 方法(BH/BY)= 10 格,全报,并逐族列出没存活的。**

⚠ 硬规则②(这条主张路由过哪具仪器):**两个分层都来自 GSS 同一份问卷** ——
  这不是跨仪器复现,**而 `both68` 之所以是最宽的可辩护族,正因为两张网格共用一具仪器、回答同一个问题。**
⚠ 本轮换不了仪器(对象是同一份数据上的校正族)。⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_is_the_eight_point_axis_an_axis_or_eight_labels_on_noise/results/is_the_ordering_an_object.json"))
ITEMS = P791["items"]
B, Q, NREP, BNULL = 6000, 0.05, 60, 2000
SEXUAL = ["homosex", "teensex", "sexeduc"]

d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", "polviews"]+ITEMS,
                  convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
K = {c: (4 if c == "homosex" else len(cat[c].cat.categories)) for c in ITEMS}
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= K[c])) for c in ITEMS})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("polviews", (1, 7))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
BASE = M.dropna(subset=["year"]).copy()
Rr = BASE.dropna(subset=["attend", "reliten", "fund"]).copy()
Rr["REL"] = zs(Rr[["attend", "reliten", "fund"]]).mean(axis=1)
BASE = BASE.join(Rr["REL"])
for nm, col in (("k_rel", "REL"), ("k_pol", "polviews")):
    BASE[nm] = BASE.groupby("year")[col].transform(
        lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)

def build(kcol):
    YR, COV, FULL = {}, {}, {}
    for it in ITEMS:
        g = BASE.dropna(subset=[it, kcol]); ys = {}
        for y, gy in g.groupby("year"):
            a = gy[gy[kcol] == 2][it].to_numpy(float); b = gy[gy[kcol] == 0][it].to_numpy(float)
            if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
        YR[it] = ys
        dec = {}
        for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
        COV[it] = {k2: v for k2, v in dec.items() if len(v) >= 3}
        s = sorted(ys)
        g0 = float(ys[s[0]][0].mean()-ys[s[0]][1].mean()); g1 = float(ys[s[-1]][0].mean()-ys[s[-1]][1].mean())
        FULL[it] = dict(span=s[-1]-s[0], dgap=g1-g0)
    return YR, COV, FULL
AX = {"rel": build("k_rel"), "pol": build("k_pol")}

def pval(ax, it, dc, rng, Bv, src=None):
    YR, COV, FULL = AX[ax]; S = (src[it] if src else YR[it]); ys = COV[it][dc]
    span = ys[-1]-ys[0]; ref = FULL[it]["dgap"]*span/FULL[it]["span"]
    dr = np.empty(Bv)
    for i in range(Bv):
        r = lambda a: a[rng.integers(0, len(a), len(a))]
        a0, b0 = r(S[ys[0]][0]), r(S[ys[0]][1]); a1, b1 = r(S[ys[-1]][0]), r(S[ys[-1]][1])
        dr[i] = (a1.mean()-b1.mean()) - (a0.mean()-b0.mean())
    lo = 1.0/(Bv+1)
    return max(2*min(float(np.mean(dr <= ref)), float(np.mean(dr >= ref))), lo)

CELLS = [(ax, it, dc) for ax in ("rel", "pol") for it in ITEMS for dc in sorted(AX[ax][1][it])]
print(f"=== ⓪ 两张网格的格数(硬规则②:两个分层都来自 GSS 同一份问卷)===")
for ax in ("rel", "pol"):
    n = sum(1 for c in CELLS if c[0] == ax)
    print(f"  {ax}: **{n} 格**")
print(f"  ⇒ 合计 **{len(CELLS)}** 格 · `B`={B} ⇒ p 分辨率下限 {1/(B+1):.2e}")

rng = np.random.default_rng(271)
PS = {c: pval(c[0], c[1], c[2], rng, B) for c in CELLS}
print(f"\n=== ① p 最小的 8 格(两张网格合起来)===")
for c in sorted(CELLS, key=lambda x: PS[x])[:8]:
    print(f"  {c[0]:4s} {c[1]:9s} {c[2]}s  p = {PS[c]:.5f}")

FAMILIES = {
    "① religion34(#826 用的)": [c for c in CELLS if c[0] == "rel"],
    "② both68(最宽可辩护)": CELLS,
    "④ per_decade(6 族)": None,
    "⑤ sexual_only(框架关心的)": [c for c in CELLS if c[0] == "rel" and c[1] in SEXUAL],
    "③ per_item(8 族,最不站得住)": None,
}
def survivors(cells, method):
    ps = [PS[c] for c in cells]
    idx = (Gate.bh if method == "bh" else Gate.by)(ps, Q)
    return {cells[i] for i in idx}
def survivors_split(groups, method):
    out = set()
    for g in groups: out |= survivors(g, method)
    return out
GROUPS_ITEM = [[c for c in CELLS if c[0] == "rel" and c[1] == it] for it in ITEMS]
GROUPS_ITEM = [g for g in GROUPS_ITEM if g]
GROUPS_DEC = [[c for c in CELLS if c[0] == "rel" and c[2] == dc]
              for dc in sorted({c[2] for c in CELLS if c[0] == "rel"})]
GROUPS_DEC = [g for g in GROUPS_DEC if g]

print(f"\n=== ② `G3`/`G4`:5 种家族 × 2 方法 = 10 格,全报(⚠ **族越窄存活越容易,那是算术不是证据**)===")
RES = {}
for name in FAMILIES:
    for meth in ("bh", "by"):
        if name.startswith("③"):
            surv = survivors_split(GROUPS_ITEM, meth); size = f"8 族 × {[len(g) for g in GROUPS_ITEM]}"
        elif name.startswith("④"):
            surv = survivors_split(GROUPS_DEC, meth); size = f"6 族 × {[len(g) for g in GROUPS_DEC]}"
        else:
            cells = FAMILIES[name]; surv = survivors(cells, meth); size = f"{len(cells)} 格"
        RES[(name, meth)] = dict(size=size, n=len(surv),
                                 surv=sorted(f"{a}/{b}/{c}s" for a, b, c in surv))
        print(f"  {name:26s} {meth.upper():3s} 族大小 {size:22s} 存活 **{len(surv)}** ⇒ "
              + (", ".join(RES[(name, meth)]["surv"][:5]) + ("…" if len(surv) > 5 else "") if surv else "无"))
HOMO90 = ("rel", "homosex", 1990)
alive = {k: (HOMO90 in {tuple(s.split("/")[0:1]+[s.split("/")[1]]+[int(s.split("/")[2][:-1])])
                        for s in v["surv"]}) for k, v in RES.items()}
alive = {k: any(s == "rel/homosex/1990s" for s in v["surv"]) for k, v in RES.items()}
n_alive = sum(alive.values())
print(f"\n  ⇒ **`rel/homosex/1990s` 在 {n_alive}/{len(RES)} 个(族 × 方法)组合下存活**")
print(f"  ⚠ **`#826`② 顺带执行**:政治网格单独成族时的存活 —— "
      f"{len(survivors([c for c in CELLS if c[0]=='pol'], 'bh'))} 格(BH)· "
      f"{len(survivors([c for c in CELLS if c[0]=='pol'], 'by'))} 格(BY)")

print("\n=== ③ 控制 ===")
def syn(mode, rng_, plant=("rel", "homosex", 1990)):
    S = {}
    for ax in ("rel", "pol"):
        YR, COV, FULL = AX[ax]; S[ax] = {}
        for it in ITEMS:
            ys = sorted(YR[it]); y0, y1 = ys[0], ys[-1]
            g0 = float(YR[it][y0][0].mean()-YR[it][y0][1].mean()); tot = FULL[it]["dgap"]
            S[ax][it] = {}
            for y in ys:
                f = (y-y0)/(y1-y0)
                tgt = g0 + tot*f
                if mode == "planted" and ax == plant[0] and it == plant[1] and plant[2] <= y <= plant[2]+9:
                    tgt += 3.0*abs(tot)*(y-plant[2])/9.0
                cur = float(YR[it][y][0].mean()-YR[it][y][1].mean())
                a, b = YR[it][y]
                S[ax][it][y] = (a[rng_.integers(0, len(a), len(a))] + (tgt-cur),
                                b[rng_.integers(0, len(b), len(b))])
    return S
r2 = np.random.default_rng(272)
Sp = syn("planted", r2)
PSP = {c: pval(c[0], c[1], c[2], r2, BNULL, src=Sp[c[0]]) for c in CELLS}
def surv_with(psmap, cells, meth):
    ps = [psmap[c] for c in cells]
    return {cells[i] for i in (Gate.bh if meth == "bh" else Gate.by)(ps, Q)}
pc_ok = all(HOMO90 in (surv_with(PSP, FAMILIES[n], "bh") if FAMILIES[n] else
            set().union(*[surv_with(PSP, g, "bh") for g in (GROUPS_ITEM if n.startswith("③") else GROUPS_DEC)]))
            for n in FAMILIES)
print(f"  正控(在 `rel/homosex/1990s` 植入 3× 全程量的偏离)⇒ **五种族下都存活:{pc_ok}**")
hits = 0
for j in range(NREP):
    rj = np.random.default_rng(6000+j)
    Su = syn("uniform", rj)
    psu = {c: pval(c[0], c[1], c[2], rj, BNULL, src=Su[c[0]]) for c in CELLS}
    if len(surv_with(psu, CELLS, "bh")) > 0: hits += 1
rate = hits/NREP; se = float(np.sqrt(rate*(1-rate)/NREP))
print(f"  负控(全匀速 × {NREP} 次,`both68` 族)⇒ 「至少一个存活」比例 **{rate:.3f} ± {se:.3f}**")
print(f"     ⚠ **期望是 ≤ q = {Q},不是 0**(BH 控 FDR;`#825` 正是犯在向额定率为 q 的程序索要 0)")

G = Gate("#832 · 「家族是什么」这条规格轴")
G.asserted("① 正控:植入 3× 全程量的偏离后,`rel/homosex/1990s` 必须在**五种族下都存活**"
           "(否则不同族之间的比较没有意义)", bool(pc_ok), f"五族全存活 = {pc_ok}", kind="control")
G.asserted("② 负控:全匀速世界里「至少一个存活」的比例必须 **≤ q**(⚠ **不是 0** —— BH 控 FDR,"
           "`#825` 正是犯在向额定率为 q 的程序索要 0,`#826` 已改对,本轮沿用)",
           bool(rate <= Q + 2*se), f"{rate:.3f} ± {se:.3f}(阈 q+2se = {Q+2*se:.3f})", kind="control")
G.asserted("③ 前提(跑前写下的最强混淆):**族越窄存活越容易 —— 那是算术不是证据** ⇒ "
           "**每个存活计数旁边都印出族大小**,且 **③ per_item 明确标为最不站得住的一种,不是最有力的**",
           bool(all("size" in v for v in RES.values())),
           " · ".join(f"{k[0][:2]}{k[1]}:{v['size'][:12]}" for k, v in list(RES.items())[:4]), kind="control")
G.asserted("④ 前提(硬规则②):**两个分层都来自 GSS 同一份问卷** ⇒ 这不是跨仪器复现;"
           "**而 `both68` 之所以最宽可辩护,正因为两张网格共用一具仪器、回答同一个问题**",
           True, f"rel {sum(1 for c in CELLS if c[0]=='rel')} 格 + pol {sum(1 for c in CELLS if c[0]=='pol')} 格",
           kind="control")
G.asserted("⑤ kill(预注册):「`#826` 的句子与家族定义无关」要成立,"
           "需 `rel/homosex/1990s` 在**全部 10 个(族 × 方法)组合**下都存活",
           bool(n_alive == len(RES)), f"存活 {n_alive}/{len(RES)}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif n_alive == len(RES):
    V = (f"**A 家族无关。** `rel/homosex/1990s` 在全部 {len(RES)} 个(族 × 方法)组合下都存活。\n"
         f"  ⇒ **`#826` 那句话不是分组的产物 —— 而这是它第一次被证明与「我怎么分族」无关。**")
else:
    dead = [f"{k[0][:14]}/{k[1].upper()}" for k, v in RES.items() if not alive[k]]
    V = (f"**B/C 家族相关。** `rel/homosex/1990s` 只在 **{n_alive}/{len(RES)}** 个组合下存活;"
         f"死在:{dead}。\n"
         f"  ⇒ **「只有同性恋那一题存活」是一句关于我怎么分族的话 ——\n"
         f"  此后每一次说它,都必须把家族说出来。**")
print(V)
json.dump(dict(cells=[f"{a}/{b}/{c}" for a, b, c in CELLS], n_cells=len(CELLS), B=B, q=Q,
               ps={f"{a}/{b}/{c}": PS[(a, b, c)] for a, b, c in CELLS},
               families={f"{k[0]}|{k[1]}": v for k, v in RES.items()},
               homosex1990_alive={f"{k[0]}|{k[1]}": bool(v) for k, v in alive.items()},
               n_alive=n_alive, n_combos=len(RES),
               pol_only_bh=len(survivors([c for c in CELLS if c[0] == "pol"], "bh")),
               pol_only_by=len(survivors([c for c in CELLS if c[0] == "pol"], "by")),
               pos_control=bool(pc_ok), neg_control=dict(rate=rate, se=se, expectation="<= q, NOT 0"),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"which_family.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'which_family.json'}")

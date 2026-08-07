"""#846 · E03·A84·R285 —— 「是世俗那边走了」是关于人的话,还是只关于这一道题?

`#838` 的 NEXT 自己写下了本轮:**「不要把『世俗那边走了』外推到别的题 —— 本轮只跑了 `homosex`。」**
⇒ 这是**一句关于人的主张的作用域检验**,而作用域正是 `§2` 说的、
**十二条撤回里十一条的死因:一个正确的数字,报的时候没带它成立的范围。**

**两个作用域,而它们是两种完全不同的心理学:**
 **S1 关于人**:世俗那一层在九十年代**普遍**比自己的步子快 ——
   那是「一群人整体加速」,而同性恋只是他们加速时经过的一道题。
 **S2 关于题**:只有在同性恋这道题上是世俗层在动 ——
   那么九十年代发生的事**属于这个题目**,不属于那群人。

G1 估计量:**八道题各自的「九十年代,哪一层的偏离排除零」** ——
   与 `#838` 完全同一构造(每层对**它自己**的全程匀速参照),只换题目。
   ⚠ **参照由每题每层自己的全程位移定 ⇒ 题目极性自动抵消**(`#789` 的规范检查)。

三个世界:
   A **S1 关于人**:世俗层在**多数**题上都排零,虔诚层在多数题上不排 ⇒ **一群人在加速。**
   B **S2 关于题**:只有 `homosex`(或极少数题)是这个形状 ⇒ **那件事属于题目。**
   C **两层在不同题上各自排零,没有稳定的「哪一层」** ⇒ **「谁在动」这个问法在八题上不成立**
     —— **元分离器:它说的是我的分解方式只在一道题上有意义。**

预测矩阵:
   | 世界 | 现在 | 世俗多数排零 | 只 homosex | 无稳定模式 |
   | A 关于人 | 0.35 | **0.85** | 0.05 | 0.10 |
   | B 关于题 | 0.40 | 0.05 | **0.85** | 0.10 |
   | C 问法不成立 | 0.25 | 0.10 | 0.10 | **0.80** |

预注册判词(条件式):
  if 正控开火(**只往某一题某一层植入位移,只有那一格该动**)
     and 负控开火(**全匀速的世界里所有格都必须落在零上**):
      世俗层排零的题数 ≥ 5/8 且 > 虔诚层 -> A
      世俗层排零的题数 ≤ 2/8            -> B
      其他                                -> C
  else: UNVERIFIED
⚠ `G3`:整族 = 8 题 × 2 层 = **16 格**,BH 与 BY 都做,**族大小印在旁边**
  (`#832`:族越窄存活越易,是算术不是证据)。

⚠ 跑之前写下的最强混淆:**八道题的天花板余量差别很大** ——
  一层若已贴边,它「没动」会被天花板伪造出来(`#838` 量过 `homosex`,别的题没量过)。
  ⇒ 控制:**逐题逐层印出九十年代首尾均值与到端点的余量**;
  **余量占用 >60% 的格,它的「排零」标注为 `CEILING?`,不计入判词。**
⚠ 本轮换不了仪器(GSS);而它不需要 —— 本轮换的是**题目**,这正是作用域检验要的那根轴。
"""
import json, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
ITEMS = {"homosex": 4, "premarsx": 4, "teensex": 4, "xmarsex": 4,
         "sexeduc": 2, "racmar": 2, "prayer": 2, "suicide2": 2}
B, Q, DEC = 4000, 0.05, 1990

cols = ["year", "attend", "reliten", "fund"] + list(ITEMS)
d = pd.read_stata(gp, columns=cols, convert_categoricals=False)
M = pd.DataFrame({"year": d.year})
for c, kk in ITEMS.items():
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, kk=kk: (v >= 1) & (v <= kk))
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
R = M.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = zs(R[["attend", "reliten", "fund"]]).mean(axis=1); M = M.join(R["REL"])
t = M.groupby("year")["REL"].transform(
    lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)
HI, LO = (t == 2), (t == 0)

print("=== ⓪ 硬规则①:每道题真的问了哪些年、每层多少人(变量名不是测量)===")
CELLS = {}
for it, kk in ITEMS.items():
    ok = M[it].notna() & (HI | LO); ys = {}
    for y, g in M[ok].groupby("year"):
        a = g[HI.loc[g.index]][it].to_numpy(float); b = g[LO.loc[g.index]][it].to_numpy(float)
        if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
    S = sorted(ys); dec = {}
    for y in S: dec.setdefault((y//10)*10, []).append(y)
    dec = {k: v for k, v in dec.items() if len(v) >= 3}
    okd = DEC in dec and len(S) >= 8
    CELLS[it] = dict(ys=ys, S=S, dec=dec, ok=okd, kk=kk)
    print(f"  {it:9s} 1–{kk} · 合格年 **{len(S):>2}** ({S[0] if S else '-'}–{S[-1] if S else '-'}) · "
          f"n={int(M[ok][it].notna().sum()):>6,} · 九十年代{'可用' if okd else '**不可用**'}")
USE = [i for i in ITEMS if CELLS[i]["ok"]]
print(f"  ⇒ 进入分析的题:**{len(USE)}/8** ⇒ {USE}")

def dep(it, lay, Bv=B, seed=846, src=None):
    C = CELLS[it]; ys = src if src else C["ys"]; S = C["S"]; span = S[-1]-S[0]; yy = C["dec"][DEC]
    ref = (ys[S[-1]][lay].mean()-ys[S[0]][lay].mean())*(yy[-1]-yy[0])/span
    obs = float(ys[yy[-1]][lay].mean()-ys[yy[0]][lay].mean()) - ref
    rg = np.random.default_rng(seed); out = np.empty(Bv)
    r = lambda a: a[rg.integers(0, len(a), len(a))]
    for i in range(Bv):
        rf = (r(ys[S[-1]][lay]).mean()-r(ys[S[0]][lay]).mean())*(yy[-1]-yy[0])/span
        out[i] = r(ys[yy[-1]][lay]).mean() - r(ys[yy[0]][lay]).mean() - rf
    lo, hi = np.quantile(out, [.025, .975])
    p = max(2*min(float(np.mean(out <= 0)), float(np.mean(out >= 0))), 1.0/(Bv+1))
    return dict(obs=obs, lo=float(lo), hi=float(hi), p=p, sd=float(np.std(out)),
                excl=bool(lo > 0 or hi < 0))

print(f"\n=== ① 八题 × 两层的九十年代偏离(B={B})· ⚠ 天花板余量同列,>60% 标 `CEILING?` ===")
G16, keys = {}, []
for it in USE:
    C = CELLS[it]; yy = C["dec"][DEC]
    line = []
    for lay, nm in ((0, "虔诚"), (1, "世俗")):
        r = dep(it, lay)
        m0 = float(C["ys"][yy[0]][lay].mean()); m1 = float(C["ys"][yy[-1]][lay].mean())
        mv = m1-m0; reach = (C["kk"]-m0) if mv > 0 else (m0-1)
        hr = abs(mv)/reach if reach > 0 else float("inf")
        r.update(head=hr, m0=m0, m1=m1, ceiling=bool(hr > 0.60))
        G16[(it, nm)] = r; keys.append((it, nm))
        line.append(f"{nm} **{r['obs']:+.4f}**[{r['lo']:+.4f},{r['hi']:+.4f}]"
                    f"{'**排零**' if r['excl'] else '     '}{' `CEILING?`' if r['ceiling'] else ''}"
                    f"(余占 {hr:.0%})")
    print(f"  {it:9s} " + "  |  ".join(line))
ps = [G16[k]["p"] for k in keys]
bh = {keys[i] for i in Gate.bh(ps, Q)}; by = {keys[i] for i in Gate.by(ps, Q)}
sec = [it for it in USE if G16[(it, "世俗")]["excl"] and not G16[(it, "世俗")]["ceiling"]]
dev = [it for it in USE if G16[(it, "虔诚")]["excl"] and not G16[(it, "虔诚")]["ceiling"]]
print(f"  `G3` 整族 **{len(keys)} 格**(8 题 × 2 层;⚠ `#832`:族越窄存活越易)⇒ "
      f"BH **{len(bh)}** · BY **{len(by)}**")
print(f"  ⇒ **世俗层排零(且非天花板)的题:{len(sec)}/{len(USE)}** ⇒ {sec}")
print(f"  ⇒ **虔诚层排零(且非天花板)的题:{len(dev)}/{len(USE)}** ⇒ {dev}")

print("\n=== ② 控制 ===")
tgt = USE[0]
C = CELLS[tgt]; yy = C["dec"][DEC]
Yp = {y: (a.copy(), b.copy()) for y, (a, b) in C["ys"].items()}
Yp[yy[-1]] = (Yp[yy[-1]][0] + 0.30, Yp[yy[-1]][1])
pc_h = dep(tgt, 0, 400, src=Yp)["obs"] - G16[(tgt, "虔诚")]["obs"]
pc_l = dep(tgt, 1, 400, src=Yp)["obs"] - G16[(tgt, "世俗")]["obs"]
print(f"  正控:只往 `{tgt}` 的**虔诚层**九十年代尾年植入 +0.30 ⇒ 虔诚层动 **{pc_h:+.4f}**"
      f"(预期 +0.3000;该年非全程端点 ⇒ 参照不变)· **世俗层动 {pc_l:+.4f}(该为 0)**")
S = C["S"]; span = S[-1]-S[0]
full = {lay: float(C["ys"][S[-1]][lay].mean()-C["ys"][S[0]][lay].mean()) for lay in (0, 1)}
Yu = {y: tuple(C["ys"][S[0]][lay] + full[lay]*(y-S[0])/span for lay in (0, 1)) for y in S}
nc_h = dep(tgt, 0, 400, src=Yu)["obs"]; nc_l = dep(tgt, 1, 400, src=Yu)["obs"]
print(f"  负控:`{tgt}` 两层都严格匀速的世界 ⇒ 虔诚 **{nc_h:+.2e}** · 世俗 **{nc_l:+.2e}** —— "
      f"⚠ **「这个零该不该是零?」该**(匀速按定义无偏离,且无离散化 ⇒ **解析零**)")

G = Gate("#846 · 「是世俗那边走了」是关于人的,还是只关于这一道题")
G.asserted("① 硬规则①:逐题印出合格年数、样本量、九十年代是否可用 —— **变量名不是测量**",
           bool(len(USE) >= 4), f"可用 {len(USE)}/8 ⇒ {USE}", kind="control")
G.asserted("② 正控:只往**一题一层**植入 +0.30 ⇒ 该格按预期动,**另一层不动**",
           bool(abs(pc_h-0.30) < 0.02 and abs(pc_l) < 0.02),
           f"虔诚 {pc_h:+.4f}/预期 +0.3000 · 世俗 {pc_l:+.4f}", kind="control")
G.asserted("③ 负控:两层都严格匀速的世界里偏离必须为 **0**(⚠ **这个零该是零**;无离散化 ⇒ 解析零)",
           bool(abs(nc_h) < 1e-9 and abs(nc_l) < 1e-9), f"{nc_h:+.2e} / {nc_l:+.2e}", kind="control")
G.asserted("④ 前提(跑前写下的最强混淆):**八题的天花板余量差别很大** ⇒ 逐格印出余量占用,"
           "**>60% 的格标 `CEILING?` 且不计入判词** —— `#838` 只量过 `homosex`,别的题没量过",
           bool(all(np.isfinite(G16[k]["head"]) for k in keys)),
           f"标记为 CEILING? 的格 {sum(1 for k in keys if G16[k]['ceiling'])}/{len(keys)}",
           kind="control")
G.asserted("⑤ kill(预注册):「这是关于人的、不是关于题的」要成立,需**世俗层在 ≥5/8 题上排零"
           "且多于虔诚层**",
           bool(len(sec) >= 5 and len(sec) > len(dev)),
           f"世俗 {len(sec)}/{len(USE)} vs 虔诚 {len(dev)}/{len(USE)}", kind="kill",
           yardstick="每题每层自己的九十年代偏离,对照它自己的 95% 自助区间",
           yardstick_noise=float(np.mean([G16[k]["sd"] for k in keys])))
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif len(sec) >= 5 and len(sec) > len(dev):
    VERD = (f"**A 这是关于人的。** 世俗层在 **{len(sec)}/{len(USE)}** 题上偏离了自己的步子"
            f"({sec}),虔诚层只有 {len(dev)} 题。\n"
            f"  ⇒ **一句关于人的话:九十年代不是同性恋这一道题上出了事 ——\n"
            f"  是不信教的那一群人在那十年里整体加速了,而同性恋只是他们经过的其中一道题。**")
elif len(sec) <= 2:
    VERD = (f"**B 这是关于题的 ⇒ `#838` 的作用域要收窄到 `homosex`。** "
            f"世俗层只在 **{len(sec)}/{len(USE)}** 题上排零({sec});虔诚层 {len(dev)} 题({dev})。\n"
            f"  ⇒ **一句关于人的话:「不信教的人走得更快」不是一句关于那群人的话 ——\n"
            f"  在八道态度题里,他们只在同性恋这一道上明显偏离了自己的步子。\n"
            f"  九十年代发生的事,属于这个题目,不属于那群人。**")
else:
    VERD = (f"**C 「谁在动」这个问法在八题上不成立。** 世俗层排零 {len(sec)} 题({sec})、"
            f"虔诚层 {len(dev)} 题({dev}),**没有稳定的「哪一层」**。\n"
            f"  ⇒ **元分离器:我的分解方式(「是哪一边在动」)只在一道题上有意义,\n"
            f"  而不是一条可以跨题目搬运的性质。**")
print(VERD)
json.dump(dict(items=list(ITEMS), used=USE, decade=DEC, B=B, q=Q, family=len(keys),
               grid={f"{i}|{l}": v for (i, l), v in G16.items()},
               bh=sorted(f"{i}|{l}" for i, l in bh), by=sorted(f"{i}|{l}" for i, l in by),
               secular_excl=sec, devout_excl=dev,
               pos_control=dict(devout=pc_h, secular=pc_l), neg_control=dict(d=nc_h, s=nc_l),
               admissible=adm, verdict=VERD, gate_ok=G.verdict()),
          open(OUT/"whose_move_on_every_item.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'whose_move_on_every_item.json'}")

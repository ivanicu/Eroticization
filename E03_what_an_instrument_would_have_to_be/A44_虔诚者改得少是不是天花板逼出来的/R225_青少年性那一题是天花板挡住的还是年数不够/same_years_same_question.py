"""#786 · E03·A44·R225 —— `teensex` 读不出来,是天花板挡住的,还是年数不够?

`#785` 判了 C,而它自己在「改不了的」里写下了一条:
   **`teensex` 的跨度只有 38 年(其余 50 年)⇒ 它的不可读里有一部分只是年数少,本轮没有拆开。**
⚠ 这正是 `#777` 那条「墙从没被查过」的形状:**一个被归因给机制的结果,可能只是功效。**
   而 `#785`① 已经写下了那个便宜的先验检查,本轮就是它。

G1 估计量:**同一道题在两种跨度下的可读性** —— 把 `homosex` 与 `premarsx` 截到
   **`teensex` 的那一组年份**(不是同样的年数,是**同一组年份**),重算比值的自助区间,
   再问 `#782` 的那个问题:**区间含不含 1.0。**

识别:这是一个**配对**设计 —— 同一道题、同一批人、同一具机器,**只动年份集合**。
   ⇒ 可读性的变化只能来自年份集合,**这比任何「控制变量」都干净**,因为其余一切逐字相同。

两个世界:
   span    **年数**:`homosex` 截到 `teensex` 的年份后**也变得不可读**
           ⇒ `teensex` 的不可读与天花板无关,`#785` 关于那一题的那半句必须撤。
   ceiling **天花板/题目本身**:`homosex` 截短后**仍然可读**
           ⇒ 年数不足解释不了 `teensex`,`#785` 的读法站得住。

预测矩阵:
   | 世界 | 现在 | 若截短后不可读 | 若截短后仍可读 |
   | span    | 0.5 | **0.90** | 0.05 |
   | ceiling | 0.5 | 0.10 | **0.90** |
   ⚠ 而我**希望**输的是 ceiling —— 因为 `#785` 刚把它写上页面。这是 `frontier §3` 的 basin 逃逸。

预注册判词(条件式):
  if 正控开火(全跨度下 `homosex` 必须可读 —— 已知答案 8/8)and 截短样本非空:
      if 截短后 `homosex` 的多数格变得**不可读** -> span:撤回 `#785` 关于 `teensex` 的归因
      else                                       -> ceiling:`#785` 那半句站得住
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**截短同时改变了年数与年代区间**(`teensex` 起始更晚)——
   所以「截短后不可读」既可能是年数少,也可能是**那段年代里本来就没什么在动**。
   ⇒ 同一轮里放对照:**再跑一个「等年数但取最早 38 年」的臂**,把「年数」与「哪一段年代」分开。

本轮换不了仪器,理由由 `R223/instrument_search.py` 跑出来:第二具仪器的三条规格本机六具全部落选。
"""
import numpy as np, pandas as pd, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate

RNG = np.random.default_rng(225)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
ITEMS = ("homosex", "premarsx", "teensex")
VALID = {**{c: (1, 4) for c in ITEMS}, "attend": (0, 8), "reliten": (1, 4), "fund": (1, 3)}
for c, rg in VALID.items():
    dr, _ = check_kept_codes(gp, c, rg)
    if dr: print(f"  #766 前瞻:{c} 删 " + " · ".join(f"码{int(a)} {b!r}({sh*100:.2f}%)" for a, b, n, sh in dr[:2]))
d = pd.read_stata(gp, columns=["year"]+list(VALID), convert_categoricals=False)
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, lo=VALID[c][0], hi=VALID[c][1]: (v >= lo) & (v <= hi)) for c in VALID})
M["year"] = d.year
cat = pd.read_stata(gp, columns=list(ITEMS), convert_categoricals=True)
for c in aligned({k: list(cat[k].cat.categories)[:4] for k in ITEMS}, "strict"): M[c] = -M[c]+5
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
z = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = z(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["a"] = pd.cut(REL.attend, [-1, 1, 5, 8], labels=[0, 1, 2]).astype(float)
REL["b"] = REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def rows_of(item, col, k, nmin=120, keep=None):
    g = REL.dropna(subset=[item, col]); out = []
    for y, gy in g[g[col] == k].groupby("year"):
        if len(gy) < nmin: continue
        if keep is not None and int(y) not in keep: continue
        out.append((int(y), float(gy[item].mean()), float((gy[item] == 4).mean())))
    return out
def readable(item, col, st, keep=None, B=4000):
    j = 1 if st == "水平" else 2
    rA, rB = rows_of(item, col, 2, keep=keep), rows_of(item, col, 0, keep=keep)
    if len(rA) < 8 or len(rB) < 8: return None
    yA = np.array([r[0] for r in rA], float); yB = np.array([r[0] for r in rB], float)
    vA = np.array([r[j] for r in rA]); vB = np.array([r[j] for r in rB])
    f = lambda ia, ib: slope(yA[ia], vA[ia])/slope(yB[ib], vB[ib])
    bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(B)])
    bs = bs[np.isfinite(bs)]
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    return dict(r=float(f(np.arange(len(yA)), np.arange(len(yB)))), lo=lo, hi=hi,
                covers1=bool(lo <= 1.0 <= hi), nyr=[len(rA), len(rB)],
                span=[float(yA[-1]-yA[0]), float(yB[-1]-yB[0])])

TY = sorted({r[0] for c in ("a", "b") for k in (0, 2) for r in rows_of("teensex", c, k)})
print(f"\n=== ① `teensex` 的年份集合:{len(TY)} 年 · {TY[0]}–{TY[-1]} ===")
CUTS = {"(a) attend三档": "a", "(b) 三题三分位": "b"}
ARMS = {"全跨度(原样)": None, f"截到 teensex 的年份({TY[0]}–{TY[-1]})": set(TY)}
# ⚠ 第三臂把「年数」与「哪一段年代」分开(跑之前写下的混淆的控制)
EARLY = None

print("\n=== ② 同一道题 · 只动年份集合 ===")
res = []
for item in ("homosex", "premarsx"):
    for arm, keep in ARMS.items():
        for cn, col in CUTS.items():
            for st in ("水平", "端点"):
                r = readable(item, col, st, keep)
                if r is None:
                    print(f"  {item:9s} {arm:32s} {cn:14s} {st:4s}  年数不足,跳过"); continue
                if EARLY is None and arm.startswith("全"):
                    ys = sorted({x[0] for x in rows_of(item, col, 2)})
                    EARLY = set(ys[:len(TY)])
                res.append(dict(item=item, arm=arm, cut=cn, stat=st, **r))
                print(f"  {item:9s} {arm:32s} {cn:14s} {st:4s}  比 {r['r']:6.3f} "
                      f"[{r['lo']:7.3f}, {r['hi']:7.3f}] · 年 {r['nyr']}"
                      f"  {'**不可读(含1.0)**' if r['covers1'] else '可读'}")

print(f"\n=== ③ 第三臂:等年数但取**最早** {len(TY)} 年 —— 把「年数」与「哪一段年代」分开 ===")
early = []
for item in ("homosex", "premarsx"):
    for cn, col in CUTS.items():
        ys = sorted({x[0] for x in rows_of(item, col, 2)})
        keep = set(ys[:len(TY)])
        for st in ("水平", "端点"):
            r = readable(item, col, st, keep)
            if r is None: continue
            early.append(dict(item=item, cut=cn, stat=st, arm=f"最早{len(TY)}年", **r))
            print(f"  {item:9s} {cn:14s} {st:4s}  比 {r['r']:6.3f} [{r['lo']:7.3f}, {r['hi']:7.3f}]"
                  f"  {'**不可读(含1.0)**' if r['covers1'] else '可读'}")

full = [r for r in res if r["arm"].startswith("全")]
trunc = [r for r in res if not r["arm"].startswith("全")]
f_ok = sum(1 for r in full if not r["covers1"]); t_ok = sum(1 for r in trunc if not r["covers1"])
e_ok = sum(1 for r in early if not r["covers1"])
print(f"\n  全跨度可读 **{f_ok}/{len(full)}** · 截到 teensex 年份可读 **{t_ok}/{len(trunc)}** · "
      f"最早{len(TY)}年可读 **{e_ok}/{len(early)}**")

G = Gate("#786 · teensex 读不出来,是天花板还是年数")
hs_full = [r for r in full if r["item"] == "homosex"]
G.asserted("① 正控:全跨度下 `homosex` 必须可读(已知答案 —— `#785` 是 8/8)",
           bool(hs_full and all(not r["covers1"] for r in hs_full)),
           f"homosex 全跨度 {sum(1 for r in hs_full if not r['covers1'])}/{len(hs_full)} 可读", kind="control")
G.asserted("② 前提:截短臂必须非空(否则 FAIL 读起来像「不可读」,而真相是没有样本 —— `#777` 那条)",
           bool(trunc), f"截短臂 {len(trunc)} 格", kind="control")
hs_tr = [r for r in trunc if r["item"] == "homosex"]
span_world = bool(hs_tr and sum(1 for r in hs_tr if r["covers1"]) > len(hs_tr)/2)
G.asserted("③ kill(预注册):`#785` 对 `teensex` 的天花板归因要站住,"
           "需 `homosex` 截到同一组年份后**多数格仍可读**",
           not span_world,
           f"homosex 截短后不可读 {sum(1 for r in hs_tr if r['covers1'])}/{len(hs_tr)}", kind="kill")
print(); print(G)

print("\n"+"="*92)
ok = bool(hs_full and all(not r["covers1"] for r in hs_full) and trunc)
if not ok:
    v = "**UNVERIFIED:正控没过或截短臂为空,本轮不下判。**"
elif span_world:
    v = (f"**span:年数就够解释。** `homosex` 截到 `teensex` 的 {len(TY)} 个年份后,"
         f"**{sum(1 for r in hs_tr if r['covers1'])}/{len(hs_tr)} 格变得不可读** —— "
         f"而它在全跨度下是 {len(hs_full)}/{len(hs_full)} 可读。"
         f"⇒ **`#785` 里「`teensex` 上天花板确实在挡」那半句撤回**:同样的年份集合会让"
         f"**任何**一题变得读不出来,那不是关于 `teensex` 的事实,是关于 38 个年份点的事实。")
else:
    v = (f"**ceiling:年数解释不了。** `homosex` 截到 `teensex` 的 {len(TY)} 个年份后仍然"
         f"**{sum(1 for r in hs_tr if not r['covers1'])}/{len(hs_tr)} 格可读** ⇒ "
         f"`teensex` 的不可读不是年数造成的,`#785` 那半句站得住。")
print(v)
json.dump(dict(teensex_years=TY, full=full, truncated=trunc, early=early,
               f_ok=f_ok, t_ok=t_ok, e_ok=e_ok, span_world=span_world,
               verdict=v, gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"same_years_same_question.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'same_years_same_question.json'}")

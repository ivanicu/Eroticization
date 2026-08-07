"""#790 · E03·A45·R229 —— 那十二道「动了但不稳」的题,是分母太小,还是年份太少?

`#789` 找到的那条轴只有 **8 道题**撑着,而它自己写下:**`n=8` 是本轮一切结论的天花板。**
另有 **12 道题「社会确实动了,但比值不稳」** —— **它们是这条轴唯一可能的增量。**
`#789`② 预注册:**先问不稳来自哪里,再谈能不能救。**

G1 估计量:**比值区间宽度的来源分解** —— 两个事先可测的候选,谁预测得了它:
   ① **分母的信噪比** `|slope_非虔诚| / se(slope_非虔诚)`(自助 se)—— 分母越接近零,比值越发散
   ② **年份点数** `n_yr`
   ⇒ 把 20 道「动了」的题(8 可读 + 12 不稳)的 `log(区间宽度)` 同时对这两个量回归,看谁在说话。

⚠⚠ **两个世界,而它们的后果完全不同 —— 这是本轮值得跑的原因:**
   A **分母驱动**:不稳由 `|slope|/se` 预测 ⇒ **年份窗口手术救不了它们**;
     那不是功效问题,是**这个估计量在社会没怎么动的题上本来就没有定义** ⇒
     **`n=8` 是一个硬天花板,而 `#789` 那条轴要永远带着「只在社会动得多的题上成立」这个 scope。**
   B **年份驱动**:不稳由 `n_yr` 预测 ⇒ **`#786` 那套等窗口/合并手术可以搬过来**,n 有希望长大。

预测矩阵:
   | 世界 | 现在 | 若 ①的系数显著而②不 | 若 ②显著而①不 | 若都不显著 |
   | A 分母 | 0.55 | **0.90** | 0.05 | 0.20 |
   | B 年份 | 0.30 | 0.05 | **0.90** | 0.20 |
   | 都不是 | 0.15 | 0.05 | 0.05 | **0.60**(⇒ 我连不稳的来源都没找对,得换问法)

预注册判词(条件式):
  if 正控开火(合成检验:人为把分母缩小,区间宽度必须涨)and 可用题 >= 12:
      比较两个标准化回归系数 |b1| 与 |b2|,并各自带自助区间:
      if |b1| 的区间排除 0 且 |b2| 的不排除 -> A 分母驱动,n=8 是硬天花板
      elif 反过来                            -> B 年份驱动,可以救
      else                                   -> 都不是/分不开,如实报,不选边
  else: UNVERIFIED
⚠ **不用 p 值,用「系数自己的自助区间排不排除 0」** —— `n≈20`,而 `#782` 那一课就是
  「两个量各自显著 ≠ 它们的关系稳」;这里直接对**要下判的那个量**取区间。

⚠ 跑之前写下的最强混淆:**两个候选量本身可能高度相关**(年份点少的题往往也是问得少、
  斜率估得糙的题)⇒ 那样任何一个都能「预测」,而回归会把功劳随机分给一个。
  ⇒ **同一轮里先报 `corr(①, ②)`**;若 |corr| > 0.7,**声明两者不可分,不许下 A/B 的判。**

⚠ `#787` 的规矩:每一个比例同时报「可读 / 尝试」。
本轮换不了仪器,理由同 `R223/instrument_search.py`(对象是世界,第二具仪器本机六具全部落选)。
"""
import numpy as np, pandas as pd, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

RNG = np.random.default_rng(229)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
PREV = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R228_那条边界不约束我真正需要的那个量/results/ratio_needs_no_polarity.json"))
ITEMS = [r["item"] for r in PREV["rows"] if r.get("moved")]     # 社会确实动了的那些题
print(f"=== ⓪ 对象:`#789` 里**社会确实动了**的 {len(ITEMS)} 道题(可读 "
      f"{sum(1 for r in PREV['rows'] if r.get('readable'))} + 不稳 "
      f"{sum(1 for r in PREV['rows'] if r.get('moved') and not r.get('readable'))})===")

d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+ITEMS, convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
KMAX = {c: (4 if c == "homosex" else len(cat[c].cat.categories)) for c in ITEMS}
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= KMAX[c])) for c in ITEMS})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
z = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = z(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["b"] = REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def series(item, k, nmin=120):
    g = REL.dropna(subset=[item]); out = []
    for y, gy in g[g["b"] == k].groupby("year"):
        if len(gy) < nmin: continue
        out.append((int(y), float(gy[item].mean())))
    return out
def measure(item, shrink=1.0):
    """shrink<1 人为把**分母层**的斜率压小(正控用),其余一切不变。"""
    rA, rB = series(item, 2), series(item, 0)
    if len(rA) < 10 or len(rB) < 10: return None
    yA = np.array([r[0] for r in rA], float); yB = np.array([r[0] for r in rB], float)
    vA = np.array([r[1] for r in rA]); vB = np.array([r[1] for r in rB])
    if shrink != 1.0:
        vB = vB.mean() + shrink*(vB - vB.mean())      # 压斜率,保留噪声结构
    sB = slope(yB, vB)
    bsB = np.array([slope(yB[i], vB[i]) for i in (RNG.integers(0, len(yB), len(yB)) for _ in range(1500))])
    seB = float(np.std(bsB, ddof=1))
    f = lambda ia, ib: slope(yA[ia], vA[ia])/slope(yB[ib], vB[ib])
    bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(3000)])
    bs = bs[np.isfinite(bs)]
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    return dict(item=item, width=hi-lo, snr=abs(sB)/seB if seB > 0 else np.nan,
                nyr=min(len(yA), len(yB)), readable=bool(not (lo <= 1.0 <= hi)))

print("\n=== ① 逐题:区间宽度 · 分母信噪比 · 年份点数 ===")
rows = [m for m in (measure(c) for c in ITEMS) if m]
for m in sorted(rows, key=lambda x: x["width"]):
    print(f"  {m['item']:9s} 宽 {m['width']:8.3f} · 分母 |slope|/se {m['snr']:6.2f} · 年 {m['nyr']:>3} "
          f"{'可读' if m['readable'] else '不稳'}")

W = np.log(np.array([m["width"] for m in rows]))
S = np.array([m["snr"] for m in rows]); Y = np.array([float(m["nyr"]) for m in rows])
rho = float(np.corrcoef(S, Y)[0, 1])
print(f"\n=== ② 跑前写下的混淆先查:两个候选量自己相关吗 ===")
print(f"  corr(分母信噪比, 年份点数) = **{rho:+.3f}** —— "
      f"{'⛔ |corr|>0.7,两者不可分,不许下 A/B 的判' if abs(rho) > 0.7 else '≤0.7,可以分开读'}")

def std_betas(y, X):
    Xz = np.column_stack([(c-c.mean())/c.std(ddof=1) for c in X.T])
    yz = (y-y.mean())/y.std(ddof=1)
    A = np.column_stack([np.ones(len(yz)), Xz])
    return np.linalg.lstsq(A, yz, rcond=None)[0][1:]

X = np.column_stack([S, Y])
b = std_betas(W, X)
BB = np.array([std_betas(W[i], X[i]) for i in (RNG.integers(0, len(W), len(W)) for _ in range(4000))])
ci = [(float(np.percentile(BB[:, j], 2.5)), float(np.percentile(BB[:, j], 97.5))) for j in (0, 1)]
print(f"\n=== ③ `log(区间宽度)` 同时对两者回归(标准化系数,各带自助区间)===")
for j, nm in enumerate(("分母信噪比", "年份点数")):
    print(f"  {nm:8s} b = {b[j]:+.3f} · 95% [{ci[j][0]:+.3f}, {ci[j][1]:+.3f}]"
          f"  {'**排除 0**' if (ci[j][0] > 0) == (ci[j][1] > 0) else '含 0'}")

print(f"\n=== ④ 正控:人为把分母压小,区间宽度必须涨(否则我量的根本不是这件事)===")
probe = [c for c in ("homosex", "spanking", "suicide2") if c in ITEMS][:2]
pc_rows = []
for c in probe:
    base = measure(c, 1.0); sh = measure(c, 0.25)
    if not base or not sh: continue
    pc_rows.append((c, base["width"], sh["width"]))
    print(f"  {c:9s} 原宽 {base['width']:.3f} → 分母压到 0.25× 后 {sh['width']:.3f} "
          f"{'✅ 涨了' if sh['width'] > base['width'] else '⛔ 没涨'}")
pc_ok = bool(pc_rows and all(s > b0 for _, b0, s in pc_rows))

G = Gate("#790 · 十二道不稳的题:分母太小,还是年份太少")
G.asserted("① 正控:人为把分母压到 0.25×,区间宽度必须涨(否则我量的不是这件事)",
           pc_ok, f"{len(pc_rows)} 道探针:" + " · ".join(f"{c} {b0:.2f}→{s:.2f}" for c, b0, s in pc_rows), kind="control")
G.asserted("② 前提(跑前写下的混淆):两个候选量必须可分,|corr| ≤ 0.7",
           bool(abs(rho) <= 0.7), f"corr = {rho:+.3f}", kind="control")
G.asserted("③ 前提:题数 >= 12", bool(len(rows) >= 12), f"题数 {len(rows)} / 尝试 {len(ITEMS)}", kind="control")
excl = [((c[0] > 0) == (c[1] > 0)) for c in ci]
A_world = bool(excl[0] and not excl[1]); B_world = bool(excl[1] and not excl[0])
G.asserted("④ kill(预注册):要判「年份驱动、可以救」,需**年份**系数的区间排除 0 而分母的不排除",
           B_world, f"分母 {'排除0' if excl[0] else '含0'} · 年份 {'排除0' if excl[1] else '含0'}", kind="kill")
print(); print(G)

print("\n"+"="*92)
if not (pc_ok and abs(rho) <= 0.7 and len(rows) >= 12):
    v = "**UNVERIFIED:正控或前提没过,本轮不下判。**"
elif A_world:
    v = (f"**A 分母驱动 —— 而这意味着 `n=8` 是一个硬天花板。** `log(区间宽度)` 上,"
         f"**分母信噪比的标准化系数 {b[0]:+.3f} [{ci[0][0]:+.3f}, {ci[0][1]:+.3f}] 排除 0**,"
         f"而年份点数 {b[1]:+.3f} [{ci[1][0]:+.3f}, {ci[1][1]:+.3f}] 含 0。\n"
         f"  ⇒ **那十二道题不稳,不是因为年份少,是因为社会在它们上面动得不够多。**\n"
         f"  ⇒ **年份窗口手术救不了它们** —— 这不是功效问题,是**这个估计量在社会没怎么动的题上本来就没定义**。\n"
         f"  ⇒ **`#789` 那条轴要永远带着它的 scope:只在社会动得多的那些题上成立,而那是 8 道。**")
elif B_world:
    v = (f"**B 年份驱动 —— 可以救。** 年份点数的系数 {b[1]:+.3f} [{ci[1][0]:+.3f}, {ci[1][1]:+.3f}] 排除 0,"
         f"而分母信噪比 {b[0]:+.3f} 含 0 ⇒ **`#786` 那套等窗口/合并手术可以搬过来,n 有希望长大。**")
else:
    v = (f"**分不开,如实报。** 分母信噪比 {b[0]:+.3f} [{ci[0][0]:+.3f}, {ci[0][1]:+.3f}] · "
         f"年份点数 {b[1]:+.3f} [{ci[1][0]:+.3f}, {ci[1][1]:+.3f}] —— "
         f"{'两者的区间都排除 0' if all(excl) else '两者的区间都含 0'};"
         f"⇒ **在 {len(rows)} 道题上,这两个来源分不开,而分不开本身要写下来,不许挑一个说。**")
print(v)
json.dump(dict(rows=rows, corr_snr_nyr=rho, betas=list(map(float, b)), ci=ci,
               positive_control=pc_rows, A_world=A_world, B_world=B_world,
               verdict=v, gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"why_unstable.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'why_unstable.json'}")

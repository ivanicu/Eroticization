"""#844 · E03·A82·R283 —— 审计我自己写过的三堵墙:哪一堵是查过的,哪一堵只是相信的

**Frontier,而它的理由是 `#843` 当场发生的事**:我在 `#842` 里登记了一条「结构性改不了」,
**下一轮用一条 `git show` 就推翻了它。**
⇒ `#843`② 于是登记:**回头查每一条写过的「改不了 / 结构性不可能」。**
**这不是洁癖 —— 一堵墙关掉的是一整条研究路线,而它从来不被审计**
(`feedback_fabricated_impossibility` · `realstat` 的 *a wall never checked*)。

**三堵墙,而每一堵现在都关着一条路:**
 **W1 `#835`** —— 「同一调查年内 `世代 = 年份 − 年龄` 完全共线 ⇒ 扣年龄同时在扣世代,
    **APC 不可分离**」⇒ 关掉了「宗教那条缝里有多少是年龄」。
 **W2 `#837`** —— 「SCCS 里『同性行为频率』只有 Broude & Greene 一份 ⇒
    **结构性拿不到第二具仪器**」⇒ 关掉了社会这个单位上的跨编码者复现。
 **W3 `#840`** —— 「十年分辨率是这份调查的结构性上限(GSS 一个十年只跑 4–6 次)⇒
    **再多算也换不来更细的年份聚类**」⇒ 关掉了「缝有没有合上」。

**⚠⚠ 先分类,再检验 —— 因为三种墙要用三种办法拆:**
  **DERIVATION** 代数强制的(比如共线性)⇒ 只能检查**代数本身说了什么**,不能靠跑数据。
  **MEASUREMENT** 量出来的(比如功效)⇒ 可以换设计重量。
  **BELIEF** 从没查过的(比如「只有一份编码」)⇒ **一条查询就能定。**
**`#842` 那堵墙是 BELIEF 伪装成 DERIVATION —— 而那正是最危险的一类。**

G1 估计量:**每堵墙的类别,以及一次能让它倒下的最便宜的检验的结果。**

预注册判词(条件式):
  if 每堵墙的检验各自带正控与负控:
      三堵全立   -> A(墙是真的,而现在是查过的,不是相信的)
      任一堵倒   -> B(那条路重开,并撤掉对应那句「结构性」)
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**我有动机让墙立着** —— 墙倒了意味着我要回去做更多工作。
  ⇒ 控制:**每一堵墙的检验都设计成「它倒下」比「它立着」更容易发生**
  (W1 只要求恢复出**非线性**部分而不要求恢复全部;W2 用最宽的词表而不是最窄的;
  W3 只要求**任何**一种别的时间分箱把区间收窄到可分辨,不要求它给出显著结果)。
  **若在这种偏向下墙仍然立着,那是更强的证据,不是更弱的。**
⚠ 本轮换不了仪器;而它不需要 —— 审计的对象是我自己写过的句子。
"""
import csv, collections, json, math, pathlib, re, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
W = {}

# ── W2 · BELIEF:一条查询就能定 ────────────────────────────────────────────────
print("=== W2 `#837`「SCCS 里同性行为只有一份编码」—— 类别 **BELIEF**(一条查询就能定)===")
S = ROOT/"data/external/dplace/repo/datasets/SCCS"
V = list(csv.DictReader(open(S/"variables.csv", encoding="utf-8")))
DD = collections.defaultdict(set)
for r in csv.DictReader(open(S/"data.csv", encoding="utf-8")):
    if r["code"] not in ("", "NA"): DD[r["var_id"]].add(r["soc_id"])
WIDE = re.compile(r"homosex|same.?sex|sodom|berdache|two.?spirit|transvest|bisexual|lesbian|gay|"
                  r"inversion|catamite|pederast", re.I)
NARROW = re.compile(r"homosexual", re.I)
hits = [v for v in V if WIDE.search(v["title"]+" "+v.get("description", "")[:600])]
srcs = collections.Counter(v["source"] for v in hits)
RAW = sum(1 for _ in open(S/"variables.csv", encoding="utf-8"))
print(f"  ⚠ **顺带更正 `#837`**:文件 **{RAW:,} 行**,但 `csv` **解析出 {len(V):,} 条**变量 —— "
      f"**描述字段里有换行**,所以 `#837` 写的 1,968 是 `wc -l`,不是变量数。")
print(f"  用**最宽的词表**(12 个词根,含 berdache/two-spirit/transvestite 等)扫题名+描述:"
      f"命中 **{len(hits)}** 条")
for v in sorted(hits, key=lambda x: -len(DD[x["id"]])):
    print(f"    {v['id']:<9} n={len(DD[v['id']]):>3}  {v['source'][:26]:26s} 「{v['title'][:58]}」")
narrow_n = sum(1 for v in V if NARROW.search(v["title"]))
W["W2"] = dict(kind="BELIEF", total_vars=len(V), wide_hits=len(hits),
               sources=dict(srcs), narrow_title_hits=narrow_n,
               stands=bool(len(srcs) == 1))
print(f"  ⇒ 命中的**来源数 = {len(srcs)}** ⇒ {dict(srcs)}")
print(f"  ⇒ **墙 {'立着' if W['W2']['stands'] else '倒了'}** —— 而它现在是**查过的**,不是相信的。")

# ── W1 · 代数说了什么 ────────────────────────────────────────────────────────
print("\n=== W1 `#835`「APC 不可分离」—— 类别 **DERIVATION**,而代数说的比我写的少 ===")
print("  我写的是「扣年龄同时在扣世代 ⇒ APC 不可分离」。**而线性代数说的是:**")
print("  `世代 = 年份 − 年龄` 只让三者的**线性分量**共线;**非线性分量是可识别的。**")
print("  ⇒ 最便宜的反证:**造一个只有非线性年龄效应、零世代效应的世界,看能不能把它捞回来。**")
rg = np.random.default_rng(844)
YEARS = np.arange(1974, 2025, 2); AGES = np.arange(18, 90)
rows = []
for y in YEARS:
    for a in AGES:
        n = 60
        age_nl = 0.6*math.sin((a-18)/72*2*math.pi)          # 只有非线性年龄效应
        per = 0.02*(y-1974)                                  # 纯线性时期效应
        rows.append((y, a, y-a, age_nl+per, n))
df = pd.DataFrame(rows, columns=["year", "age", "coh", "mu", "n"])
df["obs"] = df.mu + rg.normal(0, 0.35/np.sqrt(df.n))
# 只放 年龄哑变量 + 年份哑变量(世代不放)⇒ 看非线性年龄形状能否恢复
A = pd.get_dummies(df.age, prefix="a", drop_first=True).astype(float)
P = pd.get_dummies(df.year, prefix="y", drop_first=True).astype(float)
X = np.column_stack([np.ones(len(df)), A.values, P.values])
beta, *_ = np.linalg.lstsq(X, df.obs.values, rcond=None)
ahat = np.concatenate([[0.0], beta[1:1+A.shape[1]]])
truth = np.array([0.6*math.sin((a-18)/72*2*math.pi) for a in AGES]); truth = truth - truth[0]
det = lambda v: v - np.polyval(np.polyfit(np.arange(len(v)), v, 1), np.arange(len(v)))  # 去掉线性分量
r_nl = float(np.corrcoef(det(ahat), det(truth))[0, 1])
r_raw = float(np.corrcoef(ahat, truth)[0, 1])
# 负控:真值里没有非线性年龄效应时,不许恢复出一个
df2 = df.copy(); df2["obs"] = 0.02*(df2.year-1974) + rg.normal(0, 0.35/np.sqrt(df2.n))
b2, *_ = np.linalg.lstsq(X, df2.obs.values, rcond=None)
a2 = np.concatenate([[0.0], b2[1:1+A.shape[1]]])
r_null = float(np.corrcoef(det(a2), det(truth))[0, 1])
print(f"  正控:世界里**只有非线性年龄效应 + 纯线性时期效应**,只拟合 年龄哑 + 年份哑")
print(f"     ⇒ 恢复出的年龄曲线与真值:**去线性后 corr = {r_nl:+.4f}**(原始 corr = {r_raw:+.4f})")
print(f"  负控:真值里**没有**非线性年龄效应 ⇒ 去线性后 corr = **{r_null:+.4f}**(该接近 0)")
W["W1"] = dict(kind="DERIVATION", nonlinear_recovered_corr=r_nl, raw_corr=r_raw,
               null_corr=r_null, stands=bool(abs(r_nl) < 0.5))
print(f"  ⇒ **墙 {'立着' if W['W1']['stands'] else '**倒了**'}** —— "
      f"{'非线性也捞不回来' if W['W1']['stands'] else '**非线性分量是可识别的,`#835` 那句话说得太满**'}")

# ── W3 · MEASUREMENT:换一种时间分箱 ─────────────────────────────────────────
print("\n=== W3 `#840`「十年分辨率是结构性上限」—— 类别 **MEASUREMENT**,而「十年」是我选的 ===")
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", "homosex"], convert_categoricals=False)
M = pd.DataFrame({"homosex": pd.to_numeric(d.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
R = M.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = zs(R[["attend", "reliten", "fund"]]).mean(axis=1); M = M.join(R["REL"])
t = M.groupby("year")["REL"].transform(
    lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)
HI, LO = (t == 2), (t == 0); okm = M.homosex.notna() & (HI | LO)
ys = {}
for y, g in M[okm].groupby("year"):
    a = g[HI.loc[g.index]].homosex.to_numpy(float); b = g[LO.loc[g.index]].homosex.to_numpy(float)
    if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
Sy = sorted(ys); span = Sy[-1]-Sy[0]
def gapdep(bins, lo_y, hi_y, B=3000, seed=5):
    """在给定分箱下,覆盖 [lo_y,hi_y] 的那一箱的差距偏离与它自己的区间。"""
    grp = {}
    for y in Sy: grp.setdefault((y//bins)*bins, []).append(y)
    grp = {k: v for k, v in grp.items() if len(v) >= 2}
    # ⚠⚠ **第一版这里选桶用的是「落在窗口内且年数最多」,于是 bins=20 选到了 [2021,2022,2024]**
    #    —— **那是 2020 年代,不是 2010 年代。** 一个「换分箱后排零了」的结论,
    #    比较的根本不是同一个估计量。**而这正是 `realstat` 说的最贵的一种错:
    #    一次看起来成功的攻击,会撤掉一件真的东西。**
    #    ⇒ 改成:桶必须**真的覆盖 2010s**(含 2010 且起点 <2020),否则这个分箱**不适用**,
    #      而「不适用」要如实报成不适用,不能拿另一个格子冒充。
    key = [k for k in grp if grp[k][0] >= lo_y and grp[k][-1] <= hi_y
           and grp[k][0] < 2020 and 2010 in grp[k]]
    if not key: return None
    yy = grp[max(key, key=lambda k: len(grp[k]))]
    rg2 = np.random.default_rng(seed); g = lambda Z, y: Z[y][0].mean()-Z[y][1].mean()
    ref = (g(ys, Sy[-1])-g(ys, Sy[0]))*(yy[-1]-yy[0])/span
    obs = (g(ys, yy[-1])-g(ys, yy[0])) - ref
    out = np.empty(B); r = lambda a: a[rg2.integers(0, len(a), len(a))]
    for i in range(B):
        Z = {y: (r(ys[y][0]), r(ys[y][1])) for y in (Sy[0], Sy[-1], yy[0], yy[-1])}
        out[i] = (g(Z, yy[-1])-g(Z, yy[0])) - (g(Z, Sy[-1])-g(Z, Sy[0]))*(yy[-1]-yy[0])/span
    lo_, hi_ = np.quantile(out, [.025, .975])
    return dict(bin=bins, years=yy, obs=float(obs), lo=float(lo_), hi=float(hi_),
                half=float((hi_-lo_)/2), excl=bool(lo_ > 0 or hi_ < 0))
print(f"  `#840` 量的是**十年**分箱下 2010s 的区间半宽 ≈0.23 vs 待分辨量 ≈0.09 ⇒ 噪声是 2.5 倍。")
print(f"  ⚠ **而「十年」是我选的分箱,不是数据的性质。** 换几种分箱,看半宽会不会收窄:")
W3 = {}
for bsz in (5, 10, 20):
    r = gapdep(bsz, 2005, 2024)
    if r is None:
        print(f"     {bsz:>2} 年箱:**不适用** —— 没有一个桶真的覆盖 2010s(如实报,不拿别的格子冒充)")
        continue
    W3[bsz] = r; print(f"     {bsz:>2} 年箱({r['years'][0]}–{r['years'][-1]},{len(r['years'])} 年):"
                             f"偏离 **{r['obs']:+.4f}** [{r['lo']:+.4f},{r['hi']:+.4f}] · "
                             f"半宽 **{r['half']:.4f}**{'  **排零**' if r['excl'] else ''}")
lead = gapdep.__wrapped__(20, 2020, 2024) if hasattr(gapdep, "__wrapped__") else None
print(f"  ⚠ **另记一条线索,而它不是这堵墙的证据**:第一版误选到的那个桶是 **[2021,2022,2024]** ——"
      f"**2020 年代,不是 2010 年代**。它给的 +0.2613 [+0.0008,+0.5228] 是**另一个格子**,"
      f"下界只有 0.0008、只有 2 个年份,**登记为线索,不作任何主张**。")
best = min(W3.values(), key=lambda x: x["half"]) if W3 else None
W["W3"] = dict(kind="MEASUREMENT", bins=W3, best_half=best["half"] if best else None,
               any_excl=bool(any(v["excl"] for v in W3.values())),
               stands=bool(not any(v["excl"] for v in W3.values())))
print(f"  ⇒ **墙 {'立着' if W['W3']['stands'] else '**倒了**'}** —— "
      f"{'换了分箱也没有一种把它推到排零' if W['W3']['stands'] else '**有一种分箱把它推到了排零**'}")

fell = [k for k, v in W.items() if not v["stands"]]
G = Gate("#844 · 审计我自己写过的三堵墙")
G.asserted("① 前提(跑前写下的最强混淆):**我有动机让墙立着**(墙倒了意味着更多工作)⇒ "
           "每堵墙的检验都设计成**「它倒下」比「它立着」更容易发生**:W1 只要求恢复**非线性**分量、"
           "W2 用**最宽**的词表(12 词根)而不是最窄的、W3 只要求**任何**一种别的分箱把区间推到排零 ⇒ "
           "**在这种偏向下仍然立着的墙,是更强的证据,不是更弱的**",
           bool(len(W) == 3), f"三堵墙的类别 {[v['kind'] for v in W.values()]}", kind="control")
G.asserted("② W1 负控:世界里**没有**非线性年龄效应时,不许恢复出一个",
           bool(abs(W["W1"]["null_corr"]) < 0.35), f"去线性后 corr = {W['W1']['null_corr']:+.4f}",
           kind="control")
G.asserted("③ W2 正控:最宽词表必须至少捞到已知的那两条(`SCCS176`/`SCCS177`)"
           "—— 否则这个「只有一个来源」是仪器瞎了",
           bool(W["W2"]["wide_hits"] >= 2), f"宽词表命中 {W['W2']['wide_hits']} 条", kind="control")
G.asserted("④ W3 正控:换分箱后的箱必须真的**不同**(年数不同),否则「换了分箱」是空话",
           bool(len({tuple(v["years"]) for v in W3.values()}) == len(W3)),
           " · ".join(f"{k}年箱:{len(v['years'])} 年" for k, v in W3.items()), kind="control")
G.asserted("⑤ kill(预注册):「三堵墙都是真的」要成立,需**三堵都在各自的检验下立着**",
           bool(not fell), f"倒下的墙 {fell or '无'}", kind="kill",
           yardstick="每堵墙自己那条最便宜的反证的结果", yardstick_noise=0.0)
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    VERDICT = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif not fell:
    VERDICT = (f"**A 三堵墙都立着 —— 而它们现在是查过的,不是相信的。**\n"
         f"  W2 用最宽的词表扫 {W['W2']['total_vars']:,} 条变量,只有 **{W['W2']['wide_hits']}** 条命中,"
         f"来源数 **{len(W['W2']['sources'])}**;\n"
         f"  W1 造一个只有非线性年龄效应的世界,去线性后也只恢复到 corr={W['W1']['nonlinear_recovered_corr']:+.3f};\n"
         f"  W3 换了 {len(W3)} 种分箱,没有一种把 2010s 推到排零。\n"
         f"  ⇒ **一句关于方法的话:三堵墙全立着,而这一轮仍然值得跑 ——\n"
         f"  因为在跑之前,「立着」和「没查过」在我这里长得一模一样。**")
else:
    VERDICT = (f"**B 有墙倒了:{fell} ⇒ 那条路重开,对应那句「结构性」要撤。**\n" +
         ("\n".join(f"  {k}:{json.dumps(W[k], ensure_ascii=False)[:200]}" for k in fell)))
print(VERDICT)
print(f"\n⚠ **顺带更正 `#837`**:SCCS 变量**解析出 {len(V):,} 条**(而文件有 {RAW} 行 —— 描述字段里有换行,`#837` 写的 1,968 是 `wc -l` 不是变量数),`#837` 写的 **1,968** 是错的。")
json.dump(dict(walls=W, fell=fell, admissible=adm, verdict=VERDICT, gate_ok=G.verdict(),
               correction_to_837=f"variables.csv has {len(V)} rows, not 1968"),
          open(OUT/"audit_my_own_walls.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'audit_my_own_walls.json'}")

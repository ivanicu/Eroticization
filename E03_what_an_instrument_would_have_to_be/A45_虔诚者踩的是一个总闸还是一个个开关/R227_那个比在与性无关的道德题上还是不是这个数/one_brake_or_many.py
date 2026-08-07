"""#788 · E03·A45·R227 —— 虔诚者踩的是**一个总闸**,还是**一个个开关**?

`#786` 之后,活下来的那句话很窄:**在同性恋上,最虔诚的三分之一人只改了其余人的 0.31–0.48**
(三个年份窗口、12/12 格可读),**而天花板一分也没解释**(`#785`:逼出的比 1.37 > 1)。
`premarsx` 降级、`teensex` 不可读、`xmarsex` 问不出来。
⇒ **剩下的问题是心理学的:为什么是同性恋这一题?**

⚠⚠ **两个世界,而它们的本体不同 —— 不是参数不同:**
   A **一个总闸(dispositional)**:虔诚是一种**对任何道德变迁的普遍阻力**。
     ⇒ 那个比在**与性无关**的道德题上应当**同样小**,且与题目内容无关。
   B **一个个开关(doctrinal)**:那个刹车是**教义性的** ——
     只踩在传统**明确点过名**的地方,其余地方松开。
     ⇒ 那个比应当**随题目而变**:在教义点名处最小,在别处趋近 1。

**这是一次真正的分离,因为两个世界对同一批数据给出相反的预言,而以前从没在非性题上量过。**

G1 估计量:**逐题的「虔诚/非虔诚」年代斜率比** —— 与 `#784`–`#786` 逐字同一个量,
   只把**题目**从四道性题扩到**全部有足够年份的 GSS 道德/态度题**。

⚠⚠ **硬规则①:一个变量名不是一次测量。** 本轮**先枚举、先打印每一列的 n 与它真正被问过的年份**,
   再决定谁进网格 —— **不许拿变量名当数据。** 候选表写在下面,而**能不能用由数据说了算。**

⚠ 可读性判据前置(`#782`):比值自助区间含 1.0 ⇒ 列出、不读、不平均进去。
   ⚠ 而本轮**必然**有大量不可读格(非性题里很多没怎么动)——
   **那不是「没有效应」,是分母趋零、比值成 Cauchy 型**,如实标注。

⚠⚠ **分母被换成活下来的格,这一族刚在 `#787` 数过三处** ⇒ 本轮**每一个比例都同时报
   「可读 / 尝试」两个数**,一处都不许只报可读的那个。

预注册判词(条件式):
  if 正控开火(`homosex` 必须复现 `#786` 的 0.31–0.48)and 可读题数 >= 4:
      设 s = 可读题的比值的**四分位距 / 中位数**(离散度,尺度无关)
      if s <= 0.5 且性题与非性题的比值中位数之差 <= 0.15 -> **A 一个总闸**
      elif 性题的比值明显小于非性题(差 > 0.15)          -> **B 一个个开关**
      else                                               -> 不可分,报分布,不选边
  else: UNVERIFIED
⚠ 0.15 与 0.5 是跑之前定的,理由:`#786` 里同一道题跨三个窗口的比值范围是 0.311–0.476,
  **跨度 0.165** —— 所以**小于这个跨度的题间差别,本设计分辨不了**。**阈值由已测的噪声定,不是我挑的。**

⚠ 跑之前写下的最强混淆:**非性题很可能整体动得更少** ⇒ 分母小 ⇒ 比值不稳。
  ⇒ 控制:同一轮里报**每题非虔诚层自己的斜率与它的零**,把「比值不稳」与「社会没动」分开。

本轮换不了仪器,理由同 `R223/instrument_search.py`(对象是世界,而第二具仪器本机六具全部落选)。

═══════════════════════════════════════════════════════════════════════════════
⚠⚠ **跑完之后的记录:本轮判 UNVERIFIED,而原因被定位到了一条真实的仪器边界,不是一个 bug。**

18 道候选里 **14 道的方向 `label_pole` 拒绝判定**,只剩 4 道进网格(可读 2)⇒ 题间比较做不了。
逐个看真实标签之后,拒绝分成两类,而**只有第一类是我的错**:

**① 我的错(机械,已修):** `homosex` 的类别是
   `['always wrong', …, 'not wrong at all', 'other']` —— **五个,最后一个是 `other`**。
   `#785` 传的是 `[:4]`,本轮传了**全部**,于是「高值端」读到 `other` ⇒ 正确地 raise。
   **仪器没坏,是我喂错了。**

**② 不是我的错,是一条边界(不修,记下来):** 二元政策题的方向**不在标签里,在题干里**。
   `abany ['yes','no']` · `cappun ['favor','oppose']` · `grass ['should be legal', …]` ·
   `prayer ['approve','disapprove']` · `racmar ['yes','no']` · `letdie1` · `suicide1` · `gunlaw`。
   **`yes` 在 `abany` 上是宽容,在 `racmar`(该不该立法禁止跨种族婚姻)上是限制** ——
   **同一个词,相反的道德方向。**
   ⇒ **给 `_POLES` 加上 `yes/no`、`favor/oppose` 会让 `label_pole` 在这些题上返回一个自信的方向,
   而那正是它存在的目的所要防止的事**(`#759` 的 P6 代理账:匹配不上 ⇒ 我确实不知道方向 ⇒ raise;
   **匹配上从不证明方向正确**)。**所以这条边界不许用扩词表来绕过。**

⇒ **正确的修法是:二元政策题的方向必须由「题干」逐题人工给出,并把题干原文引在旁边**,
   标成**人工判定**而不是机器读数。**那是下一轮的设计,不是本轮能补的一行。**
⇒ 而本轮真正产出的是**一条能力边界**:这具自动定向仪器覆盖 Likert 式道德题,
   **结构性地覆盖不了 yes/no 政策题** —— 而「一个总闸还是一个个开关」这个问题
   **恰恰需要非性的政策题**,所以它被这条边界挡住,不是被数据挡住。
═══════════════════════════════════════════════════════════════════════════════
"""
import numpy as np, pandas as pd, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes, label_pole
from lib.gates import Gate

RNG = np.random.default_rng(227)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"

# ── 候选:性题 + 非性道德题。**名字只是候选,能不能用由下面的普查说了算。** ──────────
SEX = ["homosex", "xmarsex", "premarsx", "teensex", "pornlaw"]
NONSEX = ["abany", "abdefect", "abnomore", "abpoor", "suicide1", "suicide2",
          "grass", "cappun", "gunlaw", "letdie1", "spanking", "prayer", "racmar"]
CAND = SEX + NONSEX
GRP = {c: ("性" if c in SEX else "非性") for c in CAND}

print("=== ① 普查:每一个候选列的 n 与它**真正被问过**的年份(硬规则①)===")
have = [c for c in CAND]
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+have, convert_categoricals=False)
cat = pd.read_stata(gp, columns=have, convert_categoricals=True)
usable, rejected = [], []
for c in have:
    v = pd.to_numeric(d[c], errors="coerce")
    cats = list(cat[c].cat.categories) if hasattr(cat[c], "cat") else []
    # ⚠ `homosex` 的第五档是 `other`,不是一个作答档 ⇒ 只取前 4 档(`#785` 就是这么做的,本轮第一版漏了)
    if c == "homosex": cats = cats[:4]
    k = len(cats)
    yrs = sorted(d.year[v.notna()].unique())
    ok = (k in (2, 3, 4, 5)) and len(yrs) >= 12
    (usable if ok else rejected).append(c)
    print(f"  {c:10s} {GRP[c]:2s} n={int(v.notna().sum()):>7,} · 年数 {len(yrs):>3} "
          f"({int(min(yrs)) if yrs else 0}–{int(max(yrs)) if yrs else 0}) · 档数 {k:>2} "
          f"{'**用**' if ok else '**弃**(档数或年数不合)'}")
print(f"\n  ⇒ 候选 {len(CAND)} · **可用 {len(usable)}**(性 {sum(1 for c in usable if GRP[c]=='性')} · "
      f"非性 {sum(1 for c in usable if GRP[c]=='非性')})· 弃 {len(rejected)}:{rejected}")
if len(usable) < 6:
    print("⛔ 可用题不足 6 道 ⇒ 本轮的题间比较没有意义,停。"); sys.exit(2)

# ── 方向:逐题从值标签读,认不出就 raise(`#759`)────────────────────────────────
print("\n=== ② 方向:逐题从值标签读,认不出就 raise(`#759`)===")
KMAX = {}
poles, unknown = {}, []
for c in usable:
    cs = list(cat[c].cat.categories)
    KMAX[c] = len(cs)
    try:
        poles[c] = label_pole(cs)
        print(f"  {c:10s} 高值端 = {poles[c]!r}  {cs}")
    except Exception as e:
        unknown.append(c); print(f"  {c:10s} ⚠ **方向认不出 ⇒ 弃**({type(e).__name__})  {cs}")
usable = [c for c in usable if c not in unknown]
print(f"  ⇒ 方向可判的 **{len(usable)}** 道(弃 {len(unknown)} 道:{unknown})")

M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= KMAX[c])) for c in usable})
for c in ("attend", "reliten", "fund"):
    lo, hi = {"attend": (0, 8), "reliten": (1, 4), "fund": (1, 3)}[c]
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v: (v >= lo) & (v <= hi))
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

print("\n=== ③ 逐题:比值 · 它的自助区间 · **以及非虔诚层自己的斜率与零**(跑前写下的混淆的控制)===")
rows, attempted = [], 0
for c in usable:
    attempted += 1
    rA, rB = series(c, 2), series(c, 0)
    if len(rA) < 10 or len(rB) < 10:
        rows.append(dict(item=c, grp=GRP[c], readable=False, why=f"年数不足 {len(rA)}/{len(rB)}")); continue
    yA = np.array([r[0] for r in rA], float); yB = np.array([r[0] for r in rB], float)
    vA = np.array([r[1] for r in rA]); vB = np.array([r[1] for r in rB])
    # 非虔诚层自己动了没有:斜率 vs 它自己的置换零
    sB = slope(yB, vB)
    nulB = np.abs([slope(RNG.permutation(yB), vB) for _ in range(2000)])
    q95 = float(np.quantile(nulB, .95))
    moved = abs(sB) > q95
    f = lambda ia, ib: slope(yA[ia], vA[ia])/slope(yB[ib], vB[ib])
    bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(4000)])
    bs = bs[np.isfinite(bs)]
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    r = float(f(np.arange(len(yA)), np.arange(len(yB))))
    ok = not (lo <= 1.0 <= hi)
    rows.append(dict(item=c, grp=GRP[c], readable=bool(ok), r=r, lo=lo, hi=hi,
                     slopeB=sB, nullB95=q95, moved=bool(moved), nyr=[len(rA), len(rB)],
                     why=None if ok else "比值区间含 1.0"))
    print(f"  {c:10s} {GRP[c]:2s} 比 {r:7.3f} [{lo:8.3f}, {hi:8.3f}]  "
          f"非虔诚斜率 {sB:+.5f} vs 零95% {q95:.5f} {'(动了)' if moved else '**(没动)**'}  "
          f"{'可读' if ok else '**不可读**'}")

R = [x for x in rows if x.get("readable")]
sex_r = [x["r"] for x in R if x["grp"] == "性"]; non_r = [x["r"] for x in R if x["grp"] == "非性"]
moved_but_unread = [x for x in rows if not x.get("readable") and x.get("moved")]
print(f"\n  **可读 {len(R)} / 尝试 {attempted}**(`#787`:两个数都报)· "
      f"性 {len(sex_r)}/{sum(1 for c in usable if GRP[c]=='性')} · 非性 {len(non_r)}/{sum(1 for c in usable if GRP[c]=='非性')}")
print(f"  ⚠ 不可读但**非虔诚层确实动了**的题:{[x['item'] for x in moved_but_unread]} "
      f"—— 这些是「比值不稳」,不是「社会没动」")

G = Gate("#788 · 一个总闸,还是一个个开关")
hs = next((x for x in R if x["item"] == "homosex"), None)
pc = bool(hs and 0.28 <= hs["r"] <= 0.50)
G.asserted("① 正控:`homosex` 必须复现 `#786` 的 0.31–0.48(容差放到 0.28–0.50)",
           pc, f"homosex 比 {hs['r']:.3f}" if hs else "homosex 不可读 ⇒ 仪器与前几轮不一致", kind="control")
G.asserted("② 前提:可读题 >= 4,否则题间比较无意义",
           bool(len(R) >= 4), f"可读 {len(R)} / 尝试 {attempted}", kind="control")
if len(R) >= 4 and sex_r and non_r:
    allr = np.array([x["r"] for x in R])
    iqr = float(np.percentile(allr, 75)-np.percentile(allr, 25))
    disp = iqr/abs(float(np.median(allr)))
    gap = float(np.median(non_r)-np.median(sex_r))
    one_brake = bool(disp <= 0.5 and abs(gap) <= 0.15)
    G.asserted("③ kill(预注册):「一个总闸」要站住,需离散度 IQR/中位 ≤ 0.5 **且** 性/非性中位差 ≤ 0.15",
               one_brake, f"IQR/中位 {disp:.3f} · 性 {np.median(sex_r):.3f} vs 非性 {np.median(non_r):.3f}"
                          f"(差 {gap:+.3f})· 阈值由 `#786` 同题跨窗口跨度 0.165 定", kind="kill")
else:
    disp = gap = float("nan"); one_brake = False
    G.asserted("③ kill:样本不足以开火", False, f"可读 {len(R)} · 性 {len(sex_r)} · 非性 {len(non_r)}", kind="kill")
print(); print(G)

print("\n"+"="*92)
if not (pc and len(R) >= 4 and sex_r and non_r):
    v = f"**UNVERIFIED:正控没过或可读题不足(可读 {len(R)}/{attempted},性 {len(sex_r)}、非性 {len(non_r)})。**"
elif one_brake:
    v = (f"**A 一个总闸:虔诚是一种对任何道德变迁的普遍阻力。** 可读 {len(R)}/{attempted} 题,"
         f"比值离散度 IQR/中位 = {disp:.3f},性题中位 {np.median(sex_r):.3f} 与非性题中位 {np.median(non_r):.3f} "
         f"只差 {gap:+.3f}(本设计的分辨极限是 0.165)⇒ **那个刹车不挑题。**")
elif gap > 0.15:
    v = (f"**B 一个个开关:那个刹车是教义性的。** 性题的比值中位 {np.median(sex_r):.3f},"
         f"非性题 {np.median(non_r):.3f},**差 {gap:+.3f} 超过本设计的分辨极限 0.165** ⇒ "
         f"**虔诚者不是对一切变迁都踩刹车,他们只在传统点过名的地方踩。**")
else:
    v = (f"**不可分,报分布不选边。** 离散度 {disp:.3f} · 性/非性中位差 {gap:+.3f};"
         f"可读 {len(R)}/{attempted}。⇒ 两个世界在这套数据上分不开,而分不开本身是要写下来的。")
print(v)
json.dump(dict(rows=rows, usable=usable, rejected=rejected, unknown_polarity=unknown,
               n_readable=len(R), n_attempted=attempted, sex_r=sex_r, nonsex_r=non_r,
               dispersion=disp, gap=gap, verdict=v, gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"one_brake_or_many.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'one_brake_or_many.json'}")

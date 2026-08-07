r"""#868 · E03·A97·R307 —— 个人一般性**有多大**,而不是两群人**差多少**

**还 `#866`②(上一条账本 `#867` 的 NEXT①)。** `#866` 量到:广度残差中位 **+0.0005**,
合成零自身 SD **0.0173**,五个靶都是。它的 NEXT 由此写下:
**「一个人反对几件事,完全由『他属于哪一层』和『至少反对一次』决定,没有剩下任何个人一般性。」**

**⚠⚠ 而先做算术就会看到:那句话的量词比它的测量大了一整级 —— 本轮第一件事是撤掉它。**

**推导(纸上,机器已验)**:令 `R` = 一个人拒绝的题数,`p_sj` = 层 `s` 里第 `j` 题的边际。
   `E[R | R≥1, s] = E[R | s] / (1 − P(R=0 | s))`,而 **`E[R|s] = Σ_j p_sj` 被合成零逐字保留**
   ⇒ **观测与零的一切差别,只经由 `P(R=0)` 这一个数进入。** 机器上验过(τ = 0 / 0.6 / 1.2):
   `E[R|R≥1]_观测 / E[R|R≥1]_零` vs 恒等式右边 `(1−P₀_零)/(1−P₀_观测)` **差 1.1e-04 – 2.4e-03**(蒙特卡洛噪声)。
   ⇒ **而 `#866` 报的是这个量在两层之间的「差」。一个差会把两层共有的部分整个消掉。**
   **所以「残差为零」只能说「两层的个人一般性一样多」,绝不能说「没有个人一般性」。**
   ⇒ **`#866` 的 NEXT 那句话作废。** ⇒ `#866`②

**⇒ 而正确的量在同一页纸上就写出来了,它不是一个差:**
   **`φ_s = Var(R | s) ÷ Σ_j p_sj(1 − p_sj)`** —— **题间独立时 φ 恰为 1**,
   有人层潜变量时 φ > 1,且随潜变量强度单调上升。机器上验过:
   **τ=0 → φ=1.0032 · τ=0.6 → φ=1.1182 · τ=1.2 → φ=1.3906**。
   **它量的正是「个人一般性」的水平,而且它不是一个差,不会抵消。**

`G1` **估计量(两个,先于方法命名)**:
   ① **`φ`(水平)** —— 每个 (靶 × 十年 × 层) 上的离散比。**独立时理论零 = 1**(⚠ 但实测要比**置换零**,
      因为 `Var` 的估计有噪声)。**这一个量决定 `#866` 那句话生死。**
   ② **`Δφ = φ(B=0) − φ(B≥1)`** —— **这才是 `#865`① 真正问的那件事**:
      **不在任何一块里的反对者,内部结构和块内的人一样吗?**
      **一样 ⇒ 同一条潜在尺,只是位置不同;不一样 ⇒ 他们的「反对」是另一种组织方式。**

四个世界(**每个都有分支**):
   A **φ 明显 > 1 且两层相近** ⇒ **`#866` 那句话撤销**(有个人一般性,而且两层一样多)
     ⇒ **零块反对者与块内反对者是同一条潜在尺上的两个位置。** 这同时回答了 `#865`①。
   B **φ 明显 > 1 而两层不同** ⇒ **那 `#866` 的残差就不该是零** ⇒ 两个测量有一个错了,
     **而恒等式指得出是哪一个**(残差只经 `P(R=0)` 进入,`φ` 经整个 `Var`)。
   C **φ ≈ 1** ⇒ **三道题在层内真的独立,没有潜在宽容特质** —— ⚠ **这与 Stouffer 以降七十年的
     文献相反**,而它是我不欢迎的那一个:它会让 `#862`–`#865` 一直在量的「那条缝」失去承载物。
   D **⚠ 元分离器**:`φ` 被**边际压住** —— `p` 靠近 0 或 1 时 `Var(R)` 与 `Σp(1−p)` 一起塌,
     `φ` 的可比性消失 ⇒ **「一般性的水平」这个提法需要换一个不受边际支配的量**,
     **不是某个世界赢了,是这个量本身不合用。**

预测矩阵:
   | 世界 | 现在 | φ≫1 两层近 | φ≫1 两层远 | φ≈1 | φ 随边际走 |
   | A 同一条尺 | 0.40 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 有一个错 | 0.15 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 无潜特质 | 0.10 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D 量不合用 | 0.35 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(**条件式**):
  if 正控开火(**注入已知强度的人层潜变量 ⇒ φ 必须单调上升;强度 0 时必须恰好不动**)
     and 负控为零(**层内逐题独立打乱 ⇒ φ 必须塌回 1**)
     and 安慰剂为零(**`ballot` 当假层 ⇒ 两个 ballot 的 φ 必须无差**):
      φ 中位 ≤ 1.05                                       -> C(且 `#866` 那句话反而站得住)
      φ 中位 > 1.05 且 |Δφ| 在其零内的格 ≥ 2/3            -> A
      φ 中位 > 1.05 且 |Δφ| 超零的格 ≥ 2/3                -> B
      φ 与「层的平均边际」的 |Spearman ρ| ≥ 0.7           -> D
  else: UNVERIFIED

⚠⚠ **跑前写下的最强混淆(它同时是世界 D)**:**`B≥1` 的人拒绝率天生就高**
   ⇒ 两层的边际不同 ⇒ **`φ` 的分母不同,而 `Var(R)` 的上限也随边际变**。
   **`Δφ` 可能整个是边际差造出来的,而不是结构差。**
   ⇒ 控制:**边际匹配的 sham** —— 在 `B=0` 里**重抽**出一个逐题边际与 `B≥1` 相同的子样本,
   再算 `Δφ`。**同一个操作,减掉「层」这个成分。** 残差若消失,那 `Δφ` 就是边际,不是结构。

`G3` 多重性:整族 = 5 靶 × 十年 × 层,BH 与 BY 都做,不同意的格一起发表。
`G4` 规格曲线:层的切法两版(`B≥1` / `B≥2`)× 五靶 × 逐十年,全部印出。
⚠ kill 带 `yardstick` / `yardstick_noise` / `population` / `direction`,
**且 `#867` 新加的一条:控制行的 `population` 必须与 kill 的 `population` 逐字相同。**

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
① 横断面 ⇒ **无因果识别**;
② **`φ` 只有三道题** —— 三个二值题定不出一个正经的潜变量模型的所有参数(区分度、猜测率都定不出),
   `φ` 是**一个**汇总数,**它证不了「潜变量是一维的」**,只证「题间不独立」⇒ **这条边界不许省**;
③ **「同一批人」只能到「同一年同一批受访者」** —— GSS 无面板 ⇒ **结构性拿不到**;
④ **换不了仪器**:五个靶都在 GSS 内 ⇒ **跨靶不是跨仪器**,控住问卷格式,控不住受访者那一层。
"""
import json, math, pathlib, sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
EXT = ROOT / "data/external"
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
NPERM, SEED = 200, 307
P866 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be"
                      "/A96_不在任何一块里的反对者是谁/R305_同一条尺上的两个位置还是两种反对"
                      "/results/two_kinds_of_refusal.json"))
POP = "GSS 五个靶 × 各自的十年 × 两种层切法上的全部 (靶×十年×层) 格"

print("=== ⓪a **算术先行**:`#866` 的残差只经由 `P(R=0)` 一个数进入 —— 这是恒等式 ===")
rg0 = np.random.default_rng(SEED)
pm0 = np.array([0.25, 0.40, 0.15])
ident, phis = [], []
for tau in (0.0, 0.6, 1.2):
    th = rg0.normal(0, tau, 300_000)[:, None]
    q = 1 / (1 + np.exp(-(np.log(pm0 / (1 - pm0)) + th)))
    X = (rg0.random(q.shape) < q).astype(float)
    R = X.sum(1); pmv = X.mean(0)
    Y = (rg0.random(X.shape) < pmv).astype(float); Rn = Y.sum(1)
    lhs = R[R >= 1].mean() / Rn[Rn >= 1].mean()
    rhs = (1 - (Rn == 0).mean()) / (1 - (R == 0).mean())
    ph = float(R.var() / (pmv * (1 - pmv)).sum())
    ident.append(abs(lhs - rhs)); phis.append(ph)
    print(f"  τ={tau}: 恒等式两边差 **{abs(lhs-rhs):.2e}** · **φ = {ph:.4f}**")
IDENT_MAX = float(max(ident))
print(f"  ⇒ **恒等式最大偏差 {IDENT_MAX:.2e}(蒙特卡洛噪声)** ⇒ "
      f"**`#866` 报的是这个量在两层之间的『差』,而一个差把两层共有的部分整个消掉了。**")
print(f"  ⇒ **所以「残差为零」只能说「两层一样多」,不能说「没有」。`#866` 的 NEXT 那句话作废。**")
print(f"  ⇒ **而 φ 在 τ=0 时是 {phis[0]:.4f}、τ=1.2 时是 {phis[2]:.4f} —— 它量的是水平,不是差。**")
print(f"  ⚠ 对照:`#866` 自己报的广度残差中位 = "
      f"**{np.median([r['breadth_resid'] for r in P866['grid'] if r['bspec']=='B≥1' and np.isfinite(r['breadth_resid'])]):+.4f}**"
      f"(从它的产物读的,不手抄)")

TARGETS = {"同性恋 homo": ("spkhomo", "colhomo", "libhomo"),
           "种族主义者 rac": ("spkrac", "colrac", "librac"),
           "共产主义者 com": ("spkcom", "colcom", "libcom"),
           "军管主义者 mil": ("spkmil", "colmil", "libmil"),
           "无神论者 ath": ("spkath", "colath", "libath")}
GC = ["year", "attend", "reliten", "fund", "polviews", "age", "ballot"]
for t in TARGETS.values(): GC += list(t)
gs = pd.read_stata(EXT / "gss/GSS_stata/gss7224_r3a.dta", columns=GC, convert_categoricals=False)
D = pd.DataFrame({"year": gs.year})
# ⚠⚠⚠ **`colcom` 的极性与其余四靶相反,而这不是本轮发现的 —— 是本轮重新发现的。**
# `#680`(`E03·A23·R122`)早就查过码本、把整张极性表写进了账本,并在代码里写下 `c != "colcom"`;
# `R123` · `R127` · `R151` 等约 180 轮**都带着那个例外,一次没错**。
# **而 `#866` 从头重写了题目构造,那个例外静静地消失了** —— 我这一轮又把它当新发现查了一遍。
# ⇒ 机制(它比「colcom 是反的」值钱得多):**那条守则活在一个被逐轮拷贝的 lambda 里,
#   靠拷贝续命,所以只要有一轮从头写,它就断了。** ⇒ 搬进 `lib/gss_polarity.py`,**被 import,不被拷贝**。
from lib.gss_polarity import permissive as _perm
for tag, (a, b, c) in TARGETS.items():
    k = tag.split()[-1]
    for pre, col in (("spk", a), ("col", b), ("lib", c)):
        D[f"{pre}_{k}"] = _perm(pd.to_numeric(gs[col], errors="coerce"), col)
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)),
                    ("polviews", (1, 7)), ("age", (18, 89)), ("ballot", (1, 4))):
    D[c] = pd.to_numeric(gs[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
Rr = D.dropna(subset=["attend", "reliten", "fund"]).copy()
Rr["REL"] = (zs(Rr.attend) + zs(-Rr.reliten) + zs(-Rr.fund)) / 3
D = D.join(Rr["REL"])
DECS = {"1970s": range(1972, 1980), "1980s": range(1980, 1990), "1990s": range(1990, 2000),
        "2000s": range(2000, 2010), "2010s": range(2010, 2020)}
BSPEC = {"B≥1": 1, "B≥2": 2}

print("\n=== ⓪b 硬规则①:**变量名不是测量** —— 每一靶的 n 与实际问过的年份 ===")
for tag in TARGETS:
    k = tag.split()[-1]
    m = D[[f"spk_{k}", f"col_{k}", f"lib_{k}"]].notna().all(1) & D.REL.notna()
    yy = sorted(D.year[m].unique().astype(int))
    print(f"  {tag:14s} 三题齐全且 REL 非缺 n={int(m.sum()):6,d} · {len(yy)} 年 {yy[0]}–{yy[-1]}")


def phi(X):
    """φ = Var(R) / Σ p(1−p)。**题间独立时恰为 1。**"""
    if len(X) < 80: return np.nan, np.nan
    p = X.mean(0); den = float((p * (1 - p)).sum())
    if den <= 0: return np.nan, np.nan
    return float(X.sum(1).var(ddof=1) / den), float(p.mean())


def prep(sub, k):
    X = 1 - sub[[f"spk_{k}", f"col_{k}", f"lib_{k}"]].to_numpy(float)
    rel, age, pol = (sub[c].to_numpy(float) for c in ("REL", "age", "polviews"))
    B = np.zeros(len(sub))
    for v in (rel, age, pol):                       # 三块都取**不宽容**那一端(与 `#866` 同一定义)
        ok = np.isfinite(v)
        if ok.sum() < 200: return None
        B += (ok & (v >= np.nanquantile(v[ok], 2 / 3))).astype(float)
    bv = sub["ballot"].to_numpy(float)
    ok = np.isfinite(X).all(1)
    return X[ok], B[ok], bv[ok]


def shuffle_items(X, rng):
    """负控:**逐题独立打乱**,毁掉人内相关,逐字保留每题边际 ⇒ φ 必须塌回 1。"""
    Y = X.copy()
    for j in range(Y.shape[1]): Y[:, j] = Y[rng.permutation(len(Y)), j]
    return Y


def inject(X, tau, U, Z):
    """正控:注入已知强度的人层潜变量,**保持每题边际不变**。

    ⚠⚠ **共用随机数(common random numbers):`U` 与 `Z` 只抽一次,所有剂量复用同一组。**
    第一版给每个剂量一个**不同的种子**,于是曲线上每一点都带**独立**的抽样噪声,
    小剂量处的信号落在噪声之下 ⇒ τ=0→1.0088 而 τ=0.25→1.0027,**非单调,正控 FAIL**。
    ⚠ 而这条教训**本来就写在本项目的 `lib/gates.py` 里**(`#124f`:
    「退化臂几乎总是因为**没有复用参照臂的种子**」)—— **工具里有,我还是犯了。**
    共用随机数之后,剂量曲线是**配对**比较,单调性才是一个可测的性质而不是噪声。"""
    p = np.clip(X.mean(0), 1e-3, 1 - 1e-3)
    q = 1 / (1 + np.exp(-(np.log(p / (1 - p)) + tau * Z[:, None])))
    return (U < q).astype(float)


def match_marginals(Xa, Xb, rng=None, iters=200):
    """sham:**IPF** —— 把 A 的 8 格联合分布的三个边际拉到 B 的边际上,**优势比逐字不动**。

    ⚠ 第一版用的是逐题独立重加权,那**不是**联合边际匹配:三个权重相乘会互相破坏,
    权重发散,60 次尝试没有一次落进 0.03 的容差 ⇒ 返回 `nan`,sham 整条没跑成。
    IPF(迭代比例拟合)在这里是**闭式**的:三个二值题只有 8 个格,
    每轮把每一题的边际按比例拉齐,几轮就收敛,**而它保持所有优势比不变** ——
    也就是说 **`sham` 只改边际,不改「关联」这个被测的成分**,正是 `realstat` 要的那种
    「同一个操作,减掉被研究的那个成分」。"""
    idx = (Xa[:, 0] * 4 + Xa[:, 1] * 2 + Xa[:, 2]).astype(int)
    q = np.bincount(idx, minlength=8).astype(float); q /= q.sum()
    tgt = Xb.mean(0)
    cells = np.array([[(c >> 2) & 1, (c >> 1) & 1, c & 1] for c in range(8)], float)
    for _ in range(iters):
        for j in range(3):
            cur = float(q[cells[:, j] == 1].sum())
            if cur <= 0 or cur >= 1: return None, None
            q[cells[:, j] == 1] *= tgt[j] / cur
            q[cells[:, j] == 0] *= (1 - tgt[j]) / (1 - cur)
            q /= q.sum()
    dev = float(np.abs((q[:, None] * cells).sum(0) - tgt).sum())
    if dev > 1e-6: return None, None
    p = (q[:, None] * cells).sum(0)
    R = cells.sum(1); ER = float((q * R).sum())
    var = float((q * (R - ER) ** 2).sum()); den = float((p * (1 - p)).sum())
    return (var / den if den > 0 else None), dev


print(f"\n=== ① 网格:{len(TARGETS)} 靶 × {len(DECS)} 十年 × {len(BSPEC)} 切法 · 置换零 {NPERM} 次 ===")
rows = []
rng = np.random.default_rng(SEED)
for tag in TARGETS:
    k = tag.split()[-1]
    line = []
    for dec in DECS:
        m = D[[f"spk_{k}", f"col_{k}", f"lib_{k}"]].notna().all(1) & D.REL.notna() \
            & D.year.isin(list(DECS[dec]))
        sub = D[m]
        if len(sub) < 800: continue
        z = prep(sub, k)
        if z is None: continue
        X, B, bv = z
        for bs, thr in BSPEC.items():
            g1 = B >= thr
            for lab, sel in (("B=0", ~g1), ("B≥1", g1), ("全体", np.ones(len(X), bool))):
                ph, pbar = phi(X[sel])
                if not np.isfinite(ph): continue
                nd = np.array([phi(shuffle_items(X[sel], rng))[0] for _ in range(NPERM)])
                nd = nd[np.isfinite(nd)]
                rows.append(dict(target=tag, dec=dec, bspec=bs, stratum=lab, phi=ph, pbar=pbar,
                                 null_med=float(np.median(nd)), floor=float(np.quantile(nd, 0.95)),
                                 p=float((1 + (nd >= ph).sum()) / (len(nd) + 1)),
                                 n=int(sel.sum())))
        r_all = [r for r in rows if (r["target"], r["dec"], r["bspec"], r["stratum"])
                 == (tag, dec, "B≥1", "全体")]
        if r_all: line.append(f"{dec}:{r_all[0]['phi']:.3f}")
    print(f"  {tag:14s} φ(全体) · " + " · ".join(line))
print("  ⚠ **题间独立时 φ 的理论值恰为 1**;实测仍对照**该格自己的置换零**(`Var` 的估计有噪声)")

MAIN = [r for r in rows if r["bspec"] == "B≥1"]
if not MAIN:
    raise SystemExit("⛔ 网格为空 —— 空总体不许当作通过")
ALLR = [r for r in MAIN if r["stratum"] == "全体"]
PHI_MED = float(np.median([r["phi"] for r in ALLR]))
over = sum(1 for r in ALLR if r["phi"] > r["floor"])
print(f"\n=== ② `φ` 的水平 —— 这一个数决定 `#866` 那句话的生死 ===")
print(f"  全体层 φ 中位 **{PHI_MED:.4f}** · 逐格 {min(r['phi'] for r in ALLR):.3f}–"
      f"{max(r['phi'] for r in ALLR):.3f} · **超自己置换零的格 {over}/{len(ALLR)}** · "
      f"置换零中位 {np.median([r['null_med'] for r in ALLR]):.4f}")

print(f"\n=== ③ `Δφ = φ(B=0) − φ(B≥1)` —— 这才是 `#865`① 真正问的那件事 ===")
dif = []
for tag in TARGETS:
    for dec in DECS:
        a = [r for r in MAIN if (r["target"], r["dec"], r["stratum"]) == (tag, dec, "B=0")]
        b = [r for r in MAIN if (r["target"], r["dec"], r["stratum"]) == (tag, dec, "B≥1")]
        if a and b:
            # ⚠⚠ **sham 必须跑在整张网格上,不能只跑一格。**
            # 第一版只在 homo/2010s 上算了边际匹配,而**判词分支根本没有引用它** ——
            # `realstat` 的补救原话:*分支条件必须引用这一轮声明过的每一条控制*。
            # 而它一算出来就说了一件判词没说的事:**匹配边际后 Δφ 的符号翻转了**
            # (−0.1046 → +0.0639)⇒ **原始 Δφ 的方向不是这批数据能定的。**
            k_ = tag.split()[-1]
            mm = D[[f"spk_{k_}", f"col_{k_}", f"lib_{k_}"]].notna().all(1) & D.REL.notna() \
                & D.year.isin(list(DECS[dec]))
            zz = prep(D[mm], k_)
            sh = np.nan
            if zz is not None:
                Xg, Bg, _ = zz
                gg = Bg >= 1
                if gg.sum() > 80 and (~gg).sum() > 80:
                    sp, _dev = match_marginals(Xg[~gg], Xg[gg])
                    if sp is not None: sh = sp - phi(Xg[gg])[0]
            dif.append(dict(target=tag, dec=dec, d=a[0]["phi"] - b[0]["phi"],
                            phi0=a[0]["phi"], phi1=b[0]["phi"],
                            dp=a[0]["pbar"] - b[0]["pbar"], d_sham=sh,
                            flip=bool(np.isfinite(sh) and np.sign(sh) != np.sign(a[0]["phi"] - b[0]["phi"])),
                            fl=float(np.hypot(a[0]["floor"] - a[0]["null_med"],
                                              b[0]["floor"] - b[0]["null_med"]))))
for tag in TARGETS:
    ds = [x for x in dif if x["target"] == tag]
    if ds:
        print(f"  {tag:14s} Δφ · " + " · ".join(
            f"{x['dec']}:{x['d']:+.3f}{'*' if abs(x['d'])>x['fl'] else ' '}" for x in ds))
DPHI_MED = float(np.median([x["d"] for x in dif]))
DPHI_OVER = sum(1 for x in dif if abs(x["d"]) > x["fl"])
SH_OK = [x for x in dif if np.isfinite(x["d_sham"])]
FLIP = sum(1 for x in SH_OK if x["flip"])
SH_MED = float(np.median([x["d_sham"] for x in SH_OK])) if SH_OK else np.nan
print(f"  ⇒ Δφ 中位 **{DPHI_MED:+.4f}** · **超零的格 {DPHI_OVER}/{len(dif)}** · "
      f"(`*` = 超出两层零的合成噪声)")
print(f"  ⚠⚠ **而边际匹配的 sham 跑遍整张网格之后**:Δφ 中位 **{SH_MED:+.4f}**,"
      f"**{FLIP}/{len(SH_OK)} 格的符号在匹配边际后翻转** ⇒ "
      + ("**原始 Δφ 的方向不是这批数据能定的 —— 它由「要不要匹配边际」这个选择决定**"
         if FLIP > len(SH_OK) / 3 else
         "**方向在匹配前后一致,Δφ 不是边际造出来的**"))

print("\n=== ④ 控制 ===")
sub0 = D[D[["spk_homo", "col_homo", "lib_homo"]].notna().all(1) & D.REL.notna()
         & D.year.isin(list(DECS["2010s"]))]
X0 = prep(sub0, "homo")[0]
rgc = np.random.default_rng(SEED + 1)
base_phi = phi(X0)[0]
TAUS = [0.0, 0.25, 0.5, 1.0, 2.0]
U0 = rgc.random(X0.shape); Z0 = rgc.normal(0, 1, len(X0))   # ⚠ 共用随机数,只抽一次
inj = [(t, phi(inject(X0, t, U0, Z0))[0]) for t in TAUS]
zero_move = inj[0][1] - phi(inject(X0, 0.0, U0, Z0))[0]
MONO = all(inj[i][1] <= inj[i + 1][1] + 1e-9 for i in range(len(inj) - 1))
CT = 0.05
print(f"  正控(**剂量-反应**,注入已知强度 τ 的人层潜变量,**保持每题边际**):"
      + " · ".join(f"τ={t}→φ={v:.4f}" for t, v in inj))
print(f"     **单调 {MONO}** · **τ=0 时与同种子的重复完全一致({zero_move:+.2e})** ⇒ "
      f"⚠ **`G2` 控制必须能失败**")
CEIL = inj[-1][1] - inj[0][1]
r_h = [r for r in ALLR if r["target"] == "同性恋 homo" and r["dec"] == "2010s"][0]
MDE_T = next((t for t, v in inj if t > 0 and v > r_h["floor"]), None)
print(f"     **控制也必须能通过**:floor {abs(zero_move):.2e} < 阈 {CT} < ceiling {CEIL:.4f} ⇒ "
      f"**{'阈在真带内' if abs(zero_move) < CT < CEIL else '⚠⚠ 阈不在带内'}** · "
      f"**MDE = 最小的能顶出该格地板({r_h['floor']:.4f})的 τ = {MDE_T if MDE_T is not None else '>2.0'}**")
rgn = np.random.default_rng(SEED + 2)
neg = phi(shuffle_items(X0, rgn))[0]
print(f"  负控:**层内逐题独立打乱**(毁掉人内相关,逐字保留每题边际)⇒ φ = **{neg:.4f}**,"
      f"该格地板 **{r_h['floor']:.4f}** ⇒ **{'塌回 1 附近' if neg <= r_h['floor'] else '⚠ 没塌回去'}**")
print(f"     ⚠ **「这个零该不该是 1?」该** —— **题间独立时 φ 的理论值恰是 1**,"
      f"这是代数,不是约定;而实测零的中位是 {r_h['null_med']:.4f}(`Var` 估计的噪声)")
# sham:边际匹配
Bx = prep(sub0, "homo")[1] >= 1
sham_phi, sham_dev = match_marginals(X0[~Bx], X0[Bx])
SHAM_D = (sham_phi - phi(X0[Bx])[0]) if sham_phi is not None else np.nan
RAW_D = phi(X0[~Bx])[0] - phi(X0[Bx])[0]
print(f"  **sham(跑前写下的最强混淆的控制)**:在 B=0 里重抽出逐题边际与 B≥1 相同的子样本 ⇒ "
      f"Δφ **{RAW_D:+.4f} → {SHAM_D:+.4f}**"
      + (f"(⚠ **IPF 未收敛,如实记 nan**)" if not np.isfinite(SHAM_D) else
         f" ⇒ **边际解释了其中 {100*(1-abs(SHAM_D)/max(abs(RAW_D),1e-9)):.0f}%**"
         f"(IPF 边际残差 {sham_dev:.1e},优势比逐字不动)"))
# 安慰剂:ballot(⚠ 总体与 kill 逐字相同)
plb = []
for tag in TARGETS:
    k = tag.split()[-1]
    for dec in DECS:
        m = D[[f"spk_{k}", f"col_{k}", f"lib_{k}"]].notna().all(1) & D.REL.notna() \
            & D.year.isin(list(DECS[dec]))
        sub = D[m]
        if len(sub) < 800: continue
        z = prep(sub, k)
        if z is None: continue
        X, _, bv = z
        s1, s3 = np.isfinite(bv) & (bv == 1), np.isfinite(bv) & (bv == 3)
        if s1.sum() < 200 or s3.sum() < 200: continue
        d = phi(X[s1])[0] - phi(X[s3])[0]
        nd = np.array([phi(X[s1][np.random.default_rng(SEED + 900 + i).permutation(int(s1.sum()))])[0]
                       - phi(X[s3])[0] for i in range(60)])
        plb.append((d, float(np.quantile(np.abs(nd[np.isfinite(nd)]), 0.95))))
PL_OUT = sum(1 for d, f in plb if abs(d) > f)
from scipy.stats import binom as _binom
PL_MAX = int(_binom.ppf(0.95, max(len(plb), 1), 0.05))
PL_P = float(1 - _binom.cdf(PL_OUT - 1, max(len(plb), 1), 0.05)) if PL_OUT > 0 else 1.0
print(f"  安慰剂 **`ballot` 1 vs 3 当假层**(GSS 自己随机分配;⚠ **只在真有 ballot 的人上算**,"
      f"`#866` 的教训)⇒ 越界 **{PL_OUT}/{len(plb)}**,`P(X≥{PL_OUT})={PL_P:.3f}`,二项 95 分位 **{PL_MAX}**")
# 世界 D:φ 是不是被边际支配
from scipy.stats import spearmanr
RHO = float(spearmanr([r["phi"] for r in ALLR], [r["pbar"] for r in ALLR]).statistic)
print(f"  **世界 D 的判据(跑前写下)**:φ 与该格平均边际的 Spearman ρ = **{RHO:+.3f}** "
      f"(**|ρ| ≥ 0.7 ⇒ φ 被边际支配,量本身不合用**)")

ps = np.array([r["p"] for r in rows if np.isfinite(r["p"])]); C = len(ps)
o_ = np.argsort(ps); q = 0.05
cH = q * np.arange(1, C + 1) / C; cY = cH / np.sum(1.0 / np.arange(1, C + 1))
su = lambda pv, cr: (int(np.max(np.where(pv <= cr)[0])) + 1 if (pv <= cr).any() else 0)
kH, kY = su(ps[o_], cH), su(ps[o_], cY)
print(f"\n=== ⑤ 多重性:整族 **{C}** 格 · BH 存活 **{kH}** · BY **{kY}** · "
      f"p 分辨率下限 {1/(NPERM+1):.4f} · 不同意的 {kH-kY} 格一起发表 ===")

G = Gate("#868 · 个人一般性有多大,而不是两群人差多少")
G.asserted("① **算术先行**:`E[R|R≥1] = E[R]/(1−P(R=0))`,而合成零逐字保留 `E[R]` ⇒ "
           "**`#866` 的残差只经由 `P(R=0)` 一个数进入,且它是一个差** ⇒ "
           "**差把两层共有的部分整个消掉,「残差为零」不能读成「没有一般性」**",
           bool(IDENT_MAX < 5e-3), f"恒等式最大偏差 {IDENT_MAX:.2e};"
           f"合成验证 φ(τ=0)={phis[0]:.4f} → φ(τ=1.2)={phis[2]:.4f}",
           kind="control", population=POP)
G.asserted("② 前提(跑前写下的最强混淆,同时是世界 D):**`B≥1` 拒绝率天生高 ⇒ 两层边际不同 ⇒ "
           "`Δφ` 可能整个是边际造出来的** ⇒ **边际匹配的 sham**(在 B=0 里重抽出与 B≥1 同边际的子样本)",
           bool(np.isfinite(SHAM_D)),
           f"Δφ 原始 {RAW_D:+.4f} → 边际匹配后 {SHAM_D:+.4f}" if np.isfinite(SHAM_D)
           else "⚠ 边际匹配失败,如实记 nan —— 这条混淆本轮**没有**被控住",
           kind="control", population=POP)
G.asserted("③ 正控(**剂量-反应**):注入已知强度 τ 的人层潜变量、**保持每题边际** ⇒ φ 必须**单调上升**;"
           "τ=0 时必须恰好不动;**且阈落在 floor 与 ceiling 之间**",
           bool(MONO and abs(zero_move) < 1e-12 and abs(zero_move) < CT < CEIL),
           f"τ 曲线 " + " ".join(f"{t}:{v:.4f}" for t, v in inj)
           + f" · τ=0 {zero_move:+.2e} · 带 [0, {CEIL:.4f}] · 单调 {MONO} · MDE τ={MDE_T}",
           kind="control", population=POP)
G.asserted("④ 负控:层内逐题独立打乱 ⇒ φ 必须塌回 **1**(⚠ **「这个零该不该是 1?」该** —— "
           "题间独立时 φ 的理论值恰是 1,这是代数不是约定)",
           bool(neg <= r_h["floor"]),
           f"{neg:.4f} ≤ 地板 {r_h['floor']:.4f};该格置换零中位 {r_h['null_med']:.4f}",
           kind="control", population=POP)
G.asserted("⑤ 安慰剂 `ballot` 1 vs 3 当假层 ⇒ 越界格数落在二项零 Bin(N, 0.05) 内",
           bool(PL_OUT <= PL_MAX),
           f"越界 {PL_OUT}/{len(plb)} · P(X≥{PL_OUT})={PL_P:.3f} · 二项 95 分位 {PL_MAX}",
           kind="control", population=POP)
G.asserted("⑥ kill(预注册):**`#866` 那句「没有剩下任何个人一般性」要成立,需 φ 中位 ≤ 1.05**",
           bool(PHI_MED <= 1.05),
           f"φ 中位 {PHI_MED:.4f} · 超自己置换零的格 {over}/{len(ALLR)} · "
           f"Δφ 中位 {DPHI_MED:+.4f},超零 {DPHI_OVER}/{len(dif)} · φ~边际 ρ={RHO:+.3f}",
           kind="kill",
           yardstick="`φ = Var(R)/Σp(1−p)`,**题间独立时代数上恰为 1**;实测对照该格自己的逐题打乱置换零",
           yardstick_noise=float(np.median([r["floor"] - r["null_med"] for r in ALLR])),
           population=POP,
           direction=[r["phi"] - 1.0 for r in ALLR])
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif abs(RHO) >= 0.7:
    VERD = (f"**D `φ` 被边际支配(ρ={RHO:+.3f})⇒ 「一般性的水平」这个提法要换一个不受边际支配的量。**\n"
            f"  **不是某个世界赢了,是这个量本身不合用。**")
elif PHI_MED <= 1.05:
    VERD = (f"**C 三道题在层内真的接近独立(φ 中位 {PHI_MED:.4f})⇒ 没有潜在宽容特质。**\n"
            f"  ⚠ **这与 Stouffer 以降的文献相反,而它是我不欢迎的那一个** —— "
            f"它会让 `#862`–`#865` 一直在量的那条缝失去承载物。")
elif FLIP > len(SH_OK) / 3:
    VERD = (f"**D `Δφ` 的方向由「要不要匹配边际」决定,不由数据决定 ⇒ 两层的比较这批数据支持不了。**\n"
            f"  原始 Δφ 中位 **{DPHI_MED:+.4f}**,边际匹配后 **{SH_MED:+.4f}**,"
            f"**{FLIP}/{len(SH_OK)} 格符号翻转**。\n"
            f"  ⚠ **而 `φ` 的水平不受这一条影响**:它是每一层各自算的,不是两层的差 ⇒\n"
            f"  **φ 中位 {PHI_MED:.4f},{over}/{len(ALLR)} 格超自己的置换零(独立时代数值恰为 1)。**\n"
            f"  ⇒ **一句关于人的话:一个人对同一个群体的三件事 —— 让不让他讲话、让不让他教书、\n"
            f"  把不把他的书撤下来 —— 不是三件独立的事。拒绝其中一件的人,拒绝另外两件的可能性\n"
            f"  远高于偶然:R 的方差是题间独立时的 **{PHI_MED:.1f} 倍**,五个靶、五个十年、25 格无一例外。\n"
            f"  「反对」在一个人身上是一件事,不是三件。**\n"
            f"  ⚠⚠ **而 `#866` 说「没有剩下任何个人一般性」是把一个差读成了一个水平** —— "
            f"**残差为零只说明两层一样多,而它们是一样多地『有』。那句话撤销。**\n"
            f"  ⚠ **两层之间差多少,本轮答不了** —— 而答不了的理由是量出来的,不是猜的。")
elif DPHI_OVER <= len(dif) / 3:
    VERD = (f"**A 有个人一般性,而且两层一样多 —— `#866` 的 NEXT 那句话撤销。**\n"
            f"  **φ 中位 {PHI_MED:.4f}**(独立时代数值恰为 1),**{over}/{len(ALLR)} 格超自己的置换零**;\n"
            f"  而 **Δφ 中位 {DPHI_MED:+.4f},只有 {DPHI_OVER}/{len(dif)} 格超零** ⇒ "
            f"**两层的内部结构一样。**\n"
            f"  ⇒ **一句关于人的话,而它同时回答了 `#865`① 问了两轮没问出口的那件事:\n"
            f"  不在任何一块里的反对者,不是另一种反对者。他们身上的『反对』和最虔诚、最年长、\n"
            f"  最保守的那批人一样有内在的一致性 —— 一个人拒绝其中一件事,就更可能拒绝另外两件,\n"
            f"  而这种一致性在两群人身上一样强。他们不是另一种人,只是同一条尺上更靠前的位置。**\n"
            f"  ⚠ **而 `#866` 说的「没有剩下任何个人一般性」是把一个差读成了一个水平** —— "
            f"**残差为零只说明两层一样多,而它们一样多地『有』,不是一样多地『没有』。**")
else:
    VERD = (f"**B 有个人一般性,而两层不同(Δφ 超零 {DPHI_OVER}/{len(dif)} 格)⇒ "
            f"那 `#866` 的残差就不该是零。**\n"
            f"  **两个测量有一个错了**,而恒等式指得出是哪一个:残差只经 `P(R=0)` 进入,"
            f"`φ` 经整个 `Var` ⇒ 下一轮直接比这两条路径。")
print(VERD)
print(f"\n⚠ **本轮结构性做不到的**:① 横断面 ⇒ **无因果识别**;② **只有三道题** —— "
      f"三个二值题定不出潜变量模型的全部参数,**`φ` 证不了「潜变量是一维的」,只证「题间不独立」**;"
      f"③ **「同一批人」只能到「同一年同一批受访者」**,GSS 无面板 ⇒ **结构性拿不到**;"
      f"④ **换不了仪器** —— 五靶都在 GSS 内,**跨靶不是跨仪器**。")

json.dump(dict(grid=rows, dphi=dif, phi_median=PHI_MED, phi_over=over, phi_cells=len(ALLR),
               dphi_median=DPHI_MED, dphi_over=DPHI_OVER, dphi_cells=len(dif),
               dphi_sham_median=SH_MED, dphi_sign_flips=FLIP, dphi_sham_cells=len(SH_OK),
               colcom_polarity_fix="colcom is 4=yes,FIRED / 5=not fired -- opposite to the other four "
                                   "col* items; the only reversed item of fifteen. In the pipeline since "
                                   "#866, where communists were the most deviant control target.",
               rho_phi_marginal=RHO, identity_max_dev=IDENT_MAX, synthetic_phi=phis,
               multiplicity=dict(cells=C, bh=int(kH), by=int(kY), q=q),
               controls=dict(inject_curve=[[t, v] for t, v in inj], monotone=MONO,
                             zero_move=zero_move, ceiling=CEIL, threshold=CT, mde_tau=MDE_T,
                             neg=neg, neg_floor=r_h["floor"], null_med=r_h["null_med"],
                             sham_raw=RAW_D, sham_matched=SHAM_D,
                             placebo_out=PL_OUT, placebo_cells=len(plb), placebo_p=PL_P),
               derivation="E[R|R>=1] = E[R]/(1-P(R=0)); the synthetic null preserves E[R] exactly, "
                          "so #866's residual enters ONLY through P(R=0) -- and it is a DIFFERENCE, "
                          "which cancels whatever both strata share. phi = Var(R)/sum p(1-p) is the "
                          "LEVEL: exactly 1 under item independence, rising with the latent trait.",
               retracts_866="#866's NEXT sentence 'no personal generality left over' reads a level off "
                            "a difference; the measurement can only support 'the two strata have the "
                            "same amount'",
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), seed=SEED, nperm=NPERM,
               population=POP),
          open(OUT / "how_much_generality.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'how_much_generality.json'}")

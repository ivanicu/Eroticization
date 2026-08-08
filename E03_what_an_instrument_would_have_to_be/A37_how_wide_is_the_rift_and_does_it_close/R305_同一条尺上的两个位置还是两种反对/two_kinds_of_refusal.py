r"""#866 · E03·A96·R305 —— 同一条尺上的两个位置,还是两种不同的反对?

**还 `#865`①(上一条账本记录的 NEXT)。** `#865` 量出:三块「不宽容侧」的少数加起来,
只把「不放行」的人抬到基率的 **1.28 倍** ⇒ **大部分不放行的人不在任何一块里。他们是谁?**
`#865` 自己写下:**不能再靠加轴解决**(`#862` 已排过七根轴)⇒ **换问法。**

**⚠⚠ 算术先行 —— 这一轮的severity全在这里,因为最直觉的那个问法是被强制的:**
「落在零块里的反对者,和落在块里的反对者,反对得一样多吗?」
**这个问题的答案被两件事强制:**① 块本来就预测「反不反对」(那正是 `#862`–`#865` 量的东西);
② 只看「至少反对一次」的人,**是在结局上做条件**(`realstat` 点名的 collider / Oldham 陷阱)。
**两者叠加会强制出一个非零差,不管有没有第二种反对。**
⇒ 所以**广度差必须减掉它自己的合成零**:一个「每题独立抽、但保留每个块的逐题概率」的世界,
   **同样施加 R≥1 的条件**。**只有残差是测量。**

**⇒ 而真正干净的那个估计量来自一条定理,不是来自我的直觉:**
**Rasch 条件似然定理** —— 单维 1PL 模型下,**给定总分,哪一道题被拒与被试的严厉程度 θ 无关**
(θ 是条件分布的充分统计量已消去)。
本轮已在机器上验过(n=4e5,θ 三档跨 3.5 logit):
`P(哪一题 | 恰好拒一题)` = **[0.5666, 0.3088, 0.1246] / [0.5662, 0.3082, 0.1255] / [0.5666, 0.3085, 0.1250]**
—— **三档到小数点后三位一致。**
⇒ **所以:如果「零块反对者」与「块内反对者」只是同一条宽容度尺上的两个位置,
   那么在「恰好拒一题」的人里,他们拒的是哪一题,分布必须相同。**
   **这个零不是我假设的,是定理推出来的;偏离它就是单维性被证伪。**

`G1` **估计量(三个,先于方法命名)**:
   ① **`TVD_profile`**(主):在**恰好拒一题**的人里,`B=0` 与 `B≥1` 两群人「拒的是哪一题」
      分布之间的总变差距离。**Rasch 下应为 0(至抽样噪声)。**
      ⚠ TVD 有正偏 ⇒ **对照它自己的置换零(在这群人里打乱 B 标签),不对照 0。**
   ② **`Δbreadth`**:`E[R | R≥1, B=0] − E[R | R≥1, B≥1]`,**减掉合成独立零下的同一个量**。
   ③ **`Δmoral`**(最像人的那个):在**反对者内部**,`B=0` 与 `B≥1` 的道德谴责差
      (`homosex`,高=宽容)。**⚠ 它部分被块的定义强制**(虔诚块本来就更谴责)⇒
      **同样减掉合成零**,并**只在 homo 靶上可算**(别的靶没有对应的道德题)。

**块的定义(⚠ 与 `#865` 不同,且这个不同必须说出来)**:
   `#865` 问的是**哪一块是孤立的**(政治那块是**自由派**);**本轮问的是哪一块是不宽容的**,
   所以政治取的是**保守**那一端。**同一根轴、不同的三分之一 —— 不是矛盾,是两个不同的问题。**
   ⇒ `B` = 一个人落在 {虔诚最高三分位 · 年龄最高三分位 · polviews 最保守三分位} 里的个数(0–3)。
   **块的定义不看任何结局** ⇒ 跨靶比较才不循环。

四个世界(**每个都有分支**,`#856`):
   A **一条尺两个位置**:`TVD` 落在自己的置换零内,`Δbreadth` 残差 ≈ 0
     ⇒ **「块」只是标签,零块反对者只是没那么严厉的同一种人。**
   B **两种反对**:`TVD` 显著超零 ⇒ **单维性被证伪,零块反对者是另一种东西。**
   C **零块的反对是噪声**:他们的模式最贴近「每题独立抽」的合成零,而块内反对者偏离它
     ⇒ **大部分「反对」根本不是立场,是作答噪声/非态度。**
     ⚠ **这是我不欢迎的那一个** —— 它把 `#862`–`#865` 量到的一切降级为「量了一小撮人」。
   D **⚠ 元分离器**:五个靶(同性恋 · 种族主义者 · 共产主义者 · 军管主义者 · 无神论者)
     **给出同样大小的 TVD** ⇒ 那不是关于这个对象的发现,**是问卷格式的性质(作答定势)**
     ⇒ **「块 / 零块」这个分解本身没有碰到内容。**

预测矩阵:
   | 世界 | 现在 | TVD 在零内 | TVD 超零且只在 homo | TVD 超零且五靶相同 | 零块最贴合成零 |
   | A 一条尺   | 0.30 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 两种反对 | 0.25 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 噪声     | 0.20 | 0.05 | 0.05 | 0.05 | **0.85** |
   | D 作答定势 | 0.25 | 0.05 | 0.05 | **0.85** | 0.05 |

预注册判词(**条件式,不是阈值**):
  if 正控开火(**把一部分 B=0 的单拒者强制改拒某一题 ⇒ TVD 必须涨;plant=0 时必须恰为 0;
     且阈落在 floor 与 ceiling 之间**)
     and 负控为零(**在单拒者里打乱 B 标签 ⇒ TVD 塌回置换零**)
     and 安慰剂为零(**`ballot` 当假块 ⇒ TVD ≈ 置换零**):
      homo 靶上 TVD 超零的十年 ≤ 1/3,且 |Δbreadth 残差| 在零内      -> A
      homo 靶上 TVD 超零 ≥ 2/3 的十年,而其余四靶 < 1/3               -> B
      五靶 TVD 超零的比例都 ≥ 2/3                                     -> D
      零块群与合成零的距离 < 块内群的一半                              -> C
  else: UNVERIFIED

⚠⚠ **跑前写下的最强混淆:作答定势(acquiescence / response set)。**
   三道题同一个题干、连着问 ⇒ **一个只会一路说「不」的人会制造出「广度」,而里面没有任何态度内容。**
   ⇒ 控制:**同一具仪器、同一个格式、换靶** —— GSS 把这三道题原封不动问过
   **种族主义者 · 共产主义者 · 军管主义者 · 无神论者**。
   **如果是作答定势,五个靶必须给出同样的 TVD;如果是关于这个对象的,homo 必须与其余四靶分开。**
   **这是本轮唯一能把「关于世界」和「关于问卷」分开的东西**(硬规则②:说出claim路由经过的仪器)。

`G3` 多重性:整族 = 5 靶 × 若干十年 × 2 统计量,BH 与 BY 都做,不同意的格一起发表。
`G4` 规格曲线:`B≥1` vs `B≥2` 两种切法 × 五靶 × 逐十年,全部印出。
⚠ kill 带 `yardstick` / `yardstick_noise` / `population` / `direction`。

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
① 横断面 ⇒ **无因果识别**;
② **「同一种人」只能到「同一年同一批受访者的作答模式」** —— GSS 无面板 ⇒ **结构性拿不到**;
③ **Rasch 只是一个模型** —— 偏离它证伪的是**单维 1PL**,不是「一定有两种反对」;
   2PL(题目区分度不同)也能造出同样的偏离,**而三道题定不出 2PL**(自由度不够)
   ⇒ **本轮报的是「单维 1PL 被证伪」,不是「因此是两种人」**,这条边界不许省;
④ **`Δmoral` 只有 homo 靶有** —— 其余四靶没有对应的道德题 ⇒ 那一格跨靶复制**不存在**;
⑤ **换不了仪器**:五个靶都在 GSS 内(`#854` 已点名盘上七具)⇒ 本轮是**跨靶**不是**跨仪器**,
   **它控住格式,控不住「GSS 的受访者」这一层。**
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
NPERM, SEED = 200, 305
P865 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be"
                      "/A95_那几块孤立的少数是不是同一批人/R304_一个人能同时落在几块里"
                      "/results/blocs_per_person.json"))

print("=== ⓪a **定理先行**:Rasch 条件似然 —— 给定总分,哪一题被拒与 θ 无关(机器验) ===")
rg0 = np.random.default_rng(SEED)
bdiff = np.array([-0.6, 0.0, 0.9])
prof = []
for lo, hi, tag in ((-3, -0.5, "低 θ"), (-0.5, 0.5, "中 θ"), (0.5, 3, "高 θ")):
    th = rg0.uniform(lo, hi, 300_000)[:, None]
    x = (rg0.random((300_000, 3)) < 1 / (1 + np.exp(-(th - bdiff)))).astype(int)
    one = x.sum(1) == 1
    d = x[one].mean(0); prof.append(d)
    print(f"  {tag} n(R=1)={one.sum():7,d} · P(哪一题 | R=1) = {np.round(d,4)}")
RASCH_TVD = float(max(0.5 * np.abs(prof[i] - prof[j]).sum() for i in range(3) for j in range(i + 1, 3)))
print(f"  ⇒ 三档 θ 两两 TVD 最大 **{RASCH_TVD:.5f}** ⇒ **定理成立,零是 0** —— "
      f"**偏离它就是单维 1PL 被证伪**(⚠ 而不是「因此有两种人」,2PL 也能造出偏离)")
print(f"  ⚠ 而实测的 TVD 仍要对照**它自己的置换零**,不对照 0 —— **TVD 有正偏**")

print("\n=== ⓪b **块的定义与 `#865` 不同,而这个不同必须说出来** ===")
print(f"  `#865` 问「哪一块是孤立的」⇒ 政治那块是**自由派**(λ_perm {P865['grid'][0]['n']and''}"
      f"见其产物);**本轮问「哪一块是不宽容的」⇒ 政治取保守端。**")
print("  **同一根轴、不同的三分之一 —— 不是矛盾,是两个不同的问题。** "
      "块的定义**不看任何结局**,跨靶比较才不循环。")

TARGETS = {"同性恋 homo": ("spkhomo", "colhomo", "libhomo"),
           "种族主义者 rac": ("spkrac", "colrac", "librac"),
           "共产主义者 com": ("spkcom", "colcom", "libcom"),
           "军管主义者 mil": ("spkmil", "colmil", "libmil"),
           "无神论者 ath": ("spkath", "colath", "libath")}
GC = ["year", "homosex", "attend", "reliten", "fund", "polviews", "age", "ballot"]
for t in TARGETS.values(): GC += list(t)
gs = pd.read_stata(EXT / "gss/GSS_stata/gss7224_r3a.dta", columns=GC, convert_categoricals=False)
D = pd.DataFrame({"year": gs.year})
D["moral"] = pd.to_numeric(gs.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
for tag, (a, b, c) in TARGETS.items():
    k = tag.split()[-1]
    D[f"spk_{k}"] = 2 - pd.to_numeric(gs[a], errors="coerce").where(lambda v: v.isin([1, 2]))
    D[f"col_{k}"] = 5 - pd.to_numeric(gs[b], errors="coerce").where(lambda v: v.isin([4, 5]))
    D[f"lib_{k}"] = pd.to_numeric(gs[c], errors="coerce").where(lambda v: v.isin([1, 2])) - 1
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)),
                    ("polviews", (1, 7)), ("age", (18, 89)), ("ballot", (1, 4))):
    D[c] = pd.to_numeric(gs[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
Rr = D.dropna(subset=["attend", "reliten", "fund"]).copy()
Rr["REL"] = (zs(Rr.attend) + zs(-Rr.reliten) + zs(-Rr.fund)) / 3
D = D.join(Rr["REL"])

DECS = {"1970s": range(1972, 1980), "1980s": range(1980, 1990), "1990s": range(1990, 2000),
        "2000s": range(2000, 2010), "2010s": range(2010, 2020), "2020s": range(2020, 2025)}
BSPEC = {"B≥1": 1, "B≥2": 2}

print("\n=== ⓪c 硬规则①:**变量名不是测量** —— 每一靶的 n 与实际问过的年份 ===")
for tag in TARGETS:
    k = tag.split()[-1]
    m = D[[f"spk_{k}", f"col_{k}", f"lib_{k}"]].notna().all(1)
    yy = sorted(D.year[m].unique().astype(int))
    print(f"  {tag:14s} 三题齐全 n={int(m.sum()):6,d} · {len(yy)} 个年份 {yy[0]}–{yy[-1]}")


def prep(sub, k):
    """返回 X(n×3 拒绝矩阵,1=不放行)· B(落在几块里)· pb(安慰剂假块)· moral。"""
    cols = [f"spk_{k}", f"col_{k}", f"lib_{k}"]
    X = 1 - sub[cols].to_numpy(float)                 # 1 = 不放行
    rel, age, pol = (sub[c].to_numpy(float) for c in ("REL", "age", "polviews"))
    B = np.zeros(len(sub))
    for v in (rel, age, pol):                          # ⚠ 三块都取**不宽容**那一端:
        ok = np.isfinite(v)                            #    虔诚高 · 年龄高 · polviews 高(保守)
        if ok.sum() < 200: return None
        q = np.nanquantile(v[ok], 2 / 3)
        B += (ok & (v >= q)).astype(float)
    bv = sub["ballot"].to_numpy(float)
    # ⚠⚠ **安慰剂的补集里不许有「这个变量根本不存在」的人。**
    # 第一版写的是 `ballot == 1` vs **其余一切**,而 `ballot` 1988 年才有 ⇒ 1980s 那一格的
    # 「对照组」其实是 **1980–87 的受访者**,于是安慰剂变成了一个**时期对比**,不是随机对比。
    # 实测:同一年内 ballot 1 vs 3 的剖面只差 0.05–0.11(n≈80,抽样量级);
    # 而 ballot 1 与 3 **拿到的靶完全一样**(2010s 五靶齐全率都是 0.92–0.97)⇒
    # **不是问卷内容的差,是我把缺失值当成了对照组。**
    # ⇒ 安慰剂只在 `ballot ∈ {1,3}` 上算,`NaN` 与 ballot 2(它一个靶都没被问)一律排除。
    pb = np.isfinite(bv) & (bv == 1)
    pb_ok = np.isfinite(bv) & np.isin(bv, [1.0, 3.0])
    return X, B, pb, sub["moral"].to_numpy(float), pb_ok


def tvd_profile(X, grp):
    """恰好拒一题的人里,两群「拒哪一题」分布的 TVD。"""
    one = X.sum(1) == 1
    a, b = one & grp, one & ~grp
    if a.sum() < 40 or b.sum() < 40: return np.nan, int(a.sum()), int(b.sum())
    return float(0.5 * np.abs(X[a].mean(0) - X[b].mean(0)).sum()), int(a.sum()), int(b.sum())


def synth(X, B, thr, rng):
    """合成零:每题独立抽,**保留每个 B 层的逐题边际** ⇒ 层水平效应原样在,只毁掉人层的一般性。"""
    Y = np.zeros_like(X)
    for lev in np.unique(B):
        m = B == lev
        if m.sum() == 0: continue
        p = X[m].mean(0)
        Y[m] = (rng.random((int(m.sum()), 3)) < p).astype(float)
    return Y


def cell(sub, k, thr, rng, plant=0.0, plant_item=0, perm=False, placebo=False):
    z = prep(sub, k)
    if z is None: return None
    X, B, pb, moral, pb_ok = z
    ok = np.isfinite(X).all(1)
    if placebo: ok = ok & pb_ok              # ⚠ 安慰剂总体 = 只有真的有 ballot 的人
    X, B, pb, moral = X[ok], B[ok], pb[ok], moral[ok]
    if len(X) < 400: return None
    grp = pb if placebo else (B >= thr)
    if perm: grp = grp[rng.permutation(len(grp))]
    if plant > 0:                                    # 正控:把一部分 B=0 单拒者改拒指定题
        one = (X.sum(1) == 1) & ~grp
        idx = np.where(one)[0]
        if len(idx):
            take = idx[rng.random(len(idx)) < plant]
            X = X.copy(); X[take] = 0.0; X[take, plant_item] = 1.0
    tv, na, nb = tvd_profile(X, grp)
    R = X.sum(1); ref = R >= 1
    if ref.sum() < 100: return None
    br = float(R[ref & ~grp].mean() - R[ref & grp].mean()) if (ref & grp).sum() > 30 else np.nan
    mo = np.nan
    if np.isfinite(moral).sum() > 200:
        f = ref & np.isfinite(moral)
        if (f & grp).sum() > 30 and (f & ~grp).sum() > 30:
            mo = float(np.nanmean(moral[f & ~grp]) - np.nanmean(moral[f & grp]))
    return dict(tvd=tv, n_one_blocless=na, n_one_bloc=nb, breadth=br, moral=mo,
                n=int(len(X)), refuse_rate=float(ref.mean()))


def synth_cell(sub, k, thr, rng):
    """同一格在合成独立零下的 breadth / moral —— 同样施加 R≥1 的条件(collider 一并复制)。"""
    z = prep(sub, k)
    if z is None: return np.nan, np.nan
    X, B, pb, moral, _pbok = z
    ok = np.isfinite(X).all(1)
    X, B, moral = X[ok], B[ok], moral[ok]
    if len(X) < 400: return np.nan, np.nan
    grp = B >= thr
    Y = synth(X, B, thr, rng)
    R = Y.sum(1); ref = R >= 1
    if (ref & grp).sum() < 30 or (ref & ~grp).sum() < 30: return np.nan, np.nan
    br = float(R[ref & ~grp].mean() - R[ref & grp].mean())
    mo = np.nan
    f = ref & np.isfinite(moral)
    if (f & grp).sum() > 30 and (f & ~grp).sum() > 30:
        mo = float(np.nanmean(moral[f & ~grp]) - np.nanmean(moral[f & grp]))
    return br, mo


print(f"\n=== ① 网格:{len(TARGETS)} 靶 × {len(DECS)} 十年 × {len(BSPEC)} 切法 · "
      f"置换零 + 合成零各 {NPERM} 次 ===")
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
        for bs, thr in BSPEC.items():
            o = cell(sub, k, thr, rng)
            if o is None or not np.isfinite(o["tvd"]): continue
            nt = np.array([cell(sub, k, thr, rng, perm=True) for _ in range(NPERM)], dtype=object)
            nt = np.array([x["tvd"] for x in nt if x and np.isfinite(x["tvd"])])
            sb = np.array([synth_cell(sub, k, thr, rng) for _ in range(NPERM)], dtype=float)
            pl = cell(sub, k, thr, rng, placebo=True)
            # ⚠⚠ **安慰剂必须有它自己的地板。** `ballot` 是 50/50,而 `B≥1` 最偏到 0.908/0.092 ——
            # **TVD 的抽样地板取决于两组的大小**,拿 B 的地板去判 ballot 是拿一个量的零去判另一个量
            # (`#862` 犯过完全同一形状:用「单个 |g|」的地板去判「最大−次大」)。
            # 第一版正是这么写的,安慰剂 13/23 越界、P≈0 —— **失败的是判据,不是安慰剂。**
            npl = np.array([x["tvd"] for x in
                            (cell(sub, k, thr, rng, placebo=True, perm=True) for _ in range(NPERM))
                            if x and np.isfinite(x["tvd"])])
            pl_fl = float(np.quantile(npl, 0.95)) if len(npl) else np.nan
            rows.append(dict(target=tag, dec=dec, bspec=bs, tvd=o["tvd"],
                             tvd_floor=float(np.quantile(nt, 0.95)) if len(nt) else np.nan,
                             tvd_null_med=float(np.median(nt)) if len(nt) else np.nan,
                             p=float((1 + (nt >= o["tvd"]).sum()) / (len(nt) + 1)) if len(nt) else np.nan,
                             breadth=o["breadth"], breadth_synth=float(np.nanmean(sb[:, 0])),
                             breadth_resid=o["breadth"] - float(np.nanmean(sb[:, 0])),
                             breadth_synth_sd=float(np.nanstd(sb[:, 0])),
                             moral=o["moral"], moral_synth=float(np.nanmean(sb[:, 1])),
                             moral_resid=(o["moral"] - float(np.nanmean(sb[:, 1]))
                                          if np.isfinite(o["moral"]) else np.nan),
                             placebo_tvd=(pl["tvd"] if pl else np.nan), placebo_floor=pl_fl,
                             n=o["n"], n_one_blocless=o["n_one_blocless"], n_one_bloc=o["n_one_bloc"],
                             refuse_rate=o["refuse_rate"]))
        r1 = [r for r in rows if (r["target"], r["dec"], r["bspec"]) == (tag, dec, "B≥1")]
        if r1:
            r = r1[0]
            line.append(f"{dec}:{r['tvd']:.3f}{'*' if r['tvd'] > r['tvd_floor'] else ' '}"
                        f"(地板{r['tvd_floor']:.3f})")
    print(f"  {tag:14s} TVD · " + " · ".join(line))
print("  ⚠ `*` = 超出该格自己的置换地板。**Rasch 下这个量的理论零是 0;实测零是置换分布,因为 TVD 有正偏。**")

MAIN = [r for r in rows if r["bspec"] == "B≥1"]
if not MAIN:
    raise SystemExit("⛔ 网格为空 —— 空总体不许当作通过")
by_t = {}
for tag in TARGETS:
    rs = [r for r in MAIN if r["target"] == tag]
    if not rs: continue
    by_t[tag] = dict(n_cells=len(rs), n_over=sum(1 for r in rs if r["tvd"] > r["tvd_floor"]),
                     tvd_med=float(np.median([r["tvd"] for r in rs])),
                     floor_med=float(np.median([r["tvd_floor"] for r in rs])))
print(f"\n=== ② 跨靶:**这是唯一能把「关于世界」和「关于问卷」分开的东西** ===")
for tag, v in by_t.items():
    print(f"  {tag:14s} TVD 中位 **{v['tvd_med']:.4f}**(地板中位 {v['floor_med']:.4f})· "
          f"**超零 {v['n_over']}/{v['n_cells']} 个十年**")
HOMO = "同性恋 homo"
homo_r = by_t.get(HOMO, dict(n_over=0, n_cells=1))["n_over"] / max(by_t.get(HOMO, dict(n_cells=1))["n_cells"], 1)
oth_r = [v["n_over"] / v["n_cells"] for t, v in by_t.items() if t != HOMO]
print(f"  ⇒ homo 超零比例 **{homo_r:.2f}** · 其余四靶 " + " · ".join(f"{x:.2f}" for x in oth_r)
      + f"(中位 {np.median(oth_r):.2f})")

print(f"\n=== ③ 广度与道德:**减掉合成零之后还剩什么** ===")
for tag in by_t:
    rs = [r for r in MAIN if r["target"] == tag and np.isfinite(r["breadth_resid"])]
    if not rs: continue
    br = [r["breadth_resid"] for r in rs]; bs_ = [r["breadth_synth_sd"] for r in rs]
    print(f"  {tag:14s} Δbreadth 实测中位 {np.median([r['breadth'] for r in rs]):+.4f} · "
          f"合成零中位 {np.median([r['breadth_synth'] for r in rs]):+.4f} ⇒ "
          f"**残差 {np.median(br):+.4f}**(合成零自身 SD 中位 {np.median(bs_):.4f})")
mrs = [r for r in MAIN if r["target"] == HOMO and np.isfinite(r["moral_resid"])]
if mrs:
    print(f"  {HOMO:14s} Δmoral 实测中位 {np.median([r['moral'] for r in mrs]):+.4f} · "
          f"合成零中位 {np.median([r['moral_synth'] for r in mrs]):+.4f} ⇒ "
          f"**残差 {np.median([r['moral_resid'] for r in mrs]):+.4f}** "
          f"(⚠ **它部分被块的定义强制,所以只有残差可读**;**其余四靶没有对应的道德题,"
          f"这一格跨靶复制不存在**)")

print("\n=== ④ 控制 ===")
sub0 = D[D[["spk_homo", "col_homo", "lib_homo"]].notna().all(1) & D.REL.notna()
         & D.year.isin(list(DECS["2010s"]))]
rgc = np.random.default_rng(SEED + 1)
b0 = cell(sub0, "homo", 1, rgc)
z0 = cell(sub0, "homo", 1, rgc, plant=0.0)
# ⚠⚠ **植入的方向必须由数据算,不能由我选。**
# 第一版把 B=0 的单拒者一律推向「图书馆」,而 B≥1 本来就最常只拒图书馆 ⇒ **两群被推到了一起**,
# TVD 从 0.1896 掉到 0.0274(−0.16),正控 FAIL。**它证明了仪器灵敏,却证反了方向。**
# `#863` 刚犯过同一个错(往高组加 +0.30 是把两组拉近),**这是四轮里第三次。**
# ⇒ 推向 **B≥1 最少拒的那一题**,那才是把两个剖面推开。
_z = prep(sub0, "homo"); _X, _B = _z[0], _z[1]
_ok = np.isfinite(_X).all(1); _X, _B = _X[_ok], _B[_ok]
_one = (_X.sum(1) == 1) & (_B >= 1)
PROF_B1 = _X[_one].mean(0)
PITEM = int(np.argmin(PROF_B1))
ITEMNM = ["发言", "教书", "图书馆"]
print(f"  ⚠ **植入方向由数据算**:B≥1 单拒者的剖面 = "
      + " · ".join(f"{ITEMNM[i]} {PROF_B1[i]:.3f}" for i in range(3))
      + f" ⇒ **推向它最少拒的「{ITEMNM[PITEM]}」才是把两群推开**")
DOSE = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00]
curve = [(f, cell(sub0, "homo", 1, np.random.default_rng(SEED + 100 + i),
                  plant=f, plant_item=PITEM)["tvd"]) for i, f in enumerate(DOSE)]
p3 = dict(tvd=[v for f, v in curve if f == 0.30][0])
c1 = dict(tvd=[v for f, v in curve if f == 1.00][0])
CT = 0.03
print(f"  正控(**剂量-反应**,把 B=0 的单拒者按比例改成只拒「{ITEMNM[PITEM]}」):"
      + " · ".join(f"{f:.0%}→{v:.4f}" for f, v in curve))
MONO = all(curve[i][1] <= curve[i + 1][1] + 1e-9 for i in range(len(curve) - 1))
r0_pre = [r for r in MAIN if r["target"] == HOMO and r["dec"] == "2010s"][0]
MDE_P = next((f for f, v in curve if f > 0 and v > r0_pre["tvd_floor"]), None)
print(f"     **单调 {MONO}** · **plant=0 时 {z0['tvd']-b0['tvd']:+.2e}**(⚠ `G2` 控制必须能失败)· "
      f"**MDE = 最小的能顶出该格地板({r0_pre['tvd_floor']:.4f})的植入比例 = "
      f"{f'{MDE_P:.0%}' if MDE_P is not None else '>100%'}**")
CEIL = abs(c1["tvd"] - b0["tvd"])
print(f"     **控制也必须能通过**:floor {abs(z0['tvd']-b0['tvd']):.2e} < 阈 {CT} < "
      f"ceiling(100% 植入){CEIL:.4f} ⇒ "
      f"**{'阈在真带内' if abs(z0['tvd']-b0['tvd']) < CT < CEIL else '⚠⚠ 阈不在带内'}**")
rgn = np.random.default_rng(SEED + 2)
nn = cell(sub0, "homo", 1, rgn, perm=True)
r0 = [r for r in MAIN if r["target"] == HOMO and r["dec"] == "2010s"][0]
print(f"  负控:在单拒者里**打乱 B 标签** ⇒ TVD = **{nn['tvd']:.4f}**,该格地板 **{r0['tvd_floor']:.4f}** ⇒ "
      f"**{'塌回地板内' if nn['tvd'] <= r0['tvd_floor'] else '⚠ 没塌回去'}**")
print(f"     ⚠ **「这个零该不该是零?」不该** —— **TVD 是非负的,有正偏**:"
      f"该格置换零的中位是 **{r0['tvd_null_med']:.4f}** ≠ 0,"
      f"**Rasch 定理给的 0 是理论零,可比的是置换零。**")
plv = [r["placebo_tvd"] for r in MAIN if np.isfinite(r["placebo_tvd"]) and np.isfinite(r["placebo_floor"])]
plfl = [r["placebo_floor"] for r in MAIN if np.isfinite(r["placebo_tvd"]) and np.isfinite(r["placebo_floor"])]
pl_out = sum(1 for a, b in zip(plv, plfl) if a > b)
from scipy.stats import binom as _binom
PL_MAX = int(_binom.ppf(0.95, len(plv), 0.05))
PL_P = float(1 - _binom.cdf(pl_out - 1, len(plv), 0.05)) if pl_out > 0 else 1.0
print(f"  安慰剂 **`ballot` 当假块**,⚠ **对照它自己的地板**(50/50 的几何,地板比 B 的小)⇒ 越界 **{pl_out}/{len(plv)}**,"
      f"`P(X≥{pl_out})={PL_P:.3f}`,二项 95 分位 **{PL_MAX}** ⇒ "
      f"**{'落在零里' if pl_out <= PL_MAX else '⚠ 超出零'}**"
      f"(⚠ **判据是二项零,不是「几乎全过」** —— `#864` 补的那一条)")

ps = np.array([r["p"] for r in rows if np.isfinite(r["p"])]); C = len(ps)
o_ = np.argsort(ps); q = 0.05
cH = q * np.arange(1, C + 1) / C; cY = cH / np.sum(1.0 / np.arange(1, C + 1))
su = lambda pv, cr: (int(np.max(np.where(pv <= cr)[0])) + 1 if (pv <= cr).any() else 0)
kH, kY = su(ps[o_], cH), su(ps[o_], cY)
print(f"\n=== ⑤ 多重性:整族 **{C}** 格(5 靶 × 十年 × 2 切法)· BH 存活 **{kH}** · BY **{kY}** · "
      f"p 分辨率下限 {1/(NPERM+1):.4f} · 不同意的 {kH-kY} 格一起发表 ===")

G = Gate("#866 · 同一条尺上的两个位置,还是两种反对")
G.asserted("① **定理先行**:Rasch 条件似然 —— 给定总分,哪一题被拒与 θ 无关 ⇒ "
           "**理论零是 0**;机器上验过(θ 三档跨 3.5 logit)",
           bool(RASCH_TVD < 0.01), f"三档 θ 两两 TVD 最大 {RASCH_TVD:.5f}", kind="control")
G.asserted("② 前提(跑前写下的最强混淆):**作答定势** —— 同题干连问会制造无内容的「广度」 ⇒ "
           "**同一仪器同一格式换靶**,五个靶一起跑;**是格式的话五靶同,是对象的话 homo 分开**",
           bool(len(by_t) >= 4),
           " · ".join(f"{t.split()[0]}:{v['n_over']}/{v['n_cells']}" for t, v in by_t.items()),
           kind="control")
G.asserted("③ **广度差被强制,只有残差可读**:块预测「反不反对」+ 在 R≥1 上做条件(collider) ⇒ "
           "**合成零同样施加 R≥1**,报残差",
           bool(all(np.isfinite(r["breadth_synth"]) for r in MAIN if np.isfinite(r["breadth"]))),
           f"残差中位 {np.median([r['breadth_resid'] for r in MAIN if np.isfinite(r['breadth_resid'])]):+.4f} · "
           f"合成零 SD 中位 {np.median([r['breadth_synth_sd'] for r in MAIN if np.isfinite(r['breadth_synth_sd'])]):.4f}",
           kind="control")
G.asserted("④ 正控(**剂量-反应**):把 B=0 的单拒者按比例改拒 **B≥1 最少拒的那一题**(方向由数据算,"
           "不由我选)⇒ TVD 必须**单调上升**;plant=0 恰为 0;**且阈落在 floor 与 ceiling 之间**",
           bool(MONO and p3["tvd"] - b0["tvd"] > CT and abs(z0["tvd"] - b0["tvd"]) < 1e-12
                and abs(z0["tvd"] - b0["tvd"]) < CT < CEIL),
           f"{b0['tvd']:.4f} → {p3['tvd']:.4f}({p3['tvd']-b0['tvd']:+.4f})· "
           f"plant=0 {z0['tvd']-b0['tvd']:+.2e} · 带 [0, {CEIL:.4f}] · 单调 {MONO} · "
           f"剂量曲线 " + " ".join(f"{f:.0%}:{v:.4f}" for f, v in curve), kind="control")
G.asserted("⑤ 负控:在单拒者里打乱 B 标签 ⇒ TVD 塌回**该格自己的置换地板**"
           "(⚠ **不是塌回 0** —— TVD 非负有正偏)",
           bool(nn["tvd"] <= r0["tvd_floor"]),
           f"{nn['tvd']:.4f} ≤ {r0['tvd_floor']:.4f};该格置换零中位 {r0['tvd_null_med']:.4f} ≠ 0",
           kind="control")
G.asserted("⑥ 安慰剂 `ballot` 当假块 ⇒ 越界格数落在二项零 Bin(N, 0.05) 内",
           bool(pl_out <= PL_MAX),
           f"越界 {pl_out}/{len(plv)} · P(X≥{pl_out})={PL_P:.3f} · 二项 95 分位 {PL_MAX}", kind="control")
G.asserted("⑦ kill(预注册):「零块反对者只是同一条尺上更轻的位置」要成立,需 "
           "**homo 靶上 TVD 超零的十年 ≤ 1/3**",
           bool(homo_r <= 1 / 3),
           f"homo 超零 {by_t.get(HOMO,{}).get('n_over','?')}/{by_t.get(HOMO,{}).get('n_cells','?')} "
           f"= {homo_r:.2f} · 其余四靶中位 {np.median(oth_r):.2f}", kind="kill",
           yardstick="`TVD_profile` 相对**它自己的年内置换零**(在单拒者里打乱 B 标签);"
                     "Rasch 条件似然给的理论零是 0,但可比的是置换零(TVD 非负有正偏)",
           yardstick_noise=float(np.median([r["tvd_floor"] for r in MAIN])),
           population=f"GSS 五个靶 × 各自的十年,共 {len(MAIN)} 个 B≥1 格;"
                      f"**kill 只判 homo 靶那 {by_t.get(HOMO,{}).get('n_cells','?')} 格**,"
                      f"其余四靶是**作答定势的控制**,不是被判的对象(`#854`①:总体里不许有在位者)",
           direction=[r["tvd"] - r["tvd_floor"] for r in MAIN if r["target"] == HOMO])
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
ITEMNAME = ["发言", "教书", "图书馆"]
if adm and by_t.get(HOMO):
    hs = [r for r in MAIN if r["target"] == HOMO]
    pr_b0, pr_b1 = [], []
    for dec in DECS:
        m = D[["spk_homo", "col_homo", "lib_homo"]].notna().all(1) & D.REL.notna() \
            & D.year.isin(list(DECS[dec]))
        sub = D[m]
        if len(sub) < 800: continue
        z = prep(sub, "homo")
        if z is None: continue
        X, B = z[0], z[1]
        ok = np.isfinite(X).all(1); X, B = X[ok], B[ok]
        one = X.sum(1) == 1; g1 = B >= 1
        if (one & ~g1).sum() > 40: pr_b0.append(X[one & ~g1].mean(0))
        if (one & g1).sum() > 40: pr_b1.append(X[one & g1].mean(0))
    PB0 = np.mean(pr_b0, 0) if pr_b0 else np.full(3, np.nan)
    PB1 = np.mean(pr_b1, 0) if pr_b1 else np.full(3, np.nan)
else:
    PB0 = PB1 = np.full(3, np.nan)

if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif np.median(oth_r) >= 2 / 3 and homo_r >= 2 / 3:
    VERD = (f"**D 五个靶给出同样的偏离 ⇒ 那不是关于这个对象的发现,是问卷格式的性质(作答定势)。**\n"
            f"  homo 超零 {homo_r:.2f},其余四靶中位 {np.median(oth_r):.2f} —— "
            f"**「块 / 零块」这个分解没有碰到内容。**")
elif homo_r >= 2 / 3 and np.median(oth_r) < 1 / 3:
    VERD = (f"**B 两种反对 —— 而且只在这个对象上。**\n"
            f"  homo 靶 TVD 超零 {homo_r:.2f} 的十年,其余四靶中位只有 {np.median(oth_r):.2f} ⇒ "
            f"**不是作答定势。**\n"
            f"  恰好拒一题的人里,拒的是哪一题:**零块 " +
            " / ".join(f"{ITEMNAME[i]} {PB0[i]:.2f}" for i in range(3)) + "** vs **块内 " +
            " / ".join(f"{ITEMNAME[i]} {PB1[i]:.2f}" for i in range(3)) + "**。\n"
            f"  ⇒ **Rasch 单维性被证伪** —— ⚠ **而这只说明 1PL 不成立,不说明「因此是两种人」**:"
            f"2PL(题目区分度不同)也能造出同样的偏离,**而三道题定不出 2PL**。")
elif homo_r <= 1 / 3:
    VERD = (f"**A 一条尺上的两个位置 —— 「块」只是标签。**\n"
            f"  homo 靶 TVD 超零只有 {homo_r:.2f} 的十年,广度残差中位 "
            f"{np.median([r['breadth_resid'] for r in MAIN if r['target']==HOMO and np.isfinite(r['breadth_resid'])]):+.4f}。\n"
            f"  恰好拒一题的人里,拒的是哪一题:**零块 " +
            " / ".join(f"{ITEMNAME[i]} {PB0[i]:.2f}" for i in range(3)) + "** vs **块内 " +
            " / ".join(f"{ITEMNAME[i]} {PB1[i]:.2f}" for i in range(3)) + "** —— **几乎同一个分布。**\n"
            f"  ⇒ **一句关于人的话:那些不在任何一块里的反对者,不是第二种反对者。\n"
            f"  他们和最虔诚、最年长、最保守的那批人反对的是同一批东西,按同样的先后顺序 ——\n"
            f"  只是没走那么远。反对不是几个阵营,是一条连续的尺,而块只是尺上人多的地方。**")
else:
    VERD = (f"**都不是**:homo 超零 {homo_r:.2f},其余四靶中位 {np.median(oth_r):.2f} —— "
            f"**四个预注册世界都没被满足,如实登记。**")
print(VERD)
print(f"\n⚠ **本轮结构性做不到的**:① 横断面 ⇒ **无因果识别**;② **「同一种人」只能到「同一年同一批"
      f"受访者的作答模式」**,GSS 无面板 ⇒ **结构性拿不到**;③ **Rasch 只是一个模型** —— "
      f"偏离它证伪的是**单维 1PL**,不是「一定有两种反对」,而**三道题定不出 2PL**;"
      f"④ **`Δmoral` 只有 homo 靶有**,其余四靶没有对应的道德题 ⇒ **那一格跨靶复制不存在**;"
      f"⑤ **换不了仪器** —— 五个靶都在 GSS 内,本轮是**跨靶**不是**跨仪器**,"
      f"**它控住格式,控不住「GSS 的受访者」这一层。**")

json.dump(dict(grid=rows, by_target=by_t, homo_over_rate=homo_r,
               other_over_rates=oth_r, rasch_theory_tvd=RASCH_TVD,
               profile_blocless=PB0.tolist(), profile_bloc=PB1.tolist(),
               multiplicity=dict(cells=C, bh=int(kH), by=int(kY), q=q),
               controls=dict(pos_from=b0["tvd"], pos_to=p3["tvd"], ceiling=CEIL, threshold=CT,
                             dose_curve=[[f, v] for f, v in curve], monotone=MONO,
                             plant_item=PITEM, profile_b1=PROF_B1.tolist(), mde_plant=MDE_P,
                             placebo_floor_note="placebo compared to ITS OWN floor; ballot is 50/50 "
                                                "while B>=1 reaches 0.908/0.092",
                             zero=z0["tvd"] - b0["tvd"], neg=nn["tvd"], neg_floor=r0["tvd_floor"],
                             null_median_not_zero=r0["tvd_null_med"],
                             placebo_out=pl_out, placebo_cells=len(plv), placebo_p=PL_P),
               derivation="Rasch conditional-likelihood: given the total score, WHICH item is refused "
                          "is independent of theta; verified numerically across 3.5 logits",
               bloc_def_differs_from_865="#865 asked which third is ISOLATED (politics = liberals); "
                                         "this round asks which third is RESTRICTIVE (politics = "
                                         "conservatives). Same axis, different third, different question.",
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), seed=SEED, nperm=NPERM),
          open(OUT / "two_kinds_of_refusal.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'two_kinds_of_refusal.json'}")

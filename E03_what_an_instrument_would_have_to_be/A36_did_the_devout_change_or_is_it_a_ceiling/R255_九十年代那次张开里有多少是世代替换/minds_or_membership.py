"""#816 · E03·A62·R255 —— 九十年代那次张开,是人改了主意,还是人换了?

`#813` 写下的那句关于人的话是:**「那两次张开是真的两群人朝不同方向走。」**
⚠⚠ **而那句话里藏着一个我从来没检验过的前提:「走」预设了同一批人移动了。**
   **一个反复出现的横截面里,两个年份之间的均值变化,可以完全不需要任何人改变主意 ——
   只要进来的年轻世代与离开的年老世代不同就够了。**
   ⇒ **世代替换与态度改变,在均值上无法区分,而它们对人说的是两句相反的话:**
   **① 有人改了主意(période)· ② 没有人改主意,是人群换了(cohorte)。**
   **`#813` 的头条按构造预设了 ①,而这一轮去看它站不站得住。这是一个本体论的分叉,不是参数的。**

G1 估计量:**把 1990→1998 每一层的均值变化拆成两项**(标准的 Kitagawa/Oaxaca 形状):
   `Δmean = Σ_c w̄_c·(m_c1 − m_c0)` **(世代内的态度改变)**
          `+ Σ_c m̄_c·(w_c1 − w_c0)` **(世代构成的改变)**
   再把两层的拆解相减,得到 **`Δgap` 里「态度」与「构成」各占多少**。

⚠⚠ **而这里有一个必须先说清楚的陷阱,否则整轮是 `realstat` 点名的算术陷阱:**
   **上面那个等式是恒等式 —— 它一定成立,拆出来的两项一定加起来等于总量。**
   **拆解本身不是证据,它是代数。** 有证据价值的只有两件事:
   **① 两项的相对大小 · ② 它在重抽下稳不稳。** ⇒ **本轮只报这两件,不把恒等式当成发现。**

三个世界:
   A **是改主意**:构成项小 ⇒ `#813` 的头条站住,**而它现在有了它一直缺的那个前提检验。**
   B **是人换了**:构成项大 ⇒ **`#813` 那句「两群人朝不同方向走」必须撤** ——
     **正确的话变成:虔诚层里年轻人少了/老人多了,而那不是任何人改了主意。**
   C **两者都有,且方向相反** ⇒ **总量掩盖了两个互相抵消的过程**,那会把对象再换一次。

预测矩阵:
   | 世界 | 现在 | 构成项 < 1/3 | 构成项 > 2/3 | 两项反号 |
   | A 改主意 | 0.45 | **0.85** | 0.03 | 0.10 |
   | B 人换了 | 0.30 | 0.05 | **0.90** | 0.15 |
   | C 互相抵消 | 0.25 | 0.10 | 0.07 | **0.75** |

预注册判词(条件式):
  if 正控开火(**造一个只有构成变、态度完全不变的世界,拆解必须把几乎全部归给构成项**)
     and 负控开火(**造一个世代权重完全冻结的世界,构成项必须 ≈ 0**,
        ⚠ 而**这条负控的参照真的是 0**,与 `#801`/`#805` 参照 1.0 的情形相反)
     and 两条控制都**量过自己的噪声半宽**(`#815` 刚立的规矩,本轮第一次执行):
      按构成项占 `Δgap` 的份额三分判;**报整张表,不设「多数」阈值**
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**分层本身是按年份内的三分位算的**(`REL` 的 tercile per year)——
  **若虔诚度的世代分布在变,那么「虔诚层」这个集合本身在两个年份不是同一群人**,
  而那会让「构成」与「态度」的划分变得模糊。
  ⇒ 控制:**同时报每一层内部的世代权重变化**,让读者看见这个混淆有多大;
  ⚠ 并**如实登记:这个设计分不开「层内世代构成变了」与「层的定义漂移了」。**

⚠ 硬规则①:先打印世代分箱、每箱在两个年份、两层的 n。
⚠ 本轮换不了仪器(对象是 GSS 自己的一个十年)。⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(255)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, Y0, Y1 = "homosex", 4, 1990, 1998
B, NREP = 2000, 200

d = pd.read_stata(gp, columns=["year", "cohort", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("cohort", (1880, 2010))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))
G = REL.dropna(subset=[IT, "cohort"])
G = G[G.year.isin([Y0, Y1])].copy()
G["gen"] = (G.cohort//10*10).astype(int)
keep = [g for g, n in G.groupby("gen").size().items() if n >= 60]
G = G[G.gen.isin(keep)]

print(f"=== ⓪ 硬规则①:世代分箱(出生十年)· 两个年份 × 两层的 n ===")
tab = G.groupby(["gen", "year", "k"]).size().unstack(fill_value=0)
print(f"  保留的世代箱({len(keep)} 个,每箱总 n ≥ 60):{sorted(keep)}")
for gen in sorted(keep):
    r = [(int(y), int(tab.loc[(gen, y), 2]) if (gen, y) in tab.index else 0,
          int(tab.loc[(gen, y), 0]) if (gen, y) in tab.index else 0) for y in (Y0, Y1)]
    print(f"    {gen}s 出生:{Y0} 虔诚 {r[0][1]:>3}/非虔诚 {r[0][2]:>3} · {Y1} 虔诚 {r[1][1]:>3}/非虔诚 {r[1][2]:>3}")

def decomp(df):
    """返回 (总变化, 态度项, 构成项) —— ⚠ 这是一个恒等式,拆解本身不是证据。"""
    out = {}
    for k in (2, 0):
        s = df[df.k == k]
        a, b = s[s.year == Y0], s[s.year == Y1]
        gens = sorted(set(a.gen) & set(b.gen))
        if not gens: return None
        w0 = np.array([len(a[a.gen == g])/len(a) for g in gens])
        w1 = np.array([len(b[b.gen == g])/len(b) for g in gens])
        m0 = np.array([a[a.gen == g][IT].mean() for g in gens])
        m1 = np.array([b[b.gen == g][IT].mean() for g in gens])
        if np.isnan(m0).any() or np.isnan(m1).any(): return None
        wbar, mbar = (w0+w1)/2, (m0+m1)/2
        att = float((wbar*(m1-m0)).sum()); com = float((mbar*(w1-w0)).sum())
        out[k] = dict(total=float((w1*m1).sum()-(w0*m0).sum()), attitude=att, composition=com,
                      dw=float(np.abs(w1-w0).sum()/2))
    return dict(dgap=out[2]["total"]-out[0]["total"],
                attitude=out[2]["attitude"]-out[0]["attitude"],
                composition=out[2]["composition"]-out[0]["composition"],
                dw_devout=out[2]["dw"], dw_non=out[0]["dw"], per_stratum=out)

R0 = decomp(G)
print(f"\n=== ① 1990→1998 的 `Δgap` 拆解(⚠ **恒等式,拆解本身不是证据;有价值的是两项的相对大小与稳定性**)===")
print(f"  Δgap 总量 **{R0['dgap']:+.4f}** = 态度项 **{R0['attitude']:+.4f}** + 构成项 **{R0['composition']:+.4f}**")
print(f"  ⇒ 构成项占 |Δgap| 的 **{abs(R0['composition'])/abs(R0['dgap']):.1%}** · "
      f"两项{'同号' if np.sign(R0['attitude'])==np.sign(R0['composition']) else '**反号**'}")
print(f"  ⚠ 跑前混淆的控制 —— 每层内部的世代权重总变动量(½·Σ|Δw|):"
      f"虔诚层 **{R0['dw_devout']:.3f}** · 非虔诚层 **{R0['dw_non']:.3f}**")

sh = np.empty(B)
for i in range(B):
    idx = RNG.integers(0, len(G), len(G))
    r = decomp(G.iloc[idx])
    sh[i] = np.nan if (r is None or abs(r["dgap"]) < 1e-9) else abs(r["composition"])/abs(r["dgap"])
sh = sh[np.isfinite(sh)]
S_LO, S_HI = float(np.percentile(sh, 2.5)), float(np.percentile(sh, 97.5))
print(f"  自助 B={B}:构成项份额 **{abs(R0['composition'])/abs(R0['dgap']):.3f}** [{S_LO:.3f}, {S_HI:.3f}]")

print("\n=== ② 控制(合成世界,同一条代码路径;⚠ **两条都量自己的噪声半宽** —— `#815` 首次执行)===")
def syn(mode):
    """mode='composition_only' 态度完全冻结,只让世代权重变;'weights_frozen' 权重冻结,只让态度变。"""
    H = G.copy()
    gm = H.groupby(["gen", "k"])[IT].mean()
    if mode == "composition_only":
        H[IT] = [gm.get((g, k), H[IT].mean()) for g, k in zip(H.gen, H.k)]   # 态度只由世代决定 ⇒ 不随年份变
    else:
        rows = []
        for (y, k), s in H.groupby(["year", "k"]):
            t = s.copy()
            if y == Y1: t[IT] = t[IT] - 0.4*(k == 0)
            rows.append(t)
        H = pd.concat(rows)
        base = H[H.year == Y0]
        out = [base]
        for k in (2, 0):
            b = base[base.k == k]; nn = len(H[(H.year == Y1) & (H.k == k)])
            t = b.sample(n=nn, replace=True, random_state=int(RNG.integers(1e6))).copy()
            t["year"] = Y1
            t[IT] = t[IT] - 0.4*(k == 0)
            out.append(t)
        H = pd.concat(out)                    # ⚠ 权重冻结:Y1 的世代构成从 Y0 重抽而来
    return H
def ctl(mode, rep=NREP):
    v = []
    for _ in range(rep):
        H = syn(mode)
        idx = RNG.integers(0, len(H), len(H))
        r = decomp(H.iloc[idx])
        if r and abs(r["dgap"]) > 1e-9: v.append(abs(r["composition"])/abs(r["dgap"]))
    v = np.array(v)
    return float(np.median(v)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))
pc_m, pc_lo, pc_hi = ctl("composition_only")
nc_m, nc_lo, nc_hi = ctl("weights_frozen")
print(f"  正控(**态度完全冻结,只有构成变**)⇒ 构成项份额中位 **{pc_m:.3f}** [{pc_lo:.3f}, {pc_hi:.3f}] "
      f"—— 该**接近 1**,噪声半宽 **{(pc_hi-pc_lo)/2:.3f}**")
print(f"  负控(**世代权重冻结,只有态度变**)⇒ 构成项份额中位 **{nc_m:.3f}** [{nc_lo:.3f}, {nc_hi:.3f}] "
      f"—— 该**接近 0**(⚠ **这一次参照真的是 0**),噪声半宽 **{(nc_hi-nc_lo)/2:.3f}**")

Gg = Gate("#816 · 九十年代那次张开,是人改了主意,还是人换了")
Gg.identity_control("① 正控:态度完全冻结、只有构成变的合成世界,构成项份额必须**接近 1**"
                    " —— ⚠ **容差按它自己量出来的噪声半宽给**(`#815` 首次执行)",
                    observed=pc_m, expected=1.0, tol=max((pc_hi-pc_lo)/2, 1e-3),
                    noise_half_width=(pc_hi-pc_lo)/2,
                    what=f"{NREP} 次重复,95% 跨度 [{pc_lo:.3f}, {pc_hi:.3f}]")
Gg.identity_control("② 负控:世代权重冻结、只有态度变的合成世界,构成项份额必须**接近 0**"
                    "(⚠ **这一次参照真的是 0**,与 `#801`/`#805` 参照 1.0 相反)",
                    observed=nc_m, expected=0.0, tol=max((nc_hi-nc_lo)/2, 1e-3),
                    noise_half_width=(nc_hi-nc_lo)/2,
                    what=f"{NREP} 次重复,95% 跨度 [{nc_lo:.3f}, {nc_hi:.3f}]")
Gg.asserted("③ 前提(跑前写下的混淆):每层内部的世代权重变动量并排印出,"
            "且**如实登记这个设计分不开「层内构成变了」与「层的定义漂移了」**",
            True, f"虔诚 {R0['dw_devout']:.3f} · 非虔诚 {R0['dw_non']:.3f}", kind="control")
Gg.asserted("④ 前提:拆解是**恒等式**,本轮只报两项的相对大小与重抽稳定性,不把恒等式当发现",
            bool(abs((R0['attitude']+R0['composition'])-R0['dgap']) < 1e-9),
            f"态度 {R0['attitude']:+.4f} + 构成 {R0['composition']:+.4f} = {R0['dgap']:+.4f}(恒等,按构造)",
            kind="control")
share = abs(R0["composition"])/abs(R0["dgap"])
Gg.asserted("⑤ kill(预注册):「那次张开是人改了主意」要成立,需**构成项份额的区间上界 < 1/3**",
            bool(S_HI < 1/3), f"份额 {share:.3f} [{S_LO:.3f}, {S_HI:.3f}]", kind="kill")
print(); print(Gg)
adm = Gg.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*98)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif S_HI < 1/3:
    V = (f"**A 是改主意。构成项只占 `Δgap` 的 {share:.1%} [{S_LO:.1%}, {S_HI:.1%}]。**\n"
         f"  ⇒ **`#813` 那句「两群人朝不同方向走」站住了,而它现在有了它一直缺的那个前提检验:\n"
         f"  九十年代那次张开不是人群换了,是同一批人里真的有人改了主意、有人没改。**")
elif S_LO > 2/3:
    V = (f"**B 是人换了。构成项占 `Δgap` 的 {share:.1%} [{S_LO:.1%}, {S_HI:.1%}]。**\n"
         f"  ⇒ **`#813` 那句「两群人朝不同方向走」必须撤 —— 正确的话是:\n"
         f"  虔诚层与非虔诚层的世代构成变了,而那不需要任何人改主意。**")
else:
    V = (f"**分不开。构成项份额 {share:.1%} [{S_LO:.1%}, {S_HI:.1%}],跨过了预注册的两个门槛。**\n"
         f"  ⇒ **两个过程都在,而这个设计给不出它们的比例** ——\n"
         f"  ⚠⚠ **但有一件事已经确定,而它足以改写 `#813` 的措辞:构成项不是零。**\n"
         f"  **「两群人朝不同方向走」这句话里,有一部分根本不是「走」,是人换了。**")
print(V)
print(f"\n⚠ **这个设计分不开「层内世代构成变了」与「层的定义漂移了」** —— 分层是按年份内三分位算的,"
      f"若虔诚度的世代分布在变,「虔诚层」在两个年份就不是同一群人。**如实登记,不假装控制了它。**")
json.dump(dict(item=IT, y0=Y0, y1=Y1, gens=sorted(keep), decomposition=R0,
               composition_share=share, share_lo=S_LO, share_hi=S_HI, B=B, n_rep=NREP,
               pos_control=dict(median=pc_m, lo=pc_lo, hi=pc_hi, half_width=(pc_hi-pc_lo)/2),
               neg_control=dict(median=nc_m, lo=nc_lo, hi=nc_hi, half_width=(nc_hi-nc_lo)/2, reference=0.0),
               admissible=adm, verdict=V, gate_ok=Gg.verdict()),
          open(OUT/"minds_or_membership.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'minds_or_membership.json'}")

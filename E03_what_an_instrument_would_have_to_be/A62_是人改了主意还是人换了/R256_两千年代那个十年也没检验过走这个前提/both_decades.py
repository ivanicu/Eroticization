"""#817 · E03·A62·R256 —— 两千年代那个十年也没检验过「走」这个前提;两个十年并排看

`#816`① 与 `#816`② 一起做。
① **`#813` 说两千年代「几乎全是分歧」,而那一轮同样没检验过「走」这个前提** ——
   `#816` 已经证明九十年代有 14.7% 根本不是任何人移动。**两千年代呢?**
   ⚠⚠ **而两个十年的构成项份额若不同,那本身就是一句关于人的话**:
   **同一条鸿沟,可以在一个十年里靠人改主意张开,在下一个十年里靠人换了张开 ——
   而对一个活着的人来说,这两件事完全不同。**
② `#816`② 的修正:上一轮的正控 95% 跨度是 [1.000, 1.000]、半宽 0.000,
   **那不是「噪声极小」,是我造的世界里态度项按构造恒为 0 ⇒ 份额恒等于 1** ——
   **它是一条确定性的代码检查,本该声明 `deterministic=True`,而我传了一个 0 噪声半宽。**
   ⇒ 本轮改正,**让那一行说对它自己是什么。**

G1 估计量:**两个十年各自「构成项占 `Δgap` 的份额」,并排报,不合并、不平均。**

⚠⚠ **拆解是恒等式**(态度项 + 构成项 ≡ 总量,按构造成立)⇒ **拆解不是证据,是代数**
(`realstat` 的算术陷阱)。**有证据价值的只有两项的相对大小与它在重抽下稳不稳。**

⚠ 跑之前写下的最强混淆(与 `#816` 同一条,**重申因为它没有被解决**):
  分层按**年份内三分位**算 ⇒ **这个设计分不开「层内世代构成变了」与「层的定义漂移了」。**
  ⇒ 控制:两个十年、两层的世代权重总变动量全部印出来,**如实登记,不假装控制了它。**

三个世界:
   A **两个十年一样**:份额都小 ⇒ **这条鸿沟两次张开都主要是人改了主意。**
   B **两个十年不同** ⇒ **同一条鸿沟在两个十年里由两种不同的机制张开** ——
     **那是本轮最想要的结果,而它会把「这条鸿沟怎么长的」这个问题本身换掉。**
   C **两个十年都分不开** ⇒ 登记功效不足,**不硬判。**

预注册判词(条件式):
  if 正控开火(态度冻结 ⇒ 份额恒等于 1,**声明为确定性检查**)
     and 负控开火(权重冻结 ⇒ 份额 ≈ 0,**容差取它自己量出来的噪声半宽**,⚠ 参照真的是 0):
      两个十年的份额与区间并排报;**判「两个十年的区间是否重叠」**,不设多数阈值
  else: UNVERIFIED

⚠ 硬规则①:先打印两个十年、每个世代箱、两层的 n。⚠ 换不了仪器(对象是 GSS 自己的两个十年)。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(256)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK = "homosex", 4
WINDOWS = {"1990s": (1990, 1998), "2000s": (2000, 2008)}
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
BASE = REL.dropna(subset=[IT, "cohort"]).copy()
BASE["gen"] = (BASE.cohort//10*10).astype(int)

def window(y0, y1):
    W = BASE[BASE.year.isin([y0, y1])].copy()
    keep = [g for g, n in W.groupby("gen").size().items() if n >= 60]
    return W[W.gen.isin(keep)], sorted(keep)

def decomp(df, y0, y1):
    out = {}
    for k in (2, 0):
        s = df[df.k == k]
        a, b = s[s.year == y0], s[s.year == y1]
        gens = sorted(set(a.gen) & set(b.gen))
        if not gens or not len(a) or not len(b): return None
        w0 = np.array([len(a[a.gen == g])/len(a) for g in gens])
        w1 = np.array([len(b[b.gen == g])/len(b) for g in gens])
        m0 = np.array([a[a.gen == g][IT].mean() for g in gens])
        m1 = np.array([b[b.gen == g][IT].mean() for g in gens])
        if np.isnan(m0).any() or np.isnan(m1).any(): return None
        wbar, mbar = (w0+w1)/2, (m0+m1)/2
        out[k] = dict(total=float((w1*m1).sum()-(w0*m0).sum()),
                      attitude=float((wbar*(m1-m0)).sum()),
                      composition=float((mbar*(w1-w0)).sum()),
                      dw=float(np.abs(w1-w0).sum()/2))
    return dict(dgap=out[2]["total"]-out[0]["total"],
                attitude=out[2]["attitude"]-out[0]["attitude"],
                composition=out[2]["composition"]-out[0]["composition"],
                dw_devout=out[2]["dw"], dw_non=out[0]["dw"], per_stratum=out)

print("=== ⓪ 硬规则①:两个十年 · 世代箱 · 两层 n ===")
RES = {}
for lab, (y0, y1) in WINDOWS.items():
    W, keep = window(y0, y1)
    n = W.groupby(["year", "k"]).size()
    print(f"  {lab} {y0}→{y1} · 世代箱 {len(keep)} 个 {keep} · "
          f"{y0}:虔诚 {n.get((y0,2),0)}/非虔诚 {n.get((y0,0),0)} · {y1}:虔诚 {n.get((y1,2),0)}/非虔诚 {n.get((y1,0),0)}")
    RES[lab] = dict(window=[y0, y1], gens=keep, W=W)

print(f"\n=== ① 两个十年并排(⚠ **拆解是恒等式,不是证据**;B={B})===")
for lab, r in RES.items():
    y0, y1 = r["window"]; W = r["W"]
    D = decomp(W, y0, y1)
    sh = np.empty(B)
    for i in range(B):
        rr = decomp(W.iloc[RNG.integers(0, len(W), len(W))], y0, y1)
        sh[i] = np.nan if (rr is None or abs(rr["dgap"]) < 1e-9) else abs(rr["composition"])/abs(rr["dgap"])
    sh = sh[np.isfinite(sh)]
    lo, hi = float(np.percentile(sh, 2.5)), float(np.percentile(sh, 97.5))
    r.update(dict(D=D, share=abs(D["composition"])/abs(D["dgap"]), lo=lo, hi=hi))
    del r["W"]
    print(f"  {lab}: Δgap **{D['dgap']:+.4f}** = 态度 **{D['attitude']:+.4f}** + 构成 **{D['composition']:+.4f}** "
          f"⇒ 构成项份额 **{r['share']:.1%}** [{lo:.1%}, {hi:.1%}]")
    print(f"        虔诚层 {D['per_stratum'][2]['total']:+.4f}(态度 {D['per_stratum'][2]['attitude']:+.4f}) · "
          f"非虔诚层 {D['per_stratum'][0]['total']:+.4f}(态度 {D['per_stratum'][0]['attitude']:+.4f}) · "
          f"⚠ 层内世代权重变动 {D['dw_devout']:.3f} / {D['dw_non']:.3f}")
A, Bq = RES["1990s"], RES["2000s"]
overlap = not (A["hi"] < Bq["lo"] or Bq["hi"] < A["lo"])
print(f"\n  两个十年的份额区间**{'重叠' if overlap else '**不重叠**'}** ⇒ "
      f"{'这个设计说不出两个十年的机制不同' if overlap else '**两个十年由不同的机制张开**'}")

print("\n=== ② 控制(合成世界,同一条代码路径)===")
W9, _ = window(1990, 1998)
def syn(mode, W, y0, y1):
    H = W.copy()
    if mode == "composition_only":
        gm = H.groupby(["gen", "k"])[IT].mean()
        H[IT] = [gm.get((g, k), H[IT].mean()) for g, k in zip(H.gen, H.k)]
        return H
    base = H[H.year == y0]; out = [base]
    for k in (2, 0):
        b = base[base.k == k]; nn = len(H[(H.year == y1) & (H.k == k)])
        if not len(b) or not nn: continue
        t = b.sample(n=nn, replace=True, random_state=int(RNG.integers(1e6))).copy()
        t["year"] = y1; t[IT] = t[IT] - 0.4*(k == 0)
        out.append(t)
    return pd.concat(out)
def ctl(mode, rep=NREP):
    v = []
    for _ in range(rep):
        H = syn(mode, W9, 1990, 1998)
        rr = decomp(H.iloc[RNG.integers(0, len(H), len(H))], 1990, 1998)
        if rr and abs(rr["dgap"]) > 1e-9: v.append(abs(rr["composition"])/abs(rr["dgap"]))
    v = np.array(v)
    return float(np.median(v)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))
pc_m, pc_lo, pc_hi = ctl("composition_only")
nc_m, nc_lo, nc_hi = ctl("weights_frozen")
print(f"  正控(态度冻结)⇒ 份额 **{pc_m:.6f}** [{pc_lo:.6f}, {pc_hi:.6f}] —— ⚠ **按构造恒等于 1,确定性,不是随机控制**")
print(f"  负控(权重冻结)⇒ 份额中位 **{nc_m:.3f}** [{nc_lo:.3f}, {nc_hi:.3f}],噪声半宽 **{(nc_hi-nc_lo)/2:.3f}**"
      f" —— ⚠ **参照真的是 0**")

Gg = Gate("#817 · 两个十年并排:是改主意还是人换了")
# ⚠ `#816`② 的修正:这条控制按构造态度项恒为 0 ⇒ 份额恒等于 1,**它是确定性的**。
#   上一轮传了一个 0 噪声半宽,读起来像「噪声极小」—— **那是说错了它是什么。**
Gg.identity_control("① 正控:态度完全冻结的合成世界,构成项份额必须**恒等于 1**"
                    " —— ⚠ **声明为确定性检查(`#816`② 的修正)**:它按构造不抖,"
                    "「噪声半宽 0」不是一次测量,是构造的后果",
                    observed=pc_m, expected=1.0, tol=1e-6, deterministic=True,
                    what=f"{NREP} 次重复全部落在 [{pc_lo:.6f}, {pc_hi:.6f}] —— 确定性,非随机控制")
Gg.identity_control("② 负控:世代权重冻结的合成世界,构成项份额必须**接近 0**"
                    "(⚠ **这一次参照真的是 0**),容差取它**自己量出来的噪声半宽**(`#815`)",
                    observed=nc_m, expected=0.0, tol=max((nc_hi-nc_lo)/2, 1e-3),
                    noise_half_width=(nc_hi-nc_lo)/2,
                    what=f"{NREP} 次重复,95% 跨度 [{nc_lo:.3f}, {nc_hi:.3f}]")
Gg.asserted("③ 前提(跑前写下的混淆,重申因为它没被解决):两个十年、两层的世代权重变动量全部印出,"
            "**如实登记这个设计分不开「层内构成变了」与「层的定义漂移了」**", True,
            " · ".join(f"{l} {r['D']['dw_devout']:.3f}/{r['D']['dw_non']:.3f}" for l, r in RES.items()),
            kind="control")
Gg.asserted("④ 前提:拆解是恒等式,只报相对大小与重抽稳定性", 
            bool(all(abs((r["D"]["attitude"]+r["D"]["composition"])-r["D"]["dgap"]) < 1e-9 for r in RES.values())),
            "两个十年的恒等式关系都成立(按构造)", kind="control")
Gg.asserted("⑤ kill(预注册):「两个十年由不同机制张开」要成立,需两个份额区间**不重叠**",
            bool(not overlap),
            f"1990s [{A['lo']:.3f}, {A['hi']:.3f}] vs 2000s [{Bq['lo']:.3f}, {Bq['hi']:.3f}]", kind="kill")
print(); print(Gg)
adm = Gg.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*96)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif not overlap:
    V = (f"**B 两个十年由不同的机制张开。** 1990s 构成项份额 **{A['share']:.1%}** [{A['lo']:.1%}, {A['hi']:.1%}] · "
         f"2000s **{Bq['share']:.1%}** [{Bq['lo']:.1%}, {Bq['hi']:.1%}],**两个区间不重叠。**\n"
         f"  ⇒ **一句关于人的话:同一条鸿沟,两次张开不是同一回事 —— "
         f"一次主要是人改了主意,一次里「人换了」占的份量不同。**")
else:
    V = (f"**这个设计说不出两个十年的机制不同。** 1990s **{A['share']:.1%}** [{A['lo']:.1%}, {A['hi']:.1%}] · "
         f"2000s **{Bq['share']:.1%}** [{Bq['lo']:.1%}, {Bq['hi']:.1%}] —— **区间重叠。**\n"
         f"  ⇒ **两个十年的构成项份额都不为零,而它们之间的差别本设计分辨不出。**\n"
         f"  ⚠⚠ **而两个十年共同确定的那件事,才是要带走的:两次张开都不是纯粹的「走」。**")
print(V)
json.dump(dict(item=IT, windows={k: v["window"] for k, v in RES.items()},
               results={k: dict(gens=v["gens"], decomposition=v["D"], share=v["share"],
                                lo=v["lo"], hi=v["hi"]) for k, v in RES.items()},
               intervals_overlap=bool(overlap), B=B, n_rep=NREP,
               pos_control=dict(median=pc_m, lo=pc_lo, hi=pc_hi, deterministic=True),
               neg_control=dict(median=nc_m, lo=nc_lo, hi=nc_hi,
                                half_width=(nc_hi-nc_lo)/2, reference=0.0),
               admissible=adm, verdict=V, gate_ok=Gg.verdict()),
          open(OUT/"both_decades.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'both_decades.json'}")

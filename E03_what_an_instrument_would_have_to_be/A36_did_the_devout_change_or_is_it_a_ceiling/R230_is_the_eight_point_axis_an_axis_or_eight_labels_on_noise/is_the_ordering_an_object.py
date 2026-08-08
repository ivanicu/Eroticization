"""#791 · E03·A45·R230 —— `#789` 那条八点的轴,是一条轴,还是八个标签贴在噪声上?

`#790`① 的信号是:「先算 8 个点上的 **MDE**,再挑候选量;**候选量在看到次序之前写下来**。」

**⚠⚠ 而那句话的后半截我已经做不到了,所以先说清楚,不许假装:我已经看过那个次序。**
`sexeduc 2.19 · racmar 1.84 · spanking 0.46 · suicide2 0.41 · homosex 0.41 · teensex 0.36 ·
prayer 0.13 · helpblk 0.11` 就在上一轮的产物和页面上。
⇒ **任何我现在提出的「解释这个次序」的候选量都是被污染的** —— 它是对着答案编的,
   而 `realstat §2` 明写:**把一个做不到的判准标成「计划中」是禁止的。**
⇒ 本轮**不提候选量**,只做那件**不需要候选量、因而不可能被污染**的事:
   **那个次序本身,在这套数据上分辨得出来吗?**
⚠ 而这不是退而求其次 —— **如果次序分辨不出来,那么「解释它」这个任务根本不存在**,
   `#789`① 那一整条路线就不必走了。**先问对象在不在,再问它为什么是那样。**

G1 估计量(两个,方法之前先命名):
   (a) **可分辨的题对数** —— 8 题 28 对里,`r_i − r_j` 的自助区间在 BH(q=0.05)之后仍排除 0 的对数
   (b) **8 点上的 MDE** —— 一个候选量要在 n=8 上被检出,需要多大的 |Spearman ρ|(80% 功效,α=0.05)
   ⚠ (b) **只依赖设计,不依赖数据,更不依赖任何候选量** ⇒ **它是本轮唯一完全免疫污染的数字。**

识别:(a) 需要**题间相关的**自助 —— 八题问的是**同一批年份里的同一批人**,
   逐题独立自助会把题间相关抹掉,**从而高估可分辨性**。
   ⇒ **联合年份自助**:每一抽把年份从并集里有放回地抽一次,八题各自取它有的那些年。

⚠⚠ 三个世界,而第三个我会很不想要(`frontier §3` 的 basin 逃逸 —— 我上一轮刚把这条轴写上页面):
   A **梯度**:多数对可分辨 ⇒ 那是一条真的连续轴,值得去解释。
   B **两堆**:只有跨越某个界的对可分辨 ⇒ 它不是轴,是一个二分,而「次序」是过度读数。
   C **不是对象**:几乎没有对可分辨 ⇒ **`#789` 的「轴」是八个标签贴在噪声上**,
     而我上一轮把它写进了页面的结论句里。

预测矩阵:
   | 世界 | 现在 | 若 ≥14 对可分辨 | 若 7–13 对 | 若 <7 对 |
   | A 梯度   | 0.35 | **0.85** | 0.10 | 0.02 |
   | B 两堆   | 0.40 | 0.10 | **0.80** | 0.10 |
   | C 非对象 | 0.25 | 0.05 | 0.10 | **0.88** |

预注册判词(条件式,不是阈值):
  if 正控开火(合成的「真有梯度」世界 >= 14 对可分辨)
     and 负控在 g=0 上**不**开火(合成的「八题同一个真值」世界 < 3 对可分辨):
      可分辨对数 >= 14 -> A · 7..13 -> B · < 7 -> C
  else: UNVERIFIED
⚠ **「这个零该不该是零?」** —— 该。零假设是**八题共享同一个真比值**,那么 `r_i − r_j` 的期望
  **恰好是 0**,不是某个偏移量 ⇒ 用 **`negative_control`,不是 `offset_control`**。
  ⇒ 而这个零**必须被合成地造出来检验**(`realstat §1`:先命名它排除的世界,再把那个世界造出来)。

⚠ 跑之前写下的最强混淆:**八题的精度差得很远**(`#789` 里区间宽度从 0.23 到 3.6)。
  ⇒ 「可分辨的对」可能只是「最窄的对最宽的」,而不是「位置真的不同」。
  ⇒ 同一轮里放控制:**报每一对的 `|Δr|` 与两端区间宽度的关系** ——
     若可分辨性由宽度而非 `|Δr|` 预测,**那 A 的读法作废。**

`G3` 多重性:28 对是一个族,**BH 全族校正,非幸存者全部列出。**
`G4` 规格曲线:自助抽数 {1000, 4000} × BH 的 q {0.05, 0.10} × 差的定义 {r, log r} = 8 格全报。

本轮换不了仪器,理由同 `R223/instrument_search.py`(对象是世界,第二具仪器本机六具全部落选)。
⚠ 而本轮还有一条**结构性做不到**的:**我自己是被污染的读者** ——
   要真正预注册一个候选量,需要一个**没看过次序的**头脑,而本会话不许派 agent。
   **登记为不可能,并写下它需要什么:一次由未看过产物的人/进程写下的候选清单。**
"""
import numpy as np, pandas as pd, json, pathlib, sys
from itertools import combinations
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

RNG = np.random.default_rng(230)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
PREV = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R228_that_boundary_does_not_bind_the_quantity_i_need/results/ratio_needs_no_polarity.json"))
ITEMS = [r["item"] for r in PREV["rows"] if r.get("readable")]
print(f"=== ⓪ 对象:`#789` 的 {len(ITEMS)} 道可读题 —— {ITEMS} ===")

# ── (b) MDE:只依赖设计,先算,因为它可能直接关掉整条路线 ────────────────────────
print("\n=== ① 8 个点上的 MDE(只依赖设计 ⇒ 本轮唯一完全免疫污染的数字)===")
def spearman(a, b):
    ra, rb = pd.Series(a).rank().to_numpy(), pd.Series(b).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])
n8 = len(ITEMS)
null_rho = np.array([abs(spearman(np.arange(n8), RNG.permutation(n8))) for _ in range(20000)])
crit = float(np.quantile(null_rho, 0.95))
mde = None
for rho_true in np.arange(0.05, 1.0, 0.01):
    hit = 0
    for _ in range(1500):
        x = RNG.normal(size=n8)
        y = rho_true*x + np.sqrt(max(1e-9, 1-rho_true**2))*RNG.normal(size=n8)
        if abs(spearman(x, y)) > crit: hit += 1
    if hit/1500 >= 0.80: mde = float(rho_true); break
print(f"  n = {n8} · 置换零的 95% 分位 |ρ| = **{crit:.3f}**")
print(f"  ⇒ **MDE(80% 功效)= |ρ| ≈ {mde:.2f}** —— 一个候选量要在这 8 个点上被检出,"
      f"它与次序的相关必须至少这么大。")
print(f"  ⚠ 那是一个**极强**的要求:社会科学里跨题的解释变量很少到 {mde:.2f}。"
      f"**所以就算我能干净地预注册一个候选量,这 8 个点也几乎注定检不出它。**")

# ── 数据 ─────────────────────────────────────────────────────────────────────
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+ITEMS, convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
KMAX = {c: (4 if c == "homosex" else len(cat[c].cat.categories)) for c in ITEMS}
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= KMAX[c])) for c in ITEMS})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["b"] = REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
SER = {}
for c in ITEMS:
    g = REL.dropna(subset=[c])
    for k in (2, 0):
        rows = [(int(y), float(gy[c].mean())) for y, gy in g[g["b"] == k].groupby("year") if len(gy) >= 120]
        SER[(c, k)] = dict(rows)
UNION = sorted({y for c in ITEMS for k in (2, 0) for y in SER[(c, k)]})
print(f"\n=== ② 联合年份自助:年份并集 {len(UNION)} 年 · 逐题独立自助会抹掉题间相关 ⇒ 高估可分辨性 ===")

def ratios(draw_years):
    out = {}
    for c in ITEMS:
        yA = [y for y in draw_years if y in SER[(c, 2)]]
        yB = [y for y in draw_years if y in SER[(c, 0)]]
        if len(set(yA)) < 8 or len(set(yB)) < 8: return None
        sA = slope(yA, [SER[(c, 2)][y] for y in yA]); sB = slope(yB, [SER[(c, 0)][y] for y in yB])
        if abs(sB) < 1e-12: return None
        out[c] = sA/sB
    return out

OBS = ratios(UNION)
print("  观测比值:" + " · ".join(f"{c} {OBS[c]:+.3f}" for c in sorted(OBS, key=lambda x: -OBS[x])))

def run(B, diff_kind):
    draws = []
    while len(draws) < B:
        r = ratios(list(RNG.choice(UNION, len(UNION), replace=True)))
        if r is not None: draws.append(r)
    P = []
    for a, b in combinations(ITEMS, 2):
        if diff_kind == "log":
            if OBS[a] <= 0 or OBS[b] <= 0: continue
            v = np.array([np.log(x[a])-np.log(x[b]) for x in draws if x[a] > 0 and x[b] > 0])
            obs = np.log(OBS[a])-np.log(OBS[b])
        else:
            v = np.array([x[a]-x[b] for x in draws]); obs = OBS[a]-OBS[b]
        p = 2*min((v <= 0).mean(), (v >= 0).mean())
        P.append(dict(a=a, b=b, obs=float(obs), p=float(max(p, 1/(len(v)+1))),
                      lo=float(np.percentile(v, 2.5)), hi=float(np.percentile(v, 97.5))))
    return P, draws

def bh(P, q):
    S = sorted(P, key=lambda r: r["p"]); C = len(S)
    keep = 0
    for i, r in enumerate(S, 1):
        if r["p"] <= q*i/C: keep = i
    for i, r in enumerate(S, 1): r["survives"] = bool(i <= keep)
    return S, keep

print("\n=== ③ `G4` 规格曲线:抽数 × BH 的 q × 差的定义 = 8 格,全报 ===")
grid = []
for B in (1000, 4000):
    for kind in ("r", "log"):
        P, draws = run(B, kind)
        for q in (0.05, 0.10):
            S, keep = bh([dict(x) for x in P], q)
            grid.append(dict(B=B, kind=kind, q=q, n_pairs=len(S), n_sig=keep))
            print(f"  B={B:<5d} 差={kind:3s} q={q:.2f} → **可分辨 {keep} / {len(S)} 对**")
MAIN = [g for g in grid if g["B"] == 4000 and g["kind"] == "r" and g["q"] == 0.05][0]
P_main, draws_main = run(4000, "r")
S_main, keep_main = bh([dict(x) for x in P_main], 0.05)

print(f"\n=== ④ 主格(B=4000 · 差=r · q=0.05):**可分辨 {keep_main} / {len(S_main)} 对**,非幸存者全部列出 ===")
for r in S_main:
    print(f"  {'✅' if r['survives'] else '  '} {r['a']:9s} − {r['b']:9s} = {r['obs']:+7.3f} "
          f"[{r['lo']:+7.3f}, {r['hi']:+7.3f}]  p={r['p']:.4f}")

# ── 跑前写下的混淆:可分辨性是不是只是「窄 vs 宽」 ────────────────────────────────
WID = {c: float(np.percentile([x[c] for x in draws_main], 97.5)
                - np.percentile([x[c] for x in draws_main], 2.5)) for c in ITEMS}
dv = np.array([abs(r["obs"]) for r in S_main])
wv = np.array([WID[r["a"]]+WID[r["b"]] for r in S_main])
sig = np.array([1.0 if r["survives"] else 0.0 for r in S_main])
print(f"\n=== ⑤ 跑前写下的混淆的控制:可分辨性由 `|Δr|` 预测,还是由「两端有多宽」预测? ===")
print(f"  corr(可分辨, |Δr|)   = **{np.corrcoef(sig, dv)[0,1]:+.3f}**")
print(f"  corr(可分辨, 宽度和) = **{np.corrcoef(sig, wv)[0,1]:+.3f}**")
conf_ok = bool(abs(np.corrcoef(sig, dv)[0, 1]) > abs(np.corrcoef(sig, wv)[0, 1]))
print(f"  ⇒ {'位置差主导 ⇒ 可以按位置读' if conf_ok else '**宽度主导 ⇒ 「可分辨」量的是精度不是位置,A 的读法作废**'}")

# ── 正控 / 负控:合成两个世界 ────────────────────────────────────────────────────
print("\n=== ⑥ 控制:两个合成世界,同一条流水线 ===")
def synth(spread, B=1000):
    base = np.linspace(0.1, 2.2, n8) if spread else np.full(n8, 0.5)
    sd = float(np.median([np.std([x[c] for x in draws_main], ddof=1) for c in ITEMS]))
    D = [{c: base[i] + RNG.normal(0, sd) for i, c in enumerate(ITEMS)} for _ in range(B)]
    obs = {c: base[i] for i, c in enumerate(ITEMS)}
    P = []
    for a, b in combinations(ITEMS, 2):
        v = np.array([x[a]-x[b] for x in D])
        p = 2*min((v <= 0).mean(), (v >= 0).mean())
        P.append(dict(a=a, b=b, p=float(max(p, 1/(len(v)+1)))))
    return bh(P, 0.05)[1]
pos_n = synth(True); neg_n = synth(False)
print(f"  正控(真有梯度 0.1→2.2,噪声取自实测中位 sd):可分辨 **{pos_n}/28** 对")
print(f"  负控(g=0:八题同一个真值 0.5,其余一切不变):可分辨 **{neg_n}/28** 对")
print(f"  ⚠ 「这个零该不该是零?」—— **该**:零假设是八题共享一个真比值,`r_i − r_j` 的期望恰为 0 "
      f"⇒ 用 `negative_control` 而不是 `offset_control`。")

# ── ⑦ ⚠⚠ 跑完之后补的:预注册数了对数,**没看那些对长什么样** ────────────────────
# 16/28 触发了「A 梯度」,而把可分辨的对**列出来**之后,图样根本不是梯度:
# 可分辨的对几乎全是**跨越某个界**的,而**界内的对一个都不可分辨**。
# ⇒ 判据数了一个**计数**,而问题问的是一个**结构** —— 「判词测错问题」那一族的又一次。
# ⚠ 而分组是**看着同一批数据划的** ⇒ 下面这两个数是**描述,不是检验**,必须这么标。
CLUMP = {"sexeduc": "顶", "racmar": "顶",
         "homosex": "中", "teensex": "中", "suicide2": "中", "spanking": "中",
         "prayer": "底", "helpblk": "底"}
win = [r for r in S_main if CLUMP[r["a"]] == CLUMP[r["b"]]]
btw = [r for r in S_main if CLUMP[r["a"]] != CLUMP[r["b"]]]
win_sig = sum(1 for r in win if r["survives"]); btw_sig = sum(1 for r in btw if r["survives"])
print(f"\n=== ⑦ 那些可分辨的对长什么样(⚠ 分组读自同一批数据 ⇒ 描述,不是检验)===")
print(f"  组内对:**{win_sig} / {len(win)} 可分辨**   组间对:**{btw_sig} / {len(btw)} 可分辨**")
print(f"  顶 = {[k for k,v in CLUMP.items() if v=='顶']} (~1.8–2.2,虔诚者改得**更多**)")
print(f"  中 = {[k for k,v in CLUMP.items() if v=='中']} (~0.36–0.46)")
print(f"  底 = {[k for k,v in CLUMP.items() if v=='底']} (~0.11–0.13)")
print(f"  ⇒ **不是一条梯度,是三堆** —— 而预注册的判据只数了 {keep_main}/28,**没问那些对长什么样。**")
print(f"  ⚠ **本轮照预注册执行(它判 A),同时把这一段写下来**:判据数了计数,问题问的是结构。")

G = Gate("#791 · 那条八点的轴,是一条轴还是八个标签")
G.negative_control("① 负控:八题同一个真值时,可分辨对数必须塌到 0 附近(<3)",
                   null=float(neg_n), effect=float(pos_n), ratio=0.5,
                   null_kind="八题共享同一个真比值的合成世界(差的期望恰为 0)")
G.asserted("② 正控:真有梯度时必须 >= 14 对可分辨(否则流水线连梯度都看不见)",
           bool(pos_n >= 14), f"合成梯度世界 {pos_n}/28", kind="control")
G.asserted("③ 前提(跑前写下的混淆):可分辨性必须由 |Δr| 而非区间宽度主导",
           conf_ok, f"corr(可分辨,|Δr|)={np.corrcoef(sig,dv)[0,1]:+.3f} vs "
                    f"corr(可分辨,宽度和)={np.corrcoef(sig,wv)[0,1]:+.3f}", kind="control")
G.asserted("⑤ 描述(非检验):组内可分辨对必须接近 0,否则「三堆」这个描述连描述都不成立",
           bool(win_sig == 0), f"组内 {win_sig}/{len(win)} · 组间 {btw_sig}/{len(btw)}", kind="control")
G.asserted("④ kill(预注册):「那是一条值得解释的轴」要站住,需主格 >= 14/28 对可分辨",
           bool(keep_main >= 14), f"主格 {keep_main}/28 · 规格曲线 {[g['n_sig'] for g in grid]}", kind="kill")
print(); print(G)

print("\n"+"="*92)
ctrl = bool(neg_n < 3 and pos_n >= 14 and conf_ok)
if not ctrl:
    v = (f"**UNVERIFIED:控制没过(负控 {neg_n}/28 需 <3 · 正控 {pos_n}/28 需 ≥14 · "
         f"混淆控制 {'过' if conf_ok else '没过'})⇒ 本轮不下判。**")
elif keep_main >= 14:
    v = (f"**预注册判 A(主格 {keep_main}/28 ≥ 14),而我照它执行 —— 但把可分辨的对列出来之后,"
         f"图样不是梯度,是三堆。**\n"
         f"  **组内 {win_sig}/{len(win)} 可分辨 · 组间 {btw_sig}/{len(btw)} 可分辨** ——\n"
         f"  顶 `sexeduc`+`racmar`(~1.8–2.2,虔诚者改得**更多**)· "
         f"中 `homosex`+`teensex`+`suicide2`+`spanking`(~0.36–0.46)· 底 `prayer`+`helpblk`(~0.11–0.13)。\n"
         f"  ⚠ **分组读自同一批数据 ⇒ 这是描述,不是检验**;要检验它需要一批**留出的题**,"
         f"而 `#790` 已经证明这具仪器上没有更多题。\n"
         f"  ⚠ **预注册的判据数了一个计数,而问题问的是一个结构** —— 「判词测错问题」那一族的又一次。\n"
         f"  ⇒ **而无论三堆还是梯度,MDE 都把解释这条路关上了:需要 |ρ| ≥ {mde:.2f},8 个点给不了。**")
elif keep_main >= 7:
    v = (f"**B 两堆,不是一条轴。** 主格只有 {keep_main}/28 对可分辨(规格曲线 {[g['n_sig'] for g in grid]})"
         f" ⇒ **数据支持的是一个二分,而不是一个次序** —— "
         f"我在 `#789` 里把八个点排成一列读,**那是过度读数**。")
else:
    v = (f"**C 不是一个对象。** 主格只有 {keep_main}/28 对可分辨(规格曲线 {[g['n_sig'] for g in grid]})"
         f" ⇒ **`#789` 的那条「轴」是八个标签贴在噪声上**,而我上一轮把它写进了页面的结论句。")
print(v)
print(f"\n⚠ 而无论上面哪一支,**MDE 那一条都独立成立**:n=8 上要检出一个候选量需要 "
      f"**|ρ| ≥ {mde:.2f}**(置换零 95% 分位 {crit:.3f})—— "
      f"**`#789`① 那条路线在这 8 个点上是关着的,不是因为我想不出候选量,是因为算术。**")
json.dump(dict(items=ITEMS, obs=OBS, clumps=CLUMP,
               within_sig=win_sig, within_n=len(win), between_sig=btw_sig, between_n=len(btw),
               clump_partition_is_descriptive=True, mde=mde, perm_crit=crit, grid=grid,
               pairs=S_main, n_sig=keep_main, widths=WID,
               corr_sig_absdiff=float(np.corrcoef(sig, dv)[0, 1]),
               corr_sig_width=float(np.corrcoef(sig, wv)[0, 1]),
               pos_control=pos_n, neg_control=neg_n, verdict=v,
               contaminated_reader=True, gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"is_the_ordering_an_object.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'is_the_ordering_an_object.json'}")

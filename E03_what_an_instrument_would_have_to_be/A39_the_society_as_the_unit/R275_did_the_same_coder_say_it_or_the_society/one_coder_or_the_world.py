"""#836 · E03·A78·R275 —— 换单位到社会:那条「谴责↔稀少」的关联,是社会说的,还是同一个编码者说的?

⚠⚠ **先说为什么换仪器换单位**:`git log` 实测,**本会话 16 轮全在 GSS 一具仪器、一个单位(人/十年)上**,
   而 `E02` 的对象陈述明写**单位也包括社会(SCCS 186)**,硬规则④明写**跨仪器复现胜过同一具再跑一轮**。
   `#835`① 又指向「另一具仪器」。⇒ **本轮把单位换成社会。**

**⚠⚠ ② 当轮更正 —— 而它是被自己的负控抓出来的,记在这里因为它是本轮最硬的一条方法论收获:**
   第一版的 `rho()` 用 `np.argsort(np.argsort(·))` 取秩,那是**序数秩**。
   **而本轮的两个变量几乎全是并列**(结局二值、谴责 4 级、n=34)——
   序数秩把并列项按**原始下标顺序**摊开,于是 `x` 与 `y` 的秩都带上一个**下标位置分量**,
   **两者因此天然正相关。** 实测:自由置换零的中心是 **+0.0468**,而它自己的 SE 只有 0.0026 ——
   **偏离 0 有 18 个标准误。一个置换零必须以 0 为中心,这不是噪声,是统计量错了。**
   ⇒ 改为**平均秩**(并列取组内平均),即真正的 Spearman ρ。
   **⚠ 值得记的不是这个 bug,是它被什么抓住的:**
   **负控问的是「打乱配对后还剩什么」,而剩下的那 0.0468 不属于数据,属于我的秩函数。**
   **`G2` 那句「一个从未返回非零的仪器给的零是沉默」有一个对偶面,这轮撞上了:**
   **一个本该返回零却不返回零的零,指的是仪器,不是世界。**

**⚠ 硬规则①先做,而它当场改了设计两次(读定义与码值,不在题名上下结论):**
   ① `SCCS176`「Homosexuality」**不是一条有序量表** —— 码 2 是 **"None"**,那是一句**频率**陈述,
     坐在一条本来有序的谴责量表(1 接受 → 5 强烈反对)里面。**按题名当 5 级有序用就是错的。**
     ⇒ **剔除码 2**(n 34,而非 38)。
   ② `SCCS634`「Control of Sex Scale」**不是宽严量表** —— 码值是「对女性控制更严 → 男女相当」,
     **它是性别双标尺度**。**拿它当「性压抑」用是构念错误,已弃用。**

**E01 的那条被降级的主张,正好在这里有它自己的两个构念:**
`corr(rarity, shame) = 0.758` 曾被降级为「**两者都是谴责的下游**」——
而 SCCS 里恰有 **稀少(`SCCS177` 频率)** 与 **谴责(`SCCS176` 态度)**,在**社会**这个单位上。

**⚠⚠⚠ 而硬规则②在这里是全项目最尖锐的一次,必须前置,不能只是命名:**
`SCCS176` 与 `SCCS177` **出自同一篇论文(Broude & Greene 1976)、同两位编码者、同一批民族志。**
**一个编码者读到「此地同性恋罕见且强烈不容」这一句,就同时填了两个格。**
**更糟:如果谴责把行为逼入隐蔽,民族志作者就记录「没有」——**
**谴责在因果上制造了「稀少」的外观,而这不是一个可控制的混淆,它就是测量本身。**
⇒ **只跑那一对,无论结果如何都不可解释。**

⇒ **所以本轮的设计是围绕「打破共享编码者」建的:**
   **同一个结局(`SCCS177` 稀少),三条臂,其中两条来自完全不同的研究团队:**
   · **A 同编码者**:`SCCS176` 谴责(Broude & Greene)· n=34 · MDE|r|₈₀ ≈ 0.464
   · **B 独立编码者**:`SCCS602`「明确认为性活动危险/污染」(**Whyte 1978**)· n=33 · MDE ≈ 0.471
   · **C 独立编码者**:`SCCS961`「婚前性行为限制」(**Frayser 1985**)· n=21 · MDE ≈ 0.578
   **B 与 C 是不同的人、读不同的材料、为不同的目的编的码。**

G1 估计量:**每条臂上,社会层面的 `Spearman ρ`(限制/谴责 越强 ↔ 记录到的同性行为 越「缺席」)。**

三个世界:
   A **独立编码者也看得见**:B 或 C 与 A 同号且 |ρ| 超过它自己的 MDE ⇒
     **这条关联经得住换编码者 ⇒ 在社会这个单位上,E01 那条被降级的主张得到跨仪器支持。**
   B **只有同编码者看得见** ⇒ **那条关联是关于 Broude & Greene 怎么读民族志的事实,不是关于社会的** ——
     **这是硬规则②预言的结果,也是我不欢迎的那个。**
   C **三条臂都看不见** ⇒ **E01 的相关不外推到社会这个单位。**

预测矩阵:
   | 世界 | 现在 | B/C 同号且超 MDE | 只有 A | 三条都无 |
   | A 经得住换编码者 | 0.35 | **0.85** | 0.05 | 0.05 |
   | B 是编码者的事实 | 0.40 | 0.05 | **0.85** | 0.10 |
   | C 不外推       | 0.25 | 0.10 | 0.10 | **0.85** |

预注册判词(条件式):
  if 正控开火(**植入一个已知的关联,三条臂都必须取回**)
     and 负控开火(**自由置换社会标签后 ρ 的分布必须以 0 为中心** ——
        ⚠ **「这个零该不该是零?」该。** 打乱配对后相关的期望就是 0,**用 `negative_control` 对 0**):
      B 或 C 与 A 同号且 |ρ| > 该臂自己的 MDE -> A
      只有 A 显著                            -> B
      三条都不显著                            -> C
  else: UNVERIFIED
⚠ **`G3` 多重性:三条臂一族 BH。⚠ 而 `#832` 已证族越窄存活越易 —— 三条是很窄的一族,如实标注。**
⚠ **`G4` 两个零都跑**:自由置换 · **按地理区块置换**(见下)。

⚠⚠ **跑之前写下的最强混淆 —— Galton 问题(跨文化统计的经典威胁):**
   **各社会不是独立单位** —— 相邻/同源的社会共享历史,任何两个变量都会因为共同祖先而相关。
   ⇒ 控制:**除自由置换外,再跑一个「只在地理区块内置换」的零**;
   区块由 `Lat/Long` 用**明写的切分**造(不靠我记 Murdock 的区号)。
   **若关联在区块内置换下消失,那它可能只是地理自相关。**

⚠⚠ **而有一条边界无论哪个世界都成立,必须写在结论里:**
   **三条臂共用同一个结局变量(`SCCS177`,Broude & Greene 编)** ——
   **所以即使世界 A 成立,能说的也只是「两种独立的『限制』编码都与同一份『稀少』编码同向」,
   而不是「限制导致稀少」。结局这一侧永远是单仪器的,这个设计改不了。**

⚠ 本轮**换不了仪器**的说法不适用 —— 本轮**正是**换仪器(GSS → SCCS)与换单位(人/十年 → 社会)。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys, csv, collections, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
S = ROOT/"data/external/dplace/repo/datasets/SCCS"
NPERM, Q, SEEDS = 20000, 0.05, [275, 276, 277]

VAR = {v["id"]: v for v in csv.DictReader(open(S/"variables.csv", encoding="utf-8"))}
SOC = {r["id"]: r for r in csv.DictReader(open(S/"societies.csv", encoding="utf-8"))}
D = collections.defaultdict(dict)
for r in csv.DictReader(open(S/"data.csv", encoding="utf-8")):
    if r["code"] not in ("", "NA"): D[r["soc_id"]][r["var_id"]] = r["code"]

print("=== ⓪ 硬规则①(已在设计前跑过一遍,这里复打以留痕)===")
for vid in ("SCCS177", "SCCS176", "SCCS602", "SCCS961"):
    n = sum(1 for s in D if vid in D[s])
    print(f"  {vid} n={n:>3d}/186 · 来源 {VAR[vid]['source'][:26]:26s} 「{VAR[vid]['title'][:52]}」")
print("  ⚠ `SCCS176` 码 2 = 'None' 是**频率**陈述,不在谴责维度 ⇒ **剔除**")
print("  ⚠ `SCCS634`「Control of Sex Scale」实为**性别双标**尺度 ⇒ **弃用**(构念错误)")

# 地理区块(切分明写,不靠记 Murdock 区号)
def region(sid):
    try: la, lo = float(SOC[sid]["Lat"]), float(SOC[sid]["Long"])
    except Exception: return "NA"
    if -35 <= la <= 37 and -20 <= lo <= 52: return "Africa"
    if 20 <= la <= 70 and -15 <= lo <= 60: return "CircumMed"
    if lo >= 60 and la >= 0: return "EastEurasia"
    if lo >= 90 and la < 0: return "InsularPacific"
    if la >= 12 and lo < -30: return "NorthAmerica"
    if la < 12 and lo < -30: return "SouthAmerica"
    return "Other"

ARMS = {
    "A 同编码者 `SCCS176` 谴责 (Broude&Greene)": ("SCCS176", True),
    "B 独立 `SCCS602` 性危险/污染 (Whyte)": ("SCCS602", False),
    "C 独立 `SCCS961` 婚前限制 (Frayser)": ("SCCS961", False),
}
def _avgrank(v):
    """⚠⚠ **平均秩,不是序数秩** —— 见文件头 `#836`② 的当轮更正。"""
    v = np.asarray(v, float); o = np.argsort(v, kind="mergesort"); r = np.empty(len(v), float)
    i = 0
    while i < len(v):
        j = i
        while j+1 < len(v) and v[o[j+1]] == v[o[i]]: j += 1
        r[o[i:j+1]] = (i+j)/2.0 + 1.0
        i = j+1
    return r

def rho(x, y):
    rx, ry = _avgrank(x), _avgrank(y)
    if rx.std() == 0 or ry.std() == 0: return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])

print(f"\n=== ① 三条臂(结局统一为 `SCCS177` 稀少:1 缺席 / 2 在场)· 置换 {NPERM} 次 ===")
RES = {}
for nm, (vid, same) in ARMS.items():
    socs = [s for s in D if vid in D[s] and "SCCS177" in D[s]]
    if vid == "SCCS176": socs = [s for s in socs if D[s][vid] != "2"]
    x = [float(D[s][vid]) for s in socs]; y = [float(D[s]["SCCS177"]) for s in socs]
    reg = [region(s) for s in socs]
    n = len(socs); r = rho(x, y)
    mde = math.tanh((1.96+0.84)/math.sqrt(n-3)) if n > 3 else np.nan
    rngs = [np.random.default_rng(sd) for sd in SEEDS]
    free = np.array([rho(x, list(rg.permutation(y))) for rg in rngs for _ in range(NPERM//len(SEEDS))])
    # 区块置换:只在同一地理区内打乱结局
    blk = []
    idx = collections.defaultdict(list)
    for i, g in enumerate(reg): idx[g].append(i)
    for rg in rngs:
        for _ in range(NPERM//len(SEEDS)):
            yy = list(y)
            for g, ii in idx.items():
                vals = [y[i] for i in ii]; rg.shuffle(vals)
                for i, v in zip(ii, vals): yy[i] = v
            blk.append(rho(x, yy))
    blk = np.array(blk)
    p_free = float((np.sum(np.abs(free) >= abs(r))+1)/(len(free)+1))
    p_blk = float((np.sum(np.abs(blk) >= abs(r))+1)/(len(blk)+1))
    RES[nm] = dict(var=vid, same_coder=same, n=n, rho=r, mde=mde, p_free=p_free, p_blk=p_blk,
                   null_free_mean=float(np.mean(free)), null_free_sd=float(np.std(free)),
                   null_blk_sd=float(np.std(blk)),
                   regions=dict(collections.Counter(reg)))
    print(f"  {nm[:44]:46s} n={n:>3d} ρ=**{r:+.3f}** (MDE {mde:.3f}) · "
          f"p_自由={p_free:.4f} · p_区块={p_blk:.4f}")
print(f"  ⚠ 区块由 Lat/Long 明写切分而来;各臂区块分布已存入产物。")

ps = [RES[nm]["p_free"] for nm in ARMS]
surv = {list(ARMS)[i] for i in Gate.bh(ps, Q)}
print(f"\n  `G3` BH(3 条臂一族,q={Q}):存活 **{len(surv)}/3** ⇒ {sorted(x[:20] for x in surv) or '无'}")
print(f"  ⚠ **三条是很窄的一族 —— `#832` 已证族越窄存活越易,如实标注,不当成加分。**")
A = RES[list(ARMS)[0]]; B = RES[list(ARMS)[1]]; C = RES[list(ARMS)[2]]
indep_ok = [nm for nm in list(ARMS)[1:]
            if np.sign(RES[nm]["rho"]) == np.sign(A["rho"]) and abs(RES[nm]["rho"]) > RES[nm]["mde"]]
print(f"  ⇒ 与 A 同号且 |ρ| 超自身 MDE 的**独立臂**:**{len(indep_ok)}** ⇒ {indep_ok or '无'}")

print("\n=== ② 控制 ===")
rg = np.random.default_rng(999)
socsA = [s for s in D if "SCCS176" in D[s] and "SCCS177" in D[s] and D[s]["SCCS176"] != "2"]
xa = np.array([float(D[s]["SCCS176"]) for s in socsA])
plant = np.where(xa >= 4, 1.0, 2.0)                    # 造一个已知强关联的结局
r_pc = rho(xa, plant)
nullp = np.array([rho(xa, list(rg.permutation(plant))) for _ in range(4000)])
p_pc = float((np.sum(np.abs(nullp) >= abs(r_pc))+1)/(len(nullp)+1))
print(f"  正控:把结局按谴责码构造成强关联 ⇒ ρ = **{r_pc:+.3f}**,p = {p_pc:.5f} —— 该**明显非零**")
ya = np.array([float(D[s]["SCCS177"]) for s in socsA])
nc = np.array([rho(xa, list(rg.permutation(ya))) for _ in range(4000)])
print(f"  负控:**打乱社会配对** ⇒ ρ 分布 中心 **{np.mean(nc):+.4f}** · SD {np.std(nc):.4f} —— "
      f"⚠ **「这个零该不该是零?」该** —— 打乱配对后相关的期望**就是 0**,所以用 `negative_control` 对 0")
print(f"     ⚠ `realstat`:置换零只回答「配对重不重要」,从不回答「为什么」——"
      f"它排除的世界是「限制码与稀少码随机配对」,而那正是本轮构造出来的世界。")

G = Gate("#836 · 那条关联是社会说的,还是同一个编码者说的")
G.asserted("① 硬规则①(在设计前跑,当场改了设计两次):`SCCS176` 码 2='None' 是频率陈述**不在谴责维度** ⇒ 剔除;"
           "`SCCS634`「Control of Sex Scale」实为**性别双标**尺度 ⇒ **弃用(构念错误)**"
           " —— **两处都是「按题名下结论」会犯的错**",
           True, f"A 臂 n={A['n']}(剔除码 2 前为 38)· SCCS634 未进入任何一条臂", kind="control")
G.asserted("② 硬规则②(本项目最尖锐的一次):结局与 A 臂**同出 Broude&Greene 一篇论文、同两位编码者、同一批民族志**"
           " ⇒ **只跑那一对无论结果如何都不可解释** ⇒ 设计围绕**打破共享编码者**建:"
           "B(Whyte)与 C(Frayser)是**不同的人读不同材料为不同目的编的码**",
           bool(not B["same_coder"] and not C["same_coder"]),
           f"B 来源 {VAR[B['var']]['source'][:20]} · C 来源 {VAR[C['var']]['source'][:20]}", kind="control")
G.asserted("③ 正控:把结局按谴责码构造成强关联,必须取回一个明显非零的 ρ",
           bool(abs(r_pc) > 0.5 and p_pc < 0.01), f"ρ={r_pc:+.3f}, p={p_pc:.5f}", kind="control")
G.identity_control("④ 负控:**打乱社会配对**后 ρ 分布中心必须 **== 0**"
                   "(⚠ **「这个零该不该是零?」该** —— 打乱配对后相关的期望**就是 0**,"
                   "所以对 0 而不是对某个观测量;⚠⚠ **而正是这一行抓出了本轮的统计量缺陷,见文件头②**)",
                   observed=float(np.mean(nc)), expected=0.0,
                   tol=3*float(np.std(nc))/math.sqrt(len(nc)),
                   what="4000 次自由置换的 ρ 分布中心",
                   noise_half_width=float(np.std(nc))/math.sqrt(len(nc)))
G.asserted("⑤ 前提(跑前写下的最强混淆 · Galton 问题):**各社会不是独立单位** ⇒ "
           "除自由置换外**再跑一个只在地理区块内置换的零**,区块由 Lat/Long **明写切分**而来",
           bool(all("p_blk" in RES[nm] for nm in ARMS)),
           " · ".join(f"{nm[0]} p区块={RES[nm]['p_blk']:.3f}" for nm in ARMS), kind="control")
G.asserted("⑥ kill(预注册):「这条关联经得住换编码者」要成立,需**至少一条独立臂**与 A 同号且 |ρ| 超它自己的 MDE",
           bool(len(indep_ok) >= 1), f"独立臂达标 {len(indep_ok)}/2 ⇒ {indep_ok or '无'}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif len(indep_ok) >= 1:
    V = (f"**A 经得住换编码者。** A(同编码者)ρ={A['rho']:+.3f};独立臂达标:{indep_ok}。\n"
         f"  ⇒ **一句关于人的话:在一百八十六个社会这个尺度上,「越是压制,越是记录不到」——\n"
         f"  而这一次说它的不只是同一批人:换一队人类学家、读不同的材料、为不同的目的编码,\n"
         f"  仍然看得见同一条关联。**")
elif abs(A["rho"]) > A["mde"]:
    V = (f"**B 只有同编码者看得见 —— 而这正是硬规则②预言的、我不欢迎的那个结果。**\n"
         f"  A(同编码者)ρ={A['rho']:+.3f}(MDE {A['mde']:.3f});"
         f"独立臂 B ρ={B['rho']:+.3f}(MDE {B['mde']:.3f})· C ρ={C['rho']:+.3f}(MDE {C['mde']:.3f}),**无一达标**。\n"
         f"  ⇒ **那条「谴责↔稀少」的关联,在这份数据上是关于 Broude & Greene 怎么读民族志的事实,\n"
         f"  而不是关于社会的事实 —— 换一队编码者就看不见了。**")
else:
    V = (f"**C 三条臂都看不见。** A ρ={A['rho']:+.3f} · B {B['rho']:+.3f} · C {C['rho']:+.3f},"
         f"皆未超各自 MDE。\n"
         f"  ⇒ **E01 那条人层面的相关,在社会这个单位上没有对应物 —— 至少这份编码看不到。**")
print(V)
print("\n⚠⚠ **一条无论哪个世界都成立的边界**:三条臂**共用同一个结局变量**(`SCCS177`,Broude&Greene 编)。")
print("   **所以即使世界 A 成立,能说的也只是「两种独立的『限制』编码与同一份『稀少』编码同向」,**")
print("   **而不是「限制导致稀少」。结局这一侧永远是单仪器的,这个设计改不了。**")
print("⚠ 且:**若谴责把行为逼入隐蔽,民族志作者就记录『没有』** —— 谴责在因果上制造『稀少』的外观,")
print("   **那不是可控制的混淆,是测量本身。本轮不能、也不打算把它分开。**")
json.dump(dict(unit="society (SCCS 186)", instrument="D-PLACE/SCCS ethnographic codes",
               arms=RES, bh_survivors=sorted(surv), indep_ok=indep_ok, q=Q, nperm=NPERM, seeds=SEEDS,
               pos_control=dict(rho=r_pc, p=p_pc),
               neg_control=dict(center=float(np.mean(nc)), sd=float(np.std(nc)), reference=0.0),
               hard_rule_1=["SCCS176 code 2 'None' is a frequency statement, excluded",
                            "SCCS634 is a gender double-standard scale, discarded"],
               boundary="outcome SCCS177 is single-instrument in every arm; condemnation may CREATE apparent absence",
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"one_coder_or_the_world.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'one_coder_or_the_world.json'}")

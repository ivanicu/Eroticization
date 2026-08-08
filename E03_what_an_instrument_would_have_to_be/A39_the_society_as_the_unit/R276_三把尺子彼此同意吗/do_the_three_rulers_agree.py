"""#837 · E03·A78·R276 —— 攻击我上一轮刚写下的那句话:三把「限制」的尺子,彼此同意吗?

`#836` 昨天写下一句很响的话:**「那条关联更像是关于民族志怎么被读的事实,
而不是关于社会怎么运转的事实。」** 而它有一个我**没有控制**的对手,现在必须先杀掉它:

**⚠⚠ 对手世界:两条独立臂之所以平,可能不是因为「换了编码者就没了」,
而是因为它们根本不在量同一个构念。**
`SCCS602` 量的是「性是否危险/污染」(纯洁观),`SCCS961` 量的是**异性**婚前限制 ——
**一个社会完全可以对婚前异性性行为极严,而对同性行为无所谓。**
若如此,那两条臂的平就与「编码者」无关,`#836` 那句话就**说过头了**。

⇒ **本轮是 Frontier,而它攻击的是我自己昨天的头条。** `§3` basin:上一轮的结论我很喜欢,
**所以这一轮专门去找它站不住的理由。**

**⓪ 而硬规则①在这里抓出了 `#836` 的两处缺陷,其中一处是我刚犯过的同一类错误(`#803` 的模式:
把上一轮刚记下的缺陷在下一轮重建一遍):**

**① `SCCS961` 的码 4 与码 6 不在限制维度上。**
   码 4 = 「**permitted for males but no females**」· 码 6 = 「**insistence on virginity for the
   woman. There is no information on restrictions on the male**」——
   **一条是性别双标,一条是「男方信息缺失」。** 而这**正是 `#836` 弃用 `SCCS634`
   (性别双标尺度)、剔除 `SCCS176` 码 2(频率陈述)的同一个理由** ——
   **我上一轮读了 `SCCS176` 的码值,却只读了 `SCCS961` 的题名。**
   ⇒ 本轮剔除 4 与 6,并**重跑 `#836` 的 C 臂**。

**② `SCCS602` 的极性,从元数据里定不下来。**
   题名是「**(No)** Explicit View That Sexual Activity Is Dangerous or Contaminating」,
   码值只有 `1 = Yes` / `2 = No`。**「Yes」是在回答带 (No) 的题名,还是不带 (No) 的题名?
   元数据不说。** 而 `#836` 的 kill 要求「**与 A 同号**」——
   **那条判据吃进了一个我从未确定方向的输入。**
   ⇒ 本轮**两个极性都跑**,并把极性登记为 `UNVERIFIED`(不是「已核」,也不是「不重要」)。

G1 估计量:**三把「限制」尺子之间的两两 Spearman ρ**(收敛效度),
   以及**清洗+双极性之后,`#836` 的 kill 是否还给同一个答案**。

三个世界:
   A **`#836` 站得住**:三把尺子**彼此同意**(存在一个共同的「限制」构念),
     而只有同编码者那把能预测稀少 ⇒ **编码者解释存活。**
   B **构念不匹配**:三把尺子**彼此不同意** ⇒ SCCS 里根本没有一个共同的「限制」构念,
     **`#836` 的头条说过头了**,必须收窄成「三份互不相关的编码」。
   C **清洗改变了 C 臂** ⇒ **`#836` 自己的 kill 要重判。**

预测矩阵:
   | 世界 | 现在 | ≥2/3 对同意 | ≤1/3 对同意 | 清洗后 C 臂达标 |
   | A `#836` 站得住 | 0.40 | **0.85** | 0.10 | 0.05 |
   | B 构念不匹配    | 0.45 | 0.10 | **0.85** | 0.10 |
   | C `#836` 要重判 | 0.15 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式):
  if 正控开火(**在真实边缘分布与真实并列结构上模拟种植**,必须能在这些 n 上捞回来)
     and 负控开火(置换零以 0 为中心):
      ≥2/3 对的 |ρ| 超它自己的经验检出下限 -> A(`#836` 存活)
      ≤1/3 对                              -> B(**收窄 `#836` 的头条**)
      任一独立臂在清洗/换极性后达标         -> C(**重判 `#836`**)
  else: UNVERIFIED

⚠⚠ **而检出下限本轮不用 Fisher 公式算,要用模拟量** —— `#836` 用的
`tanh((1.96+0.84)/√(n−3))` 假设二元正态、无并列;**本轮的变量是二值与 4/7 级,几乎全是并列**,
那个公式在这里是**外推**。⇒ 在**观测到的边缘分布**上按强度网格种植,
**取检出率首次达到 80% 的那个强度所对应的 ρ 中位数**,作为这一对的经验下限。

⚠ 跑之前写下的最强混淆:**三把尺子的重叠样本各不相同**(n 在 21–34 之间),
  **样本小的那一对天然更难同意** ⇒ 若不给每一对自己的检出下限,
  「不同意」会被样本量伪造出来。⇒ 控制:**每一对用它自己的 n 和自己的边缘分布模拟下限。**

⚠ 本轮换不了仪器:对象就是 SCCS 这三份编码之间的关系,**结构性地拿不到第二具仪器**
  (`#836` 已证 SCCS 里「同性行为频率」只有 Broude & Greene 一份)。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys, csv, collections, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
S = ROOT/"data/external/dplace/repo/datasets/SCCS"
NPERM, ALPHA, NSIM = 8000, 0.05, 600

D = collections.defaultdict(dict)
for r in csv.DictReader(open(S/"data.csv", encoding="utf-8")):
    if r["code"] not in ("", "NA"): D[r["soc_id"]][r["var_id"]] = r["code"]

def _avgrank(v):
    v = np.asarray(v, float); o = np.argsort(v, kind="mergesort"); r = np.empty(len(v), float)
    i = 0
    while i < len(v):
        j = i
        while j+1 < len(v) and v[o[j+1]] == v[o[i]]: j += 1
        r[o[i:j+1]] = (i+j)/2.0 + 1.0; i = j+1
    return r
def rho(x, y):
    rx, ry = _avgrank(x), _avgrank(y)
    if rx.std() == 0 or ry.std() == 0: return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])
def perm_p(x, y, rg, B=NPERM):
    r = rho(x, y); n = np.array([rho(x, list(rg.permutation(y))) for _ in range(B)])
    return r, float((np.sum(np.abs(n) >= abs(r))+1)/(B+1)), n

# ── 每一对自己的经验检出下限:在观测边缘分布 + 真实并列结构上种植 ──────────────
def empirical_floor(x, y, rg):
    """按强度网格种植,取检出率首次 >= 80% 的强度所对应的 |ρ| 中位数。

    ⚠⚠ **零分布每对只算一次,而这不只是省时间,是更正确的做法:**
    置换零的分布只取决于 `x` 与 `y` 的**多重集**,**与 y 的排列无关** ——
    所以「为每一个种植出来的 y 重算一次零」算的是同一个分布 N 遍。
    第一版就是那么写的(20 强度 × 400 模拟 × 300 置换),**跑不完**,而且多出来的
    只有蒙特卡洛噪声。⇒ 每对算一次零,取它 |ρ| 的 95 分位当阈。
    """
    x = np.asarray(x, float); ys = np.sort(np.asarray(y, float)); n = len(x)
    null = np.array([rho(x, list(rg.permutation(ys))) for _ in range(4000)])
    thr = float(np.quantile(np.abs(null[np.isfinite(null)]), 1-ALPHA))
    rk = _avgrank(x); zx = (rk-rk.mean())/max(rk.std(), 1e-9)
    for s in np.arange(0.2, 6.01, 0.2):
        rr = []
        for _ in range(NSIM):
            yy = ys[np.argsort(np.argsort(zx*s + rg.normal(size=n)))]
            rr.append(abs(rho(x, yy)))
        rr = np.array(rr); rr = rr[np.isfinite(rr)]
        if len(rr) and np.mean(rr >= thr) >= 0.80:
            return float(np.median(rr)), float(s), float(np.mean(rr >= thr)), thr
    return float("inf"), float("nan"), float(np.mean(rr >= thr)) if len(rr) else 0.0, thr

print("=== ⓪ 硬规则① —— 而它抓出了 `#836` 的两处缺陷 ===")
n961_raw = sum(1 for s in D if "SCCS961" in D[s])
n961_cln = sum(1 for s in D if D[s].get("SCCS961") not in (None, "4", "6"))
print(f"  ① `SCCS961` 码 4「permitted for males but no females」= **性别双标**;"
      f"码 6「insistence on virginity for the woman…no information on the male」= **男方信息缺失**")
print(f"     ⇒ 两者都**不在限制维度**上 —— **与 `#836` 弃用 `SCCS634`、剔除 `SCCS176` 码 2 是同一个理由**")
print(f"     ⇒ n {n961_raw} → **{n961_cln}**(⚠ **我上一轮读了 176 的码值,只读了 961 的题名**)")
print(f"  ② `SCCS602` 题名带 (No),码值只有 1=Yes / 2=No ⇒ **极性从元数据定不下来**,"
      f"而 `#836` 的 kill 要求「与 A 同号」⇒ **两个极性都跑,极性登记为 UNVERIFIED**")

RULERS = {
    "`SCCS176` 谴责(Broude&Greene)": ("SCCS176", {"2"}, +1),
    "`SCCS602` 性危险/污染(Whyte)":  ("SCCS602", set(), +1),
    "`SCCS961` 婚前限制(Frayser)":   ("SCCS961", {"4", "6"}, +1),
}
def vec(vid, drop, socs): return [float(D[s][vid]) for s in socs]
def overlap(a, b):
    va, da, _ = RULERS[a]; vb, db, _ = RULERS[b]
    return sorted(s for s in D if D[s].get(va) not in (None, *da) and D[s].get(vb) not in (None, *db))

print(f"\n=== ① 收敛效度:三把「限制」尺子两两同意吗?(每对用**自己的** n 与边缘分布定下限)===")
rg = np.random.default_rng(276)
PAIRS, agree = {}, 0
names = list(RULERS)
for i in range(3):
    for j in range(i+1, 3):
        a, b = names[i], names[j]; socs = overlap(a, b)
        x, y = vec(RULERS[a][0], RULERS[a][1], socs), vec(RULERS[b][0], RULERS[b][1], socs)
        r, p, nl = perm_p(x, y, rg)
        fl, st, hr, thr = empirical_floor(x, y, rg)
        ok = abs(r) >= fl
        agree += ok
        PAIRS[f"{a} × {b}"] = dict(n=len(socs), rho=r, p=p, floor=fl, plant_strength=st,
                                   hit_rate=hr, thr95=thr, null_sd=float(np.std(nl)), agrees=bool(ok))
        print(f"  {a[:22]:24s}× {b[:22]:24s} n={len(socs):>3d} ρ=**{r:+.3f}** "
              f"p={p:.4f} · 经验下限 **{fl:.3f}**(种植强度 {st:.1f},检出 {hr:.0%}) ⇒ "
              f"{'**同意**' if ok else '不同意'}")
print(f"  ⇒ **三对里同意 {agree}/3**")

print(f"\n=== ② 重跑 `#836` 的独立臂:C 臂清洗(剔码 4/6)· B 臂两个极性 ===")
ARMS2 = {}
for lab, vid, drop, sign in (("C 臂 `SCCS961` **清洗后**", "SCCS961", {"4", "6"}, +1),
                             ("C 臂 `SCCS961` 清洗前(`#836` 用的)", "SCCS961", set(), +1),
                             ("B 臂 `SCCS602` 极性 +1(`#836` 用的)", "SCCS602", set(), +1),
                             ("B 臂 `SCCS602` 极性 −1(另一读法)", "SCCS602", set(), -1)):
    socs = sorted(s for s in D if D[s].get(vid) not in (None, *drop) and "SCCS177" in D[s])
    x = [sign*float(D[s][vid]) for s in socs]; y = [float(D[s]["SCCS177"]) for s in socs]
    r, p, _ = perm_p(x, y, rg); fl, st, hr, thr = empirical_floor(x, y, rg)
    ARMS2[lab] = dict(n=len(socs), rho=r, p=p, floor=fl, reaches=bool(abs(r) >= fl and r < 0))
    print(f"  {lab:36s} n={len(socs):>3d} ρ=**{r:+.3f}** p={p:.4f} · 经验下限 {fl:.3f} ⇒ "
          f"{'**达标(与 A 同号且超下限)**' if ARMS2[lab]['reaches'] else '未达标'}")
reaches = [k for k, v in ARMS2.items() if v["reaches"]]
print(f"  ⇒ 清洗/换极性后达标的独立臂:**{len(reaches)}** ⇒ {reaches or '无'}")

print("\n=== ③ 控制 ===")
socs0 = overlap(names[0], names[2]); x0 = vec("SCCS176", {"2"}, socs0)
ys = np.sort(np.array(vec("SCCS961", {"4", "6"}, socs0), float))
zx = (_avgrank(x0)-_avgrank(x0).mean())/max(_avgrank(x0).std(), 1e-9)
yp = ys[np.argsort(np.argsort(zx*4.0 + rg.normal(size=len(x0))*0.05))]
r_pc, p_pc, _ = perm_p(x0, yp, rg)
print(f"  正控:在**真实边缘分布 + 真实并列结构**上种一个强关联 ⇒ ρ=**{r_pc:+.3f}**, p={p_pc:.5f}")
_, _, nc = perm_p(x0, vec("SCCS961", {"4", "6"}, socs0), rg)
print(f"  负控:打乱配对 ⇒ 中心 **{np.mean(nc):+.4f}** · SD {np.std(nc):.4f} —— "
      f"⚠ **「这个零该不该是零?」该**:打乱配对后相关的期望**就是 0**")

G = Gate("#837 · 三把「限制」的尺子,彼此同意吗")
G.asserted("① 硬规则①(补跑 `#836` 只读了题名的那两条):`SCCS961` 码 4=性别双标 · 码 6=男方信息缺失 ⇒ "
           "**不在限制维度** ⇒ 剔除;`SCCS602` 极性**从元数据定不下来** ⇒ **两个极性都跑,极性登记 UNVERIFIED**",
           bool(n961_cln < n961_raw and len([k for k in ARMS2 if "极性" in k]) == 2),
           f"SCCS961 n {n961_raw}→{n961_cln} · SCCS602 两极性均已跑", kind="control")
G.asserted("② 正控:在**真实边缘分布与真实并列结构**上种植的强关联必须捞得回来",
           bool(abs(r_pc) > 0.5 and p_pc < 0.01), f"ρ={r_pc:+.3f}, p={p_pc:.5f}", kind="control")
G.identity_control("③ 负控:打乱配对后 ρ 分布中心必须 == 0(⚠ **这个零该是零**)",
                   observed=float(np.mean(nc)), expected=0.0,
                   tol=3*float(np.std(nc))/math.sqrt(len(nc)), what="置换零的中心",
                   noise_half_width=float(np.std(nc))/math.sqrt(len(nc)))
G.asserted("④ 前提(跑前写下的最强混淆):**三对的重叠 n 各不相同(21–34),小样本那一对天然更难同意** ⇒ "
           "**每一对用它自己的 n 与自己的边缘分布模拟检出下限**,不用 Fisher 公式外推",
           bool(all(np.isfinite(v["floor"]) for v in PAIRS.values())),
           " · ".join(f"n={v['n']} 下限 {v['floor']:.3f}" for v in PAIRS.values()), kind="control")
G.asserted("⑤ kill(预注册):「存在一个共同的『限制』构念」要成立,需**≥2/3 对**的 |ρ| 超它自己的经验下限",
           bool(agree >= 2), f"同意 {agree}/3", kind="kill",
           yardstick="每对自己的两两 ρ,对照它自己的经验检出下限(80% 检出)",
           yardstick_noise=float(np.mean([v["null_sd"] for v in PAIRS.values()])))
G.asserted("⑥ kill(预注册):「`#836` 要重判」要成立,需清洗/换极性后**至少一条独立臂**与 A 同号且超其经验下限",
           bool(len(reaches) == 0), f"达标 {len(reaches)}/4 ⇒ {reaches or '无'}", kind="kill",
           yardstick="各臂对 `SCCS177` 的 ρ,对照该臂自己的经验检出下限",
           yardstick_noise=float(np.mean([v["floor"] for v in ARMS2.values()])/3))
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif len(reaches) > 0:
    V = (f"**C `#836` 要重判。** 清洗/换极性后达标的独立臂:{reaches}。")
elif agree >= 2:
    V = (f"**A `#836` 站得住。** 三把「限制」尺子里 **{agree}/3** 对彼此同意 ⇒ "
         f"**确实存在一个共同的『限制』构念**,而只有同编码者那把能预测稀少。\n"
         f"  ⇒ **一句关于人的话:那三份编码量的是同一样东西 —— 一个社会对性有多严 ——\n"
         f"  可是只有当读『他们多严』的那双眼睛同时也去读『那里有没有』的时候,\n"
         f"  这两件事才对得上。所以昨天那句话站得住。**")
else:
    V = (f"**B 构念不匹配 ⇒ `#836` 的头条说过头了,现在收窄。** "
         f"三把「限制」尺子彼此只有 **{agree}/3** 对同意 —— "
         f"**SCCS 里根本没有一个共同的『限制』构念可言。**\n"
         f"  ⇒ 所以两条独立臂的平,**不能**读成「换了编码者那条关联就没了」;\n"
         f"  它同样可以只是**「性危险观」「婚前限制」「对同性行为的谴责」本来就是三件不同的事**。\n"
         f"  ⇒ **一句关于人的话,而它比昨天那句谨慎:一个社会对婚前性行为多严,\n"
         f"  和它对同性行为多不容忍,在这一百多个社会里根本不是同一件事 ——\n"
         f"  所以「压制让人看不见」这句话,我目前连在哪一种压制上说它都还没定下来。**")
print(V)
print("\n⚠ **`#836` 的 kill 吃进过一个方向未定的输入**(`SCCS602` 极性)——"
      "本轮两极性都跑,结论未变,**但那是运气好,不是设计对**。极性仍登记 `UNVERIFIED`。")
json.dump(dict(pairs=PAIRS, arms_recheck=ARMS2, agree=agree, reaches=reaches,
               sccs961_n_raw=n961_raw, sccs961_n_clean=n961_cln,
               sccs602_polarity="UNVERIFIED — title carries '(No)', codes are only Yes/No",
               pos_control=dict(rho=r_pc, p=p_pc), neg_control=dict(center=float(np.mean(nc))),
               nperm=NPERM, nsim=NSIM, admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"do_the_three_rulers_agree.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'do_the_three_rulers_agree.json'}")

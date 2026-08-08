"""#810 · E03·A59·R249 —— 八轮一具仪器;而估计量早就换了,仪器搜索却没重跑过

⚠⚠ **本轮的起点是一个关于我自己的观察,不是关于世界的**:
   `#802`–`#809` 连着八轮全在 GSS 上。**而每一轮我都写「本轮换不了仪器」,依据是 `R223`。**
   **可 `R223` 关的是 `r = 斜率之比` 那个估计量的门** —— 它要求 **≥8 个时间点**,
   因为 `r_forced` 要标定一个五十年的潜在位移。
   **而 `#802` 之后对象已经换成了 `Δgap = gap1 − gap0`,它只需要两个时间点。**
   ⇒ **一句「换不了仪器」被从一个估计量搬到了另一个估计量,搬运过程没人检查。**
   **这正是 `#777` 那条「墙从没被查过」:一条限制的有效期,不长于它当初被测量时的那个问题。**

⚠ 而放松规格本身是危险的(**把规格放松到有东西能通过,就是制造第二具仪器**)。
  ⇒ **放松必须由估计量推出来,不是由我想通过推出来**:`Δgap` 字面上就是两个时间点的差。
  **而「题目必须是性道德题」这一条我不放松** —— 换了话题就换了主张。

**找到的**:NSFG 女性问卷两波(2011-2013 · 2017-2019)同时有
   ① `SAMESEX`「两个同性成年人之间的性关系是可以的」(5 档)——**与 GSS `homosex` 同一个话题**
   ② `ATTNDNOW`(现在多久参加一次宗教活动,7 档,**两波都无缺失**)
   ⇒ **三条规格全部满足,而 `R223` 当时按 ≥8 波把 NSFG 判出局是对的 —— 对那个估计量。**

G1 估计量:**同一个 `Δgap`,在两具仪器上、同一个话题、同一段时间(2011→2019)。**

⚠⚠ **三处必须匹配,否则「两具仪器一致」是假的一致:**
   ① **分层变量**:GSS 一直用 `attend+reliten+fund` 三项合成,而 NSFG 的 `RELDLIFE`/`FUNDAM1`
     两波各有 20–31% 空白 ⇒ **NSFG 只能用 `ATTNDNOW` 一项。**
     ⇒ 控制:**GSS 也跑一遍「只用 `attend`」的规格**,两个规格都报 —— **分层定义的代价要看得见。**
   ② **总体**:NSFG 是 **15–49 岁女性** ⇒ 控制:**GSS 限制到 18–49 岁女性**。
     ⚠ 下界差 3 岁(GSS 不访问 18 岁以下)—— **这一条匹配不上,如实登记,不假装匹配了。**
   ③ **极性与档数**:GSS `homosex` 1=总是错 → 4=完全没错(**高 = 宽容**);
     NSFG `SAMESEX` 是同意量表(**低 = 宽容**)⇒ **翻成 `6 − v`,让两边都「高 = 宽容」。**
     ⚠ 而**极性判断本身是个假设**(`#789` 已证明极性翻转会翻转符号)
     ⇒ **正控:对齐之后,两具仪器的总体均值在 2011→2019 之间都必须变得更宽容。**
       **若 NSFG 对齐后反而变保守,是我的编码假设错了,不是发现。**
     ⚠ 档数 4 vs 5 ⇒ **同时报 `Δgap ÷ 跨度`,而结论只用符号,不跨仪器比大小。**

三个世界:
   A **复现**:两具仪器的 `Δgap` 同号 ⇒ **这条鸿沟在 2010 年代变宽(或变窄)不是 GSS 的性质。**
   B **不复现**:异号 ⇒ **`realstat §2.5`:设计分歧本身就是发现** —— 去找两者差在哪个假设上,
     **不取平均,不挑一个喜欢的。**
   C **NSFG 分辨不出**:区间含 0 且宽 ⇒ 两波、六年,**本来就可能没有功效** ——
     那是关于第二具仪器的真收获,如实登记。

预测矩阵:
   | 世界 | 现在 | 同号且都排除 0 | 异号 | NSFG 含 0 |
   | A 复现       | 0.40 | **0.85** | 0.03 | 0.25 |
   | B 不复现     | 0.20 | 0.05 | **0.90** | 0.15 |
   | C 没功效     | 0.40 | 0.10 | 0.07 | **0.60** |

预注册判词(条件式):
  if 极性正控开火(两具仪器对齐后总体都变宽容)
     and 每具仪器的正控开火(合成一个已知 `Δgap` 必须取回)
     and 负控开火(差距恒定的世界里 `Δgap` 区间含 0):
      两具仪器 `Δgap` 同号且都排除 0 -> A
      异号                          -> B(**报两个数,不取平均**)
      NSFG 含 0                      -> C(登记功效不足)
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**`SAMESEX` 的题号在两波之间从 `IH-3` 变成了 `IH-1`** ——
  **问卷位置变了,而位置会影响作答(顺序效应)。** 这不是我能控制的,
  ⇒ **登记为 NSFG 侧一个无法排除的混淆,并且它只会影响 NSFG,不影响 GSS 侧** ——
  **所以两具仪器一致时它更有说服力,不一致时它是第一个嫌疑人。**

⚠ 硬规则①:先打印两具仪器每个变量的 n、真正被问过的年份/波次、档数。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(249)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
ND = ROOT/"data/external/nsfg"
B = 3000

# ── NSFG ─────────────────────────────────────────────────────────────────────
RXD = re.compile(r'_column\((\d+)\)\s+(\w+)\s+(\w+)\s+%(\d+)f\s+"([^"]*)"')
def nsfg_spec(f):
    out = {}
    for line in (ND/"setup"/f).read_text(errors="replace").split("\n"):
        m = RXD.search(line)
        if m: out[m.group(3).upper()] = (int(m.group(1)), int(m.group(4)), m.group(5))
    return out
def nsfg_read(dct, dat, cols):
    S = nsfg_spec(dct); rows = []
    with open(ND/dat, errors="replace") as fh:
        for line in fh:
            r = {}
            for c in cols:
                p, w, _ = S[c]; v = line[p-1:p-1+w].strip()
                r[c] = float(v) if v.isdigit() else np.nan
            rows.append(r)
    return pd.DataFrame(rows), S

print("=== ⓪ 硬规则①:NSFG 两波,每个变量的 n / 档 / 题号 ===")
NS = {}
for tag, dct, dat in (("2011-2013", "2011_2013_FemRespSetup.dct", "2011_2013_FemRespData.dat"),
                      ("2017-2019", "2017_2019_FemRespSetup.dct", "2017_2019_FemRespData.dat")):
    df, S = nsfg_read(dct, dat, ["SAMESEX", "ATTNDNOW", "AGER"])
    df["SAMESEX"] = df.SAMESEX.where(df.SAMESEX.between(1, 5))      # 8/9 = 拒答/不知道,剔除
    df["ATTNDNOW"] = df.ATTNDNOW.where(df.ATTNDNOW.between(1, 7))
    df["PERM"] = 6 - df.SAMESEX                                      # ⚠ 翻极性 ⇒ 高 = 宽容
    NS[tag] = df
    print(f"  {tag}: 行 {len(df):,} · `SAMESEX` n={int(df.SAMESEX.notna().sum()):,}(档 5,题号「{S['SAMESEX'][2][:6]}」)"
          f" · `ATTNDNOW` n={int(df.ATTNDNOW.notna().sum()):,}(档 7,**无缺失**)"
          f" · 年龄 {int(df.AGER.min())}–{int(df.AGER.max())}")
print("  ⚠ **题号 2011-13 = `IH-3`,2017-19 = `IH-1` —— 问卷位置变了,顺序效应无法排除,只影响 NSFG 侧**")

# ── GSS,匹配到同一段时间 · 同一个人群 · 两种分层规格 ──────────────────────────
print("\n=== ⓪ 硬规则①:GSS 匹配窗口 ===")
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
g = pd.read_stata(gp, columns=["year", "age", "sex", "homosex", "attend", "reliten", "fund"],
                  convert_categoricals=False)
for c, (lo, hi) in (("homosex", (1, 4)), ("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)),
                    ("age", (18, 89)), ("sex", (1, 2))):
    g[c] = pd.to_numeric(g[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
g["reliten"] = -g["reliten"]; g["fund"] = -g["fund"]
W0 = [2012, 2014]; W1 = [2018, 2021, 2022]        # 匹配 NSFG 的两个窗口
G = g[(g.sex == 2) & g.age.between(18, 49) & g.year.isin(W0+W1)].dropna(subset=["homosex", "attend"]).copy()
G["win"] = np.where(G.year.isin(W0), 0, 1)
print(f"  GSS 女性 18–49 · 年份 {W0} vs {W1} · n={len(G):,} "
      f"(窗口 0: {int((G.win==0).sum()):,} · 窗口 1: {int((G.win==1).sum()):,})")
print("  ⚠ **年龄下界差 3 岁(NSFG 从 15 岁起,GSS 从 18 岁起)—— 匹配不上,如实登记**")

def terc(s):
    q = s.quantile([1/3, 2/3]).values
    return np.where(s > q[1], 2, np.where(s > q[0], 1, 0))     # ⚠ `#794`:右开,用 `>`

def dgap(v0, k0, v1, k1):
    """两波 × 两层的均值差之差。k=2 虔诚 · k=0 非虔诚。"""
    f = lambda v, k, t: float(v[k == t].mean())
    return (f(v1, k1, 2)-f(v1, k1, 0)) - (f(v0, k0, 2)-f(v0, k0, 0))

def boot(v0, k0, v1, k1, rng):
    r = lambda v, k: (lambda i: (v[i], k[i]))(rng.integers(0, len(v), len(v)))
    a0, b0 = r(v0, k0); a1, b1 = r(v1, k1)
    try: return dgap(a0, b0, a1, b1)
    except Exception: return np.nan

def arm(v0, s0, v1, s1, label, span):
    k0, k1 = terc(pd.Series(s0)), terc(pd.Series(s1))
    v0, v1 = np.asarray(v0, float), np.asarray(v1, float)
    pt = dgap(v0, k0, v1, k1)
    dr = np.array([boot(v0, k0, v1, k1, RNG) for _ in range(B)]); dr = dr[np.isfinite(dr)]
    lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
    m0, m1 = float(v0.mean()), float(v1.mean())
    print(f"  {label:34s} 总体均值 {m0:.3f} → {m1:.3f}(**{'更宽容' if m1 > m0 else '更保守'}**) · "
          f"Δgap = **{pt:+.4f}** [{lo:+.4f}, {hi:+.4f}] · ÷跨度 {pt/span:+.4f} · 排除 0:"
          f"**{'是' if (lo > 0 or hi < 0) else '否'}**")
    return dict(label=label, dgap=float(pt), lo=lo, hi=hi, span=span, per_span=float(pt)/span,
                mean0=m0, mean1=m1, more_permissive=bool(m1 > m0), excl0=bool(lo > 0 or hi < 0))

print(f"\n=== ① 同一个 `Δgap`,两具仪器,同一个话题,2011→2019(B={B})===")
n0 = NS["2011-2013"].dropna(subset=["PERM", "ATTNDNOW"]); n1 = NS["2017-2019"].dropna(subset=["PERM", "ATTNDNOW"])
A_NSFG = arm(n0.PERM.values, n0.ATTNDNOW.values, n1.PERM.values, n1.ATTNDNOW.values,
             "NSFG `SAMESEX`(分层=`ATTNDNOW`)", 4)
G0, G1 = G[G.win == 0], G[G.win == 1]
A_GSS_att = arm(G0.homosex.values, G0.attend.values, G1.homosex.values, G1.attend.values,
                "GSS `homosex`(分层=`attend` 一项,匹配)", 3)
zs = lambda d: ((d - d.mean())/d.std(ddof=1))
GC = G.dropna(subset=["reliten", "fund"]).copy()
GC["REL"] = zs(GC[["attend", "reliten", "fund"]]).mean(axis=1)
C0, C1 = GC[GC.win == 0], GC[GC.win == 1]
A_GSS_3 = arm(C0.homosex.values, C0.REL.values, C1.homosex.values, C1.REL.values,
              "GSS `homosex`(分层=三项合成,原规格)", 3)
ARMS = [A_NSFG, A_GSS_att, A_GSS_3]

print("\n=== ② 控制 ===")
pol_ok = all(a["more_permissive"] for a in ARMS)
print(f"  ⚠ **极性正控**:对齐之后两具仪器的总体都必须变得更宽容 ⇒ "
      + " · ".join(f"{a['label'][:14]}:{'✓' if a['more_permissive'] else '✗'}" for a in ARMS)
      + f" ⇒ **{'通过' if pol_ok else '不通过 —— 我的编码假设错了,不是发现'}**")
# ⚠⚠ 第一版的合成器**坏了,而两条控制当场把它抓住了 —— 这正是它们存在的理由。**
#    我用 `k`(第 0 波的层标签)去造 `v1`,却 `return` 了一组**全新随机**的 `k1` ⇒
#    第 1 波的值与层标签被打散 ⇒ `gap1 ≈ 0` 而 `gap0 = −0.5` ⇒ **Δgap 恒等于 +0.5,
#    与我植入的 `dg` 完全无关**。实测正控 +0.4995、负控 +0.4574 —— **两条都停在 +0.5 附近,
#    而那正是「层标签被打散」的指纹,不是「植入量没取回」的指纹。**
#    ⇒ 修法:两波各自的标签各自生成,**而值必须由自己那一波的标签造**;
#      并把增量挂在 `k==2`(虔诚层)上,使 `gap = b`、`Δgap = b1 − b0 = dg` **直接对应植入量**。
# ⚠ 而这条合成臂**不经过 `terc`**(它直接给层标签)⇒ **它检验 `dgap`/`boot`,不检验分层那一段。**
#   如实说:**这两条控制没有覆盖整条路径。**
def synth(dg, n=4000, const=False):
    k0 = RNG.integers(0, 3, n); k1 = RNG.integers(0, 3, n)
    v0 = 2.0 + 0.50*(k0 == 2) + RNG.normal(0, 1, n)
    v1 = 2.0 + (0.50 + (0.0 if const else dg))*(k1 == 2) + RNG.normal(0, 1, n)
    return v0, k0.astype(float), v1, k1.astype(float)
def syn_arm(dg, const=False):
    v0, k0, v1, k1 = synth(dg, const=const)
    pt = dgap(v0, k0, v1, k1)
    dr = np.array([boot(v0, k0, v1, k1, RNG) for _ in range(1200)]); dr = dr[np.isfinite(dr)]
    return float(pt), float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
pc = syn_arm(-0.40); nc = syn_arm(0.0, const=True)
print(f"  正控:合成一个已知 Δgap = −0.400 ⇒ 取回 **{pc[0]:+.4f} [{pc[1]:+.4f}, {pc[2]:+.4f}]**,"
      f"排除 0:**{'是' if pc[2] < 0 else '否'}**")
print(f"  负控:差距按构造恒定 ⇒ **{nc[0]:+.4f} [{nc[1]:+.4f}, {nc[2]:+.4f}]**,"
      f"含 0:**{'是' if nc[1] <= 0 <= nc[2] else '否'}**(⚠ 这一次参照真的是 0)")

Gg = Gate("#810 · 同一个问题换一份调查,还在不在")
Gg.asserted("① 极性正控:两具仪器对齐之后,总体均值在 2011→2019 之间都必须变得更宽容"
            "(⚠ 若 NSFG 反而变保守,是我的编码假设错了,不是发现 —— `#789` 已证极性翻符号)",
            pol_ok, " · ".join(f"{a['label'][:16]} {a['mean0']:.3f}→{a['mean1']:.3f}" for a in ARMS),
            kind="control")
Gg.asserted("② 正控:合成一个已知 Δgap = −0.400 必须取回且区间排除 0",
            bool(pc[2] < 0 and abs(pc[0]+0.40) < 0.15), f"取回 {pc[0]:+.4f} [{pc[1]:+.4f}, {pc[2]:+.4f}]",
            kind="control")
Gg.asserted("③ 负控:差距按构造恒定的世界里,Δgap 区间必须**含 0**(⚠ 参照真的是 0)",
            bool(nc[1] <= 0 <= nc[2]), f"{nc[0]:+.4f} [{nc[1]:+.4f}, {nc[2]:+.4f}]", kind="control")
Gg.asserted("④ 前提(跑前写下的三处匹配):分层两个规格都跑 · 总体限到女性 18–49 · 极性对齐后同时报 ÷跨度",
            bool(len(ARMS) == 3 and all("per_span" in a for a in ARMS)),
            f"三条臂 · 跨度 {[a['span'] for a in ARMS]} · GSS n={len(G):,}", kind="control")
Gg.asserted("⑤ kill(预注册):「换一份调查还在」要成立,需两具仪器的 `Δgap` **同号且都排除 0**",
            bool(np.sign(A_NSFG["dgap"]) == np.sign(A_GSS_att["dgap"])
                 and A_NSFG["excl0"] and A_GSS_att["excl0"]),
            f"NSFG {A_NSFG['dgap']:+.4f}(排0={A_NSFG['excl0']}) · "
            f"GSS-attend {A_GSS_att['dgap']:+.4f}(排0={A_GSS_att['excl0']})", kind="kill")
print(); print(Gg)
adm = Gg.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*98)
same = np.sign(A_NSFG["dgap"]) == np.sign(A_GSS_att["dgap"])
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif same and A_NSFG["excl0"] and A_GSS_att["excl0"]:
    V = (f"**A 复现 —— 而这是这个项目第一次在第二具仪器上问同一个问题。**\n"
         f"  NSFG(女性 15–49,`SAMESEX`)Δgap = **{A_NSFG['dgap']:+.4f}** [{A_NSFG['lo']:+.4f}, {A_NSFG['hi']:+.4f}]\n"
         f"  GSS(女性 18–49,`homosex`,同一段时间)Δgap = **{A_GSS_att['dgap']:+.4f}** "
         f"[{A_GSS_att['lo']:+.4f}, {A_GSS_att['hi']:+.4f}]\n"
         f"  ⇒ **同号,且都排除 0 ⇒ 2010 年代这条宗教鸿沟的走向不是 GSS 的性质。**")
elif not same:
    V = (f"**B 不复现 —— 而 `realstat §2.5` 说得清楚:设计分歧本身就是发现,不取平均,不挑一个喜欢的。**\n"
         f"  NSFG **{A_NSFG['dgap']:+.4f}** [{A_NSFG['lo']:+.4f}, {A_NSFG['hi']:+.4f}] vs "
         f"GSS **{A_GSS_att['dgap']:+.4f}** [{A_GSS_att['lo']:+.4f}, {A_GSS_att['hi']:+.4f}] —— **异号。**\n"
         f"  ⇒ **去找两者差在哪个假设上:人群(15–49 女性 vs 18–49 女性)· 题目措辞 · 题号位置(`IH-3`→`IH-1`)·\n"
         f"  分层变量(一项 vs 三项)。第一个嫌疑人是题号位置,因为它只影响 NSFG 一侧。**")
else:
    V = (f"**C 第二具仪器没有分辨力。** NSFG Δgap = {A_NSFG['dgap']:+.4f} "
         f"[{A_NSFG['lo']:+.4f}, {A_NSFG['hi']:+.4f}] 含 0 —— **两波、六年,本来就可能没有功效。**\n"
         f"  ⇒ **这是关于第二具仪器的真收获,不是关于世界的判断。**\n"
         f"  ⚠ 而 GSS 侧在同一段时间上是 {A_GSS_att['dgap']:+.4f} "
         f"[{A_GSS_att['lo']:+.4f}, {A_GSS_att['hi']:+.4f}] —— **两者并排报,不合并。**")
print(V)
print("\n⚠ **无法排除的混淆,只影响 NSFG 一侧**:`SAMESEX` 的题号在两波间从 `IH-3` 变成 `IH-1`(顺序效应)。")
print("⚠ **匹配不上的一处,如实登记**:年龄下界 NSFG 15 岁 vs GSS 18 岁。")
json.dump(dict(arms=ARMS, gss_years=[W0, W1], gss_n=int(len(G)),
               nsfg_n={t: int(len(NS[t])) for t in NS},
               pos_control=dict(point=pc[0], lo=pc[1], hi=pc[2], planted=-0.40),
               neg_control=dict(point=nc[0], lo=nc[1], hi=nc[2], reference=0.0),
               polarity_control=pol_ok, same_sign=bool(same),
               admissible=adm, verdict=V, gate_ok=Gg.verdict(),
               confounds=["SAMESEX 题号 IH-3 -> IH-1(顺序效应,只影响 NSFG)",
                          "年龄下界 15 vs 18,匹配不上"]),
          open(OUT/"second_instrument.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'second_instrument.json'}")

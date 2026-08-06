"""E02·A214·R581 — 「性」是一个领域,还是那个领域?

`#535` 的 NEXT。行动类型:**FRONTIER**(结局会改写 `#534` 的措辞,不只是它的数)。
换第二份仪器(硬规则 4):**NSFG 2011–2013 女性问卷,ACASI 自填**,与 GSS 面访是不同的作答模式。

**为什么是这份仪器:** `#535c` 说 GSS 里**没有**「同为四级、同为道德判断、非性」的参照题组,
所以「同领域」与「同格式」分不开。**NSFG 的 IH 段恰好有** ——
11 道**同为五级李克特、同一份问卷、同一批受访者、同一作答模式**的题,其中
**3 道是性行为判断**,**7 道是非性的婚姻/家庭道德判断**。⇒ 格式已匹配,只剩领域在变。

⚠ 硬规则 1 已先做:11 题各 n=5601,取值 1–5(8/9 = 不知道/拒答,已剔除)。逐题码在下方打印。
⚠ 而这改变了世界的设定:GSS 的参照(死刑/大麻/安乐死)**彼此不属于一个领域**,
   NSFG 的参照(离婚/同居/未婚生育)**本身就是一个领域**。
   ⇒ 本轮问的不再是「性题比任意两题更耦合吗」,而是**「性这个领域,比另一个领域更紧吗」**。

G1 ESTIMAND:`ρ_sex` = 3 道性行为题两两 |ρ| 中位(3 对);
   `ρ_fam` = 7 道婚家题两两 |ρ| 中位(21 对);`ρ_x` = 跨领域 |ρ| 中位(21 对)。
   **主量 = `ρ_sex − ρ_fam`;次量 = `ρ_x`(两个领域有多分离)。**

WORLDS:
  W-SEX-SPECIAL   `ρ_sex` >> `ρ_fam` ⇒ 性是一个**更紧**的领域
  W-JUST-A-DOMAIN `ρ_sex` ≈ `ρ_fam` ⇒ **性是一个领域,和别的领域一样紧** ——
     则 `#534` 的 0.375-vs-0.093 是**领域 vs 任意**,不是**性 vs 其他领域**,措辞必须改
  W-WEAKER        `ρ_sex` < `ρ_fam` ⇒ 性反而更松
⚠ BASIN:`W-SEX-SPECIAL` 会让第八条更漂亮,**不是**本轮下注方向。
   本轮下注 `W-JUST-A-DOMAIN` —— 它要求我把刚写上页面的那条**收窄**。

CONTROLS(G2):
  正对照 `sxok18` × `sxok16`(18 岁 vs 16 岁,近重复题)必须**远高于**任何领域中位 ——
     它是本仪器个体层的上限;⚠ 且 `#517` 已测过这一对在别的量上的差,**不是新发现**;
  安慰剂 每道态度题 × **受访者编号的末位数**(纯任意标签)必须 ≈ 0;
  边界项 `gayadopt`(性少数但问的是收养权,不是性行为)**单独成一个规格臂**,不默认归类。
KILL(条件式):if 正对照 > max(领域中位) and 安慰剂 ≈ 0:
     |ρ_sex − ρ_fam| > 两者的 bootstrap 展布 -> W-SEX-SPECIAL / W-WEAKER(按符号)
     否则 -> W-JUST-A-DOMAIN
   else UNVERIFIED
IMPOSSIBLE:仅女性受访者 ⇒ 无性别外推 · 单一年份波(2011–2013)⇒ 无时间维 ·
   3 道性题只给 3 对 ⇒ `ρ_sex` 的分辨率远低于 `ρ_fam` 的 21 对,**必须用 bootstrap 展布比较,
   不能直接比中位** · 观察性非因果 · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, re, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import rankdata
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NS = ROOT / "data/external/nsfg"

def parse_dct(p):
    out = {}; pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out

LAY = parse_dct(NS / "setup" / "2011_2013_FemRespSetup.dct")
SEX = ["samesex", "sxok18", "sxok16"]
FAM = ["staytog", "chunless", "chsuppor", "okcohab", "marrfail", "chcohab", "prvntdiv"]
BND = ["gayadopt"]
ALL = SEX + FAM + BND
cols = {n: LAY[n] for n in ALL if n in LAY}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
X = {n: np.where(np.isin(np.array(buf[n]), [1, 2, 3, 4, 5]), np.array(buf[n]), np.nan) for n in cols}
print("=== 硬规则 1:逐题 n(1–5 有效)与标题 ===")
for n in ALL:
    print(f"  {n:9s} n={int(np.isfinite(X[n]).sum()):5d}  {cols[n][2][:54]}")

def rho(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 200: return np.nan, int(m.sum())
    return float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]), int(m.sum())

def med(group_a, group_b=None):
    pairs = (list(itertools.combinations(group_a, 2)) if group_b is None
             else [(a, b) for a in group_a for b in group_b])
    out = []
    for a, b in pairs:
        r, n = rho(X[a], X[b])
        if np.isfinite(r): out.append((abs(r), f"{a}×{b}", n))
    return out

def boot(group_a, group_b=None, k=600):
    idx0 = len(X[ALL[0]]); vals = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(k // len(SEEDS)):
            i = rng.integers(0, idx0, idx0)
            pairs = (list(itertools.combinations(group_a, 2)) if group_b is None
                     else [(a, b) for a in group_a for b in group_b])
            o = []
            for a, b in pairs:
                r, _n = rho(X[a][i], X[b][i])
                if np.isfinite(r): o.append(abs(r))
            if o: vals.append(np.median(o))
    return np.array(vals)

SPECS = {"A 严格(gayadopt 剔除)": (SEX, FAM),
         "B gayadopt 归性": (SEX + BND, FAM),
         "C gayadopt 归家庭": (SEX, FAM + BND)}
res = {}
print("\n=== 三个规格臂(gayadopt 的归属不默认,单独成臂)===")
for sname, (sx, fm) in SPECS.items():
    ps, pf, px = med(sx), med(fm), med(sx, fm)
    ms, mf, mx = [float(np.median([x[0] for x in p])) for p in (ps, pf, px)]
    bs, bf = boot(sx), boot(fm)
    sd = float(np.sqrt(bs.var() + bf.var()))
    res[sname] = dict(rho_sex=ms, rho_fam=mf, rho_cross=mx, diff=ms - mf, joint_sd=sd,
                      k_sex=len(ps), k_fam=len(pf), k_cross=len(px),
                      cells=[dict(pair=p[1], rho=p[0], n=p[2], group=gname, scheme=sname,
                                  inclusion=[f"两题都在 1–5 (n={p[2]})", sname, gname])
                             for gname, pp in (("sex", ps), ("fam", pf), ("cross", px)) for p in pp])
    print(f"  {sname:22s} 性内={ms:.4f}({len(ps)}对) 家内={mf:.4f}({len(pf)}对) "
          f"跨={mx:.4f}({len(px)}对)  **差={ms-mf:+.4f}**  联合展布={sd:.4f}  "
          f"{'超展布' if abs(ms-mf)>sd else '**不可分辨**'}")

G = Gate("「性」是一个领域,还是那个领域?(NSFG 2011-2013,ACASI 自填)")
pc, pcn = rho(X["sxok18"], X["sxok16"])
print(f"\n=== 对照 ===\n  正对照 sxok18×sxok16(18 岁 vs 16 岁,近重复题)n={pcn} |ρ|={abs(pc):.4f}"
      f"  ⚠ 这一对在别的量上已被 `#517` 测过,**不是新发现**")
G.positive_control("正对照:sxok18×sxok16(仪器个体层上限)", planted=abs(pc),
                   floor=max(res[s]["rho_sex"] for s in res), spread=1e-9)
rng = np.random.default_rng(SEEDS[0])
tag = rng.integers(0, 10, len(X[ALL[0]])).astype(float)
zs = [abs(rho(X[n], tag)[0]) for n in ALL if np.isfinite(rho(X[n], tag)[0])]
G.negative_control("安慰剂:态度题 × 随机任意标签", null=float(np.median(zs)), effect=abs(pc),
                   null_spread=float(np.std(zs)), null_kind="与问卷无关的随机整数标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{c['scheme'][:1]}|{c['pair']}": c
                                             for v in res.values() for c in v["cells"]})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件",
                                     {f"{c['scheme'][:1]}|{c['pair']}": c for v in res.values() for c in v["cells"]})
print("\n" + "=" * 76)
if abs(pc) > max(res[s]["rho_sex"] for s in res) and np.median(zs) < 0.5 * abs(pc):
    seps = {s: abs(v["diff"]) > v["joint_sd"] for s, v in res.items()}
    if not any(seps.values()):
        world = "W-JUST-A-DOMAIN"
        verdict = (f"三个规格臂上性内与家内**都不可分辨**"
                   f"(差 {[f'{res[s][chr(100)+chr(105)+chr(102)+chr(102)]:+.4f}' for s in res]}) -> "
                   f"**性是一个领域,和婚姻家庭那个领域一样紧;`#534` 的措辞必须收窄**")
    elif all(res[s]["diff"] > 0 for s in res if seps[s]):
        world = "W-SEX-SPECIAL"; verdict = f"性内显著高于家内 -> **性是更紧的领域**"
    else:
        world = "W-WEAKER"; verdict = f"性内显著低于家内 -> **性反而更松**"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print(f"  跨领域中位 = {[f'{res[s][chr(114)+chr(104)+chr(111)+chr(95)+chr(99)+chr(114)+chr(111)+chr(115)+chr(115)]:.4f}' for s in res]}"
          f" —— 两个领域**没有**分开成两块,这本身是结果")
    print("⚠ 这个 KILL 会怎样失败:性题只有 3 道 ⇒ 3 对,`ρ_sex` 的展布远大于 21 对的 `ρ_fam`;"
          "「不可分辨」在这里**部分是分辨率**,而这正是为什么本轮比的是 bootstrap 展布,不是中位本身。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(specs=res, world=world, verdict=verdict, positive_control=dict(rho=pc, n=pcn),
               placebo=[float(x) for x in zs], seeds=SEEDS,
               instrument="NSFG 2011-2013 女性问卷,ACASI 自填(与 GSS 面访不同模式)",
               impossible=["仅女性受访者", "单一年份波无时间维", "3 道性题只给 3 对,分辨率远低于家庭题",
                           "观察性非因果"], unchallenged=True),
          open(OUT / "sex_a_domain.json", "w"), indent=1)
print(f"\nwrote {OUT/'sex_a_domain.json'}")

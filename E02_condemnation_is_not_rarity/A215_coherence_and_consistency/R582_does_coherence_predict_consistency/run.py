"""E02·A215·R582 — 态度彼此更紧的人,态度与行为也更一致吗?

`#536` 的 NEXT。行动类型:**FRONTIER**。
`#536` 说性态度在**人群层**内部更紧(0.425 vs 0.140)。本轮问一个**个体层**的问题:
**那种「紧」是一个人的性质吗?** —— 如果是,态度彼此一致的人,应该也是言行更一致的人。

⚠ 硬规则 1 已先做,而它砍掉了一半候选变量:NSFG 的行为变量里
   `vry1stsx`(4867)· `sexonce`(4858)· `parts1yr`(4858)· `fsexpage`(4858)· `timesmar`(2453)
   **都条件于结局**(只问有过性行为/结过婚的人)—— `#492` 那一类。
   **全样本的只有:** `hadsex`(5601)· `lifprtnr`(5601)· `samesexany`(5575)·
   `oppsexany`(5589)· `evrmarry`(5601)· `cohever`(5601)。**本轮只用全样本变量。**
⚠ `samesexany` 码为 **1=是 / 5=否 / 7=其他**,**不是 {0,1}** —— 本轮对基础率下 assert(`#495b` 的守卫)。

G1 ESTIMAND(先于方法):
   个体的 **coherence** = −(她三道对齐后 z 分的标准差);**extremity** = |三道 z 分的均值|。
   对齐方式**预注册**:按第一主成分载荷符号翻转,不看结果调。
   **主量 = 态度→行为关联在 coherence 三分位上的差**
   (关联 = 该组内 `samesex` 态度与 `samesexany` 行为的 |lnOR|,以及 `sxok18` 与 `hadsex`)。

⚠ **写在跑之前的最强混淆:coherence 与 extremity 必然相关** ——
   三个答案全是 1 或全是 5 的人,标准差小**且**极端。
   若不控制,「一致的人言行更一致」可能只是「极端的人言行更一致」。
   ⇒ **同轮控制:在 extremity 三分位内部再看 coherence 的效应(2×3 网格全格公布)。**

WORLDS:
  W-TRAIT      控制 extremity 后,coherence 高的人关联更强 ⇒ **「紧」是一个人的性质**
  W-EXTREMITY  控制后效应消失 ⇒ **测到的是极端度,不是一致性**
  W-NEITHER    两者都无效应 ⇒ 人群层的「紧」不落到个体层的言行一致上
⚠ BASIN:`W-TRAIT` 会让 `#536` 顺理成章地延伸到个体,**不是**本轮下注方向。本轮下注 `W-EXTREMITY`。

CONTROLS(G2):
  正对照 `sxok18` 态度 → `hadsex` 行为(公认强关联)必须在全样本上显著,且**置换后消失**;
  安慰剂 coherence → **随机整数标签**的关联必须 ≈ 0;
  基础率 assert:`samesexany` 阳性率必须落在 [0.05, 0.30](`#495b` 的守卫);
  规格曲线 两条态度-行为线 × 3 个 coherence 三分位 × 3 个 extremity 三分位,**全格公布**。
KILL(条件式):if 正对照通过 and 置换后消失 and 安慰剂 ≈ 0 and 基础率合理:
     控制 extremity 后 coherence 效应仍超展布 -> W-TRAIT;消失 -> W-EXTREMITY else W-NEITHER
   else UNVERIFIED
IMPOSSIBLE:仅女性 · 单一波 · 横断面 ⇒ 态度与行为的**时序不可分**(行为多在态度之前发生)⇒
   **绝不可读成因果**,连方向都不可读 · 3 道题算 coherence ⇒ 该指标本身极不精确 · [unchallenged]
"""
import os, sys, pathlib, json, re, math, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
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
ATT = ["samesex", "sxok18", "sxok16"]
BEH = ["samesexany", "hadsex", "oppsexany", "cohever", "evrmarry"]
cols = {n: LAY[n] for n in ATT + BEH}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
A = np.column_stack([np.where(np.isin(np.array(buf[n]), [1, 2, 3, 4, 5]), buf[n], np.nan) for n in ATT])
sse = np.array(buf["samesexany"]); SSA = np.where(sse == 1, 1.0, np.where(sse == 5, 0.0, np.nan))
had = np.array(buf["hadsex"]); HAD = np.where(had == 1, 1.0, np.where(had == 2, 0.0, np.nan))
print(f"=== 规则①/#495b 守卫 ===")
br = float(np.nanmean(SSA)); print(f"  samesexany 阳性率 = {br:.4f}  (码 1=是/5=否/7=其他)")
assert 0.05 <= br <= 0.30, f"阳性率 {br:.4f} 越界 -> 码又读错了"
print(f"  hadsex 阳性率 = {np.nanmean(HAD):.4f} · 三道态度题齐全的人 n = {int(np.isfinite(A).all(1).sum())}")

ok = np.isfinite(A).all(1)
Z = (A[ok] - np.nanmean(A[ok], 0)) / np.nanstd(A[ok], 0)
# 预注册的对齐:第一主成分载荷符号
u, s, vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
sign = np.sign(vt[0]); sign[sign == 0] = 1
Za = Z * sign
COH = -Za.std(1); EXT = np.abs(Za.mean(1))
print(f"  对齐符号(第一主成分)= {sign}  · corr(coherence, extremity) = "
      f"{np.corrcoef(COH, EXT)[0,1]:+.4f}  <- 这就是那个混淆,已量出")

def lnor(att, beh, mask):
    a = att[mask]; b = beh[mask]
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 150: return np.nan, int(m.sum())
    hi = a[m] >= np.median(a[m])
    if hi.sum() < 50 or (~hi).sum() < 50: return np.nan, int(m.sum())
    p1, p0 = b[m][hi].mean(), b[m][~hi].mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, int(m.sum())
    return float(math.log((p1 / (1 - p1)) / (p0 / (1 - p0)))), int(m.sum())

SS = SSA[ok]; HD = HAD[ok]
LINES = {"samesex→samesexany": (Za[:, 0], SS), "sxok18→hadsex": (Za[:, 1], HD)}
ct = np.quantile(COH, [1/3, 2/3]); et = np.quantile(EXT, [1/3, 2/3])
cb = np.digitize(COH, ct); eb = np.digitize(EXT, et)
rows = []
print("\n=== 2×3×3 全格:态度→行为 |lnOR|,按 coherence × extremity 三分位 ===")
for lname, (att, beh) in LINES.items():
    for c in range(3):
        for e in range(3):
            m = (cb == c) & (eb == e)
            v, n = lnor(att, beh, m)
            rows.append(dict(line=lname, coh_tercile=c, ext_tercile=e, lnor=v if np.isfinite(v) else None,
                             n=n, inclusion=[f"coherence 三分位 {c}", f"extremity 三分位 {e}",
                                             f"该格 n={n}", "两列都非缺失", "每臂 >=50"]))
        vs = [r["lnor"] for r in rows if r["line"] == lname and r["coh_tercile"] == c and r["lnor"] is not None]
        print(f"  {lname:22s} coh={c}  各 ext 格 |lnOR| = "
              f"{[f'{abs(x):.3f}' for x in vs] or '全部不可算'}")

def coh_effect(lname):
    """控制 extremity:在每个 ext 格内比 coh 高低,再平均。"""
    d = []
    for e in range(3):
        hi = [r["lnor"] for r in rows if r["line"] == lname and r["ext_tercile"] == e
              and r["coh_tercile"] == 2 and r["lnor"] is not None]
        lo = [r["lnor"] for r in rows if r["line"] == lname and r["ext_tercile"] == e
              and r["coh_tercile"] == 0 and r["lnor"] is not None]
        if hi and lo: d.append(abs(hi[0]) - abs(lo[0]))
    return (float(np.mean(d)), len(d)) if d else (np.nan, 0)

G = Gate("态度彼此更紧的人,态度与行为也更一致吗?(NSFG,控制极端度)")
pcv, pcn = lnor(Za[:, 1], HD, np.ones(len(HD), bool))
rng = np.random.default_rng(SEEDS[0])
perm = [abs(lnor(Za[:, 1], HD[rng.permutation(len(HD))], np.ones(len(HD), bool))[0]) for _ in range(200)]
print(f"\n=== 对照 ===\n  正对照 sxok18→hadsex 全样本 |lnOR|={abs(pcv):.4f} n={pcn} · 置换 q95={np.quantile(perm,.95):.4f}")
G.positive_control("正对照:sxok18→hadsex", planted=abs(pcv), floor=float(np.quantile(perm, .95)), spread=1e-9)
G.negative_control("g=0:置换行为后必须消失", null=float(np.median(perm)), effect=abs(pcv),
                   null_spread=float(np.std(perm)), null_kind="个体层行为标签置换")
tagr = rng.integers(0, 2, len(COH)).astype(float)
G.negative_control("安慰剂:coherence → 随机标签",
                   null=abs(lnor(COH, tagr, np.ones(len(COH), bool))[0]), effect=abs(pcv),
                   null_spread=1e-9, null_kind="与问卷无关的随机二值标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['line'][:8]}|c{r['coh_tercile']}e{r['ext_tercile']}": r for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件",
                                     {f"{r['line'][:8]}|c{r['coh_tercile']}e{r['ext_tercile']}": r for r in rows})
effs = {l: coh_effect(l) for l in LINES}
print(f"\n  控制 extremity 后的 coherence 效应:{ {l: (f'{v:+.4f}' if np.isfinite(v) else 'NA', k) for l,(v,k) in effs.items()} }")
print("\n" + "=" * 76)
if abs(pcv) > np.quantile(perm, .95) and np.median(perm) < 0.5 * abs(pcv):
    vals = [v for v, k in effs.values() if np.isfinite(v)]
    spread = float(np.std([abs(r["lnor"]) for r in rows if r["lnor"] is not None]))
    if vals and all(abs(v) > spread for v in vals) and all(v > 0 for v in vals):
        world = "W-TRAIT"; verdict = f"控制极端度后 coherence 效应 {vals} 仍超全格展布 {spread:.4f} -> **「紧」是一个人的性质**"
    elif vals and all(abs(v) <= spread for v in vals):
        world = "W-EXTREMITY"; verdict = (f"控制极端度后 coherence 效应 {[f'{v:+.4f}' for v in vals]} "
            f"落在全格展布 {spread:.4f} 内 -> **测到的是极端度,不是一致性**")
    else:
        world = "W-NEITHER"; verdict = f"效应符号不一致或不可算 -> **人群层的「紧」不落到个体的言行一致上**"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print(f"  ⚠ corr(coherence, extremity) = {np.corrcoef(COH,EXT)[0,1]:+.4f} —— 混淆是真实的,已同轮控制")
    print("⚠ 这个 KILL 会怎样失败:横断面 ⇒ 行为多在态度之前发生,**时序不可分**;"
          "而且 coherence 由三道题算出,该指标本身极不精确,**低估**任何真实的调节效应。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rows=rows, coh_effects={l: (v if np.isfinite(v) else None, k) for l, (v, k) in effs.items()},
               corr_coh_ext=float(np.corrcoef(COH, EXT)[0, 1]), world=world, verdict=verdict,
               base_rates=dict(samesexany=br, hadsex=float(np.nanmean(HAD))), align_sign=[float(x) for x in sign],
               positive_control=dict(lnor=pcv, n=pcn), seeds=SEEDS,
               instrument="NSFG 2011-2013 女性 ACASI,横断面",
               impossible=["仅女性", "单一波", "横断面时序不可分,不可读成因果",
                           "3 道题算 coherence,指标极不精确,低估调节"], unchallenged=True),
          open(OUT / "coherence_consistency.json", "w"), indent=1)
print(f"\nwrote {OUT/'coherence_consistency.json'}")

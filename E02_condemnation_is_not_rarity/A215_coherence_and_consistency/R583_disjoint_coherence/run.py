"""E02·A215·R583 — 用不相交的题算一致性,问题问得出来了吗?

`#537` 的 NEXT。行动类型:**FRONTIER**。
**修的是 `#537b` 的定义冲突**:上一轮用**性题自身**的内部离散当 coherence,
又用**同一批性题**当预测变量 —— 条件在高 coherence 上就抽干了预测变量的方差。
**禁令(`#537d`):算 coherence 的题,必须与作预测变量的题不相交。**
⇒ 本轮 coherence 由 **7 道家庭题**算,预测变量仍是性题。

**同时修一个逻辑错(`#537c` = `#533a` 的同型,相隔四条):**
上一轮的三分写着「符号不一致**或不可算** -> W-NEITHER」,于是把「算不出来」打印成了一个世界。
**本轮改成:任一格不可算 -> 该行直接 `UNVERIFIED`,并在三分之前就打印哪些格不可算。**

G1 ESTIMAND:`coherence_fam` = −(7 道家庭题对齐后 z 分的标准差);`extremity_fam` = |均值|。
   主量 = 性态度→性行为的 |lnOR| 在 `coherence_fam` 三分位上的差,**在 `extremity_fam` 格内比较**。
   ⚠ 这问的已经不是「她的性态度一致吗」,而是**「她在别处的道德一致性,预测她这里的言行一致吗」** ——
     **估计量变了,必须承认,不能当成上一轮那个问题的更好版本。**

WORLDS:
  W-GENERAL   家庭题的一致性预测性领域的言行一致 ⇒ **道德一致性是一个跨领域的人格性质**
  W-LOCAL     无效应 ⇒ 一致性是**领域内**的,不外溢
  W-UNCOMPUTABLE 仍有格算不出 ⇒ `UNVERIFIED`,并按 `#111c` 改方向(这是这条线的第二次)
⚠ BASIN:`W-GENERAL` 是更漂亮的故事,**不是**下注方向。本轮下注 `W-LOCAL`。
CONTROLS:正对照 `sxok18`→`hadsex` 全样本 + 置换消失 · 安慰剂 coherence_fam → 随机标签 ·
   基础率 assert(`#495b`)· 全格公布,含不可算的格
KILL(条件式):if 正对照通过 and 置换消失 and 安慰剂≈0 and **全部格可算**: 按三分判
   else UNVERIFIED(**不可算就是 UNVERIFIED,不再跳过**)
IMPOSSIBLE:横断面时序不可分 · 仅女性 · 单一波 · [unchallenged]
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
FAM = ["staytog", "chunless", "chsuppor", "okcohab", "marrfail", "chcohab", "prvntdiv"]
SEXQ = ["samesex", "sxok18"]; BEH = ["samesexany", "hadsex"]
cols = {n: LAY[n] for n in FAM + SEXQ + BEH}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip(); buf[n].append(float(v) if v not in ("", ".") else np.nan)
F = np.column_stack([np.where(np.isin(np.array(buf[n]), [1, 2, 3, 4, 5]), buf[n], np.nan) for n in FAM])
S = np.column_stack([np.where(np.isin(np.array(buf[n]), [1, 2, 3, 4, 5]), buf[n], np.nan) for n in SEXQ])
sse = np.array(buf["samesexany"]); SSA = np.where(sse == 1, 1.0, np.where(sse == 5, 0.0, np.nan))
had = np.array(buf["hadsex"]); HAD = np.where(had == 1, 1.0, np.where(had == 2, 0.0, np.nan))
br = float(np.nanmean(SSA)); assert 0.05 <= br <= 0.30, f"基础率 {br:.4f} 越界"
ok = np.isfinite(F).all(1) & np.isfinite(S).all(1)
print(f"=== 规则① ===\n  samesexany 阳性率={br:.4f} · 七题+两题齐全 n={int(ok.sum())}")
Z = (F[ok] - F[ok].mean(0)) / F[ok].std(0)
u, s, vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
sign = np.sign(vt[0]); sign[sign == 0] = 1
Za = Z * sign
COH = -Za.std(1); EXT = np.abs(Za.mean(1))
SX = S[ok]; SS = SSA[ok]; HD = HAD[ok]
print(f"  对齐符号={sign.astype(int)} · corr(coh_fam, ext_fam)={np.corrcoef(COH,EXT)[0,1]:+.4f}"
      f" · corr(coh_fam, 性题离散)={np.corrcoef(COH, -np.nanstd((SX-np.nanmean(SX,0))/np.nanstd(SX,0),1))[0,1]:+.4f}"
      f"  <- 不相交是否真的不相交,量出来")
def lnor(att, beh, m):
    a, b = att[m], beh[m]
    k = np.isfinite(a) & np.isfinite(b)
    if k.sum() < 150: return None, int(k.sum()), "n<150"
    hi = a[k] >= np.median(a[k])
    if hi.sum() < 50 or (~hi).sum() < 50: return None, int(k.sum()), "劈开后某臂<50"
    p1, p0 = b[k][hi].mean(), b[k][~hi].mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return None, int(k.sum()), "某臂比例为 0 或 1"
    return float(math.log((p1 / (1 - p1)) / (p0 / (1 - p0)))), int(k.sum()), "ok"
ct = np.quantile(COH, [1/3, 2/3]); et = np.quantile(EXT, [1/3, 2/3])
cb, eb = np.digitize(COH, ct), np.digitize(EXT, et)
LINES = {"samesex→samesexany": (SX[:, 0], SS), "sxok18→hadsex": (SX[:, 1], HD)}
rows, bad = [], []
print("\n=== 全格,含不可算的格(**在三分之前打印**)===")
for lname, (att, beh) in LINES.items():
    for c in range(3):
        for e in range(3):
            m = (cb == c) & (eb == e)
            v, n, why = lnor(att, beh, m)
            rows.append(dict(line=lname, coh=c, ext=e, lnor=v, n=n, status=why,
                             inclusion=[f"coh_fam 三分位 {c}", f"ext_fam 三分位 {e}", f"n={n}", why]))
            if v is None: bad.append(f"{lname}|c{c}e{e}({why},n={n})")
            print(f"  {lname:22s} c={c} e={e}  n={n:4d}  " +
                  (f"|lnOR|={abs(v):.4f}" if v is not None else f"**不可算:{why}**"))
print(f"\n  不可算的格:{len(bad)}/{len(rows)}  {bad if bad else '无'}")
G = Gate("用不相交的题算一致性,问题问得出来了吗?")
pcv, pcn, _ = lnor(SX[:, 1], HD, np.ones(len(HD), bool))
rng = np.random.default_rng(SEEDS[0])
perm = [abs(lnor(SX[:, 1], HD[rng.permutation(len(HD))], np.ones(len(HD), bool))[0] or 0) for _ in range(200)]
G.positive_control("正对照:sxok18→hadsex", planted=abs(pcv), floor=float(np.quantile(perm, .95)), spread=1e-9)
G.negative_control("g=0:置换行为后必须消失", null=float(np.median(perm)), effect=abs(pcv),
                   null_spread=float(np.std(perm)), null_kind="个体层行为标签置换")
tg = rng.integers(0, 2, len(COH)).astype(float)
G.negative_control("安慰剂:coh_fam → 随机标签", null=abs(lnor(COH, tg, np.ones(len(COH), bool))[0] or 0),
                   effect=abs(pcv), null_spread=1e-9, null_kind="与问卷无关的随机二值标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['line'][:8]}|c{r['coh']}e{r['ext']}": r for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件",
                                     {f"{r['line'][:8]}|c{r['coh']}e{r['ext']}": r for r in rows})
print("\n" + "=" * 76)
if bad:
    world = "UNVERIFIED"
    verdict = (f"**{len(bad)}/{len(rows)} 格不可算 -> 直接 UNVERIFIED,不进三分逻辑** "
               f"(修 `#533a`/`#537c` 的同型:不可算不是一个世界)")
    print(f"⚠ {verdict}")
elif abs(pcv) > np.quantile(perm, .95) and np.median(perm) < 0.5 * abs(pcv):
    effs = {}
    for lname in LINES:
        d = [abs([r for r in rows if r["line"] == lname and r["coh"] == 2 and r["ext"] == e][0]["lnor"]) -
             abs([r for r in rows if r["line"] == lname and r["coh"] == 0 and r["ext"] == e][0]["lnor"])
             for e in range(3)]
        effs[lname] = float(np.mean(d))
    spread = float(np.std([abs(r["lnor"]) for r in rows]))
    vals = list(effs.values())
    if all(abs(v) > spread for v in vals) and all(v > 0 for v in vals):
        world = "W-GENERAL"; verdict = f"控制 ext_fam 后 coh_fam 效应 {[f'{v:+.4f}' for v in vals]} 超展布 {spread:.4f} -> **道德一致性是跨领域的人格性质**"
    elif all(abs(v) <= spread for v in vals):
        world = "W-LOCAL"; verdict = (f"控制 ext_fam 后 coh_fam 效应 {[f'{v:+.4f}' for v in vals]} "
            f"落在全格展布 {spread:.4f} 内 -> **一致性是领域内的,不外溢**")
    else:
        world = "W-LOCAL"; verdict = f"效应符号不一致 {[f'{v:+.4f}' for v in vals]} -> **不外溢**"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 估计量已经变了:本轮问的是「她在**别处**的道德一致性,预测她这里的言行一致吗」,"
          "不是上一轮那个问题的更好版本。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rows=rows, uncomputable=bad, world=world, verdict=verdict,
               corr_coh_ext=float(np.corrcoef(COH, EXT)[0, 1]), n=int(ok.sum()), base_rate=br,
               positive_control=dict(lnor=pcv, n=pcn), seeds=SEEDS,
               estimand_changed="coherence 来自家庭题,预测变量是性题 —— 与 #537 不是同一个问题",
               instrument="NSFG 2011-2013 女性 ACASI,横断面",
               impossible=["横断面时序不可分", "仅女性", "单一波"], unchallenged=True),
          open(OUT / "disjoint_coherence.json", "w"), indent=1)
print(f"\nwrote {OUT/'disjoint_coherence.json'}")

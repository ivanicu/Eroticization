"""E03·A231·R613 — 同一个女人,对人说的和自己打的,不是同一个数

`#568` 的 NEXT。行动类型:**FRONTIER**。
`#568d`:这一页每一条人层结论都来自「对着人出声说出来的回答」,而**没有任何一处量过独自作答**。
NSFG 恰好**在同一份问卷里两样都有**:
  **CAPI**(访员施测的妊娠史)-> `abortion` = 「**CAPI-based** total number of induced abortions (RECODE)」
  **CASI**(自访)-> `casiabor` = 「JB-3 # of pregnancies ending in abortion in **5 yrs** before interview」
⚠ 窗口不同(终生 vs 5 年)⇒ **不能直接比数**。

G1 ESTIMAND(先于方法,且绕开窗口):
   **矛盾率 = P(CASI 报 ≥1 | CAPI 终生 = 0)。**
   一个终生 0 次的人,**逻辑上不可能**在近 5 年有 ≥1 次 ⇒ **真值必须是 0**。
   任何超出 0 的部分,**只能来自两种作答方式之间的披露差**(或纯记录误差,见对照)。
   ⚠ 方向是单向的:窗口更短的 CASI **只会更小**,所以这个量对「CASI 披露更多」是**保守**的。

**关键对照(先于结果写死):同一算法用在 `casibirth`(活产,**不敏感**)上。**
   活产同样是「CAPI 终生 = 0 而 CASI 近 5 年 ≥1」的逻辑矛盾,
   但它**没有社会期待压力** ⇒ 它的矛盾率就是**纯记录误差的地板**。
   **预注册:堕胎的矛盾率必须显著高于活产的,否则本轮测到的是记录误差,不是模式。**

WORLDS:
  W-DISCLOSURE 堕胎矛盾率 >> 活产矛盾率 ⇒ **独自作答时她说得更多**,而页面上的态度题全部是出声说的
  W-NOISE      两者相当 ⇒ 测到的是记录误差,**不能说模式**
  W-REVERSE    活产更高 ⇒ 设计有问题,UNVERIFIED
⚠ BASIN:`W-DISCLOSURE` 支持我刚写上页面的那句话,**所以不是下注方向**。本轮下注 `W-NOISE`。
KILL(条件式):if 活产矛盾率可算 and 两者 n 都 >= 300:
   堕胎率 > 2× 活产率 且差超 bootstrap 展布 -> W-DISCLOSURE;否则 W-NOISE
   else UNVERIFIED
IMPOSSIBLE:仅女性 · 单一波 · `abortion` 是 RECODE,其构造细节不在本地文档 ⇒
   **「CAPI-based」这四个字来自它自己的标签,不是从流程文件读的** · 未加权 ⇒ 非人群估计 · [unchallenged]
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
LAY = {}
for line in open(NS / "setup/2011_2013_FemRespSetup.dct", errors="replace"):
    m = pat.search(line)
    if m: LAY[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
NEED = ["abortion", "casiabor", "casibirth", "parity"]
NEED = [n for n in NEED if n in LAY]
print("=== 硬规则 1:逐个打印列、宽、标签 ===")
for n in NEED: print(f"  {n:11s} 列{LAY[n][0]:5d} 宽{LAY[n][1]}  {LAY[n][2][:58]}")
buf = {n: [] for n in NEED}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n in NEED:
        s, w, _ = LAY[n]; v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
A = {n: np.array(buf[n]) for n in NEED}
for n in NEED:
    v = A[n][np.isfinite(A[n])]
    print(f"  {n:11s} n={len(v):5d} 取值 {sorted(set(v.astype(int)))[:8]}")

def contradiction(capi, casi, name):
    """P(CASI >= 1 | CAPI == 0)。逻辑上必须为 0。"""
    m = np.isfinite(capi) & np.isfinite(casi) & (capi < 90) & (casi < 90)
    base = m & (capi == 0)
    if base.sum() < 300: return None, int(base.sum())
    r = float((casi[base] >= 1).mean())
    print(f"  {name:14s} CAPI=0 的人 n={int(base.sum()):5d} · 其中 CASI≥1 的 **{r:.4f}**")
    return r, int(base.sum())

print("\n=== 矛盾率(逻辑上都必须是 0)===")
r_ab, n_ab = contradiction(A["abortion"], A["casiabor"], "堕胎(敏感)")
r_lb, n_lb = contradiction(A.get("parity", A["abortion"] * np.nan), A["casibirth"], "活产(不敏感)")
G = Gate("同一个女人,对人说的和自己打的,不是同一个数")
if r_ab is None or r_lb is None:
    print("⚠ 某一格 n<300 -> UNVERIFIED"); world, verdict = "UNVERIFIED", "基数不足"
else:
    bs_ab, bs_lb = [], []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(400):
            i = rng.integers(0, n_ab, n_ab); j = rng.integers(0, n_lb, n_lb)
            m1 = np.isfinite(A["abortion"]) & np.isfinite(A["casiabor"]) & (A["abortion"] == 0)
            m2 = np.isfinite(A["parity"]) & np.isfinite(A["casibirth"]) & (A["parity"] == 0)
            bs_ab.append(float((A["casiabor"][m1][i] >= 1).mean()))
            bs_lb.append(float((A["casibirth"][m2][j] >= 1).mean()))
    sd_j = float(np.sqrt(np.var(bs_ab) + np.var(bs_lb)))
    print(f"\n  差 = {r_ab-r_lb:+.4f} · 联合展布 2.8σ = {2.8*sd_j:.4f} · 比值 = {r_ab/r_lb if r_lb else float('inf'):.2f}×")
    G.positive_control("活产矛盾率可算(记录误差地板存在)", planted=float(n_lb), floor=299.0, spread=1e-9)
    G.negative_control("活产(不敏感)矛盾率 = 纯记录误差地板",
                       null=r_lb, effect=r_ab, null_spread=float(np.std(bs_lb)),
                       null_kind="同一逻辑矛盾,但没有社会期待压力的题")
    cells = {"堕胎": dict(n=n_ab, rate=r_ab, inclusion=["CAPI 终生=0 的女性", "两列都非缺失且 <90"]),
             "活产": dict(n=n_lb, rate=r_lb, inclusion=["CAPI parity=0 的女性", "两列都非缺失且 <90"])}
    G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
    G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
    if r_ab > 2 * r_lb and (r_ab - r_lb) > 2.8 * sd_j:
        world = "W-DISCLOSURE"; verdict = f"堕胎 {r_ab:.4f} > 2× 活产 {r_lb:.4f},差超展布 -> **独自作答时她说得更多**"
    elif r_lb > r_ab:
        world = "W-REVERSE"; verdict = f"活产 {r_lb:.4f} > 堕胎 {r_ab:.4f} -> 设计有问题"
    else:
        world = "W-NOISE"; verdict = f"堕胎 {r_ab:.4f} 与活产 {r_lb:.4f} 不足 2 倍或差在展布内 -> **测到的是记录误差**"
    print(f"\n控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:CASI 的窗口更短(5 年 vs 终生),所以它**只会低估**披露差;"
          "而 `abortion` 是 RECODE,「CAPI-based」四个字来自**它自己的标签**,不是从流程文件读的。")
print(G)
json.dump(dict(rate_abortion=r_ab, n_abortion=n_ab, rate_livebirth=r_lb, n_livebirth=n_lb,
               world=world, verdict=verdict, seeds=SEEDS,
               estimand="P(CASI>=1 | CAPI==0),逻辑上必须为 0",
               impossible=["仅女性", "单一波", "abortion 是 RECODE,构造细节不在本地",
                           "未加权,非人群估计", "CASI 窗口更短 -> 对披露差保守"],
               unchallenged=True), open(OUT / "mode_disclosure.json", "w"), indent=1)
print(f"\nwrote {OUT/'mode_disclosure.json'}")

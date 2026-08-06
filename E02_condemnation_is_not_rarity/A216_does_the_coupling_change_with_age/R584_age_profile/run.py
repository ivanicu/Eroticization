"""E02·A216·R584 — 那个 0.425 的耦合,在人的一生里变吗?

`#538` 的 NEXT。行动类型:**FRONTIER**。改方向后的第一轮。
`#536` 测出性态度题内部耦合 0.4246 vs 家庭题 0.1401 —— 但那是**一个横断面的单一数字**。
本轮问:**它随年龄变吗?**

⚠ **这个问题的答案在结构上是二义的,而这是结论的一部分,不是限制条款:**
   同一波调查里的年龄差 **既是发展(同一个人会变),也是世代(不同年代长大的人不同)**。
   NSFG 2011–2013 是**单一波** ⇒ **年龄与出生队列完全共线,不可分。**
   ⇒ 任何走向都只能写成「**随年龄或队列**」,**永远不能写成「人会变」**。

G1 ESTIMAND(先于方法):在每个年龄组内,
   `ρ_sex(age)` = 3 道性题两两 |ρ| 中位;`ρ_fam(age)` = 7 道家庭题两两 |ρ| 中位;
   **主量 = `ρ_sex(age) − ρ_fam(age)` 随年龄的走向**(用差,不用 `ρ_sex` 本身,
   因为整体作答一致性可能随年龄变,差把它抵消掉)。

**预注册(`#538c`,写在看到任何计算之前):判决按行。**
   每一行(`ρ_sex` · `ρ_fam` · 差)独立报 CONFIRMED/OVERTURNED/UNVERIFIED;
   **任一行有不可算的格,该行 UNVERIFIED,其余行照常判。**

WORLDS:
  W-STABLE    差不随年龄变 ⇒ 「性是更紧的领域」是一个**跨年龄稳定**的结构
  W-TIGHTENS  差随年龄上升 ⇒ 年长/更早出生的人那里,性更是一个整块
  W-LOOSENS   差随年龄下降
⚠ BASIN:`W-STABLE` 让 `#536` 更省事,**不是**下注方向。本轮下注 `W-TIGHTENS` ——
   它一旦为真,`#536` 的 0.425 就是**一个平均**,而不是一个结构常数。
CONTROLS:正对照 `sxok18`×`sxok16` 在每个年龄组必须最高(仪器上限,逐组算)·
   安慰剂 每组的态度题 × 随机标签 ≈ 0 · 组内 n 全部打印 · 全格公布
KILL(条件式,按行):if 该行全部格可算 and 正对照在该组通过:
   斜率 |slope| > bootstrap MDE -> W-TIGHTENS/W-LOOSENS(按符号) else W-STABLE
   否则该行 UNVERIFIED
IMPOSSIBLE:**年龄与队列完全共线,不可分** · 仅女性 · 15–44 岁 ⇒ **不覆盖中老年** ·
   横断面 ⇒ 非因果 · 3 道性题只给 3 对 ⇒ `ρ_sex` 分辨率远低于 `ρ_fam` · [unchallenged]
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
AGEV = next((v for v in ["ager", "age_r", "ager_i"] if v in LAY), None)
assert AGEV, f"没有找到年龄变量 —— 候选都不在 dct 里"
cols = {n: LAY[n] for n in SEX + FAM + [AGEV]}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip(); buf[n].append(float(v) if v not in ("", ".") else np.nan)
AGE = np.array(buf[AGEV])
print(f"=== 硬规则 1 ===\n  年龄变量 `{AGEV}`:{cols[AGEV][2][:44]}")
print(f"  n={int(np.isfinite(AGE).sum())} 范围 {np.nanmin(AGE):.0f}-{np.nanmax(AGE):.0f} 中位 {np.nanmedian(AGE):.0f}")
assert 14 <= np.nanmin(AGE) <= 16 and 43 <= np.nanmax(AGE) <= 45, "年龄范围不是 15-44,码读错了"
X = {n: np.where(np.isin(np.array(buf[n]), [1, 2, 3, 4, 5]), buf[n], np.nan) for n in SEX + FAM}
BINS = [(15, 20), (21, 25), (26, 30), (31, 35), (36, 40), (41, 44)]
def rho(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 150: return None, int(m.sum())
    return float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]), int(m.sum())
def med(items, sel):
    out, nmin = [], 10 ** 9
    for a, b in itertools.combinations(items, 2):
        r, n = rho(X[a][sel], X[b][sel])
        if r is None: return None, n
        out.append(abs(r)); nmin = min(nmin, n)
    return float(np.median(out)), nmin
rows = []
print("\n=== 逐年龄组:先打 n,再看 ρ(差 = 性内 − 家内,抵消整体作答一致性)===")
for lo, hi in BINS:
    sel = (AGE >= lo) & (AGE <= hi)
    ms, ns = med(SEX, sel); mf, nf = med(FAM, sel)
    pc, npc = rho(X["sxok18"][sel], X["sxok16"][sel])
    rows.append(dict(band=f"{lo}-{hi}", n=int(sel.sum()), rho_sex=ms, rho_fam=mf,
                     diff=(ms - mf) if (ms is not None and mf is not None) else None,
                     pos_ctrl=abs(pc) if pc is not None else None, n_min_sex=ns, n_min_fam=nf,
                     inclusion=[f"年龄 {lo}-{hi},n={int(sel.sum())}", "每对 n>=150",
                                "3 道性题 / 7 道家庭题", "单一波:年龄与队列共线"]))
    print(f"  {lo}-{hi}  n={int(sel.sum()):4d}  性内={ms if ms is None else f'{ms:.4f}'}  "
          f"家内={mf if mf is None else f'{mf:.4f}'}  "
          f"**差={'不可算' if rows[-1]['diff'] is None else f'{rows[-1][chr(100)+chr(105)+chr(102)+chr(102)]:+.4f}'}**  "
          f"上限(18×16)={pc if pc is None else f'{abs(pc):.4f}'}")
mid = np.array([(lo + hi) / 2 for lo, hi in BINS])
G = Gate("那个 0.425 的耦合,在人的一生里变吗?(NSFG 年龄剖面)")
LINES = {"rho_sex": "性内耦合", "rho_fam": "家内耦合", "diff": "差(主量)"}
out = {}
for key, label in LINES.items():
    vals = [r[key] for r in rows]
    if any(v is None for v in vals):
        out[key] = dict(verdict="UNVERIFIED", why="有不可算的格(预注册:该行 UNVERIFIED)")
        print(f"\n  [{label}] **UNVERIFIED —— 有不可算的格**"); continue
    v = np.array(vals); slope = float(np.polyfit(mid, v, 1)[0])
    bs = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(2000):
            i = rng.integers(0, len(mid), len(mid))
            if np.ptp(mid[i]) > 0: bs.append(np.polyfit(mid[i], v[i], 1)[0])
    MDE = 2.8 * float(np.std(bs))
    w = ("W-STABLE" if abs(slope) <= MDE else ("W-TIGHTENS" if slope > 0 else "W-LOOSENS"))
    out[key] = dict(verdict=w, slope=slope, MDE=MDE, values=[float(x) for x in v])
    print(f"\n  [{label}] 斜率={slope:+.6f}/岁  MDE={MDE:.6f}  -> **{w}**  值={[f'{x:.4f}' for x in v]}")
pcs = [r["pos_ctrl"] for r in rows if r["pos_ctrl"] is not None]
G.positive_control("正对照:每组 sxok18×sxok16 为上限", planted=float(np.median(pcs)),
                   floor=float(np.median([r["rho_sex"] for r in rows if r["rho_sex"] is not None])), spread=1e-9)
rng = np.random.default_rng(SEEDS[0]); tg = rng.integers(0, 5, len(AGE)).astype(float)
zs = [abs(rho(X[n], tg)[0] or 0) for n in SEX + FAM]
G.negative_control("安慰剂:态度题 × 随机标签", null=float(np.median(zs)), effect=float(np.median(pcs)),
                   null_spread=float(np.std(zs)), null_kind="与问卷无关的随机整数标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {r["band"]: r for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {r["band"]: r for r in rows})
print("\n" + "=" * 76)
print(f"逐行判决(预注册按行):{ {LINES[k]: v['verdict'] for k, v in out.items()} }")
print("⚠ 无论走向如何,**只能写成「随年龄或队列」** —— 单一波里年龄与出生队列完全共线,"
      "**永远不能写成「人会变」**。这是结论的一部分,不是限制条款。")
print(G)
json.dump(dict(rows=rows, lines=out, age_var=AGEV, seeds=SEEDS,
               prereg="判决按行;任一行有不可算的格,该行 UNVERIFIED(#538c,写于本轮之前)",
               instrument="NSFG 2011-2013 女性 ACASI,单一波",
               impossible=["年龄与出生队列完全共线,不可分", "仅女性", "15-44 岁不覆盖中老年",
                           "横断面非因果", "3 道性题只给 3 对"], unchallenged=True),
          open(OUT / "age_profile.json", "w"), indent=1)
print(f"\nwrote {OUT/'age_profile.json'}")

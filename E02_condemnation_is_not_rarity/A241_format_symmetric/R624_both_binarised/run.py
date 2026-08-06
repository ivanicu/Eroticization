"""E02·A241·R624 — 把格式完全对称之后,那条差还在吗?

`#579` 的 NEXT。行动类型:**FRONTIER**(结局决定 `#534`/`#579` 要不要带上「只在格式不对称时成立」)。

**要解决的老限制(`#535c`,自那以来一直开着):**
GSS 上非性道德题**原生二值**、性题**四级** ⇒ 两组的信息量不同,
`#535` 只量了「格式解释掉 13%–45%」,**没有消掉它**。
⚠ 而 **NSFG 做得到 GSS 做不到的事**:那里性题与家庭题**同为五级**(`#536a`)——
所以可以**把两边一起二值化**,得到一个**格式完全对称**的对比。

G1 ESTIMAND(先于方法):`性内 − 家内`,四个格式规格并排:
  F1 **两边都五级**(`#536` 原版)· F2 **两边都二值·严切点** · F3 **两边都二值·宽切点** ·
  F4 **性二值 / 家五级**(**故意做成不对称**,当作那个混淆的正对照)
预注册:
  **F2/F3(对称)仍为正,且量级落在 F1 的 ±50% 内** -> **格式不是这条结论的解释**;
  **F2/F3 塌到 0 附近或翻号** -> **`#534`/`#579` 必须带上「只在格式不对称时成立」**;
  **F4 明显偏离 F2/F3** -> 说明不对称本身确实会移动这个量(那正是 `#535` 量到的东西)。
CONTROLS:正对照 `sxok18`×`sxok16`(近重复题)在每个规格下必须最高 ·
  安慰剂 性题 × 随机整数标签 ≈ 0 · 逐格 n 全部打印
IMPOSSIBLE:仅女性 · 单一波 · 三道性题只给三对 ⇒ `性内` 的分辨率远低于 `家内` 的 21 对 ·
  二值化丢信息 ⇒ 对称版可能整体偏低 · [unchallenged]
"""
import os, sys, pathlib, json, re, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import rankdata
from lib.gates import Gate
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
LAY = {}
for line in open(NS / "setup/2011_2013_FemRespSetup.dct", errors="replace"):
    m = pat.search(line)
    if m: LAY[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
SEX = ["samesex", "sxok18", "sxok16"]
FAM = ["staytog", "chunless", "chsuppor", "okcohab", "marrfail", "chcohab", "prvntdiv"]
ALL = SEX + FAM
buf = {n: [] for n in ALL}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n in ALL:
        s, w, _ = LAY[n]; v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
X5 = {n: np.where(np.isin(np.array(buf[n]), [1, 2, 3, 4, 5]), buf[n], np.nan) for n in ALL}
def binar(v, cut): return np.where(np.isnan(v), np.nan, np.isin(v, cut).astype(float))
XA = {n: binar(X5[n], [1]) for n in ALL}          # 严:仅码 1
XB = {n: binar(X5[n], [1, 2]) for n in ALL}       # 宽:码 1–2
print("=== 硬规则 1:逐题 n 与两种二值化后的阳性率 ===")
for n in ALL:
    print(f"  {n:10s} n={int(np.isfinite(X5[n]).sum()):5d}  严={np.nanmean(XA[n]):.4f} 宽={np.nanmean(XB[n]):.4f}")
def rho(D, a, b):
    m = np.isfinite(D[a]) & np.isfinite(D[b])
    if m.sum() < 200 or np.std(D[a][m]) == 0 or np.std(D[b][m]) == 0: return np.nan
    return abs(float(np.corrcoef(rankdata(D[a][m]), rankdata(D[b][m]))[0, 1]))
def med(D, items):
    o = [rho(D, a, b) for a, b in itertools.combinations(items, 2)]
    o = [x for x in o if np.isfinite(x)]
    return float(np.median(o)) if o else np.nan
def mixed(Dsex, Dfam, a, b, insex):
    D = {**{k: Dsex[k] for k in SEX}, **{k: Dfam[k] for k in FAM}}
    return rho(D, a, b)
SPECS = {}
SPECS["F1 两边都五级(#536 原版)"] = (med(X5, SEX), med(X5, FAM))
SPECS["F2 两边都二值·严"] = (med(XA, SEX), med(XA, FAM))
SPECS["F3 两边都二值·宽"] = (med(XB, SEX), med(XB, FAM))
MIX = {**{k: XA[k] for k in SEX}, **{k: X5[k] for k in FAM}}
SPECS["F4 性二值/家五级(故意不对称)"] = (med(MIX, SEX), med(MIX, FAM))
print("\n=== 四个格式规格(全格公布)===")
res = {}
for k, (s, f) in SPECS.items():
    res[k] = dict(sex=s, fam=f, gap=s - f)
    print(f"  {k:26s} 性内={s:.4f} 家内={f:.4f}  **差={s-f:+.4f}**")
F1, F2, F3, F4 = (res[k]["gap"] for k in SPECS)
print(f"\n  对称版 F2={F2:+.4f} · F3={F3:+.4f}  vs 原版 F1={F1:+.4f}  "
      f"(保留 {F2/F1:.0%} / {F3/F1:.0%})· 不对称 F4={F4:+.4f}")
G = Gate("把格式完全对称之后,那条差还在吗?")
# ⚠ 第一版把「上限对在三个规格上的**中位**」去比「性内中位在三个规格上的**最大**」——
#   两个不同的聚合,`#576` 那一类错(拿一个配对方式不同的不确定度去比)。
#   改:**按规格逐个比** —— 每个规格里,上限对必须高于该规格自己的性内中位。
_SPEC_D = {"F1 两边都五级(#536 原版)": X5, "F2 两边都二值·严": XA, "F3 两边都二值·宽": XB}
for _k, _D in _SPEC_D.items():
    _ceil = rho(_D, "sxok18", "sxok16")
    G.positive_control(f"正对照[{_k[:2]}]:sxok18×sxok16 高于该规格的性内中位",
                       planted=float(_ceil), floor=float(med(_D, SEX)), spread=1e-9)
rng = np.random.default_rng(SEEDS[0]); tg = rng.integers(0, 5, len(X5["samesex"])).astype(float)
zs = []
for n in SEX:
    m = np.isfinite(XA[n])
    zs.append(abs(float(np.corrcoef(rankdata(XA[n][m]), rankdata(tg[m]))[0, 1])))
G.negative_control("安慰剂:二值性题 × 随机标签", null=float(np.median(zs)),
                   effect=abs(F2), null_spread=float(np.std(zs)), null_kind="与问卷无关的随机整数标签")
cells = {k: dict(n=int(np.isfinite(X5["samesex"]).sum()), **res[k],
                 inclusion=[k, "每对 n>=200", "同一批 5,600 名女性"]) for k in SPECS}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 74)
sym_ok = (F2 > 0 and F3 > 0 and 0.5 * F1 <= F2 <= 1.5 * F1 and 0.5 * F1 <= F3 <= 1.5 * F1)
if sym_ok:
    world = "FORMAT-NOT-THE-EXPLANATION"
    verdict = f"两个对称规格都为正且落在原版 ±50% 内(F2 {F2:+.4f}·F3 {F3:+.4f} vs F1 {F1:+.4f}) -> **格式不是这条结论的解释**"
elif F2 <= 0 or F3 <= 0:
    world = "FORMAT-IS-THE-EXPLANATION"; verdict = f"对称之后塌掉(F2 {F2:+.4f}·F3 {F3:+.4f}) -> **必须带上「只在格式不对称时成立」**"
else:
    world = "PARTIAL"; verdict = (f"对称版仍为正但量级移出 ±50%(F2 {F2:+.4f}·F3 {F3:+.4f} vs F1 {F1:+.4f})"
        f" -> **方向存活,数字必须报成跨格式区间**")
print(f"评判:**{world}** —— {verdict}")
print(f"⚠ 而 F4(故意不对称)= {F4:+.4f},与对称版差 {abs(F4-np.mean([F2,F3])):+.4f} —— "
      f"这就是 `#535` 量到的那个格式效应本身,现在它被单独拿出来了。")
print("⚠ 这个 KILL 会怎样失败:三道性题只给三对,`性内` 的分辨率远低于 `家内` 的 21 对;"
      "二值化丢信息,所以对称版整体偏低是**预期**的,判据用的是**比值**不是差值。")
print(G)
json.dump(dict(specs=res, F1=F1, F2=F2, F3=F3, F4=F4, world=world, verdict=verdict, seeds=SEEDS,
               impossible=["仅女性", "单一波", "三道性题只给三对", "二值化丢信息,对称版整体偏低是预期"],
               unchallenged=True), open(OUT / "format_symmetric.json", "w"), indent=1)
print(f"\nwrote {OUT/'format_symmetric.json'}")

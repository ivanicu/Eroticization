"""E02·A225·R600 — 那个 0.985,是性专有的,还是整体风险行为在下降?

`#554` 的 NEXT,**预测矩阵已写在 `#554d`,先于本轮任何计算**。行动类型:**FRONTIER**。

**要分离的三个世界(预注册,`#554` NEXT ②):**
  W-GENERAL-RISK   `q55`(注射毒品,**非性**)也与五题同步 ⇒ 同步是**整体青少年风险行为下降**,
                   与性无关 ⇒ `#554b` 降级
  W-SEX-SPECIFIC   `q55` **不**同步,而 `q19`–`q21`(性暴力,**性相关但不下游于性活跃**)同步
                   ⇒ 同步是**性专有**的
  W-ACTIVITY-ONLY  `q19`–`q21` 也**不**同步 ⇒ 同步**只存在于下游于性活跃的题之间**
                   ⇒ `#554c` 的推导成立,那个 0.985 是同一件事被数了五遍
**发布规则(预注册):** 仅在 `W-SEX-SPECIFIC` 或 `W-GENERAL-RISK` 下才上页面;
   `W-ACTIVITY-ONLY` 下**页面只写一句「这五道题不是五种行为」**。

方法与 `#599` 完全相同(噪声配平 + 受访者不相交 + 打乱年份安慰剂),**只加变量,不改设计**。
IMPOSSIBLE:同 `#554e` · 且 `q19`–`q21` 的年份覆盖可能短于 17 年 ⇒ **逐题打印年份**(硬规则 1)
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, re
Y = ROOT / "data/external/yrbs"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
sas = (Y / "2023-SADC-SAS-Input-Program.sas").read_text(errors="replace")
lab = dict(re.findall(r'^\s*(\w+)\s*=\s*"([^"]*)"', sas, re.M))
pos = {m.group(1): (int(m.group(3)) - 1, int(m.group(4)))
       for m in re.finditer(r'^\s*(\w+)\s+(\$\s+)?(\d+)-(\d+)\s*$', sas, re.M)}
ACT = ["q56", "q57", "q58", "q59", "q61"]          # 下游于「是否性活跃」
VIO = ["q19", "q20", "q21"]                         # 性相关,**不**下游于性活跃
RISK = ["q55"]                                      # 非性风险行为
NEED = ["year", "sex"] + ACT + VIO + RISK
assert all(n in pos for n in NEED)
rec = {n: [] for n in NEED}
for f in sorted(Y.glob("sadc_2023_state_*.dat")):
    for line in open(f, errors="replace"):
        for n in NEED:
            s, e = pos[n]; v = line[s:e].strip()
            rec[n].append(float(v) if v not in ("", ".") else np.nan)
A = {n: np.array(rec[n]) for n in NEED}
print(f"=== 硬规则 1:{len(A['year']):,} 行,逐题年份与 n ===")
for n in ACT + VIO + RISK:
    ok = np.isfinite(A[n]); ys = sorted(set(A["year"][ok].astype(int)))
    print(f"  {n} [{'活跃下游' if n in ACT else ('性暴力' if n in VIO else '非性风险')}] "
          f"n={int(ok.sum()):8,} 年 {len(ys):2d} 个 {ys[0]}–{ys[-1]}  {lab.get(n,'')[:30]}")

def series(q, sx):
    m = np.isfinite(A[q]) & (A["sex"] == sx)
    return {y: float((A[q][m & (A["year"].astype(int) == y)] == 1).mean())
            for y in sorted(set(A["year"][m].astype(int)))
            if (m & (A["year"].astype(int) == y)).sum() >= 300}

ALL = ACT + VIO + RISK
S = {(q, sx): series(q, sx) for q in ALL for sx in (1, 2)}
def dc(a, b, perm=None):
    yr = sorted(set(a) & set(b))
    if len(yr) < 8: return np.nan
    va = np.array([a[y] for y in yr]); vb = np.array([b[y] for y in yr])
    if perm is not None: vb = vb[perm.permutation(len(yr))]
    da, db = np.diff(va), np.diff(vb)
    if np.std(da) == 0 or np.std(db) == 0: return np.nan
    return float(np.corrcoef(da, db)[0, 1])

def blockmed(G1, G2=None):
    out = []
    pairs = list(itertools.combinations(G1, 2)) if G2 is None else [(a, b) for a in G1 for b in G2]
    for a, b in pairs:
        v = [dc(S[(a, i)], S[(b, j)]) for i, j in ((1, 2), (2, 1))]
        v = [x for x in v if np.isfinite(x)]
        if v: out.append(float(np.mean(v)))
    return (float(np.median(out)), len(out)) if out else (np.nan, 0)

same = [dc(S[(q, 1)], S[(q, 2)]) for q in ALL]
same = [x for x in same if np.isfinite(x)]
SM = float(np.median(same))
shuf = []
for sd in SEEDS:
    r = np.random.default_rng(sd)
    for _ in range(60):
        for a, b in itertools.combinations(ALL, 2):
            x = dc(S[(a, 1)], S[(b, 2)], perm=r)
            if np.isfinite(x): shuf.append(abs(x))
Q95 = float(np.quantile(shuf, .95))
cells = {
    "活跃下游 × 活跃下游": blockmed(ACT),
    "性暴力 × 性暴力": blockmed(VIO),
    "活跃下游 × 性暴力": blockmed(ACT, VIO),
    "活跃下游 × 非性风险": blockmed(ACT, RISK),
    "性暴力 × 非性风险": blockmed(VIO, RISK),
}
print(f"\n=== 全网格(同题上限 {SM:+.4f} · 安慰剂 q95 {Q95:.4f})===")
for k, (m, n) in cells.items():
    tag = "**同步**" if np.isfinite(m) and m >= 0.5 * SM else ("**不同步**" if np.isfinite(m) and abs(m) < Q95 else "不判")
    print(f"  {k:22s} 中位 {m:+.4f}({n:2d} 对)  {tag}")
print("\n" + "=" * 74)
def sync(k): 
    m, n = cells[k]
    return (np.isfinite(m) and m >= 0.5 * SM), (np.isfinite(m) and abs(m) < Q95)
risk_s, risk_n = sync("活跃下游 × 非性风险")
vio_s, vio_n = sync("活跃下游 × 性暴力")
if abs(SM) > Q95:
    if risk_s:
        world = "W-GENERAL-RISK"; verdict = "非性风险也同步 -> **同步是整体风险行为下降,与性无关;`#554b` 降级**"
    elif vio_s and risk_n:
        world = "W-SEX-SPECIFIC"; verdict = "性暴力同步而非性风险不同步 -> **同步是性专有的**"
    elif vio_n:
        world = "W-ACTIVITY-ONLY"; verdict = ("性暴力**不**同步 -> **同步只存在于下游于性活跃的题之间;"
                                             "`#554c` 的推导成立**")
    else:
        world = "UNVERIFIED"; verdict = "格与格之间不一致 -> 不判"
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print(f"  发布规则(预注册):{'上页面' if world in ('W-SEX-SPECIFIC','W-GENERAL-RISK') else '**页面只写「这五道题不是五种行为」**'}")
else:
    world, verdict = "UNVERIFIED", f"同题上限 {SM:.4f} 未越过安慰剂 q95 {Q95:.4f}"
    print(f"⚠ {verdict}")
json.dump(dict(cells={k: [None if not np.isfinite(v[0]) else v[0], v[1]] for k, v in cells.items()},
               same_median=SM, placebo_q95=Q95, world=world, verdict=verdict, seeds=SEEDS,
               groups=dict(activity=ACT, violence=VIO, nonsex_risk=RISK),
               prereg="预测矩阵与发布规则写于 #554 NEXT,先于本轮任何计算",
               inclusion=["YRBS SADC 州文件", "每年每性别 n>=300", "码 1 = 是", "未加权"],
               impossible=["青少年不可外推", "校内抽样辍学者缺席", "未加权", "行为≠态度"],
               unchallenged=True), open(OUT / "sex_or_risk.json", "w"), indent=1)
print(f"\nwrote {OUT/'sex_or_risk.json'}")

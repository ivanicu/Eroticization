"""E02·A225·R599 — 年代那一层的第二份仪器:不同性行为是否一起变?

`#553` 的 NEXT。行动类型:**FRONTIER**。
E02 的十二条断言里,**「年代」那一层只有 GSS 一个来源**(`#532`/`#533`)。
YRBS SADC(1991–2023)是**从未被这条线用过**的第二份:
**不同的仪器 · 不同的人群(青少年)· 不同的模式(校内自填)· 不同的机构。**
⚠ 但它问的是**行为**,不是态度 ⇒ **这不是 `#532` 的复制,是它的行为版**。
   **估计量变了,必须明说**,不能当成「同一个结论的第二个证据」。

G1 ESTIMAND(先于方法):对每一年 `t`,每道行为题算一个流行率 `p_i(t)`;
   **主量 = `corr(Δp_i, Δp_j)` 跨相邻调查年**(与 `#532` 同形),取非对角中位。
   ⚠ 与 `#577` 同样的噪声配平:**同题按性别分半**(`sex` 1=女 2=男),
   目标 = 男的题 A × 女的题 B,**受访者不相交**;上限 = 同题男×女。

WORLDS:
  W-SYNC   跨题 Δ 相关 ≥ 同题上限的一半 ⇒ 行为**一起动**
  W-ONE-BY-ONE 跨题 ≈ 安慰剂 ⇒ **一件一件动**(`#532` 在行为上、在另一份仪器上复制)
  W-BLIND  介于两者 ⇒ UNVERIFIED-by-power
⚠ BASIN:`W-ONE-BY-ONE` 会让 `#532` 漂亮地跨仪器复制,**所以不是下注方向**。本轮下注 `W-SYNC`。
CONTROLS:上限 同题男Δ×女Δ · 安慰剂 打乱年份顺序 · 硬规则 1 逐题打印**实际被问的年份与 n**
KILL(条件式):if 上限 > 安慰剂 q95: 按三分判 else UNVERIFIED
IMPOSSIBLE:青少年 ⇒ 不可外推成人 · 校内抽样 ⇒ 辍学者缺席 · 未加权(SADC 有权重,本轮不加权)
  ⇒ **任何比例都不是人群估计** · 行为≠态度 ⇒ 与 `#532` 不可直接并列 · [unchallenged]
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
NEED = ["year", "sex", "q56", "q57", "q58", "q59", "q61"]
assert all(n in pos for n in NEED), [n for n in NEED if n not in pos]
files = sorted(Y.glob("sadc_2023_state_*.dat"))          # 州文件,避开 405 MB 的学区文件(P2)
print(f"=== 硬规则 1:读 {len(files)} 个州文件,只取 {len(NEED)} 个字段 ===")
rec = {n: [] for n in NEED}
for f in files:
    for line in open(f, errors="replace"):
        for n in NEED:
            s, e = pos[n]; v = line[s:e].strip()
            rec[n].append(float(v) if v not in ("", ".") else np.nan)
A = {n: np.array(rec[n]) for n in NEED}
print(f"  行 {len(A['year']):,}")
for n in NEED[2:]:
    ok = np.isfinite(A[n]); ys = sorted(set(A["year"][ok].astype(int)))
    print(f"  {n} n={int(ok.sum()):7,} 年份 {len(ys)} 个 {ys[0]}–{ys[-1]}  {lab.get(n,'')[:34]}")

def series(q, sx):
    m = np.isfinite(A[q]) & (A["sex"] == sx)
    out = {}
    for y in sorted(set(A["year"][m].astype(int))):
        k = m & (A["year"].astype(int) == y)
        if k.sum() >= 300: out[y] = float((A[q][k] == 1).mean())    # 码 1 = 「是」
    return out

Q = ["q56", "q57", "q58", "q59", "q61"]
S = {(q, sx): series(q, sx) for q in Q for sx in (1, 2)}
def dcorr(a, b, perm=None):
    yr = sorted(set(a) & set(b))
    if len(yr) < 8: return np.nan
    va = np.array([a[y] for y in yr]); vb = np.array([b[y] for y in yr])
    if perm is not None: vb = vb[perm.permutation(len(yr))]
    da, db = np.diff(va), np.diff(vb)
    if np.std(da) == 0 or np.std(db) == 0: return np.nan
    return float(np.corrcoef(da, db)[0, 1])
same = [dcorr(S[(q, 1)], S[(q, 2)]) for q in Q]
same = [x for x in same if np.isfinite(x)]
cross = []
for a, b in itertools.combinations(Q, 2):
    v = [dcorr(S[(a, i)], S[(b, j)]) for i, j in ((1, 2), (2, 1))]
    v = [x for x in v if np.isfinite(x)]
    if v: cross.append(float(np.mean(v)))
shuf = []
for sd in SEEDS:
    r = np.random.default_rng(sd)
    for _ in range(80):
        for a, b in itertools.combinations(Q, 2):
            x = dcorr(S[(a, 1)], S[(b, 2)], perm=r)
            if np.isfinite(x): shuf.append(abs(x))
SM, CR = float(np.median(same)), float(np.median(cross))
Q95 = float(np.quantile(shuf, .95))
print(f"\n  **同题上限(男Δ×女Δ)中位 = {SM:+.4f}({len(same)} 题)**")
print(f"  **跨题目标(受访者不相交)中位 = {CR:+.4f}({len(cross)} 对)**  比值 {CR/SM:.3f}")
print(f"  安慰剂(打乱年份)中位 {np.median(shuf):.4f} · q95 {Q95:.4f}")
print("\n" + "=" * 74)
if abs(SM) > Q95:
    world = ("W-SYNC" if CR >= 0.5 * SM else ("W-ONE-BY-ONE" if abs(CR) < Q95 else "W-BLIND"))
    verdict = f"{world}: 跨题 {CR:+.4f} vs 上限 {SM:+.4f} vs 安慰剂 q95 {Q95:.4f}"
    print(f"控制齐备 ⇒ 评判。**{world}** —— {verdict}")
    print("⚠ 这个 KILL 会怎样失败:YRBS 问的是**行为**,GSS 问的是**态度**,"
          "**估计量不同** —— 本轮不能写成「`#532` 被复制了」,只能写成「行为上也如此/不如此」。")
else:
    world, verdict = "UNVERIFIED", f"上限 {SM:.4f} 未越过安慰剂 q95 {Q95:.4f}"
    print(f"⚠ {verdict}")
json.dump(dict(same_median=SM, cross_median=CR, ratio=CR / SM if SM else None,
               placebo_median=float(np.median(shuf)), placebo_q95=Q95, k_same=len(same),
               k_cross=len(cross), world=world, verdict=verdict, seeds=SEEDS,
               items={q: lab.get(q, "") for q in Q},
               inclusion=["YRBS SADC 州文件(避开 405MB 学区文件)", "每年每性别 n>=300 才计入",
                          "码 1 = 是", "未加权"],
               impossible=["青少年不可外推成人", "校内抽样,辍学者缺席", "未加权,非人群估计",
                           "行为≠态度,与 #532 不可直接并列"], unchallenged=True),
          open(OUT / "yrbs_synchrony.json", "w"), indent=1)
print(f"\nwrote {OUT/'yrbs_synchrony.json'}")

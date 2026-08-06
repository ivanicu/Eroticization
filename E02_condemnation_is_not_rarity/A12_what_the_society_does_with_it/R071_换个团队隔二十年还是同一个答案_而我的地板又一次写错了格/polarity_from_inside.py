"""E02·A12·R675 —— 从数据内部把那条 11 级量表的方向定出来

`#638` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。
`#638d`:`SCCS453–456` 的码只写着 `(1 of 11)`,**方向从发布里读不出来**。

⚠ **两处必须写在前面的偏离/更正:**
① **预注册的门槛单位错了。** `#638` 的 NEXT 要求把三组中位数之差与 `#530` 的 **0.187 ± 0.152** 比 ——
   **那是一个相关系数**,而我的量是**中位数之差**。**单位对不上**,`#653` 抓过同一类错。
   ⇒ 改用**该估计量自身的块 bootstrap 地板**,并标注。
② **而预注册的那个锚根本用不了,理由是本会话第二次同型。**
   `SCCS619` 的三组是 **2 / 59 / 7** —— 判别性对照落在 **2 对 7** 个社会上。
   联合 `n = 68` **越过了我写的地板 30**,而**杀死它的是锚变量几乎是常数**。
   **`#635` 已经记过一次这个缺陷,而我这一次仍然把地板写在了联合 n 上。**
   > **停止规则的地板要写在「对照里最小的那一格」上,不是写在联合 n 上。**

⇒ 改用一个**不依赖组间平衡**的锚,而它恰好是**跨编码团队**的(硬规则④):
   `SCCS1766`(Lang 1998)量的是**同一个构念** ——「晚童年男孩的体罚」——
   而它的码**完全自述且单调**:`10` 不体罚 · `21` 少罚 · `22` 常罚(`20` = 有罚但频率不明,**剔除**)。
   与 `SCCS455`(Barry 1977)联合 n = 69。

G1 ESTIMAND(先于方法):`ρ(SCCS455, SCCS1766_ordered)`,社会作单位,`1766 ∈ {10,21,22}`。
  **符号即答案**:`ρ > 0` ⇒ 11 级量表「数字越大越重」;`ρ < 0` ⇒ 「越大越轻」。
CONTROLS:
  正对照:`SCCS455 × SCCS453`(**同源、同构念、不同年龄**)必须强正 —— 仪器能看见一致。
  **g=0 / 安慰剂**:`SCCS455 × |纬度|` 必须 ≈ 0。
  ⚠ **offset 用不了同队基线**(`#530` 是**同一编码团队**的耦合,而本轮是**跨团队**)——
    如实记:**跨团队耦合的基线本项目没有测过**,所以只用「CI 是否排除零 + 安慰剂」两道。
KILL(条件式,预注册):
  if 正对照强正 and 安慰剂 ≈0:
      块 bootstrap 的 95% CI **排除零** -> **极性由符号确定**
      CI 含零 -> **极性仍记未定**
  else: UNVERIFIED
G3:三个量 + 每个的 n。G4:`1766` 是否保留 `20` 档两档对比。
IMPOSSIBLE(不写 planned):它定的是**方向**,不是**刻度** —— 11 级之间的间距仍未知 ·
  两个变量都出自民族志,**共享的是同一批民族志,不是同一批观察** · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from lib.gates import Gate

SEEDS = [20260806, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
S = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(S/"data.csv"); SOC = pd.read_csv(S/"societies.csv")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
BLK = {r.id: (int(np.floor(r.Lat/10)), int(np.floor(r.Long/10)))
       for r in SOC.dropna(subset=["Lat", "Long"]).itertuples()}
LAT = {r.id: abs(r.Lat) for r in SOC.dropna(subset=["Lat"]).itertuples()}

ORD = {10.0: 0, 21.0: 1, 22.0: 2}


def block_boot(df, fn, n=600):
    blocks = sorted({BLK.get(s, ("na", "na")) for s in df.index})
    bymap = {b: [s for s in df.index if BLK.get(s, ("na", "na")) == b] for b in blocks}
    out = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            pick = [bymap[blocks[i]] for i in rng.integers(0, len(blocks), len(blocks))]
            socs = [s for g in pick for s in g]
            v = fn(df.loc[socs])
            if np.isfinite(v): out.append(v)
    return np.array(out)


rows = []
def cell(name, a, b, mask=None, ordered=False):
    m = W[[a, b]].dropna()
    if ordered: m = m[m[b].isin(ORD)]; m = m.assign(**{b: m[b].map(ORD)})
    if mask is not None: m = m[mask(m)]
    if len(m) < 10: rows.append(dict(cell=name, n=len(m), rho=np.nan)); return None, m
    f = lambda d: float(spearmanr(d[a], d[b]).statistic)
    r = f(m); bs = block_boot(m, f)
    lo, hi = np.quantile(bs, [.025, .975])
    rows.append(dict(cell=name, n=int(len(m)), rho=r, lo=float(lo), hi=float(hi), sd=float(bs.std())))
    print(f"  {name:34s} n={len(m):3d}  ρ = **{r:+.4f}**  95%CI [{lo:+.4f}, {hi:+.4f}]")
    return r, m


print("=== G3 ===")
main, m_main = cell("主:455(Barry) × 1766(Lang,有序)", "SCCS455", "SCCS1766", ordered=True)
pos, _ = cell("正对照:455 × 453(同源同构念)", "SCCS455", "SCCS453")
lat = W[["SCCS455"]].dropna(); lat["_lat"] = [LAT.get(s, np.nan) for s in lat.index]
lat = lat.dropna()
f = lambda d: float(spearmanr(d["SCCS455"], d["_lat"]).statistic)
pl = f(lat)
print(f"  {'安慰剂:455 × |纬度|':34s} n={len(lat):3d}  ρ = **{pl:+.4f}**")
rows.append(dict(cell="安慰剂:455 × |纬度|", n=int(len(lat)), rho=pl))

G = Gate("从数据内部把那条 11 级量表的方向定出来")
pos_ok = G.positive_control("正对照:455 × 453 必须强正", planted=float(pos or 0), floor=0.30, spread=0.05)
pla_ok = G.negative_control("安慰剂:455 × |纬度| 必须 ≈0", null=float(abs(pl)), effect=float(abs(main or 0)),
                            null_spread=0.05, null_kind="与体罚无关的地理属性")
r0 = [x for x in rows if x["cell"].startswith("主")][0]
excl0 = np.isfinite(r0.get("lo", np.nan)) and (r0["lo"] * r0["hi"] > 0)
if pos_ok and pla_ok:
    verdict = ("**极性确定:11 级量表「数字越大 = 体罚越重」**" if excl0 and r0["rho"] > 0 else
               "**极性确定:数字越大 = 体罚越轻**" if excl0 and r0["rho"] < 0 else
               "**极性仍记未定 —— CI 含零**")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · 安慰剂 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4:1766 是否保留 `20`(有罚但频率不明)===")
ORD2 = {10.0: 0, 20.0: 1, 21.0: 1, 22.0: 2}
m2 = W[["SCCS455", "SCCS1766"]].dropna(); m2 = m2[m2.SCCS1766.isin(ORD2)]
m2 = m2.assign(SCCS1766=m2.SCCS1766.map(ORD2))
r2 = float(spearmanr(m2.SCCS455, m2.SCCS1766).statistic)
print(f"  剔除 20:n={r0['n']} ρ={r0['rho']:+.4f} · 保留 20(并入中档):n={len(m2)} ρ={r2:+.4f}")
json.dump(dict(cells=rows, spec_keep20=dict(n=int(len(m2)), rho=r2), verdict=verdict,
               anchor_note="SCCS619 三组 2/59/7 ⇒ 判别对照落在 2 对 7,弃用;改用跨团队的 SCCS1766",
               unchallenged=True),
          open(OUT/"polarity.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'polarity.json'}")

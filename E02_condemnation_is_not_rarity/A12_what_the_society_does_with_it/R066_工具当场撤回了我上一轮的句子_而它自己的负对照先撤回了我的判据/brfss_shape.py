"""E02·A224·R596 — 不需要语义就能知道的事:这份数据的形状

`#550` 的 NEXT ③。行动类型:**PRODUCTION**。
`#550b`:四级闸对 BRFSS **完全失明**,因为 `.XPT` 不携带变量标签。
**但「有多少行、多少列、缺失长什么样」不需要语义** —— 而它决定了这份数据**值不值得去配一份码本**。

⚠ `P2`:1.2 GB,**分块单遍**(每块 5 万行),只累计计数,不把整表放进内存。
⚠ 本轮**不解释任何一列的含义**(没有码本),**不做任何统计推断** —— 只报形状。
IMPOSSIBLE:无码本 ⇒ 不知道任何一列问的是什么 · 不知道缺失是「没问」还是「拒答」 ·
  BRFSS 有复杂抽样权重,**本轮不加权**,因此**任何比例都不是人群估计**,只是文件的性质 · [unchallenged]
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

F = ROOT / "data/external/brfss/LLCP2023.XPT"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
print(f"=== {F.name} {F.stat().st_size/1e9:.2f} GB —— 分块单遍,每块 5 万行(P2:不整读)===")
n_rows, cols, miss = 0, None, None
for i, ch in enumerate(pd.read_sas(F, format="xport", chunksize=50000, encoding="latin-1")):
    if cols is None:
        cols = list(ch.columns); miss = np.zeros(len(cols), dtype=np.int64)
    n_rows += len(ch)
    miss += ch.isna().sum().values.astype(np.int64)
    if i % 3 == 0: print(f"  …{n_rows:,} 行", flush=True)
rate = miss / n_rows
print(f"\n**行 {n_rows:,} · 列 {len(cols)}**")
q = np.quantile(rate, [0, .25, .5, .75, .9, 1])
print(f"逐列缺失率分位:min={q[0]:.4f} q25={q[1]:.4f} **中位={q[2]:.4f}** q75={q[3]:.4f} "
      f"q90={q[4]:.4f} max={q[5]:.4f}")
band = {"全满 (=0)": int((rate == 0).sum()), "<10%": int(((rate > 0) & (rate < .1)).sum()),
        "10–50%": int(((rate >= .1) & (rate < .5)).sum()), "50–90%": int(((rate >= .5) & (rate < .9)).sum()),
        "≥90%": int((rate >= .9).sum())}
print("逐列缺失率分档:" + " · ".join(f"{k} {v}" for k, v in band.items()))
usable = int((rate < .5).sum())
print(f"\n⇒ **{usable}/{len(cols)} 列的缺失率 <50%**;而 **{band['≥90%']} 列缺失 ≥90%** —— "
      f"后者多半是**分州模块**(只在部分州问),不是坏数据,但**不配码本就分不清这两件事**。")
top = sorted(zip(cols, rate), key=lambda x: -x[1])[:8]
print("缺失最多的 8 列(**列名不是含义**,本轮不解释):" + " · ".join(f"{c}={r:.3f}" for c, r in top))
json.dump(dict(n_rows=int(n_rows), n_cols=len(cols), miss_quantiles=[float(x) for x in q],
               bands=band, usable_lt50=usable,
               top_missing=[[str(c), float(r)] for c, r in top],
               inclusion=[f"{F.name} 全文件单遍", "每块 5 万行", "未加权"],
               impossible=["无码本,不知任何一列问什么", "不知缺失是没问还是拒答",
                           "未加权,任何比例都不是人群估计"], unchallenged=True),
          open(OUT / "brfss_shape.json", "w"), indent=1)
print(f"\nwrote {OUT/'brfss_shape.json'}")

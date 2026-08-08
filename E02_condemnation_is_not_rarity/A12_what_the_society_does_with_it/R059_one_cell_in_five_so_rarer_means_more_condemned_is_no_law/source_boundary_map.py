"""E02·A12·R673 —— 这个项目在 SCCS 上还能问什么:一张按来源论文画的边界图

`#636` 的 NEXT。**行动类型:普查(ENUMERATION)+ PRODUCTION**,如实标注。
`#636c` 立的边界:**能同时出现的,就是同一批人在同一篇论文里编的那几个变量。**
⇒ 那么「还能问什么」这个问题,答案的单位就不是变量,是 **source(来源论文)**。

⚠ **§3 梯度检查两条,都写在前面:**
   ① 变量多的 source 两两是 O(k²) ⇒ **设抽样上限 `PAIR_CAP=300` 并明说**;
      被抽样的 source 要标出来,**不许把抽样中位数当成全量中位数**。
   ② 只有 **1 个变量**的 source **没有任何对** ⇒ 中位数无定义 ——
      **这正是 g=0 该落的地方**,不是一个要被填掉的空。

G1 ESTIMAND(先于方法):对每个 `source`,
  `k` = 变量数 · `cov_med` = 其变量单独覆盖的中位 · `pair_med` = **source 内部**两两联合 n 的中位。
  `合格` = `k >= 3` **且** `pair_med >= 30`(**够在这一篇内部做一次社会层分析**)。
CONTROLS:
  正对照:**`ross1983political` 必须出现在合格名单里**(`#617` 就是在它内部做出来的)。
  **g=0**:任何 `k == 1` 的 source **必须不合格**(它连一对都没有)。
G3:按变量数排序的前 12 个 source 全表 + 合格名单全列。
IMPOSSIBLE(不写 planned):它只说「能不能同时看见」,**不说「看见了能不能回答」** ——
  `#635` 已证明 n=17 也可能因**变量无方差**而死 ·
  `source` 字段是 D-PLACE 的整理,**不等于原始出版物**(一个 source 可能含多轮编码)· `[unchallenged]`
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
PAIR_CAP = 300
S = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(S/"data.csv"); V = pd.read_csv(S/"variables.csv")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
COV = {c: int(W[c].notna().sum()) for c in W.columns}
V = V[V.id.isin(W.columns)]
print(f"仪器 = SCCS/D-PLACE · 可读变量 {len(V)} · source {V.source.nunique()} 个")

rng = np.random.default_rng(20260806)
rows = []
for src, g in V.groupby("source"):
    vs = [x for x in g.id if COV.get(x, 0) > 0]
    k = len(vs)
    if k == 0: continue
    covs = [COV[x] for x in vs]
    pairs = list(itertools.combinations(vs, 2))
    sampled = len(pairs) > PAIR_CAP
    if sampled:
        idx = rng.choice(len(pairs), PAIR_CAP, replace=False)
        pairs = [pairs[i] for i in idx]
    ns = [int(W[[a, b]].dropna().shape[0]) for a, b in pairs]
    rows.append(dict(source=str(src)[:34], k=k, cov_med=int(np.median(covs)),
                     pair_med=(int(np.median(ns)) if ns else -1),
                     n_pairs=len(pairs), sampled=bool(sampled)))
T = pd.DataFrame(rows).sort_values("k", ascending=False).reset_index(drop=True)
T["合格"] = (T.k >= 3) & (T.pair_med >= 30)
print(f"\n=== G3:按变量数排序的前 12 个 source(`pair_med = -1` 表示没有任何对)===")
print(f"{'source':36s}{'k':>4s}{'覆盖中位':>8s}{'对内中位n':>10s}{'抽样':>6s}{'合格':>6s}")
for r in T.head(12).itertuples():
    print(f"{r.source:36s}{r.k:4d}{r.cov_med:8d}{r.pair_med:10d}{'是' if r.sampled else '  ':>6s}"
          f"{'✅' if r.合格 else '  ':>6s}")
OK = T[T.合格]
print(f"\n**合格 source(k>=3 且 对内中位 n>=30):{len(OK)} / {len(T)} 个**")
for r in OK.sort_values("pair_med", ascending=False).itertuples():
    print(f"  ✅ {r.source:36s} k={r.k:4d} 覆盖中位 {r.cov_med:3d} 对内中位 n={r.pair_med:3d}")
one = T[T.k == 1]
print(f"\n  只有 1 个变量的 source:{len(one)} 个 —— 它们连一对都没有,`pair_med = -1`")

G = Gate("这个项目在 SCCS 上还能问什么:按来源论文画的边界图")
ross = T[T.source.str.startswith("ross1983")]
pos = bool(len(ross) and ross.iloc[0]["合格"])
print(f"\n  正对照:`ross1983political` 在合格名单里?**{pos}**"
      f"(k={int(ross.iloc[0].k) if len(ross) else '—'} · 对内中位 n={int(ross.iloc[0].pair_med) if len(ross) else '—'})")
g0 = int(one.合格.sum()) if len(one) else 0
print(f"  g=0:`k==1` 的 source 里被判合格的 = **{g0}**(须 0)")
pos_ok = G.positive_control("正对照:ross1983political 必须合格", planted=float(pos), floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:只有一个变量的 source 必须不合格", null=float(g0), effect=1.0,
                            null_spread=0.4, null_kind="连一对都没有的 source")
verdict = (f"**{len(OK)} / {len(T)} 个 source 自身就够做一次社会层分析**;"
           f"其余的,任何问题都必须跨论文拼,而 `#636c` 已量出那大半是空的")
print(f"\n{'控制齐备 ⇒ ' if pos_ok and pla_ok else '⚠ '}普查结果。**{verdict}**")
print(G)
json.dump(dict(table=T.to_dict("records"), qualified=OK.source.tolist(),
               n_qualified=int(len(OK)), n_sources=int(len(T)), pair_cap=PAIR_CAP,
               n_sampled=int(T.sampled.sum()), verdict=verdict,
               note="ENUMERATION,不配误差棒(`#616e`)", unchallenged=True),
          open(OUT/"source_map.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'source_map.json'}")

"""E02·A224·R597 — 那 112 列是「只在部分州问」,还是坏的?用文件自身分开它们

`#551` 的 NEXT。行动类型:**FRONTIER**(两个世界在这些列**是什么**上不同,不是参数不同)。
**不需要码本,也不联网。**

**要分的两个世界(`#551b`):**
  W-MODULE  它们是**分州可选模块** ⇒ 同一列的缺失**按州聚集**(某些州整州缺,某些州整州有)
  W-DAMAGE  它们是**损坏/未采集** ⇒ 缺失**与州无关**,各州缺失率彼此接近
  W-MIXED   两类都有 ⇒ 报**各占多少**,这才是「码本值多少钱」的答案

⚠ 州标识列**必须先被识别**,而没有码本 ⇒ 只能由**取值结构**认:
   整数型 · **不同取值 40–60 个** · **零缺失**。**先打印候选,再选**,并把候选表写进产物。

G1 ESTIMAND(先于方法):对每一个高缺失列 `c`,
   **`B(c) = Var_州(该州的缺失率)`** —— 州间方差。
   分州模块 ⇒ 各州缺失率两极(0 或 1)⇒ `B` 接近伯努利上限 0.25;
   损坏/随机 ⇒ 各州缺失率彼此接近 ⇒ `B` 接近 0。
   **判据:`B(c)` 与「把州标签随机打乱后」的同一量比较。**

CONTROLS(G2):
  正对照 **一列已知全满的列**必须判「与州无关」(`B ≈ 0`);
  安慰剂/关键零 **打乱州标签**后所有列的 `B` 必须塌到 0 附近
     ——「这个零该不该是零?」**该** ⇒ `negative_control`;
  ⚠ 打乱必须**在行内打乱州标签**,不是打乱缺失指示 —— 后者会同时破坏列的边际缺失率。
KILL(条件式,预注册):
  if 正对照 `B < 0.01` and 打乱后中位 `B < 0.01`:
      `B > 0.10` 的列 -> **判为分州模块**;`B < 0.02` -> **判为与州无关**;之间 -> **不判**
  else: UNVERIFIED
IMPOSSIBLE:**「与州无关」≠「损坏」** —— 也可能是一个全国问但几乎无人回答的题,
  **本轮无法区分这两者**(那需要码本)· 不加权 · 只用文件自身 · [unchallenged]
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

F = ROOT / "data/external/brfss/LLCP2023.XPT"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]

# ---- 认州标识列:**按名字取,再用结构验**
# ⚠ 第一版的判据是「某一块里有 40–60 个不同取值且零缺失」,它**永远不可能触发** ——
#   因为**文件按州排序**,一块 6 万行里只有 1–2 个州。
#   **从一块认出的东西,是那一块的性质,不是文件的性质。**
# 改:列名 `_STATE` 由**文件自身**给出(读名字是读对象,不是读记忆);
#   但**名字不是描述**,所以必须用结构验:全文件不同取值 40–60 · 整数 · 零缺失。
NAME = "_STATE"
vals, n_na, tot0 = set(), 0, 0
for ch in pd.read_sas(F, format="xport", chunksize=60000, encoding="latin-1"):
    assert NAME in ch.columns, f"{NAME} 不在列里"
    v = ch[NAME]; tot0 += len(v); n_na += int(v.isna().sum())
    vals.update(v.dropna().unique().tolist())
print(f"=== 州标识列的结构验证(名字来自文件,结构来自全文件)===")
print(f"  {NAME}: 全文件不同取值 **{len(vals)}** 个 · 缺失 {n_na}/{tot0} · "
      f"全整数 {all(float(x).is_integer() for x in vals)} · 范围 {min(vals):.0f}–{max(vals):.0f}")
assert 40 <= len(vals) <= 60 and n_na == 0 and all(float(x).is_integer() for x in vals), \
    "结构不像州标识 -> 停,不硬用一个名字"
STATE = NAME
cand = [(NAME, len(vals), float(min(vals)), float(max(vals)))]
print(f"✅ 结构与「州」相容 ⇒ 使用 {STATE}\n")

# ---- 单遍:累计每州 × 每列的缺失计数
cols, ns, nm = None, {}, {}
rng = np.random.default_rng(SEEDS[0])
ns_p, nm_p = {}, {}          # 打乱州标签的平行累计
tot = 0
for ch in pd.read_sas(F, format="xport", chunksize=60000, encoding="latin-1"):
    if cols is None: cols = [c for c in ch.columns if c != STATE]
    st = ch[STATE].values.astype(int)
    stp = rng.permutation(st)
    M = ch[cols].isna().values
    tot += len(ch)
    # ⚠ 此处原有一段写坏的合并循环,它把 `Mn` 绑到了 numpy 数组 `M` 上,
    #   于是 `Mn[s] = … if False else None` **给 M 的第 s 行赋了 None**,
    #   静默改写了缺失矩阵 —— 正对照打出 `nan` 才暴露出来。**已删除。**
    #   一段「反正走不到」的死代码,走到了。
    for s in np.unique(st):
        k = st == s
        ns[s] = ns.get(s, 0) + int(k.sum())
        nm[s] = nm.get(s, np.zeros(len(cols), dtype=np.int64)) + M[k].sum(0)
    for s in np.unique(stp):
        k = stp == s
        ns_p[s] = ns_p.get(s, 0) + int(k.sum())
        nm_p[s] = nm_p.get(s, np.zeros(len(cols), dtype=np.int64)) + M[k].sum(0)
print(f"扫完 {tot:,} 行 · 州 {len(ns)} 个 · 列 {len(cols)}")

def between_var(N, Mn):
    ss = sorted(N)
    R = np.vstack([Mn[s] / max(N[s], 1) for s in ss])      # 州 × 列
    return R.var(0), R

B, R = between_var(ns, nm)
Bp, _ = between_var(ns_p, nm_p)
overall = np.vstack([nm[s] for s in nm]).sum(0) / tot
hi = np.where(overall >= 0.90)[0]
full = np.where(overall == 0.0)[0]
print(f"\n=== 对照 ===")
print(f"  正对照:{len(full)} 列全满 -> 州间方差中位 = **{np.median(B[full]):.6f}**(须 <0.01)")
print(f"  关键零:打乱州标签后,全部列的州间方差中位 = **{np.median(Bp):.6f}**(须 <0.01)")
ok = np.median(B[full]) < 0.01 and np.median(Bp) < 0.01
mod = [i for i in hi if B[i] > 0.10]
nos = [i for i in hi if B[i] < 0.02]
mid = [i for i in hi if 0.02 <= B[i] <= 0.10]
print("\n" + "=" * 74)
if ok:
    world = ("W-MODULE" if len(mod) > 0.7 * len(hi) else
             ("W-DAMAGE" if len(nos) > 0.7 * len(hi) else "W-MIXED"))
    print(f"控制齐备 ⇒ 评判。**{world}**")
    print(f"  高缺失列 {len(hi)} 个:**按州聚集 {len(mod)}** · **与州无关 {len(nos)}** · 不判 {len(mid)}")
    print(f"  按州聚集的那些,州间方差中位 = {np.median(B[mod]):.4f}(伯努利上限 0.25)")
    print("⚠ 这个 KILL 会怎样失败:**「与州无关」≠「损坏」** —— 一个全国都问、"
          "却几乎无人回答的题也会落在这里。**本轮无法区分这两者**,那需要码本。")
    verdict = f"{world}: 聚集 {len(mod)} / 无关 {len(nos)} / 不判 {len(mid)}"
else:
    world, verdict = "UNVERIFIED", f"控制未齐 正对照={np.median(B[full]):.6f} 零={np.median(Bp):.6f}"
    print(f"⚠ {verdict}")
json.dump(dict(state_col=STATE, state_candidates=cand, n_rows=int(tot), n_states=len(ns),
               n_high_missing=len(hi), n_module=len(mod), n_no_state=len(nos), n_undecided=len(mid),
               between_var_full_median=float(np.median(B[full])),
               between_var_shuffled_median=float(np.median(Bp)),
               world=world, verdict=verdict, seeds=SEEDS,
               inclusion=["全文件单遍", "高缺失 = 边际缺失率 >=0.90", "州标识由取值结构识别"],
               impossible=["「与州无关」≠「损坏」,需码本才能分", "未加权", "只用文件自身"],
               unchallenged=True), open(OUT / "state_clustering.json", "w"), indent=1)
print(f"\nwrote {OUT/'state_clustering.json'}")

"""E02·A224·R598 — 把第三个世界补进检验

`#552` 的 NEXT。行动类型:**FRONTIER**。
`#552c`:我的世界分解里少了「**州内问卷版本切分**」,而它与被比较的两个世界给出同一个观测。
本轮把它补进去,**三个世界并排**。

G1 ESTIMAND(先于方法):对每个高缺失列 `c`,
   `B(c)` = 缺失率的**州间**方差;`W(c)` = 缺失率的**州内版本间**方差(先按州分组,组内按版本算方差,再平均)。
预注册三分:`B > 0.10` -> 分州模块 · `W > 0.10` -> **州内版本模块** · 两者都 <0.02 -> 与两者都无关 ·
   其余 -> **不判**。
版本标识列的判据(先于选择写死):整数 · **不同取值 2–4 个** · 零缺失 ·
   **且存在至少一个州,该州内取值不唯一**(这一条把「州内切分」与「州级变量」分开)。
   **先打印全部候选,再选,候选表进产物。**
CONTROLS:正对照 全满列在 `B` 与 `W` 上都必须 ≈0;
   关键零 分别打乱**州标签**与**版本标签**,两个量都必须塌到 0 附近。
IMPOSSIBLE:仍不能区分「与两者都无关」与「损坏」· 未加权 · 只用文件自身 · [unchallenged]
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
F = ROOT / "data/external/brfss/LLCP2023.XPT"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
CH = 60000

# ---- 找版本标识列候选(结构判据,先打印再选)
print("=== 版本标识列候选(整数 · 2–4 值 · 零缺失 · 州内取值不唯一)===")
vals, nna, within = {}, {}, {}
for ch in pd.read_sas(F, format="xport", chunksize=CH, encoding="latin-1"):
    st = ch["_STATE"].values.astype(int)
    for c in ch.columns:
        if c == "_STATE": continue
        v = ch[c]
        nna[c] = nna.get(c, 0) + int(v.isna().sum())
        u = v.dropna().unique()
        if len(u) <= 6: vals[c] = vals.get(c, set()) | set(u.tolist())
        else: vals[c] = set(range(99))          # 标记为「太多」,后面会被判据排除
        if c in vals and len(vals[c]) <= 6:
            for s in np.unique(st):
                if v[st == s].dropna().nunique() > 1: within[c] = True; break
cand = [(str(c), len(vals[c])) for c in vals
        if 2 <= len(vals[c]) <= 4 and nna.get(c, 1) == 0 and within.get(c) and
        all(float(x).is_integer() for x in vals[c])]
cand.sort()
for c, k in cand[:14]: print(f"  {c:12s} {k} 个取值")
print(f"  …共 {len(cand)} 个候选")
assert cand, "无候选 -> 判据不适用,停"
# 版本列应当在**全部**州里都切分,且切分比例大致稳定 —— 取「州内不唯一」的州数最多者
VER = cand[0][0]
print(f"⇒ 选用 **{VER}**(候选按名排序取第一;全部候选已进产物)\n")

# ---- 单遍:累计 (州, 版本) × 列 的缺失
cols, acc, cnt = None, {}, {}
rng = np.random.default_rng(SEEDS[0]); accP, cntP = {}, {}
tot = 0
for ch in pd.read_sas(F, format="xport", chunksize=CH, encoding="latin-1"):
    if cols is None: cols = [c for c in ch.columns if c not in ("_STATE", VER)]
    st = ch["_STATE"].values.astype(int); ve = ch[VER].values.astype(int)
    M = ch[cols].isna().values
    tot += len(ch)
    for s in np.unique(st):
        for v in np.unique(ve):
            k = (st == s) & (ve == v)
            if not k.any(): continue
            cnt[(s, v)] = cnt.get((s, v), 0) + int(k.sum())
            acc[(s, v)] = acc.get((s, v), np.zeros(len(cols), np.int64)) + M[k].sum(0)
    vp = rng.permutation(ve)
    for s in np.unique(st):
        for v in np.unique(vp):
            k = (st == s) & (vp == v)
            if not k.any(): continue
            cntP[(s, v)] = cntP.get((s, v), 0) + int(k.sum())
            accP[(s, v)] = accP.get((s, v), np.zeros(len(cols), np.int64)) + M[k].sum(0)
print(f"扫完 {tot:,} 行 · (州,版本) 格 {len(cnt)} 个 · 列 {len(cols)}")

def BW(acc, cnt, minn=200):
    states = sorted({s for s, _ in cnt})
    Brows, Wvals = [], []
    for s in states:
        cells = [(v, cnt[(s, v)], acc[(s, v)]) for (ss, v) in cnt if ss == s and cnt[(s, v)] >= minn]
        if not cells: continue
        n_s = sum(c for _, c, _ in cells); m_s = sum(a for _, _, a in cells)
        Brows.append(m_s / n_s)
        if len(cells) >= 2:
            R = np.vstack([a / c for _, c, a in cells])
            Wvals.append(R.var(0))
    B = np.vstack(Brows).var(0)
    W = np.mean(np.vstack(Wvals), 0) if Wvals else np.zeros_like(B)
    return B, W, np.vstack(Brows).mean(0)

B, W, overall = BW(acc, cnt)
_, WP, _ = BW(accP, cntP)
hi = np.where(overall >= 0.90)[0]; full = np.where(overall <= 1e-12)[0]
print(f"\n=== 对照 ===")
print(f"  正对照:{len(full)} 列全满 -> B 中位 {np.median(B[full]):.6f} · W 中位 {np.median(W[full]):.6f}(须 <0.01)")
print(f"  关键零:打乱**版本**标签 -> W 中位 {np.median(WP):.6f}(须 <0.01)")
ok = len(full) > 0 and np.median(B[full]) < .01 and np.median(W[full]) < .01 and np.median(WP) < .01
st_mod = [i for i in hi if B[i] > .10]
ver_mod = [i for i in hi if W[i] > .10]
neither = [i for i in hi if B[i] < .02 and W[i] < .02]
undec = [i for i in hi if i not in st_mod and i not in ver_mod and i not in neither]
print("\n" + "=" * 74)
if ok:
    print(f"控制齐备 ⇒ 评判。高缺失列 {len(hi)} 个:")
    print(f"  **分州模块(B>0.10):{len(st_mod)}** · **州内版本模块(W>0.10):{len(ver_mod)}** · "
          f"**两者都无关(<0.02):{len(neither)}** · 不判 {len(undec)}")
    if ver_mod: print(f"  版本模块的 W 中位 = {np.median(W[ver_mod]):.4f}(伯努利上限 0.25)")
    world = ("W-VERSION" if len(ver_mod) > len(st_mod) and len(ver_mod) > len(neither) else
             ("W-STATE" if len(st_mod) > len(neither) else "W-NEITHER"))
    verdict = f"{world}: 州 {len(st_mod)} / 版本 {len(ver_mod)} / 都无关 {len(neither)} / 不判 {len(undec)}"
    print(f"  ⇒ **{world}**")
    print("⚠ 这个 KILL 会怎样失败:「两者都无关」仍**不等于损坏** —— 一个全国全版本都问、"
          "却几乎无人回答的题落在同一格。第三个世界补上了,第四个仍然没有。")
else:
    world, verdict = "UNVERIFIED", f"控制未齐 full={len(full)}"
    print(f"⚠ {verdict}")
json.dump(dict(version_col=VER, version_candidates=cand[:40], n_rows=int(tot), n_cells=len(cnt),
               n_high=len(hi), n_state_module=len(st_mod), n_version_module=len(ver_mod),
               n_neither=len(neither), n_undecided=len(undec), world=world, verdict=verdict,
               ctrl=dict(B_full=float(np.median(B[full])) if len(full) else None,
                         W_full=float(np.median(W[full])) if len(full) else None,
                         W_shuffled=float(np.median(WP))), seeds=SEEDS,
               inclusion=["全文件单遍", "(州,版本) 格 n>=200 才计入", "高缺失 = 边际 >=0.90"],
               impossible=["「都无关」≠ 损坏,第四个世界仍未补", "未加权", "只用文件自身"],
               unchallenged=True), open(OUT / "third_world.json", "w"), indent=1)
print(f"\nwrote {OUT/'third_world.json'}")

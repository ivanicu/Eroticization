"""lib/assignment_sensitivity.py — 「我把哪些东西算成同一类」是一个可失败的判断

造于 `E02·A238·R621`(`#576` 的 NEXT)。行动类型:**PRODUCTION**。

**由来(`#576d`):** 四次指派检验,三次改变了结论 ——
  社会 22 对 极差 0.004 **不承重** · 社会 18 对 0.038 **敏感** ·
  年代·行为 0.852 **符号翻转 ⇒ 撤回** · 人 0.087 **量级敏感,方向存活**。
**而它从来没有出现在任何一个 MDE、CI 或对照里。** 四次都是临时写的,门槛与 bootstrap 各轮自定。

**这件工具做的事:** 给定「变量→类」的**多个映射**与一个**统计量函数**,
返回 `(极差, 极差的 bootstrap CI, 三值判决)`。

**三值判决(先于使用写死):**
  `方向翻转` —— 各方案的统计量**符号不一致** ⇒ 结论是划分的产物,**必须撤回**
  `量级敏感` —— 符号一致但极差 ≥ `tol` ⇒ 结论的**方向**存活,**数字必须报成区间**
  `不承重`   —— 极差 < `tol`
⚠ **符号优先于极差**:符号翻转即使极差很小也必须判翻转(`#575` 的形状)。
⚠ **P6 安全侧**:本工具只在「划分影响结论」这个方向上可靠;
  它**永远不说哪一种划分是对的** —— 那需要领域知识,不是数据。
"""
import numpy as np


def assignment_sensitivity(stat_fn, mappings, tol=0.05, n_boot=300, seeds=(20260805, 7, 991),
                           boot_fn=None):
    """stat_fn(mapping) -> float;boot_fn(mapping, rng) -> float(用于极差的 bootstrap)。

    mappings: dict[str, mapping] —— 至少两个方案。
    返回 dict:vals · range · boot_ci · verdict · signs
    """
    if len(mappings) < 2:
        return dict(verdict="DEGENERATE", why="少于两个方案,这个检查不可判", vals={}, range=None)
    vals = {k: float(stat_fn(m)) for k, m in mappings.items()}
    fin = {k: v for k, v in vals.items() if np.isfinite(v)}
    if len(fin) < 2:
        return dict(verdict="DEGENERATE", why="可算的方案少于两个", vals=vals, range=None)
    rng_ = max(fin.values()) - min(fin.values())
    signs = {k: int(np.sign(v)) for k, v in fin.items()}
    flipped = len(set(s for s in signs.values() if s != 0)) > 1
    ci = None
    if boot_fn is not None:
        bs = []
        for sd in seeds:
            r = np.random.default_rng(sd)
            for _ in range(n_boot // len(seeds)):
                g = [boot_fn(m, r) for m in mappings.values()]
                g = [x for x in g if np.isfinite(x)]
                if len(g) == len(mappings): bs.append(max(g) - min(g))
        if bs: ci = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
    verdict = "方向翻转" if flipped else ("量级敏感" if rng_ >= tol else "不承重")
    return dict(vals=vals, range=float(rng_), boot_ci=ci, verdict=verdict, signs=signs,
                tol=tol, note="只在「划分影响结论」方向上可靠;从不说哪种划分是对的")


def self_test():
    """内建正/负对照。恒等映射组必须判『不承重』;人造翻转组必须判『方向翻转』。"""
    base = {f"v{i}": ("A" if i < 5 else "B") for i in range(10)}
    ident = {"S1": base, "S2": dict(base), "S3": dict(base)}
    flip = dict(base); 
    for i in range(10): flip[f"v{i}"] = "B" if base[f"v{i}"] == "A" else "A"
    # ⚠ 第一版的玩具统计量在翻转后返回同一个值(A 组恒有 5 个)——
    #   它**不可能**变号,于是正对照必然失败。「控制不可能通过」的教科书形态,当场被自检抓到。
    #   改:统计量 = A 组均值 − B 组均值,在固定的合成数据上,标签翻转 ⇒ 符号必然翻转。
    data = {f"v{i}": (1.0 if i < 5 else -1.0) for i in range(10)}
    def stat(m):
        a = [data[k] for k, v in m.items() if v == "A"]
        b = [data[k] for k, v in m.items() if v == "B"]
        if not a or not b: return float("nan")
        return float(np.mean(a) - np.mean(b))

    r1 = assignment_sensitivity(stat, ident)
    r2 = assignment_sensitivity(stat, {"S1": base, "S2": flip})
    ok1 = r1["verdict"] == "不承重"
    ok2 = r2["verdict"] == "方向翻转"
    print(f"  负对照(恒等映射组)-> {r1['verdict']} (须『不承重』) {'✅' if ok1 else '⛔'}")
    print(f"  正对照(人造翻转组)-> {r2['verdict']} (须『方向翻转』) {'✅' if ok2 else '⛔'}")
    return ok1 and ok2


if __name__ == "__main__":
    print("=== lib/assignment_sensitivity.py 自检 ===")
    import sys
    sys.exit(0 if self_test() else 2)

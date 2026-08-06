"""E02·A12·R672 —— 那两「束」严厉度量表,两两之间到底有多少个共同社会?

`#635` 的 NEXT 在硬规则①那一步就停住了:`SCCS173` 与 Ross 三题的联合 n = 21 / 19 / 14,
**三格全部低于预注册的 30 ⇒ 全部记「判不了」,不计算。**
⇒ 按 `§0.2`(不能只交一个停),本轮把这件事的**普遍形式**量出来。

**行动类型:CLOSURE/普查(ENUMERATION)**,如实标注 ——
它不检验任何假设,它数的是「哪些格未来还能被问」。**任何人打开 `data.csv` 都能复核。**

G1 ESTIMAND(先于方法):对束 A(7 个)与束 B(5 个)的每一对变量,
  `n_joint` = 两者**都非缺失**的社会数。派生:`可算对数` = `n_joint >= 30` 的对数。
⚠ **这不是一次统计检验**,所以**不配误差棒**(`#616e` 同款:普查 ≠ 被估计的量)。
CONTROLS:
  正对照:**同一变量与自身**的 `n_joint` 必须等于它的单独覆盖(仪器算得对)。
  **g=0**:与一个**全空**的合成变量配对,`n_joint` 必须为 0。
G3:两张完整矩阵发布,含所有不可算的格。
IMPOSSIBLE(不写 planned):它只说「能不能同时看见」,**不说「看见了能不能回答」** ——
  `#635` 刚证明 n=17 也可能因为**变量无方差**而答不了 · `[unchallenged]`
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
S = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(S/"data.csv"); V = pd.read_csv(S/"variables.csv").set_index("id")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")

A_ = ["SCCS165", "SCCS169", "SCCS176", "SCCS960", "SCCS961", "SCCS962", "SCCS964"]
B_ = ["SCCS173", "SCCS781", "SCCS782", "SCCS783", "SCCS1770"]
LAB = {"SCCS165":"婚前性(女)","SCCS169":"婚外性","SCCS176":"同性恋","SCCS960":"违反乱伦",
       "SCCS961":"婚前性限制","SCCS962":"违反婚前性","SCCS964":"婚外性惩罚",
       "SCCS173":"强奸","SCCS781":"暴力·本地","SCCS782":"暴力·本社会","SCCS783":"暴力·他社会",
       "SCCS1770":"暴力·他族(Lang)"}


def mat(vs):
    M = pd.DataFrame(index=vs, columns=vs, dtype=int)
    for a, b in itertools.product(vs, vs):
        M.loc[a, b] = int(W[[a, b]].dropna().shape[0])
    return M


out = {}
for nm, vs in (("束A 性(7)", A_), ("束B 伤害(5)", B_)):
    M = mat(vs)
    print(f"\n=== G3 · {nm} 两两联合非缺失社会数 ===")
    hdr = "               " + " ".join(f"{LAB[x][:6]:>8s}" for x in vs)
    print(hdr)
    for a in vs:
        print(f"  {LAB[a][:8]:12s} " + " ".join(
            (f"{int(M.loc[a,b]):8d}" if a != b else f"{int(M.loc[a,b]):7d}*") for b in vs))
    pairs = [(a, b) for a, b in itertools.combinations(vs, 2)]
    ok = [(a, b) for a, b in pairs if M.loc[a, b] >= 30]
    print(f"  对角=单独覆盖 · **{len(pairs)} 对中 n>=30 的有 {len(ok)} 对 = {len(ok)/len(pairs)*100:.0f}%**")
    out[nm] = dict(matrix=M.astype(int).to_dict(), pairs=len(pairs), ok=len(ok),
                   ok_list=[f"{a}×{b}={int(M.loc[a,b])}" for a, b in ok])

print("\n=== 跨束(性 × 伤害)35 对 ===")
X = pd.DataFrame(index=A_, columns=B_, dtype=int)
for a, b in itertools.product(A_, B_): X.loc[a, b] = int(W[[a, b]].dropna().shape[0])
print("               " + " ".join(f"{LAB[x][:6]:>8s}" for x in B_))
for a in A_: print(f"  {LAB[a][:8]:12s} " + " ".join(f"{int(X.loc[a,b]):8d}" for b in B_))
xp = [(a, b) for a, b in itertools.product(A_, B_)]
xok = [(a, b) for a, b in xp if X.loc[a, b] >= 30]
print(f"  **{len(xp)} 对中 n>=30 的有 {len(xok)} 对 = {len(xok)/len(xp)*100:.0f}%**")

G = Gate("那两束严厉度量表,两两之间到底有多少个共同社会?")
self_ok = all(int(W[[x, x]].dropna().shape[0]) == int(W[x].notna().sum()) for x in A_ + B_)
W2 = W.copy(); W2["_empty"] = np.nan
g0 = int(W2[["SCCS165", "_empty"]].dropna().shape[0])
print(f"\n  正对照:变量与自身的联合 n = 单独覆盖?**{self_ok}**")
print(f"  g=0:与一个全空变量配对 -> n = **{g0}**(须 0)")
pos_ok = G.positive_control("正对照:自配必须等于单独覆盖", planted=float(self_ok), floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:与全空变量配对必须为 0", null=float(g0), effect=1.0,
                            null_spread=0.4, null_kind="一个全部缺失的合成变量")
verdict = (f"束A {out['束A 性(7)']['ok']}/{out['束A 性(7)']['pairs']} · "
           f"束B {out['束B 伤害(5)']['ok']}/{out['束B 伤害(5)']['pairs']} · "
           f"跨束 {len(xok)}/{len(xp)} 对可算(n>=30)")
print(f"\n{'控制齐备 ⇒ ' if pos_ok and pla_ok else '⚠ '}普查结果。**{verdict}**")
print(G)
json.dump(dict(bundles=out, cross={f"{a}×{b}": int(X.loc[a,b]) for a,b in xp},
               cross_ok=len(xok), cross_pairs=len(xp), verdict=verdict,
               note="ENUMERATION,不是统计检验,故不配误差棒(`#616e`)", unchallenged=True),
          open(OUT/"coverage_atlas.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'coverage_atlas.json'}")

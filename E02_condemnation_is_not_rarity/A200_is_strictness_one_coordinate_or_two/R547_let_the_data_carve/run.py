"""E02·A200·R547 — 不要我划域,让数据划

`#501` 的 NEXT,**它指向我不愿意看到的结果**:若数据的划分与我的不一致,
`#501a` 必须降级为「按我这样划时如此」,页面要改写。

G1 ESTIMAND(方法与判据**写在跑之前**):
  ① 8 道规范的严格指标(0/1)-> 8×8 相关矩阵 -> 前两个主成分。
  ② **PC1 预期是「整体严格度」**(全部同号);**分域信息若存在,应在 PC2 的符号里。**
  ③ 按 **PC2 载荷符号**分两簇,与我的划分(性 3 / 家庭 5)**逐道比对**。
  判据:一致 **≥7/8** -> 数据同意;**≤5/8** -> 数据不同意;**6/8** -> 两可。
  ④ **载荷稳定性**:随机劈半重算 PC2 载荷,3 组种子;
     **载荷符号在半样本间不稳,则它划不动任何东西** -> 直接判「数据划不出域」。

WORLDS:
  W-AGREE    数据的划分与我的一致 -> `#501a` 的域效应不再依赖我的判断
  W-DISAGREE 不一致              -> **`#501a` 降级,页面改写**
  W-NOSTRUCT PC2 载荷不稳        -> **数据划不出域**,`#501a` 的划分是外加的
  | World      | now | ≥7/8 | ≤5/8 | 不稳 |
  | W-AGREE    | 0.4 | 0.85 | 0.05 | 0.05 |
  | W-DISAGREE | 0.3 | 0.05 | 0.85 | 0.10 |
  | W-NOSTRUCT | 0.3 | 0.05 | 0.10 | 0.85 |

⛔ STRONGEST CONFOUND,写在跑之前:8 个变量的 PCA **点太少**,PC2 很容易是噪声;
  ④ 的稳定性检查就是为它准备的,**不是事后补的**。
CONTROLS:正对照 = PC1 必须是「全部同号」(若不是,这个分解本身失败);
  阴性 = 用**打乱人的行内配对**造的矩阵,其 PC2 载荷稳定性作底(应当很低)。
KILL(条件式,预注册):
  if PC1 全部同号(分解可用) and 真实 PC2 稳定性 > 打乱底:
      一致 ≥7/8 -> W-AGREE;≤5/8 -> W-DISAGREE;6/8 -> 两可
  else: **W-NOSTRUCT / UNVERIFIED**
IMPOSSIBLE:8 个变量 ⇒ 成分结构本就脆弱 · 只有女性卷 · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEXN = ["samesex", "sxok16", "sxok18"]
FAMN = ["staytog", "chsuppor", "okcohab", "chcohab", "gayadopt"]
NORMS = SEXN + FAMN
MINE = {n: ("sex" if n in SEXN else "fam") for n in NORMS}
REVERSE = {"okcohab"}

def parse_dct(p):
    out = {}; pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out

LAY = parse_dct(NS / "setup" / "2011_2013_FemRespSetup.dct")
cols = {n: LAY[n] for n in NORMS}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
D = {n: np.array(buf[n]) for n in cols}
X = []
for n in NORMS:
    ok = np.isin(D[n], [1, 2, 3, 4]); sc = [1, 2] if n in REVERSE else [3, 4]
    X.append(np.where(ok, np.isin(D[n], sc).astype(float), np.nan))
X = np.column_stack(X)
keep = np.all(np.isfinite(X), axis=1)
X = X[keep]
print(f"完整作答 8 道的人 n = {len(X)}(共 {len(keep)} 行)")

def pcs(M):
    Z = (M - M.mean(0)) / M.std(0)
    C = np.corrcoef(Z.T)
    w, v = np.linalg.eigh(C)
    o = np.argsort(w)[::-1]
    return w[o], v[:, o], C

w, V, C = pcs(X)
print(f"\n特征值 = {np.round(w,3)}   前两个解释 {100*w[:2].sum()/w.sum():.1f}%")
pc1, pc2 = V[:, 0], V[:, 1]
if pc1.sum() < 0: pc1 = -pc1
print(f"\n{'规范':10s} {'我的域':>6s} {'PC1':>8s} {'PC2':>8s} {'数据簇':>7s}")
data_cluster = {}
for i, n in enumerate(NORMS):
    dc = "A" if pc2[i] >= 0 else "B"
    data_cluster[n] = dc
    print(f"{n:10s} {MINE[n]:>6s} {pc1[i]:+8.3f} {pc2[i]:+8.3f} {dc:>7s}")
pc1_same_sign = bool(np.all(pc1 > 0) or np.all(pc1 < 0))
print(f"\n正对照:PC1 全部同号 = {pc1_same_sign}(否则这个分解本身失败)")

# 一致度:两种簇标对应(A↔sex 或 A↔fam),取更好的
def agree(mapping):
    return sum(1 for n in NORMS if mapping[data_cluster[n]] == MINE[n])
a1 = agree({"A": "sex", "B": "fam"}); a2 = agree({"A": "fam", "B": "sex"})
best = max(a1, a2)
print(f"与我的划分一致度 = {best}/8(两种簇标对应取优:{a1} / {a2})")

# 稳定性:随机劈半的 PC2 载荷符号一致度 vs 打乱底
def half_stability(M, shuffle=False, B=200, seed=0):
    rng = np.random.default_rng(seed); out = []
    for _ in range(B):
        MM = M.copy()
        if shuffle:
            for j in range(MM.shape[1]): MM[:, j] = MM[rng.permutation(len(MM)), j]
        p = rng.permutation(len(MM)); h = len(MM) // 2
        try:
            _, V1, _ = pcs(MM[p[:h]]); _, V2, _ = pcs(MM[p[h:]])
        except Exception: continue
        s = np.sign(V1[:, 1]) * np.sign(V2[:, 1])
        out.append(max(np.mean(s > 0), np.mean(s < 0)))   # 符号一致度(允许整体翻号)
    return np.array(out)

real = np.concatenate([half_stability(X, False, 200, s) for s in SEEDS])
null = np.concatenate([half_stability(X, True, 200, s) for s in SEEDS])
print(f"\nPC2 载荷符号的劈半一致度:真实 = {real.mean():.3f} ± {real.std():.3f};"
      f"打乱底 = {null.mean():.3f} ± {null.std():.3f}")
stable = bool(real.mean() > np.quantile(null, .95))
print(f"  超过打乱底 q95({np.quantile(null,.95):.3f}) = {stable}")

G = Gate("不要我划域,让数据划")
G.asserted("正对照:PC1 全部同号(分解可用)", pc1_same_sign,
           f"PC1 载荷 {np.round(pc1,3).tolist()}", kind="control")
G.negative_control("阴性:打乱行内配对后的 PC2 稳定性(这个零应当低)",
                   null=float(null.mean()), effect=float(real.mean()),
                   null_spread=float(null.std()), null_kind="列内独立置换(破坏人内相关)")
print("\n" + "=" * 70)
if pc1_same_sign and stable:
    verdict = ("W-AGREE:数据的划分与我的一致" if best >= 7 else
               "**W-DISAGREE:数据不同意我的划分 -> `#501a` 降级为「按我这样划时如此」,页面改写**"
               if best <= 5 else "两可(6/8)-> 不足以支持也不足以推翻")
    print(f"控制齐备 ⇒ 评判。一致 {best}/8 -> {verdict}")
elif pc1_same_sign and not stable:
    verdict = "**W-NOSTRUCT:PC2 载荷在半样本间不稳 -> 数据划不出域,我的划分是外加的**"
    print(f"控制齐备 ⇒ 评判。{verdict}")
else:
    verdict = "UNVERIFIED —— PC1 不同号,分解本身失败"
    print(f"⚠ {verdict}")
print("⚠ 通过的 KILL 会怎样失败:8 个变量的 PCA 点太少,PC2 极易是噪声;"
      "而「一致度」用的是簇标对应取优,这**偏向一致**,是保守方向的反面 —— 已照实报两种对应。")
print(G)
json.dump(dict(eigen=w.tolist(), pc1=pc1.tolist(), pc2=pc2.tolist(),
               mine=MINE, data_cluster=data_cluster, agree=best, agree_both=[a1, a2],
               pc1_same_sign=pc1_same_sign, stab_real=float(real.mean()),
               stab_null=float(null.mean()), stab_null_q95=float(np.quantile(null, .95)),
               stable=stable, verdict=verdict, n=int(len(X)), seeds=SEEDS, unchallenged=True),
          open(OUT / "let_the_data_carve.json", "w"), indent=1)
print(f"\nwrote {OUT/'let_the_data_carve.json'}")

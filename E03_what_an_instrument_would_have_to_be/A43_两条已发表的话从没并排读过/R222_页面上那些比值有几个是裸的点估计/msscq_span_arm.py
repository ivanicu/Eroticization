"""#783 第二臂 · E03·A43·R222 —— 「规格跨度当误差」是我的毛病,还是那个设计的毛病?

第一臂在 GSS 上判了世界 C:页面公布的「区间」全是**规格跨度**,而同一批数的**自助宽度是它的 6.3 倍**。
⚠ **但那句话现在的口气是「页面一直在这么干」——一个关于我的报告习惯的普遍主张,
   而它只在一具仪器上量过。** `readme_gate` 的 `single_instrument` 正是对着这一条开的火,
   而 `/loop` 的硬规则④ 写的就是:**跨仪器复制胜过在同一具上再来一轮。**

G1 估计量(与第一臂**同一个**,只换仪器与单位):
   **`自助宽度 ÷ 规格跨度`**,算在一处**页面已发表的、非 GSS 的、带两种口径的数**上。
   对象:MSSCQ 面内一致度 —— 页面同时写过 **0.5645** 与 **0.5856**。
   **单位换了:GSS 那边重抽的是 28 个年份点,这边重抽的是 12,789 个人。**

⚠⚠ **第一版在这里当场被自己的正控打掉,而那一击本身是本臂最有价值的东西。**
   第一版**手搓**了「逐面中位再取中位」,得 0.5468,复现不了页面的 0.5645(差 0.0177 > 容差 0.01)
   ⇒ 判 **UNVERIFIED**,正确。查下去发现两件事,第二件比第一件重:
   ① 我重写了一个**项目里已经有的函数**(`lib.blocks.median_of_group_medians`,
      它的 docstring 甚至写着「`#721` 实测:0.5645 对 0.5856」)—— `P4` 的先例阶梯 L1 一层都没走。
   ② **而调用它也不给 0.5645。** 逐格枚举 2×2 才对上:
      **`0.5645` = 归一 ρ(有符号)+ 逐面中位;`0.5856` = 生 |ρ| + 200 对合并。**
   ⇒ **页面上并排写的这两个数,同时差了两个自由度(矩阵口径 × 汇总口径),不是一个规格跨度。**
   本臂因此改成**跑满 2×2**(`G4`:规格曲线不是一格),跨度取整张网格的极差。

⚠⚠ 两个世界,而这一次我**希望**输的是第一个(`frontier §3` 的 basin 逃逸:
   连着一轮判「我的报告习惯有问题」之后,故意去撞「其实是那个设计的问题」):
   C-general **我的报告习惯**:凡我写「跨度」的地方,抽样宽度都远大于它 ⇒ 这里也应 >2×。
   C-design  **那是 28 个年份点的毛病**:GSS 的比值有效 n 只有 28 个时间点,
             而 MSSCQ 有 12,789 个人 ⇒ 这里的自助宽度应当**小于或接近**规格跨度(≤2×)。

预测矩阵:
   | 世界 | 现在 | 若 >2× | 若 ≤2× |
   | C-general | 0.5 | 0.85 | 0.10 |
   | C-design  | 0.5 | 0.15 | 0.90 |
   两个分支都强,而**较差的那一支(≤2×)反而更改变结论** —— 它会把第一臂那句话从
   「页面一直在公布小的那种不确定性」缩成「**在有效 n 只有 28 的那一类量上**才如此」。

预注册判词(条件式):
   if 正控通过(全样本复现两个已发表值,容差 0.01)and 自助真的在变(宽 > 0):
       if 自助宽度 / 规格跨度 > 2 -> C-general:第一臂那句话按原样成立
       else                      -> C-design:第一臂那句话**必须挂上「在这一类设计上」的限定**
   else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**两个口径不是独立的两次测量,是同一批相关的两种汇总**
   ⇒ 它们的跨度天然小,自助宽度也天然相关。**所以这里比的不是「谁更准」,
   而是「跨度能不能代替宽度」** —— 与第一臂问的完全同一个问题,这一条是它的可比性前提。
⚠ 本臂结构上做不到的:MSSCQ **单次采集、无年代** ⇒ 换不了「随时间变化」这一层;
   它能换的只有**单位(人 vs 年份点)与仪器**,而那正好是要换的两样。
"""
import numpy as np, pandas as pd, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import median_of_group_medians, pooled_median   # ⚠ 用项目自己的,不手搓
from lib.gates import Gate

RNG = np.random.default_rng(2220)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"; OUT.mkdir(exist_ok=True)
P = ROOT/"data/external/openpsych/MSSCQ/MSSCQ/"

Q = [f"Q{i}" for i in range(1, 101)]
D = pd.read_csv(P/"data.csv", sep="\t")
X = D[Q].replace(0, np.nan).dropna()
X = X[(X >= 1).all(axis=1) & (X <= 5).all(axis=1)]
n = len(X)
print(f"=== MSSCQ · 完整作答 n = {n:,}(`#726`:这是测量的 n,不是 wc -l 的发布行数)===")

iu = np.triu_indices(100, 1)
fac = np.array([i % 20 for i in range(100)])
same = fac[iu[0]] == fac[iu[1]]
assert (int(same.sum()), int((~same).sum())) == (200, 4750), "分面结构变了,停"
print(f"  正控①:面内 {int(same.sum())} 对 · 面间 {int((~same).sum())} —— 与 `#542`/`#163` 逐字相同")

GROUPS = [[i for i in range(100) if i % 20 == f] for f in range(20)]
CELLS = ("生|ρ| × 逐面中位", "生|ρ| × 200对合并", "归一ρ × 逐面中位", "归一ρ × 200对合并")

def grid(A):
    """2×2 规格网格:矩阵口径(生|ρ| / 归一 ρ)× 汇总口径(逐面中位 / 200 对合并)。"""
    R = pd.DataFrame(A).rank().to_numpy(float)
    C = np.corrcoef(R.T)
    CEIL = np.abs(np.corrcoef(np.sort(R, axis=0).T))
    N = np.where(CEIL > 1e-9, C/CEIL, np.nan)
    out = []
    for M in (np.abs(C), N):
        out.append(median_of_group_medians(M, GROUPS))
        out.append(float(np.median([M[a, b] for g in GROUPS
                                    for i, a in enumerate(g) for b in g[i+1:]])))
    return out

g = grid(X[Q].to_numpy())
g_facet, g_pooled = g[2], g[1]          # 页面那两个数:归一ρ×逐面 / 生|ρ|×合并
span = float(max(g) - min(g))
print(f"\n=== ① 全样本 · 2×2 规格网格(`G4`:曲线不是一格)===")
for nm, v in zip(CELLS, g):
    tag = ""
    if abs(v-0.5645) < 1e-3: tag = "  ← **页面写的 0.5645 是这一格**"
    if abs(v-0.5856) < 1e-3: tag = "  ← **页面写的 0.5856 是这一格**"
    print(f"  {nm:18s} {v:.4f}{tag}")
print(f"  ⇒ **页面并排写的两个数差了两个自由度**(矩阵口径 × 汇总口径),不是一个规格跨度")
print(f"  ⇒ **整张网格的规格跨度 {span:.4f}**(GSS 那边 `0.409–0.431` 的跨度是 0.022)")

B = 400
print(f"\n=== ② 自助(重抽的是**人**,不是年份点)· B={B} ===")
A = X[Q].to_numpy()
BS = np.array([grid(A[RNG.integers(0, n, n)]) for _ in range(B)])
widths = []
for j, nm in enumerate(CELLS):
    lo, hi = np.percentile(BS[:, j], 2.5), np.percentile(BS[:, j], 97.5)
    widths.append(float(hi-lo))
    print(f"  {nm:18s} [{lo:.4f}, {hi:.4f}] 宽 **{hi-lo:.4f}**")
w_facet, w_pooled = widths[2], widths[1]
w_med = float(np.median(widths))
print(f"  ⇒ **自助宽度中位 {w_med:.4f} ÷ 规格跨度 {span:.4f} = {w_med/span:.2f}×**(GSS 那边是 6.3×,阈值 2×)")

GSS_MED = 6.3
G = Gate("#783 第二臂 · 「规格跨度当误差」是我的毛病还是那个设计的毛病")
pc = bool(abs(g_facet-0.5645) < 0.01 and abs(g_pooled-0.5856) < 0.01)
G.asserted("① 正控:全样本必须复现页面已发表的两个值(容差 0.01),否则不是同一个对象",
           pc, f"逐面 {g_facet:.4f} vs 0.5645(差 {abs(g_facet-0.5645):.4f})· "
               f"合并 {g_pooled:.4f} vs 0.5856(差 {abs(g_pooled-0.5856):.4f})", kind="control")
G.asserted("② 正控:自助真的在变(2×2 四格的区间都非零宽)",
           bool(min(widths) > 1e-6), f"四格宽 {[round(w,4) for w in widths]} · B={B}", kind="control")
G.asserted("③ 负控:跨度不得为零(否则比值分母退化,判词无意义)",
           bool(span > 1e-6), f"规格跨度 {span:.4f}", kind="control")
gen = bool(w_med/span > 2.0)
G.asserted("④ kill(预注册):C-general 要成立,需自助/跨度 > 2×,与 GSS 那边同号",
           gen, f"MSSCQ {w_med/span:.2f}× vs GSS {GSS_MED}× · 阈值 2×", kind="kill")
print(); print(G)

print("\n" + "="*92)
ok = pc and min(widths) > 1e-6 and span > 1e-6
if not ok:
    v = "**UNVERIFIED:正控没过,本臂不下判。**"
elif gen:
    v = (f"**C-general**:换到 MSSCQ(单位=人,n={n:,})之后,自助宽度仍是规格跨度的 **{w_med/span:.2f}×** —— "
         f"与 GSS 的 {GSS_MED}× 同号。⇒ 第一臂那句「页面一直在公布小的那种不确定性」**按原样成立**,"
         f"它是我的报告习惯,不是某一个设计的毛病。")
else:
    v = (f"**C-design —— 而这是我不想要的那一支,所以它更值钱。** 换到 MSSCQ(单位=人,n={n:,})之后,"
         f"自助宽度只有规格跨度的 **{w_med/span:.2f}×**(GSS 那边是 {GSS_MED}×),**没过 2× 的阈值。**\n"
         f"  ⇒ **第一臂那句话必须挂上限定**:不是「页面一直在公布小的那种不确定性」,而是\n"
         f"  **「在有效 n 只有 28 个年份点的那一类量上,规格跨度远小于抽样宽度」** ——\n"
         f"  在 n={n:,} 个人上,规格跨度反而是**更保守**的那一个,拿它当误差条并不失真。\n"
         f"  ⚠ 差别不在我怎么写,在**被重抽的单位有多少个**:28 个年份点 vs {n:,} 个人。")
print(v)
json.dump(dict(n=n, grid_cells=dict(zip(CELLS, g)), grid_widths=dict(zip(CELLS, widths)),
               gauge_facet=g_facet, gauge_pooled=g_pooled, span=span,
               width_facet=w_facet, width_pooled=w_pooled, width_median=w_med,
               ratio=w_med/span, gss_ratio=GSS_MED, B=B,
               c_general=gen, verdict=v, gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"msscq_span_arm.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'msscq_span_arm.json'}")

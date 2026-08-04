import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A40 R241 -- Δ 是一条通则,还是几个类别的性质

**换方向**(`#195c` 的盆地信号)。回到 `#128` 的 Δ = `mean_i rho_i`,
其中 rho_i 是这个人「起始年龄残差 × 类别稀有度」的**人内相关**(尺度无关)。

⚠ **rho 对类别不可加**(分母有 √Σy²)——「逐类别分解」在数学上走不通。
第一版正是这么做的,而它的和是 +0.0590、`#128` 是 −0.0328:**分解的是 b,不是 rho。**
**改用留一类别** —— 它对任何统计量都良定义。

ESTIMAND        `Δ_(-j)` = 去掉类别 j 之后重算的 `mean_i rho_i`;判 `|Δ − Δ_(-j)|` 的分布。
KILL            **若去掉单个类别就能让 |Δ| 变化超过它自身的一半 -> Δ 不是通则,
                是几个特定类别的性质,`#128` 的措辞要从"人"改成"这些类别"。**
POSITIVE CTRL   把信号**只种在 3 个类别**上,强度调到 Δ 明显变大 -> 留一必须点出那 3 个。
NEGATIVE CTRL   **均匀**种在全部类别上 -> 留一曲线必须平。
NULL            题内跨人置换(`#173` 自己的 `perm_null`)。
NOISE FLOOR     人层 bootstrap 60 次给 Δ 的展布;留一的变化要对着它读。
IMPOSSIBLE      去掉一个类别会同时改变每个人的 `x` 中心化与 `v` —— 留一测的是
                「这个类别在场与否」的**总影响**,不是它的"贡献份额"。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
M_CAT=V.shape[1]
print(f"载入:V {V.shape} · KEEP {int(KEEP.sum()):,} 人 · 类别 {M_CAT}",flush=True)

_obs_full=obs.copy(); _rar_full=rar.copy(); _V_full=V.copy()
def delta_of(Vm, drop=None):
    """去掉类别 drop 后重算 mean_i rho_i。用轮次自己的 `betas`。"""
    global obs, rar, V
    cols=[j for j in range(M_CAT) if j!=drop] if drop is not None else list(range(M_CAT))
    obs=_obs_full[:,cols]; rar=_rar_full[cols]; V=Vm[:,cols]
    try:
        _,rho=betas(V); m=np.isfinite(rho)&KEEP
        return float(np.mean(rho[m]))
    finally:
        obs=_obs_full; rar=_rar_full; V=_V_full

D0=delta_of(_V_full)
print(f"\nΔ = mean_i rho_i = {D0:+.4f}   (`#128` 报的是 −0.0328)")
rng=np.random.default_rng(20260803)
boot=[]
for _ in range(60):
    ii=rng.choice(np.flatnonzero(KEEP),int(KEEP.sum()),replace=True)
    Kb=np.zeros(len(_V_full),bool); Kb[ii]=True
    sk=KEEP.copy(); KEEP[:]=Kb; boot.append(delta_of(_V_full)); KEEP[:]=sk
sdD=float(np.std(boot)); print(f"人层 bootstrap sd = {sdD:.4f}  -> {abs(D0)/sdD:.1f}×")
Dn=delta_of(perm_null(_V_full,rng)); print(f"置换零 Δ = {Dn:+.4f}")

def loo(Vm, base):
    return np.array([base-delta_of(Vm,drop=j) for j in range(M_CAT)])

L=loo(_V_full,D0)
print(f"\n留一:|ΔΔ| 最大 {np.abs(L).max():.4f} · 中位 {np.median(np.abs(L)):.4f} · "
      f"最大 / |Δ| = {np.abs(L).max()/abs(D0):.2f}")
top=np.argsort(np.abs(L))[::-1][:5]
for j in top: print(f"  ΔΔ {L[j]:+.5f}  稀有度 {_rar_full[j]:.2f}  n={int(_obs_full[:,j].sum()):,}  {str(ons[j])[:52]}")

# 对照:把强度调到 Δ 明显变大(`#165b`:强度必须匹配可观测量级)
def plant_cats(cats,g):
    # ⚠ #196b:第一版种的是 `x[cats]=1.0` —— 那是一个**人主效应**(在这几个类别里整体偏移),
    #   而收敛双向去均值正好把人主效应吃掉,所以 Δ 纹丝不动、g 扫完五档都没找到。
    #   要种的是与**稀有度的交互**(`#173` 自己的 `plant()` 就是这么写的)。
    # ⚠ #196c:`rho_i` 是**尺度无关**的,所以种一个 `u_i ~ N(0,1)` 的交互会把每个人的 rho
    #   推向 `sign(u_i)·1` —— **平均下来是 0**,Δ 纹丝不动。`#173` 自己的 `plant()` 也是这个形式,
    #   但它评估的是 `rho_corr_S`(rho 与 S 的相关),不是 `mean rho`。
    #   **要移动一个人层平均的统计量,种植的人载荷必须符号一致。**
    u=np.abs(rng.standard_normal(len(_V_full)))+0.5     # 全正
    x=np.zeros(M_CAT)
    x[cats]=_rar_full[cats]-_rar_full.mean()
    return _V_full+g*np.outer(u,x)*_obs_full
# ⚠ #196d:第一版把正对照种在**最罕见的 3 个类别**上 —— 而那 3 个只有约 500/9,944 人答过,
#   种进去只够到 **5% 的人**,g 扫到 64 都推不动一个人层平均。**种植的触及面也是强度的一部分。**
#   这本身就是留一为什么平的一半答案。正对照改种在**人数最多**的 3 个类别上。
PCATS=np.argsort(_obs_full.sum(0))[-3:]
print(f"正对照种植的 3 个类别:n = {_obs_full.sum(0)[PCATS]} · 稀有度 {np.round(_rar_full[PCATS],2)}")
G=None
for g_ in [0.25,0.5,1.0,2.0,4.0,8.0,16.0]:
    d=delta_of(plant_cats(PCATS,g_))
    if abs(d-D0)>3*sdD: G=g_; break
print(f"\n对照用的种植强度 g = {G}(调到 Δ 移动 >3 个 bootstrap sd)")
assert G is not None, "扫完所有强度都没让 Δ 动 —— 种植没生效,对照不可用"
Vp=plant_cats(PCATS,G); Dp=delta_of(Vp); Lp=loo(Vp,Dp)
Vu=plant_cats(np.arange(M_CAT),G);          Du=delta_of(Vu); Lu=loo(Vu,Du)
pos_hit=set(np.argsort(np.abs(Lp))[::-1][:3])==set(PCATS)
print(f"正对照(只种 3 类,Δ {D0:+.4f}->{Dp:+.4f}):留一前 3 名{'命中' if pos_hit else '未命中'}那 3 个类别")
print(f"负对照(均匀种,Δ {D0:+.4f}->{Du:+.4f}):留一 |ΔΔ| 最大/中位 = {np.abs(Lu).max()/max(np.median(np.abs(Lu)),1e-12):.1f}"
      f"  vs 正对照 {np.abs(Lp).max()/max(np.median(np.abs(Lp)),1e-12):.1f}")

T=pd.DataFrame(dict(cat=[str(ons[j])[:52] for j in range(M_CAT)],rarity=_rar_full,
                    n_people=_obs_full.sum(0),loo_delta=L))
check_columns(T,'R241'); T.to_csv(pathlib.Path(__file__).parent/'results'/'loo.csv',index=False)
g=Gate('Δ 是通则还是几个类别的性质')
g.asserted('正对照:只种 3 类时留一必须点出那 3 个',pos_hit,
           f"前 3 名 = {sorted(np.argsort(np.abs(Lp))[::-1][:3])} vs 种植 {sorted(PCATS)}")
g.asserted('负对照:均匀种植时留一曲线必须比集中种植平',
           (np.abs(Lu).max()/max(np.median(np.abs(Lu)),1e-12)) < (np.abs(Lp).max()/max(np.median(np.abs(Lp)),1e-12)),
           f"均匀 {np.abs(Lu).max()/max(np.median(np.abs(Lu)),1e-12):.1f} vs 集中 "
           f"{np.abs(Lp).max()/max(np.median(np.abs(Lp)),1e-12):.1f}")
g.negative_control('题内跨人置换的 Δ',abs(Dn),abs(D0))
g.resolvable('Δ 本身',D0,sdD)
g.asserted('注册的 kill:去掉单个类别让 |Δ| 变化超过它自身的一半',
           float(np.abs(L).max()/abs(D0))>0.5,f"最大 ΔΔ / |Δ| = {np.abs(L).max()/abs(D0):.2f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

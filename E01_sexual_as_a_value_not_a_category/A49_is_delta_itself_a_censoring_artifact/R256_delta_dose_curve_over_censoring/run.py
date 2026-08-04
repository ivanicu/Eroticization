import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A49 R256 -- `#128` 的 Δ 本身,是不是截断

`#210` 撤掉了一条**关于年龄**的观察,而截断这个机制**对整条起始年龄线都成立** ——
`#128` 的 Δ(`mean rho_i = −0.0328`)也是在**全部**起始年龄上算的。
`#255` 顺带给出受限臂的 `mean rho_i = +0.0055` —— **它也反号了。**
**这是本项目最中心的一条声明,必须正面查。**

ESTIMAND        用 `R173` 自己的管道,在 cap ∈ {14, 16, 18, ∞} 上重算 `mean rho_i`;
                每档报 `KEEP` · 置换零 · bootstrap sd。
⚠ 最强混杂(跑之前写下)
                **截断会优先删掉罕见类别**(它们本来就获得得晚)-> 压缩人内 `x` 的变异
                -> **低 cap 档的零可能是范围受限,不是效应不存在。**
控制(同一迭代内)
                ① 每档报**人内稀有度 sd 的中位数**(x 的变异还剩多少);
                ② **正对照必须在同一档上开火** —— 在每个 cap 上种一个已知 Δ,
                   若种植在 cap=14 上也测不到,那一档的零**不可读**。
KILL            **若 Δ 在正对照仍开火的最低 cap 上不可分辨或反号 ->
                `#128` 的 Δ 本身要重估。**
NEGATIVE CTRL   每档各跑一次题内跨人置换。
IMPOSSIBLE      低 cap 档剩下的人不是随机子集(`#210c`)——
                **足以证伪,不足以立反向声明。**
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V_full=V.copy(); obs_full=obs.copy(); KEEP_full=KEEP.copy(); rar_full=rar.copy()
print(f"原管道:V {V.shape} · KEEP {int(KEEP.sum()):,}")

def with_cap(Vm, fn):
    global V, obs, KEEP
    V=Vm; obs=np.isfinite(Vm); KEEP=(obs.sum(1)>=8)
    try: return fn()
    finally: V=V_full; obs=obs_full; KEEP=KEEP_full

def delta_and_x(Vm):
    def f():
        _,rho=betas(V); m=np.isfinite(rho)&KEEP
        xs=[np.std(rar_full[np.flatnonzero(obs[i])]) for i in np.flatnonzero(KEEP)[:4000]]
        return float(np.mean(rho[m])), int(KEEP.sum()), int(m.sum()), float(np.median(xs)), rho, m
    return with_cap(Vm,f)

rng=np.random.default_rng(20260803)
# ⚠ #211a:第一版把种植强度写死成 0.8 —— 在 ∞ 档只把 Δ 从 −0.0328 推到 −0.0299,
#   **正对照一档都没开火**,于是低 cap 档的零全都不可读。`#196d` 的同一形状:
#   **种植强度必须先扫到"它真的推得动"为止,再用同一个 g 跑全部档。**
d0,_,_,_,_,_=delta_and_x(V_full)
sd0=0.0036
G=None
for g_ in [0.8,2,5,10,20,40,80]:
    u=np.abs(rng.standard_normal(len(V_full)))+0.5
    Vp=V_full+g_*np.outer(u,rar_full-rar_full.mean())*obs_full
    dp,_,_,_,_,_=delta_and_x(Vp)
    if abs(dp-d0)>3*sd0: G=g_; break
assert G is not None, '扫完所有强度都推不动 Δ —— 种植没生效,对照不可用'
print(f"种植强度扫描 -> g = {G}(在 ∞ 档把 Δ 推动 >3 sd)")
rows=[]
for cap in [14,16,18,None]:
    Vm=V_full if cap is None else np.where(V_full<=cap,V_full,np.nan)
    dm,nk,nu,xsd,rho,m=delta_and_x(Vm)
    bs=[float(np.mean(rho[i])) for i in
        [rng.choice(np.flatnonzero(m),int(m.sum()),replace=True) for _ in range(300)]]
    dn,_,_,_,_,_=delta_and_x(perm_null(Vm,np.random.default_rng(600+(cap or 99))))
    # 正对照:在同一档上种一个已知 Δ
    u=np.abs(rng.standard_normal(len(V_full)))+0.5
    Vp=Vm+G*np.outer(u,rar_full-rar_full.mean())*np.isfinite(Vm)
    dp,_,_,_,_,_=delta_and_x(Vp)
    rows.append(dict(cap=('∞' if cap is None else f'≤{cap}'),keep=nk,used=nu,
                     x_sd_median=xsd,delta=dm,sd=float(np.std(bs)),null=dn,plant=dp))
T=pd.DataFrame(rows); check_columns(T,'R256'); T.to_csv(pathlib.Path(__file__).parent/'results'/'dose.csv',index=False)
print(f"\n{'cap':<6}{'KEEP':>8}{'人内稀有度 sd':>14}{'Δ':>11}{'sd':>9}{'比':>7}{'置换零':>10}{'种植后 Δ':>11}")
for _,r in T.iterrows():
    print(f"{r.cap:<6}{r.keep:>8,}{r.x_sd_median:>14.3f}{r.delta:>+11.4f}{r.sd:>9.4f}"
          f"{abs(r.delta)/r.sd:>7.1f}{r.null:>+10.4f}{r.plant:>+11.4f}")

full=T[T.cap=='∞'].iloc[0]
ok_plant=T[(T.plant-T.delta).abs()>3*T.sd]      # 正对照在哪些档上开火
low=ok_plant.iloc[0] if len(ok_plant) else None
print(f"\n正对照开火的档:{list(ok_plant.cap)}")
g=Gate('#128 的 Δ 本身是不是截断')
g.asserted('可判前提:∞ 档复现 `#128` 的 −0.0328',abs(full.delta+0.0328)<0.003,f"{full.delta:+.4f}")
g.negative_control('∞ 档的置换零',abs(float(full.null)),abs(float(full.delta)))
g.asserted('可判前提二:正对照在最低档上也开火(否则那一档的零不可读)',
           '≤14' in list(ok_plant.cap),
           f"开火的档 {list(ok_plant.cap)}")
if low is not None:
    g.resolvable(f'最低可读档({low.cap})的 Δ',float(low.delta),float(low.sd))
    g.offset_control(f'{low.cap} vs ∞',float(low.delta),float(full.delta),
                     float(np.hypot(low.sd,full.sd)),
                     null_kind='同一条管道在全部起始年龄上的 Δ —— 不是零假设,是"若截断无关,受限档该落在哪"')
    g.asserted('注册的 kill:最低可读档上 Δ 不可分辨或反号 -> `#128` 要重估',
               (abs(low.delta)<=2*low.sd) or (np.sign(low.delta)!=np.sign(full.delta)),
               f"{low.cap}: Δ {low.delta:+.4f} ± {low.sd:.4f} vs ∞ {full.delta:+.4f}")
g.asserted('⚠ 低 cap 档的人不是随机子集 —— 足以证伪,不足以立反向声明(`#210c`)',True,
           f"KEEP ∞ {full.keep:,} -> ≤14 {int(T[T.cap=='≤14'].keep.iloc[0]):,}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

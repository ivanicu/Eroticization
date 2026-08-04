import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A19 R04 -- "共享时间表解释掉 80%" 这个数,依赖我用哪个汇总量吗?

#148a:「最早那批更常见」的 **80%** 是共享发育时间表的投影(减去题目**均值**起始年龄),
而方差匹配的随机题目向量只解释 36%,时间表特异的部分 44 个百分点(20.0×)。

**唯一的薄弱处**:那张时间表是 `#75` 的**题目均值序**,而 `#131c` 已经证明
**中位数序与均值序不同**(留出成对顺序 63.30% vs 66.70%),而且
**中位数序与稀有度对齐(Spearman +0.437),均值序不对齐(+0.011)**。

**预测(跑之前写死):** 既然 Δ 是一个**关于稀有度**的量,而中位数序才是与稀有度对齐的那个,
**中位数序应当解释掉比均值序更多**。若不是,那么"哪个汇总量与稀有度对齐"与
"哪个汇总量能解释 Δ"是两件事,而 `#131c`/`#130e` 那条线的读法要改。

    STABLE  各汇总量解释掉的比例彼此相差 < 15 个百分点 -> #148a 拿到规格稳健性
    MOVES   相差很大 -> #148a 的**量级**要加限定,"80%" 必须写成一个区间

ESTIMAND        z(m=1,按人内置换归一化)在各臂上的值,以及每个汇总量解释掉的比例;
                每个汇总量都配一个**方差匹配**的随机题目向量作稀释基准(#148d)。
IDENTIFICATION  所有臂算的是同一个统计量(m=1,#148c 的修正),同一批人。
SCOPE           >=8 个类别起始年龄的人。
WORLDS          STABLE / MOVES
KILL            条件式:每个汇总量的随机稀释基准必须落在合理范围(不为零也不等于该臂),
                才读解释比例。
POSITIVE CTRL   raw 臂本身必须可分辨(#148 已知 45.3×)。
NEGATIVE CTRL   方差匹配的随机题目向量,**每个汇总量各一个**(稀释基准,不是零)。
NOISE FLOOR     按人自助 200 次。
MULTIPLICITY    5 个汇总量 x {真实, 随机稀释},整格发表。
IMPOSSIBLE      "哪个汇总量是真的时间表"本身没有金标准;本轮只判**解释比例是否随它变动**。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R01_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

NPERM=200
def zstat(Vm,seed):
    rg=np.random.default_rng(seed); tie=np.random.default_rng(20260803)
    Z=np.full(len(Vm),np.nan)
    for i in np.flatnonzero(KEEP):
        j=np.flatnonzero(obs[i]); r=rar[j]; y=Vm[i,j]; k=len(j)
        if k<2: continue
        cand=np.flatnonzero(y==np.nanmin(y)); pick=cand[tie.integers(len(cand))]
        d=r[pick]-r.mean()
        idx=rg.integers(0,k,(NPERM,1)); dr=r[idx].mean(1)-r.mean()
        if dr.std()<1e-9: continue
        Z[i]=(d-dr.mean())/dr.std()
    return Z

SUM={'mean':   lambda: np.array([np.nanmean(V[obs[:,j],j]) for j in range(V.shape[1])]),
     'median': lambda: np.array([np.nanmedian(V[obs[:,j],j]) for j in range(V.shape[1])]),
     'q25':    lambda: np.array([np.nanpercentile(V[obs[:,j],j],25) for j in range(V.shape[1])]),
     'q75':    lambda: np.array([np.nanpercentile(V[obs[:,j],j],75) for j in range(V.shape[1])]),
     'trim20': lambda: np.array([np.nanmean(np.clip(V[obs[:,j],j],
                   np.nanpercentile(V[obs[:,j],j],10),np.nanpercentile(V[obs[:,j],j],90)))
                   for j in range(V.shape[1])])}
from scipy.stats import spearmanr
Zraw=zstat(V,1); mr=np.isfinite(Zraw); iraw=np.flatnonzero(mr)
rb=np.random.default_rng(7)
z0=float(Zraw[mr].mean())
b0=float(np.std([Zraw[iraw[rb.integers(0,len(iraw),len(iraw))]].mean() for _ in range(200)]))
print(f"raw 臂 z = {z0:+.4f} ± {b0:.4f}({abs(z0)/b0:.1f}×,#148 测得 −0.4462)",flush=True)

rows=[]
print(f"\n{'汇总量':<9}{'sd':>7}{'ρ(稀有度)':>11}{'真实 z':>10}{'随机稀释 z(20 抽)':>18}{'解释掉':>9}{'稀释解释':>10}{'净':>8}")
for nm,f in SUM.items():
    vec=f(); sp=spearmanr(rar,vec).statistic
    Zi=zstat(np.where(obs,V-vec[None,:],np.nan),zlib.crc32(nm.encode())%9991)
    # ⚠ #148d 的稀释基准是**一次抽样**。重抽 20 次 —— 单次抽样的实现差异极大
    #   (本轮首跑就出现了负的"稀释解释",即随机向量反而让效应变大)。
    zrs=[]
    for t_ in range(20):
        rg=np.random.default_rng(4242+len(nm)*97+t_)
        rvec=rg.normal(vec.mean(),vec.std(),len(vec))
        Zt=zstat(np.where(obs,V-rvec[None,:],np.nan),zlib.crc32((nm+str(t_)).encode())%9991)
        zrs.append(float(Zt[np.isfinite(Zt)].mean()))
    zrs=np.array(zrs)
    mi=np.isfinite(Zi)
    zi=float(Zi[mi].mean()); zr=float(zrs.mean())
    ex=1-zi/z0; exr=1-zr/z0
    rows.append(dict(summary=nm,sd=float(vec.std()),sp=float(sp),z=zi,z_rand=zr,
                     z_rand_sd=float(zrs.std()),z_rand_lo=float(zrs.min()),z_rand_hi=float(zrs.max()),
                     expl=ex,expl_rand=exr,net=ex-exr))
    print(f"{nm:<9}{vec.std():>7.2f}{sp:>+11.3f}{zi:>+10.4f}{zr:>+9.4f}±{zrs.std():.3f}"
          f"{100*ex:>8.0f}%{100*exr:>9.0f}%{100*(ex-exr):>7.0f}")

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'summaries.csv',index=False)
sp_ex=100*float(T.expl.max()-T.expl.min()); sp_net=100*float(T.net.max()-T.net.min())
mn=T[T.summary=='mean'].iloc[0]; md=T[T.summary=='median'].iloc[0]
print(f"\n  解释比例的极差 {sp_ex:.0f} 个百分点;净(扣掉稀释)极差 {sp_net:.0f} 个百分点")
print(f"  预测:中位数序与稀有度对齐({md.sp:+.3f})而均值序不({mn.sp:+.3f}),"
      f"所以中位数应解释更多 —— 实测 中位数 {100*md.expl:.0f}% vs 均值 {100*mn.expl:.0f}%")

g=Gate('"解释掉 80%"依赖汇总量吗')
g.asserted('raw 臂可分辨(正对照)',abs(z0)/b0>10,f"{z0:+.4f} ± {b0:.4f} = {abs(z0)/b0:.1f}×")
g.asserted('⚠ #148d 的稀释基准是一次抽样,而单次实现差异极大',
           bool((T.z_rand_sd>0.02).any()),
           " ".join(f"{s}:{100*(1-lo/z0):.0f}%..{100*(1-hi/z0):.0f}%"
                    for s,lo,hi in zip(T.summary,T.z_rand_lo,T.z_rand_hi))
           + " —— 20 次抽样的解释比例区间。**#148 的 36% 是其中一次**")
g.asserted('20 抽平均后,稀释基准落在中间(不为零也不等于该臂)',
           bool((T.expl_rand>0.05).all() and (T.expl_rand<T.expl).all()),
           " ".join(f"{s}:{100*e:.0f}%" for s,e in zip(T.summary,T.expl_rand)))
g.asserted('预测:中位数序解释掉的比均值序多',md.expl>mn.expl,
           f"中位数 {100*md.expl:.0f}%(ρ={md.sp:+.3f})vs 均值 {100*mn.expl:.0f}%(ρ={mn.sp:+.3f})"
           + ("" if md.expl>mn.expl else " —— **预测反了**"))
# ⚠ #148a 的「44 个百分点,20.0×」用错了分母:它拿**按人自助**的展布(0.0098)去判一个
#   由**基准自身实现**主导的差。正确的分母是稀释基准的实现展布。
mnr=T[T.summary=='mean'].iloc[0]
net_sd=float(mnr.z_rand_sd/abs(z0))
print(f"\n  正确的分母:均值臂稀释基准的实现展布 = {mnr.z_rand_sd:.3f}(z 单位)"
      f" = {100*net_sd:.0f} 个百分点")
print(f"  所以 #148a 的净额应作 **{100*mnr.net:.0f} ± {100*net_sd:.0f} 个百分点**"
      f"({abs(mnr.net)/net_sd:.1f}×),而不是 44 个百分点 20.0×")
g.require_resolvable_first('净额相对**基准自身实现展布**是否可分辨',abs(float(mnr.net)),net_sd)
g.equivalent_within('净解释比例在各汇总量之间等价(边界 15 个百分点)',
                    float(T.net.max()-T.net.min()),net_sd,0.15)
g.no_sign_crossing('每个汇总量都解释掉正的比例',list(T.expl.values))
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

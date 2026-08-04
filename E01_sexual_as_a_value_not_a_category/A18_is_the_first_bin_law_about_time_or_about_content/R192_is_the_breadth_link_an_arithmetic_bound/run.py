import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A18 R02 -- |Δ| 随广度变(+0.661),那是心理学,还是一条纯算术的界?

#146c:|Δ| 在 12 个侧格上变 27 倍,corr(|Δ|, 合格者平均类别数) = +0.661。
#146 的 NEXT 写的是"在广度上做卡钳 1:1 匹配"。**先纠正我自己的路标:那是错的工具。**
匹配比较**两个组**,而 Δ 没有组 —— 它的零是**人内**的,每个人的 E[Δ_i | null] 恰好为 0。
所以广度不可能让 Δ_i 偏离零;它只能改变**幅度**。要判的是那个幅度差是不是真的。

而在做任何心理学解释之前,有一条**纯算术**的界必须先排除:

    一个人有 k 个类别,最早一格占 m 个。Δ = (那 m 个的平均稀有度) − (全部 k 个的平均)。
    m 个的均值与 k 个的均值之差,幅度被 **(k−m)/k** 直接卡住 ——
    k=8、m=4 时最早一格是集合的一半,Δ 几乎动不了;
    k=27、m=4 时它只占 15%,Δ 能偏离得远得多。
    **类别越多 -> |Δ| 越大,这在算术上是强制的,与心理学无关。**

判据:用**每个人自己的人内置换零的展布**去归一化。那条界被精确除掉。

  ARITH  归一化后 z 在各广度层上**平坦** -> +0.661 就是那条算术界,
         「最早的更常见」是**普遍的**,只是原始幅度不可跨人比较 -> #130a 得救,但要换单位
  PSYCH  归一化后 z 仍随广度上升 -> 广度效应是真的,#130a 只对口味广的人成立

ESTIMAND        z_i = Δ_i / sd(Δ_i | 人内置换该人自己的起始年龄标签),按广度分层。
IDENTIFICATION  分母来自**这个人自己的**置换分布,所以 (k−m)/k 的界在分子分母里同样出现,
                精确抵消。这是结构免疫,不是回归控制。
SCOPE           >=8 个类别起始年龄、最早一格不是全部的人。
WORLDS          ARITH / PSYCH
KILL            条件式:归一化的正对照必须开火(合成一个只有算术界、无心理学信号的世界,
                它的 z 必须在各层平坦且为零),才读真实 z 的层间趋势。
POSITIVE CTRL   种植一个**与广度无关**的真实效应,z 必须在各层平坦地为负。
NEGATIVE CTRL   人内置换(每人 200 次,同时给出分母)。
NOISE FLOOR     按人自助 200 次/层。
MULTIPLICITY    5 个广度层 x {原始 Δ, 归一化 z},整格发表。
IMPOSSIBLE      广度与"这个人答了多少块"在本 release 里分不开(#5)。本轮只判
                **归一化后趋势是否还在**,不判广度本身是什么。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

NPERM=200
def per_person(Vm,seed=0):
    """返回 (Δ_i, 该人置换零的均值, 该人置换零的 sd, k, m)。"""
    rg=np.random.default_rng(seed)
    D=np.full(len(Vm),np.nan); M0=np.full(len(Vm),np.nan); S0=np.full(len(Vm),np.nan)
    KK=np.full(len(Vm),np.nan); MM=np.full(len(Vm),np.nan)
    for i in np.flatnonzero(KEEP):
        j=np.flatnonzero(obs[i]); r=rar[j]; y=Vm[i,j]
        lo=y.min(); m=int((y==lo).sum()); k=len(j)
        if m==0 or m>=k: continue
        D[i]=r[y==lo].mean()-r.mean(); KK[i]=k; MM[i]=m
        # 人内置换零:从这个人自己的 k 个里随机抽 m 个
        idx=rg.integers(0,k,(NPERM,m))
        draws=r[idx].mean(1)-r.mean()
        M0[i]=draws.mean(); S0[i]=draws.std()
    return D,M0,S0,KK,MM

D0,M0,S0,KK,MM=per_person(V,zlib.crc32(b'A18R02'))
ok=np.isfinite(D0)&np.isfinite(S0)&(S0>1e-9)
Z0=(D0-M0)/S0
print(f"{ok.sum():,} 人;k 范围 {np.nanmin(KK[ok]):.0f}–{np.nanmax(KK[ok]):.0f},"
      f"m 范围 {np.nanmin(MM[ok]):.0f}–{np.nanmax(MM[ok]):.0f}",flush=True)
print(f"整体:Δ = {np.nanmean(D0[ok]):+.4f}   归一化 z = {np.nanmean(Z0[ok]):+.4f}",flush=True)

qs=np.percentile(NCAT[ok],[20,40,60,80])
def strat(Dv,Zv,mask):
    rows=[]
    b=np.digitize(NCAT,qs)
    for g in range(5):
        m=mask&(b==g)
        if m.sum()<200: continue
        rb=np.random.default_rng(1000+g); ii=np.flatnonzero(m)
        bd=float(np.std([Dv[ii[rb.integers(0,len(ii),len(ii))]].mean() for _ in range(200)]))
        bz=float(np.std([Zv[ii[rb.integers(0,len(ii),len(ii))]].mean() for _ in range(200)]))
        rows.append(dict(band=g,n=int(m.sum()),ncat=float(NCAT[m].mean()),
                         delta=float(Dv[m].mean()),d_boot=bd,
                         z=float(Zv[m].mean()),z_boot=bz))
    return pd.DataFrame(rows)

T=strat(D0,Z0,ok)
print(f"\n{'层':<4}{'n':>7}{'平均类别数':>10}{'Δ':>10}{'展布':>8}{'归一化 z':>11}{'展布':>8}{'z 倍数':>8}")
for _,r in T.iterrows():
    print(f"{int(r.band):<4}{int(r.n):>7,}{r.ncat:>10.1f}{r.delta:>+10.4f}{r.d_boot:>8.4f}"
          f"{r.z:>+11.4f}{r.z_boot:>8.4f}{abs(r.z)/r.z_boot:>8.1f}x")

sp_d=float(T.delta.max()-T.delta.min()); sp_z=float(T.z.max()-T.z.min())
tr_d=float(np.polyfit(T.ncat,T.delta,1)[0]); tr_z=float(np.polyfit(T.ncat,T.z,1)[0])
print(f"\n  Δ  层间极差 {sp_d:.4f}  趋势 {tr_d:+.5f}/类别   ({sp_d/T.d_boot.mean():.1f}x 单层展布)")
print(f"  z  层间极差 {sp_z:.4f}  趋势 {tr_z:+.5f}/类别   ({sp_z/T.z_boot.mean():.1f}x 单层展布)")

# 正对照 A:纯算术世界 —— 起始年龄与稀有度**完全独立**(人内置换生成),z 必须平坦且为零
rgs=np.random.default_rng(77); Vs=V.copy()
for i in np.flatnonzero(KEEP):
    j=np.flatnonzero(obs[i]); Vs[i,j]=V[i,j][rgs.permutation(len(j))]
Ds,Ms,Ss,_,_=per_person(Vs,999); oks=np.isfinite(Ds)&np.isfinite(Ss)&(Ss>1e-9)
Zs=(Ds-Ms)/Ss; Ts=strat(Ds,Zs,oks)
print(f"\n正对照 A(纯算术世界:人内打乱起始年龄):"
      f"Δ 层间极差 {Ts.delta.max()-Ts.delta.min():.4f}(趋势 {np.polyfit(Ts.ncat,Ts.delta,1)[0]:+.5f}),"
      f"z 层间极差 {Ts.z.max()-Ts.z.min():.4f}(趋势 {np.polyfit(Ts.ncat,Ts.z,1)[0]:+.5f}),"
      f"整体 z = {np.nanmean(Zs[oks]):+.4f}")
# 正对照 B:种植一个与广度无关的真实效应 —— 每人以固定概率把最常见的一个搬进最早一格
rgp=np.random.default_rng(55); Vp=V.copy()
for i in np.flatnonzero(KEEP):
    if rgp.random()<0.5:
        j=np.flatnonzero(obs[i]); t=j[np.argmin(rar[j])]; Vp[i,t]=np.nanmin(V[i,j])
Dp,Mp,Sp,_,_=per_person(Vp,321); okp=np.isfinite(Dp)&np.isfinite(Sp)&(Sp>1e-9)
Zp=(Dp-Mp)/Sp; Tp=strat(Dp,Zp,okp)
print(f"正对照 B(种植与广度无关的真效应):整体 z = {np.nanmean(Zp[okp]):+.4f},"
      f"z 层间极差 {Tp.z.max()-Tp.z.min():.4f}(趋势 {np.polyfit(Tp.ncat,Tp.z,1)[0]:+.5f})")

T.to_csv(pathlib.Path(__file__).parent/'results'/'by_breadth.csv',index=False)
g=Gate('|Δ| 随广度变,是心理学还是算术界')
g.asserted('正对照 A:纯算术世界的 z 平坦且为零',
           abs(np.nanmean(Zs[oks]))<0.05 and abs(np.polyfit(Ts.ncat,Ts.z,1)[0])<0.01,
           f"整体 z {np.nanmean(Zs[oks]):+.4f},z 趋势 {np.polyfit(Ts.ncat,Ts.z,1)[0]:+.5f}")
g.asserted('正对照 B:种植的真效应在 z 上被检出',abs(np.nanmean(Zp[okp]))>0.1,
           f"整体 z {np.nanmean(Zp[okp]):+.4f}")
g.asserted('原始 Δ 确实随广度变',abs(tr_d)>0.005,f"趋势 {tr_d:+.5f}/类别,极差 {sp_d:.4f}")
# ⚠ 决定性的比较:正对照 B 种的是一个**与广度无关**的效应,所以它的 z 趋势就是
#   "归一化没除干净的那部分"的大小。真实趋势必须**超过**它,才谈得上真实的广度依赖。
tr_p=float(np.polyfit(Tp.ncat,Tp.z,1)[0])
g.asserted('真实的 z 趋势是否超过"与广度无关的效应"所产生的残余趋势',abs(tr_z)>abs(tr_p),
           f"真实 {tr_z:+.5f}/类别 vs 正对照 B(种的效应与广度无关){tr_p:+.5f}/类别 —— "
           + ("超过,广度依赖是真的" if abs(tr_z)>abs(tr_p) else
              "**没超过**。归一化后的残余梯度完全被一个与广度无关的效应解释掉,"
              "所以 #146c 的 +0.661 是那条算术界,不是心理学"))
g.require_resolvable_first('归一化后 z 的层间趋势是否还可分辨',sp_z,float(T.z_boot.mean()))
g.offset_control('归一化后最广层 vs 最窄层的 z',float(T.z.values[-1]),float(T.z.values[0]),
                 float(T.z_boot.mean()),null_kind='同一分层下的最窄广度层(不是零假设,是基准层)')
g.no_sign_crossing('每一层的 z 都为负',list(T.z.values))
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

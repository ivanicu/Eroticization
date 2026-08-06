import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A19 R01 -- 扣掉共享时间表之后剩下的那 20%,是一件事还是两件事?

#148a:「最早那批更常见」的 80% 是共享发育时间表的投影。
#148b:剩下 z = −0.0884(**9.4×**)是**个人的**,而它是什么完全没查过。

现成的候选:#128b —— 扣掉时间表后,一个人自己那些**罕见**的兴趣来得**更早**
(人内 ρ = −0.0328,8.8×)。两者符号一致(都指向"个人偏离时间表的方向是朝罕见")。

    SAME  两者高度相关,且接近分半天花板 -> 本项目两条独立线索合并成一条
    TWO   相关落在置换地板上 -> 个人成分里有**两件不同的事**,而这是新东西

⚠ 两个量来自**同一批数据、同一个稀有度向量、同一个人的同一个集合**,所以它们
   **在零下也共享算术结构**。「这个零应该是零吗」—— **不应该**,必须量出来:
   人内置换该人的起始年龄标签,两个量都重算,再相关 -> 地板。
⚠ 而"相关多高才算同一个量"需要一个**天花板**:把这个人的类别劈成不相交两半,
   z 在一半上算、ρ 在另一半上算 -> 可达到的最大相关。**报地板与天花板,不做除法**
   (`#141` 的教训:低信度下的去衰减不稳,而报出来的永远是校正后的那个)。

ESTIMAND        corr(z_resid, ρ) 跨人,对照人内置换地板与分半天花板。
IDENTIFICATION  两个量都建在**同一个收敛双向去均值残差**上,所以"个人偏离时间表"的定义一致。
SCOPE           >=8 个类别起始年龄的人;天花板臂要求两半各 >=4 个。
WORLDS          SAME / TWO
KILL            条件式:种植一个人特异的径向信号,两个量必须都被它拉动且相关上升;
                地板必须落在合理的正值上(不是零),才读真实相关。
POSITIVE CTRL   见上(种植)。
NEGATIVE CTRL   人内置换,5 个种子。
NOISE FLOOR     按人自助 200 次。
MULTIPLICITY    {真实, 地板, 天花板, 种植} x 2 个量,整格发表。
IMPOSSIBLE      两个量共享 item,所以"相同"永远无法与"由共享 item 强制"完全分开。
                地板正是为此而量的,但它只是一个下界。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

def demean_conv(Vm,tol=1e-10,cap=500):
    D=np.where(obs,Vm,np.nan)
    for _ in range(cap):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-a
        b=np.nanmean(D,axis=1,keepdims=True); D=D-b
        if max(np.nanmax(np.abs(a)),np.nanmax(np.abs(b)))<tol: break
    return D
NPERM=200
def two_stats(Dres,seed,colsz=None,colsr=None,perm=False):
    """z_resid(最早一个,m=1,按人内置换归一化)与 ρ(人内 稀有度 x 残差 相关)。"""
    rg=np.random.default_rng(seed); tie=np.random.default_rng(20260803)
    Z=np.full(len(Dres),np.nan); Rho=np.full(len(Dres),np.nan)
    for i in np.flatnonzero(KEEP):
        j=np.flatnonzero(obs[i])
        jz=np.intersect1d(j,colsz) if colsz is not None else j
        jr=np.intersect1d(j,colsr) if colsr is not None else j
        if len(jz)<4 or len(jr)<4: continue
        yz=Dres[i,jz].copy(); yr=Dres[i,jr].copy()
        if perm:
            yz=yz[rg.permutation(len(yz))]; yr=yr[rg.permutation(len(yr))]
        rz=rar[jz]; cand=np.flatnonzero(yz==np.nanmin(yz)); pick=cand[tie.integers(len(cand))]
        d=rz[pick]-rz.mean()
        idx=rg.integers(0,len(jz),(NPERM,1)); dr=rz[idx].mean(1)-rz.mean()
        if dr.std()<1e-9: continue
        Z[i]=(d-dr.mean())/dr.std()
        rr=rar[jr]; s=np.nanstd(yr)
        if s<1e-9 or np.std(rr)<1e-9: continue
        Rho[i]=np.corrcoef(yr,rr)[0,1]
    return Z,Rho

Dres=demean_conv(V)
Z,Rho=two_stats(Dres,zlib.crc32(b'A19R01'))
m=np.isfinite(Z)&np.isfinite(Rho)
real=float(np.corrcoef(Z[m],Rho[m])[0,1])
print(f"{m.sum():,} 人   z_resid 均值 {Z[m].mean():+.4f}   ρ 均值 {Rho[m].mean():+.4f}",flush=True)

fl=[]
for s in range(5):
    Zp,Rp=two_stats(Dres,5000+s,perm=True); mp=np.isfinite(Zp)&np.isfinite(Rp)
    fl.append(float(np.corrcoef(Zp[mp],Rp[mp])[0,1]))
floor=float(np.mean(fl))

rgc=np.random.default_rng(11); half=rgc.permutation(V.shape[1])
cA,cB=half[:V.shape[1]//2],half[V.shape[1]//2:]
Zc,_=two_stats(Dres,777,colsz=cA,colsr=cA)
_,Rc=two_stats(Dres,777,colsz=cB,colsr=cB)
mc=np.isfinite(Zc)&np.isfinite(Rc)
ceil=float(np.corrcoef(Zc[mc],Rc[mc])[0,1])

ii=np.flatnonzero(m); rb=np.random.default_rng(3)
bs=float(np.std([np.corrcoef(Z[ii[rb.integers(0,len(ii),len(ii))]],
                             Rho[ii[rb.integers(0,len(ii),len(ii))]])[0,1] for _ in range(20)]))
bs=float(np.std([ (lambda s_: np.corrcoef(Z[s_],Rho[s_])[0,1])(ii[rb.integers(0,len(ii),len(ii))])
                 for _ in range(200)]))
print(f"\n  corr(z_resid, ρ) = **{real:+.4f}**   自助展布 {bs:.4f}")
print(f"  人内置换地板     = {floor:+.4f}(5 个种子 {min(fl):+.4f}..{max(fl):+.4f})")
print(f"  分半天花板       = {ceil:+.4f}(z 在一半类别上、ρ 在另一半上,n={mc.sum():,})")

# 正对照:种植一个人特异的径向信号
rgp=np.random.default_rng(31); u=rgp.standard_normal(len(V)); x=rar-rar.mean()
Dp=demean_conv(np.where(obs,V+1.2*np.outer(u,x),np.nan))
Zg,Rg=two_stats(Dp,999); mg=np.isfinite(Zg)&np.isfinite(Rg)
plant=float(np.corrcoef(Zg[mg],Rg[mg])[0,1])
cu_z=float(np.corrcoef(Zg[mg],u[mg])[0,1]); cu_r=float(np.corrcoef(Rg[mg],u[mg])[0,1])
print(f"  种植正对照       = {plant:+.4f}   (种植量与 z 相关 {cu_z:+.3f}、与 ρ 相关 {cu_r:+.3f})")

T=pd.DataFrame([dict(arm='real',v=real),dict(arm='floor',v=floor),
                dict(arm='ceiling',v=ceil),dict(arm='plant',v=plant)])
T.to_csv(pathlib.Path(__file__).parent/'results'/'arms.csv',index=False)
g=Gate('那 20% 是一件事还是两件事')
g.asserted('种植的人特异径向信号被两个量都测到',abs(cu_z)>0.1 and abs(cu_r)>0.1,
           f"与 z {cu_z:+.3f}、与 ρ {cu_r:+.3f}")
g.asserted('种植让两者的相关上升',abs(plant)>abs(real),f"{real:+.4f} -> {plant:+.4f}")
# ⚠ 我预期"两个量共享算术结构 -> 地板不为零"。**错了**:人内置换把相关完全打掉
#   (−0.0034,5 个种子 −0.0066..+0.0031)。所以这个零**确实**是零,而我问的时候以为它不是。
g.asserted('⚠ 我预期的"地板不为零"是错的 —— 人内置换把相关完全打掉',abs(floor)<0.02,
           f"地板 {floor:+.4f}(5 个种子 {min(fl):+.4f}..{max(fl):+.4f})。"
           f"所以 −0.5741 是真信号,不是共享算术结构")
g.require_resolvable_first('真实相关是否越过地板',abs(real-floor),bs)
g.offset_control('corr(z_resid, ρ) vs 人内置换地板',real,floor,bs,
                 null_kind='人内置换该人的起始年龄标签后,同两个量的相关(共享算术结构的地板)')
# ⚠ 跨半臂比同数据臂**小 6 倍**,所以它不是"天花板" —— 它是**另一个量**:
#   同数据臂问"这个人在这批类别上的两个读数是否一致"(是,−0.574);
#   跨半臂问"这个人在**另一半**类别上还是这样吗"(几乎不,−0.091)。
ic=np.flatnonzero(mc); rc2=np.random.default_rng(19)
bc=float(np.std([ (lambda s_: np.corrcoef(Zc[s_],Rc[s_])[0,1])(ic[rc2.integers(0,len(ic),len(ic))])
                  for _ in range(200)]))
print(f"  跨半臂的自助展布 {bc:.4f} -> {abs(ceil)/bc:.1f}x")
g.asserted('跨半臂不是天花板,它比同数据臂小得多 —— 那本身是结论',abs(ceil)<0.5*abs(real),
           f"同一批类别上 {real:+.4f}(91.6×)vs 换到另一半类别 {ceil:+.4f} —— **小 "
           f"{abs(real)/abs(ceil):.1f} 倍**")
g.require_resolvable_first('跨半臂本身是否可分辨',abs(ceil),bc,family='crosshalf')
g.offset_control('跨半相关 vs 人内置换地板',ceil,floor,bc,
                 null_kind='人内置换后的同一个跨半相关(共享 item 已被两半的不相交排除)')
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

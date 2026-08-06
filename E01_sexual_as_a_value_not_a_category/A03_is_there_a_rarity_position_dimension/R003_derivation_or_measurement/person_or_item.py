import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A21 R01 -- 那 20% 是谁的?人的,还是东西的?

#148a:「最早那批更常见」的大部分是共享发育时间表(五种规格下 61–104%)。
#149d:剩下的个人成分**是真的但不是一个特质** —— 同一批兴趣上两个读数一致 −0.574,
       换到这个人**另一半**兴趣只剩 −0.091。**#149 只证明了它不稳,没问为什么。**

双向去均值之后两个边际都为零,所以"是不是关于东西的"只能有一个意思:
**人×题目的交互,是不是集中在特定的类别上?**

  PERSON  z 自身的分半信度像样 -> 它终究是个人性质,#149 的 −0.091 是跨统计量的衰减
  ITEM    z 的分半接近零,但**每类别载荷 λ_j 的跨人分半可复现** -> 效应由**哪些类别**承载,
          不由**谁**承载。那 20% 有了名字:它是题目侧的,不是人侧的
  NEITHER 两者都低 -> 它是细胞级噪声,聚合出一个真实的均值,却没有稳定的归属者

ESTIMAND        ① z 的分半信度(把一个人的类别劈成两半,各算一次 z,跨人相关,SB 校正)
                ② 每类别载荷 λ_j = corr_人(d_ij, z_i);把**人**劈成两半,两套 λ 的相关
IDENTIFICATION  ① 的两半在**类别**上不相交;② 的两半在**人**上不相交。两条边分开量。
SCOPE           >=10 个类别起始年龄的人(两半各 >=5)。
WORLDS          PERSON / ITEM / NEITHER
KILL            条件式:两个正对照必须开火 —— 种植一个**人**特异信号,①必须升;
                种植一个**类别**特异信号,②必须升。否则两条边都不可读。
POSITIVE CTRL   见上,两个。
NEGATIVE CTRL   人内置换起始年龄标签,两条边各跑一次。
NOISE FLOOR     200 次按人自助;5 个劈分种子。
MULTIPLICITY    2 条边 x {真实, 人内置换, 两种种植},整格发表。
IMPOSSIBLE      "类别载荷可复现"不等于"因果在类别上" —— 它只说方差的归属,不说机制。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

def demean_conv(Vm,ob,tol=1e-10,cap=500):
    D=np.where(ob,Vm,np.nan)
    for _ in range(cap):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-a
        b=np.nanmean(D,axis=1,keepdims=True); D=D-b
        if max(np.nanmax(np.abs(a)),np.nanmax(np.abs(b)))<tol: break
    return D
NPERM=200
KEEP2=obs.sum(1)>=10
def zvec(D,cols_mask,seed,perm=False):
    """每人在给定类别子集上的 z(最早一个,m=1,人内置换归一化)。"""
    rg=np.random.default_rng(seed); tie=np.random.default_rng(20260803)
    Z=np.full(len(D),np.nan)
    for i in np.flatnonzero(KEEP2):
        j=np.flatnonzero(obs[i]&cols_mask[i]); k=len(j)
        if k<5: continue
        y=D[i,j].copy(); r=rar[j]
        if perm: y=y[rg.permutation(k)]
        cand=np.flatnonzero(y==np.nanmin(y)); pick=cand[tie.integers(len(cand))]
        d=r[pick]-r.mean()
        idx=rg.integers(0,k,(NPERM,1)); dr=r[idx].mean(1)-r.mean()
        if dr.std()<1e-9: continue
        Z[i]=(d-dr.mean())/dr.std()
    return Z

def edge_person(D,seed,perm=False):
    """边①:把每人的类别劈成两半,各算 z,跨人相关 + Spearman-Brown。"""
    rg=np.random.default_rng(seed)
    mA=np.zeros_like(obs); mB=np.zeros_like(obs)
    for i in np.flatnonzero(KEEP2):
        j=np.flatnonzero(obs[i]); p=rg.permutation(j); h=len(j)//2
        mA[i,p[:h]]=True; mB[i,p[h:]]=True
    a=zvec(D,mA,seed+1,perm); b=zvec(D,mB,seed+2,perm)
    m=np.isfinite(a)&np.isfinite(b)
    r=float(np.corrcoef(a[m],b[m])[0,1])
    return (2*r/(1+r) if r>-0.99 else np.nan), r, int(m.sum())

def edge_item(D,seed,perm=False):
    """边②:每类别载荷 λ_j = corr_人(d_ij, z_i);把**人**劈成两半,两套 λ 相关。"""
    rg=np.random.default_rng(seed)
    Z=zvec(D,obs.copy(),seed+3,perm)
    who=np.flatnonzero(np.isfinite(Z)); p=rg.permutation(who); h=len(p)//2
    lam=[]
    for half in (p[:h],p[h:]):
        v=np.full(V.shape[1],np.nan)
        for j in range(V.shape[1]):
            m=half[obs[half,j]&np.isfinite(D[half,j])]
            if len(m)<200: continue
            v[j]=np.corrcoef(D[m,j],Z[m])[0,1]
        lam.append(v)
    m=np.isfinite(lam[0])&np.isfinite(lam[1])
    return float(np.corrcoef(lam[0][m],lam[1][m])[0,1]), int(m.sum())

D0=demean_conv(V,obs)
print(f"{KEEP2.sum():,} 人(>=10 个类别)",flush=True)
rows=[]
for tag,Dx,perm in [('real',D0,False),('人内置换',D0,True)]:
    e1=[edge_person(Dx,100+s,perm) for s in range(5)]
    e2=[edge_item(Dx,200+s,perm) for s in range(5)]
    rows.append(dict(arm=tag,sb=float(np.nanmean([x[0] for x in e1])),
                     r_person=float(np.nanmean([x[1] for x in e1])),
                     n_person=int(np.mean([x[2] for x in e1])),
                     r_item=float(np.nanmean([x[0] for x in e2])),
                     n_item=int(np.mean([x[1] for x in e2]))))
    print(f"  {tag:<8} 边①人侧 r={rows[-1]['r_person']:+.4f}(SB {rows[-1]['sb']:+.4f})  "
          f"边②题目侧 r={rows[-1]['r_item']:+.4f}",flush=True)

# 正对照一:种植**人**特异的径向信号 -> 边①必须升
rgp=np.random.default_rng(31); u=rgp.standard_normal(len(V)); x=rar-rar.mean()
Dp=demean_conv(np.where(obs,V+1.5*np.outer(u,x),np.nan),obs)
e1p=float(np.nanmean([edge_person(Dp,300+s)[1] for s in range(3)]))
e2p=float(np.nanmean([edge_item(Dp,310+s)[0] for s in range(3)]))
# 正对照二:种植**类别**特异的信号 -> 边②必须升
rgq=np.random.default_rng(77); w=rgq.standard_normal(V.shape[1]); uu=rgq.standard_normal(len(V))
Dq=demean_conv(np.where(obs,V+1.5*np.outer(uu,w),np.nan),obs)
e1q=float(np.nanmean([edge_person(Dq,400+s)[1] for s in range(3)]))
e2q=float(np.nanmean([edge_item(Dq,410+s)[0] for s in range(3)]))
print(f"\n  正对照一(种植**人**特异径向):边① {e1p:+.4f}  边② {e2p:+.4f}")
print(f"  正对照二(种植**类别**特异载荷):边① {e1q:+.4f}  边② {e2q:+.4f}")

# ---- 判据:边②的 0.96 是不是**平凡**的?λ_j 可能只是 rar_j 的确定性函数,
#      而 rar_j 在两半人里当然相同 —— 那样"可复现"什么也没说。
rgl=np.random.default_rng(555)
Zf=zvec(D0,obs.copy(),999)
whof=np.flatnonzero(np.isfinite(Zf))
lam=np.full(V.shape[1],np.nan)
for j in range(V.shape[1]):
    m=whof[obs[whof,j]&np.isfinite(D0[whof,j])]
    if len(m)<200: continue
    lam[j]=np.corrcoef(D0[m,j],Zf[m])[0,1]
ok=np.isfinite(lam)
c_rar=float(np.corrcoef(lam[ok],rar[ok])[0,1])
# 把 rar 从 λ 里回归掉,再看两半的残差 λ 还相不相关
def lam_half(seed,resid=False):
    rg=np.random.default_rng(seed); p_=rg.permutation(whof); h=len(p_)//2
    out=[]
    for half in (p_[:h],p_[h:]):
        v=np.full(V.shape[1],np.nan)
        for j in range(V.shape[1]):
            m=half[obs[half,j]&np.isfinite(D0[half,j])]
            if len(m)<200: continue
            v[j]=np.corrcoef(D0[m,j],Zf[m])[0,1]
        if resid:
            o=np.isfinite(v)
            Zr=np.c_[np.ones(o.sum()),rar[o],rar[o]**2]
            v[o]=v[o]-Zr@np.linalg.lstsq(Zr,v[o],rcond=None)[0]
        out.append(v)
    m=np.isfinite(out[0])&np.isfinite(out[1])
    return float(np.corrcoef(out[0][m],out[1][m])[0,1])
raw_r=float(np.nanmean([lam_half(700+s) for s in range(5)]))
res_r=float(np.nanmean([lam_half(700+s,resid=True) for s in range(5)]))
print(f"\n=== 判据:边②是不是平凡的 ===")
print(f"  corr(λ_j, 稀有度) = **{c_rar:+.4f}**(λ 有多像稀有度的函数)")
print(f"  两半 λ 的相关:原始 {raw_r:+.4f}  ->  **把 rar 与 rar² 回归掉之后 {res_r:+.4f}**")
g2=Gate('边②的可复现性是不是平凡的')
g2.asserted('λ 与稀有度的相关有多高',True,f"{c_rar:+.4f}")
g2.require_resolvable_first('去掉稀有度后,λ 的跨人可复现性还剩多少',abs(res_r),0.05)
g2.offset_control('去稀有度后的 λ 复现 vs 原始',res_r,raw_r,0.05,
                  null_kind='未去稀有度的同一个两半 λ 相关(它含 rar_j 在两半里恒同这一平凡成分)')
# 把它变成一句关于内容的话:去掉稀有度之后,哪些类别的载荷最极端
import re as _re
lab=[_re.sub(r'\s*\([a-z0-9]+\)$','',c) for c in ons]
lab=[_re.sub(r'^.*?(?:interest in|interested in)\s*','',l)[:30] for l in lab]
o=np.isfinite(lam); Zr=np.c_[np.ones(o.sum()),rar[o],rar[o]**2]
lr=np.full(len(lam),np.nan); lr[o]=lam[o]-Zr@np.linalg.lstsq(Zr,lam[o],rcond=None)[0]
oo=np.argsort(lr)
print("\n=== 去掉稀有度之后,载荷最极端的类别(λ 残差)===")
print("  比时间表**更早**到来的一侧:")
for j in oo[:5]:
    if np.isfinite(lr[j]): print(f"     {lr[j]:+.3f}  {lab[j]}(基率 {prev[j]:.0%},题目均值起始 "
                                 f"{np.nanmean(V[obs[:,j],j]):.1f} 岁)")
print("  比时间表**更晚**到来的一侧:")
for j in oo[::-1][:5]:
    if np.isfinite(lr[j]): print(f"     {lr[j]:+.3f}  {lab[j]}(基率 {prev[j]:.0%},题目均值起始 "
                                 f"{np.nanmean(V[obs[:,j],j]):.1f} 岁)")
g2.asserted('⚠ 边②没有开火的正对照(它已在天花板 0.96),所以它的读数靠上面那个平凡性检验',
            True,'正对照二 0.9616 -> 0.9753,升不动。**记录,不当作它通过了**')
print(g2)

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'edges.csv',index=False)
R=T[T.arm=='real'].iloc[0]; N=T[T.arm=='人内置换'].iloc[0]
rb=np.random.default_rng(9)
g=Gate('那 20% 是人的还是东西的')
g.asserted('正对照一:种植人特异信号,边①升',e1p>R.r_person+0.05,
           f"真实 {R.r_person:+.4f} -> 种植 {e1p:+.4f}")
g.asserted('正对照二:种植类别特异信号,边②升',e2q>R.r_item+0.05,
           f"真实 {R.r_item:+.4f} -> 种植 {e2q:+.4f}")
g.negative_control('边①的人内置换零',float(N.r_person),float(R.r_person))
g.negative_control('边②的人内置换零',float(N.r_item),float(R.r_item))
g.asserted('两条边的相对大小,就是答案',True,
           f"人侧 r={R.r_person:+.4f}(SB {R.sb:+.4f})  vs  题目侧 r={R.r_item:+.4f}")
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

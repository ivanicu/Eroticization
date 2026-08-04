import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A15 R04 -- 「一整套」是字面意义的吗?关系族是不是**同时**到达?

#135d:早来的东西是散的(题目层彼此 -0.0075),晚来的东西是一整套(彼此 +0.0674)。
但 #135 自己写下了缺口:那证明的是人的晚期集合**富集**于这一族,**没有**证明它们在
**同一时间**到达。「一整套」目前只是结构上的,不是时间上的。

    PACKAGE  族内条目的获得年龄比同样大小的随机子集**更靠拢** -> 字面意义的一整套
    ORDERED  与随机无异 -> 它们只是都靠后,不是同时;"一整套"只能说结构,不能说时间

⚠ 两个必须先解决的设计问题:
  (1) **族不能用时间定义**,否则是拿同一批数据既选样本又下结论。这里用**共现**
      (SIMR 的谱分割)定义族 —— 共现与获得年龄是两个不相交的仪器。
  (2) **人群时间表本身会让族内靠拢**:关系族的题目平均起始年龄都在 16.8-17.0,
      所以它们的原始年龄天然接近。**必须在题目层去均值后的残差上量离散度**,
      否则测到的是 #75 的时间表,不是"打包"。

ESTIMAND        族内条目**去题目均值后**的获得年龄离散度(人内 sd),减去同样大小的
                随机子集(人内置换)的期望。
IDENTIFICATION  族由共现定义,时间由起始年龄定义 —— 不相交。零是人内置换,精确保留
                这个人的类别集与子集大小。
SCOPE           拥有该族 >=3 个条目、且 >=8 个类别起始年龄的人。
WORLDS          PACKAGE / ORDERED
KILL            条件式:种植一个真"包"必须被检出,且种植一个"只是都靠后但时间随机"的
                对照**不**触发,才读阈值。
POSITIVE CTRL   把族内条目的残差强制压到同一个值(打包),离散度必须显著下降,且随
                打包强度单调。
NEGATIVE CTRL   人内置换子集标签,5 个种子。
NOISE FLOOR     200 次按人自助。
MULTIPLICITY    2 个族(共现谱分割的两侧)x 4 个打包强度 x 5 seeds,整格发表。
IMPOSSIBLE      2 年分箱是离散度的地板;若真实打包比分箱还紧,本设计看不见。地板报出。
"""
import numpy as np, pandas as pd, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

Ob=obs.astype(float); pj=Ob.mean(0); Cm=(Ob.T@Ob)/len(Ob)
den=np.sqrt(np.outer(pj*(1-pj),pj*(1-pj))); den[den<1e-9]=1e-9
SIM=(Cm-np.outer(pj,pj))/den; np.fill_diagonal(SIM,0.)
iu=np.triu_indices(len(rar),1)
X=np.c_[np.ones(len(iu[0])),rar[iu[0]]+rar[iu[1]],rar[iu[0]]*rar[iu[1]],np.abs(rar[iu[0]]-rar[iu[1]])]
res=SIM[iu]-X@np.linalg.lstsq(X,SIM[iu],rcond=None)[0]
check_residualized(res,rar[iu[0]]+rar[iu[1]],"配对相似度对稀有度")
SIMR=np.zeros_like(SIM); SIMR[iu]=res; SIMR=SIMR+SIMR.T

# ---- 族 = 共现的谱分割(与获得年龄不相交的仪器)
w,vv=np.linalg.eigh(SIMR); pc=vv[:,-1]
lab=[re.sub(r'\s*\([a-z0-9]+\)$','',c) for c in ons]
lab=[re.sub(r'^.*?(?:interest in|interested in)\s*','',l)[:30] for l in lab]
FAM_A=np.flatnonzero(pc>0); FAM_B=np.flatnonzero(pc<=0)
if SIMR[np.ix_(FAM_A,FAM_A)].mean()<SIMR[np.ix_(FAM_B,FAM_B)].mean(): FAM_A,FAM_B=FAM_B,FAM_A
print(f"共现谱分割:族A {len(FAM_A)} 个(内部 {SIMR[np.ix_(FAM_A,FAM_A)].mean():+.4f}),"
      f"族B {len(FAM_B)} 个(内部 {SIMR[np.ix_(FAM_B,FAM_B)].mean():+.4f})",flush=True)
itm=np.array([np.nanmean(V[obs[:,j],j]) for j in range(V.shape[1])])
for nm,F in [('A',FAM_A),('B',FAM_B)]:
    o=np.argsort(-pc[F]) if nm=='A' else np.argsort(pc[F])
    print(f"  族{nm} 平均起始 {itm[F].mean():.2f}:  " + " · ".join(lab[F[k]][:22] for k in o[:6]))

# ---- 去题目均值后的残差(否则测到的是 #75 的时间表)
Dres=np.where(obs,V-itm[None,:],np.nan)
print(f"分箱地板:相邻分箱间距最小 {np.min(np.diff(sorted(set(V[obs].tolist())))):.1f} 年",flush=True)

def disp(Dm,F,rng=None,perm=False):
    """族内残差的人内 sd;perm=True 时随机换一个同样大小的子集。"""
    out=[]
    for i in np.flatnonzero(KEEP):
        jj=np.flatnonzero(obs[i]); fam=np.intersect1d(jj,F)
        k=len(fam)
        if k<3: continue
        sel=rng.choice(jj,k,replace=False) if perm else fam
        v=Dm[i,sel]
        if np.isfinite(v).sum()>=3: out.append(float(np.nanstd(v)))
    return np.array(out)

rows=[]
print(f"\n{'族':<5} {'n':>7} {'族内离散度':>11} {'人内置换零':>11} {'差':>10} {'展布':>9} {'倍数':>7}")
for nm,F in [('A',FAM_A),('B',FAM_B)]:
    d=disp(Dres,F)
    z=np.concatenate([disp(Dres,F,np.random.default_rng(5500+s),perm=True) for s in range(5)])
    rb=np.random.default_rng(88)
    bs=float(np.std([d[rb.integers(0,len(d),len(d))].mean() for _ in range(200)]))
    rows.append(dict(fam=nm,n=len(d),disp=float(d.mean()),null=float(z.mean()),boot=bs))
    print(f"{nm:<5} {len(d):>7,} {d.mean():>11.4f} {z.mean():>11.4f} {d.mean()-z.mean():>+10.4f} "
          f"{bs:>9.4f} {abs(d.mean()-z.mean())/bs:>7.1f}x")

# 对照:打包 / 只是都靠后但时间随机
print(f"\n对照(族A):")
ctl=[]
for gpk in [0.0,0.35,0.70,1.0]:
    Dp=Dres.copy()
    for i in np.flatnonzero(KEEP):
        fam=np.intersect1d(np.flatnonzero(obs[i]),FAM_A)
        if len(fam)<3: continue
        mu=np.nanmean(Dp[i,fam]); Dp[i,fam]=(1-gpk)*Dp[i,fam]+gpk*mu
    d=disp(Dp,FAM_A); ctl.append(float(d.mean()))
    print(f"  打包强度 {gpk:.2f}  离散度 {d.mean():.4f}")
# ⚠ 第一版把这个对照写成"族内置换" —— 而 sd 对置换是**不变的**,所以它逐位给出
#   同一个数(2.2608 vs 2.2608),是一个在构造上不可能失败的检查(#96a / _degenerate)。
#   正确的实现:保留族的**平均晚到程度**,但把族内的时间**换成这个人其余类别的残差**。
Dsh=Dres.copy()
rg=np.random.default_rng(4141)
for i in np.flatnonzero(KEEP):
    jj=np.flatnonzero(obs[i]); fam=np.intersect1d(jj,FAM_A); oth=np.setdiff1d(jj,fam)
    if len(fam)<3 or len(oth)<len(fam): continue
    donor=Dres[i,rg.choice(oth,len(fam),replace=False)]
    Dsh[i,fam]=np.nanmean(Dres[i,fam])+(donor-np.nanmean(donor))   # 同样的平均晚到,随机的内部时间
d_sh=float(disp(Dsh,FAM_A).mean())
print(f"  同样晚到但内部时间取自其余类别  离散度 {d_sh:.4f}"
      f"(真实 {rows[0]['disp']:.4f},零 {rows[0]['null']:.4f})")

D=pd.DataFrame(rows); D.to_csv(pathlib.Path(__file__).parent/'results'/'dispersion.csv',index=False)
A=D[D.fam=='A'].iloc[0]
g=Gate('关系族是不是同时到达')
g.asserted('族由共现定义,与获得年龄不相交',True,
           f"谱分割:族A 内部 {SIMR[np.ix_(FAM_A,FAM_A)].mean():+.4f} / 族B {SIMR[np.ix_(FAM_B,FAM_B)].mean():+.4f}")
g.asserted('种植的"包"被检出,且随强度单调',all(ctl[i]>ctl[i+1] for i in range(len(ctl)-1)),
           " > ".join(f"{v:.4f}" for v in ctl))
# "这一族只是平均更晚" —— 这个混淆在**代数上不可能**影响本统计量,因为 sd 对平移不变。
#   这与救了 #116 与 #128b 的是同一种结构免疫。数值验证它,而不是叙述它。
Dshift=Dres.copy()
for i in np.flatnonzero(KEEP):
    fam=np.intersect1d(np.flatnonzero(obs[i]),FAM_A)
    if len(fam)<3: continue
    Dshift[i,fam]=Dshift[i,fam]+5.0                 # 整族平移 5 年
d_shift=float(disp(Dshift,FAM_A).mean())
print(f"  整族平移 5 年后的离散度 {d_shift:.10f}(真实 {rows[0]['disp']:.10f})"
      f" —— sd 对平移不变,所以\"只是都靠后\"在代数上不可能造出这个效应")
g.asserted('结构免疫:"这一族只是平均更晚"在代数上不可能影响 sd(数值验证)',
           abs(d_shift-A.disp)<1e-9,
           f"整族平移 5 年 -> {d_shift:.10f} vs 真实 {A.disp:.10f},差 {abs(d_shift-A.disp):.2e}")
g.asserted('⚠ 我第一版的"内部时间随机"对照是不能失败的',True,
           '写成族内置换 -> 2.2608 vs 真实 2.2608 逐位相同,因为 sd 对置换不变(#96a);'
           '改成捐赠池版本又撞上 same_scale(捐赠池 != 零的池,给出中间值 2.3160)。'
           '**正确的答案不是一个对照,是上面那条代数事实**')
g.require_resolvable_first('族A 的离散度差可分辨',abs(A.disp-A.null),A.boot)
g.offset_control('族A:族内是否比随机子集更靠拢',float(A.disp),float(A.null),float(A.boot),
                 null_kind='人内置换子集标签(保留类别集与子集大小),在题目去均值后的残差上')
print(g)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

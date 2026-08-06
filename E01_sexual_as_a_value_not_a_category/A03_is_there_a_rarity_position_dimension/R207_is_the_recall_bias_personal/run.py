import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A21 R04 -- #114 的回忆偏差,是人群的还是个人的?

#160c/#161b:性偏好的「什么」是个人的(0.62),「何时」不是 —— 无论问
「时间表贴合你多少」(0.16)还是「你偏离多远」(0.11),都低一个数量级。

那么 `#114`(**人把最爱的性兴趣记得更早,约九个月**,−0.2000 年/评分 sd,19.8×)
必须被重新定位:**它是一条关于人群的规律,还是关于个人的?**

    POPULATION  按人估计的斜率 b_i 分半信度 ~0.1 -> 回忆偏差和「何时」的其余部分一样,
                是人群层的。**「某些人比别人更会把最爱的记早」这句话不成立**
    PERSONAL    ~0.6 -> 「何时」里**有**一个人层成分,只是它不在**顺序**里而在**记忆**里。
                那会是本会话第一次在"何时"这一侧找到个人性质

ESTIMAND        b_i = 这个人**人内**的「起始年龄残差 ~ 该类别自己的评分」斜率;
                它的分半信度(把这个人的类别劈成不相交两半,各估一次,跨人相关 + SB),
                在与 `#160`/`#161` **同一台机器**、同一个 k 上。
IDENTIFICATION  b_i 是人内斜率,所以人层的整体早熟与整体评分水平恰好抵消;
                双向去均值已移除题目固定效应(= `#75` 的时间表)与人的整体早熟。
SCOPE           两半各 >=k 个**同时有起始年龄与匹配评分**的类别的人。
WORLDS          POPULATION / PERSONAL
KILL            条件式:人特异种植必须把 b 的信度推上去(证明仪器能测到人层的回忆偏差),
                且人内置换零必须为零,才读真实值。
POSITIVE CTRL   种植一个**人特异**的回忆偏差(每人一个自己的斜率),信度必须升。
NEGATIVE CTRL   人内置换评分标签。
NOISE FLOOR     5 个劈分种子。
MULTIPLICITY    k ∈ {4,5,6} x {真实, 置换, 种植},整格发表。
IMPOSSIBLE      斜率的人内估计天然比均值噪声大,所以低信度**部分**是估计噪声 ——
                这正是种植臂要回答的:仪器能不能把一个**真的**人层斜率测出来。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A21_who_owns_the_personal_20_percent'
          /'R205_what_versus_when_same_design'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('KS=[4,5,6,7]')[0])   # 跨轮依赖显式声明(P16)

RM=np.full_like(V,np.nan)
for j,ri in best.items(): RM[:,j]=R[:,ri]
zr=(RM-np.nanmean(RM))/np.nanstd(RM)
usable=obs&np.isfinite(zr)
print(f"同时有起始年龄与匹配评分的格子:{usable.sum():,};覆盖 {len(best)}/{V.shape[1]} 个类别",flush=True)

def slope_on(cols,i,perm,rg,plant=0.,u=0.):
    j=np.intersect1d(np.flatnonzero(usable[i]),cols)
    if len(j)<3: return np.nan
    y=D0[i,j].copy(); x=zr[i,j].copy()
    if perm: x=x[rg.permutation(len(x))]
    if plant: y=y+plant*u*x                       # 人特异的回忆偏差斜率
    x=x-x.mean(); y=y-np.nanmean(y)
    v=float((x*x).sum())
    return float((x*y).sum()/v) if v>1e-9 else np.nan

def rel_b(k,seed,perm=False,plant=0.):
    rg=np.random.default_rng(seed); u=rg.standard_normal(len(V))
    A=[];B=[]
    for i in np.flatnonzero(KEEP&okA):
        av=np.flatnonzero(usable[i])
        if len(av)<2*k: continue
        p=rg.permutation(av); h1,h2=p[:k],p[k:2*k]
        x1=slope_on(h1,i,perm,rg,plant,u[i]); x2=slope_on(h2,i,perm,rg,plant,u[i])
        if np.isfinite(x1) and np.isfinite(x2): A.append(x1); B.append(x2)
    A=np.array(A); B=np.array(B)
    if len(A)<300: return np.nan,np.nan,len(A)
    r=float(np.corrcoef(A,B)[0,1])
    return (2*r/(1+r) if r>-0.99 else np.nan), r, len(A)

# 先复现 #114 的人群层斜率,作为正对照
f=usable&np.isfinite(D0)
pop_slope=float(np.polyfit(zr[f],D0[f],1)[0])
print(f"人群层斜率(复现 #114 的 −0.2000):{pop_slope:+.4f} 年/评分 sd",flush=True)

KS=[4,5,6]; rows=[]
print(f"\n{'k':<4}{'n':>7}{'分半 r':>10}{'SB 信度':>10}{'置换零(5 抽)':>14}{'种植':>9}")
for k in KS:
    sd_=zlib.crc32(f'b{k}'.encode())%9973
    sb,r,n=rel_b(k,sd_)
    # ⚠ #162:第一版每个 k 只抽**一次**置换零,k=5 那一抽给 +0.106 而 k=4/6 给 −0.03/−0.05 ——
    #   正是 `#153` 那条教训("单次抽样的零不自动致命,但伤害 = 零的实现 sd ÷ 效应")
    #   打在我自己身上。改成 5 抽,报分布。
    nulls=np.array([rel_b(k,sd_+100+t_,perm=True)[0] for t_ in range(5)])
    sbp,_,_=rel_b(k,sd_+2,plant=1.5)
    rows.append(dict(k=k,n=n,r=r,sb=sb,sb_null=float(np.nanmean(nulls)),
                     sb_null_sd=float(np.nanstd(nulls)),sb_null_lo=float(np.nanmin(nulls)),
                     sb_null_hi=float(np.nanmax(nulls)),sb_plant=sbp))
    print(f"{k:<4}{n:>7,}{r:>+10.4f}{sb:>+10.4f}"
          f"{np.nanmean(nulls):>+8.4f}±{np.nanstd(nulls):.3f}{sbp:>9.4f}",flush=True)

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'recall_rel.csv',index=False)
mb=float(T.sb.mean()); sd=float(T.sb.std())
print(f"\n  b(回忆偏差斜率)的人侧信度均值 **{mb:+.4f}**")
print(f"  同一台机器上的参照:S(喜欢什么)0.62 · a(时间表贴合)0.16 · z(偏离多远)0.11")

g=Gate('#114 的回忆偏差是人群的还是个人的')
g.asserted('人群层斜率复现 #114',abs(pop_slope-(-0.2000))<0.06,
           f"{pop_slope:+.4f} vs #114 报的 −0.2000")
g.asserted('人特异种植把 b 的信度推上去(仪器能测到人层的回忆偏差)',
           bool((T.sb_plant>T.sb+0.15).all()),
           " ".join(f"k={int(r.k)}:{r.sb:.3f}->{r.sb_plant:.3f}" for _,r in T.iterrows()))
g.asserted('人内置换零(5 抽均值)在零附近',bool((T.sb_null.abs()<0.08).all()),
           " ".join(f"k={int(r.k)}:{r.sb_null:+.3f}±{r.sb_null_sd:.3f}" for _,r in T.iterrows()))
g.asserted('⚠ 零的实现展布本身有多大(#153 的规则)',True,
           " ".join(f"k={int(r.k)}:[{r.sb_null_lo:+.3f},{r.sb_null_hi:+.3f}]" for _,r in T.iterrows())
           + f" —— 零的 sd ÷ 效应 = {T.sb_null_sd.mean()/max(abs(mb),1e-9):.1f}x")
g.equivalent_within('b 与 z(0.114)是不是同一个量级(边界 0.10)',mb-0.114,sd,0.10)
g.require_resolvable_first('b 与 S(0.622)的差',abs(mb-0.622),sd)
g.offset_control('b vs S',mb,0.622,sd,
                 null_kind='同一台机器上 S 的人侧分半信度(#160/#161,不是零假设,是被比较的对象)')
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

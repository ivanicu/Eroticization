import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A14 R03 -- 一个人报告的**第一个**兴趣,比他自己后来的那些更罕见吗?

#128b 说罕见的兴趣来得更早(人内梯度,8.8x),#129e 说口味越罕见的人提前得越多(匹配后
3.1x)。两个都建立在**双向去均值后的连续残差**上 —— 而 #128 的第四十四个错正是去均值本身
把符号弄反过一次。所以这一轮问同一个假说的**离散版本**:

    在一个人自己的曲目库里,他最早报告的那些兴趣,是不是比**从他自己库里随机抽**的更罕见?

这个提法的全部价值在于它的**零**:零是这个人自己的曲目库。所以它
  - 对双向去均值的做法**完全免疫**(根本不去均值)
  - 对"这个人喜欢多少东西"**完全免疫**(库是他自己的)
  - 对"这个人整体早熟不早熟"**完全免疫**(只看顺序,不看年龄值)
  - 对覆盖度**完全免疫**(#5 的定律在这里无处着力)
一个不依赖任何一条我这一周踩过的坑的检验。

ESTIMAND        Delta_i = (这个人最早一格里那些类别的平均稀有度)
                        - (从他自己的类别集里随机抽同样多个的平均稀有度)。
                以及 rank_i = 最早一格的平均稀有度在他自己库内的分位(零 = 0.5)。
                第二个量:corr(Delta_i, S_i),按类别数匹配。
IDENTIFICATION  零 = **人内置换起始年龄标签**,精确保留这个人的曲目库与他的年龄分布,
                只摧毁"哪个兴趣配到哪个年龄"。这正是本轮要检验的配对。
SCOPE           报告 >=8 个类别起始年龄、且最早一格不是他全部类别的人。
                起始年龄按 2 年分箱,所以"最早一格"通常含数个类别 —— 统计量按集合定义。
WORLDS          start   Delta > 0:第一个东西就已经比后来的更罕见。世界 B 从"轨迹形状"
                        升级为"起点本身"。
                smooth  Delta = 0:最早的那些是从库里随机抽的,#128b 的梯度是平滑的,
                        没有一个被特殊标记的起点。
                common  Delta < 0:版图从常见的中心开始 —— 与 #128b 直接冲突,
                        那么两个统计量之一是伪影,而这是本轮最有价值的结果。
KILL            条件式:正对照必须开火**且**人内置换零必须为零,才读阈值。
POSITIVE CTRL   对一部分人,强制把他库里最罕见的类别搬进最早一格。Delta 与 corr 必须随
                种植比例单调。g=0 必须逐位复现真实臂。
NEGATIVE CTRL   人内置换起始年龄标签(见上),5 个种子。
CONFOUND        #114:人把最爱的记得更早。若罕见类别被这个人评分更高,回忆偏差会把它们
                拉进最早一格 -> 伪造 start。控制:把每个人的评分残差化掉后重排,报实际贡献。
                反向通路(#128f):罕见 = 近期获得 -> 推向更晚,与 start 方向相反。
NOISE FLOOR     200 次按人自助 + 5 个置换种子。
MULTIPLICITY    2 个统计量 x 4 个种植水平 x 5 seeds x {含/不含评分校正},整格发表。
IMPOSSIBLE      "第一个"最早只到 2 年分箱的第一格;真正的首次获得顺序不可见。
"""
import pandas as pd, numpy as np, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import (Gate, check_columns, check_coverage, check_residualized,
                       check_disjoint_items)

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

# ---- 每个人:最早一格 vs 他自己的库
def first_bin_stats(Vm,rng=None,perm=False):
    """perm=True -> 人内置换起始年龄标签(保留曲目库与年龄分布,只摧毁配对)。"""
    dlt=np.full(len(Vm),np.nan); rnk=np.full(len(Vm),np.nan)
    for i in np.flatnonzero(KEEP):
        m=obs[i]; y=Vm[i,m].copy(); r=rar[m]
        if perm: y=y[rng.permutation(len(y))]
        lo=y.min(); k=int((y==lo).sum())
        if k==0 or k>=len(y): continue            # 最早一格不能是他的全部类别
        obsv=r[y==lo].mean()
        # 零 = 从他自己的库里抽同样多个的期望与分位(闭式期望 + 蒙特卡洛分位)
        exp=r.mean()
        dlt[i]=obsv-exp
        draws=np.array([r[np.random.default_rng(i*7+t).choice(len(r),k,replace=False)].mean()
                        for t in range(24)])
        rnk[i]=(draws<obsv).mean()+0.5*(draws==obsv).mean()
    return dlt,rnk

def report(dlt,rnk,tag):
    m=np.isfinite(dlt)&KEEP; ms=m&np.isfinite(S)
    return dict(tag=tag,n=int(m.sum()),delta=float(np.mean(dlt[m])),
                rank=float(np.mean(rnk[m])),corr_S=float(np.corrcoef(dlt[ms],S[ms])[0,1]))

rows=[]
d0,r0=first_bin_stats(V); rows.append(report(d0,r0,'real'))
for s_ in range(5):
    d,r=first_bin_stats(V,np.random.default_rng(8800+s_),perm=True)
    rows.append(report(d,r,'perm')); rows[-1]['seed']=s_
print(f"可用 {rows[0]['n']:,} 人",flush=True)

# 正对照:把库里最罕见的类别搬进最早一格
for g_ in [0.0,0.10,0.25,0.50]:
    rg=np.random.default_rng(9100); Vp=V.copy()
    if g_>0:
        pick=rg.random(len(V))<g_
        for i in np.flatnonzero(KEEP&pick):
            m=np.flatnonzero(obs[i]); j=m[np.argmax(rar[m])]
            Vp[i,j]=np.nanmin(V[i,obs[i]])
    d,r=first_bin_stats(Vp); rows.append(report(d,r,f'plant{g_}'))
    print(f"  plant {g_}",flush=True)

# #114 校正臂:先把这个人的评分通道从年龄里剥掉,再重排
RM=np.full_like(V,np.nan)
for j,ri in best.items(): RM[:,j]=R[:,ri]
zr=(RM-np.nanmean(RM))/np.nanstd(RM)
f=obs&np.isfinite(zr)
Vc=np.where(f,V-(-0.2000)*zr,np.where(obs,V,np.nan))     # #114 独立测得的系数
dc,rc=first_bin_stats(np.nan_to_num(Vc,nan=1e9)); rows.append(report(dc,rc,'rating_corrected'))

D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'
G=D.groupby('tag')[['delta','rank','corr_S']].mean()
print("\n=== 整格 ===")
print(G.round(4).to_string())

real=G.loc['real']; perm=G.loc['perm']
mm=np.isfinite(d0)&KEEP; ii=np.flatnonzero(mm)
rb=np.random.default_rng(555)
BS=np.array([[np.mean(d0[s_]),np.mean(r0[s_])] for s_ in
             (ii[rb.integers(0,len(ii),len(ii))] for _ in range(200))])
db,rkb=BS.std(0)
mms=mm&np.isfinite(S); jj=np.flatnonzero(mms)
cb=float(np.std([np.corrcoef(d0[s_],S[s_])[0,1] for s_ in
                 (jj[rb.integers(0,len(jj),len(jj))] for _ in range(200))]))
print(f"\n按人自助(200):Delta 展布 {db:.4f}   rank 展布 {rkb:.4f}   corr(Delta,S) 展布 {cb:.4f}")

# 按类别数匹配 corr(Delta,S)
NC=NCAT.astype(float); med=np.median(S[jj])
def matched(seed):
    rg=np.random.default_rng(seed)
    hi=jj[S[jj]>med]; lo=jj[S[jj]<=med]
    c=(NC-NC[jj].mean())/NC[jj].std(); used=np.zeros(len(S),bool); P=[]
    for a in hi[rg.permutation(len(hi))]:
        d_=np.abs(c[lo]-c[a]); d_[used[lo]]=np.inf; k=int(np.argmin(d_))
        if d_[k]<0.25: used[lo[k]]=True; P.append((a,lo[k]))
    P=np.array(P); sel=np.r_[P[:,0],P[:,1]]
    return float(np.corrcoef(d0[sel],S[sel])[0,1]), abs(NC[P[:,0]].mean()-NC[P[:,1]].mean())/NC[jj].std(), len(P)
mv=[matched(700+s) for s in range(5)]
mcorr=float(np.mean([x[0] for x in mv])); mbal=float(np.mean([x[1] for x in mv])); mn=int(np.mean([x[2] for x in mv]))
print(f"按类别数匹配后 corr(Delta,S) = {mcorr:+.4f}(未匹配 {real.corr_S:+.4f},保留 "
      f"{100*mcorr/real.corr_S:.0f}%)  残差 {mbal:.3f} sd  配对 {mn:,}")

# ---- 判据:Delta 与 #128b 反号。差别只有一个:Delta 里**人群时间表还在**,
#      而 #128b 是在类别固定效应被精确扣掉之后测的。把时间表也从 Delta 里扣掉。
def demean_conv(Vm,tol=1e-10,cap=500):
    Dm=np.where(obs,Vm,np.nan)
    for k in range(cap):
        a=np.nanmean(Dm,axis=0,keepdims=True); Dm=Dm-a
        b=np.nanmean(Dm,axis=1,keepdims=True); Dm=Dm-b
        if max(np.nanmax(np.abs(a)),np.nanmax(np.abs(b)))<tol: break
    return Dm
Dres=demean_conv(V)
dd,rr_=first_bin_stats(np.where(obs,Dres,1e9))
rows.append(report(dd,rr_,'schedule_removed'))
dn,rn_=first_bin_stats(np.where(obs,Dres,1e9),np.random.default_rng(9500),perm=True)
rows.append(report(dn,rn_,'schedule_removed_perm'))
md=np.isfinite(dd)&KEEP; id_=np.flatnonzero(md)
rbd=np.random.default_rng(1234)
db2=float(np.std([np.mean(dd[id_[rbd.integers(0,len(id_),len(id_))]]) for _ in range(200)]))
from scipy.stats import spearmanr
itm=np.array([np.nanmean(V[obs[:,j],j]) for j in range(V.shape[1])])
sp=spearmanr(rar,itm)
print(f"\n=== 判据:把人群时间表也扣掉,再问同一个问题 ===")
print(f"  原始年龄上          Delta {float(G.loc['real','delta']):+.4f}   rank {float(G.loc['real','rank']):.4f}")
print(f"  时间表扣掉后        Delta {np.mean(dd[md]):+.4f}   rank {np.mean(rr_[md]):.4f}"
      f"   置换零 {np.mean(dn[np.isfinite(dn)&KEEP]):+.4f}   展布 {db2:.4f}")
print(f"  题目层 稀有度 vs 平均起始年龄:Pearson {np.corrcoef(rar,itm)[0,1]:+.3f}  "
      f"Spearman {sp.statistic:+.3f} (p={sp.pvalue:.3f}, n=31)")

g=Gate('第一个兴趣是不是就已经更罕见')
g.degenerate_matches_reference('g=0 逐位复现 real',float(G.loc['plant0.0','delta']),float(real.delta))
mono=[float(G.loc[f'plant{q}','delta']) for q in [0.0,0.10,0.25,0.50]]
g.asserted('种植的 Delta 随比例单调',all(mono[i]<mono[i+1] for i in range(len(mono)-1)),
           " < ".join(f"{v:+.4f}" for v in mono))
g.require_resolvable_first('Delta 相对人内置换零可分辨',abs(real.delta-perm.delta),db)
g.offset_control('Delta 高于人内置换零',float(real.delta),float(perm.delta),db,
                 null_kind='人内置换起始年龄标签(保留曲目库与年龄分布,只摧毁配对)')
g.offset_control('rank 高于 0.5',float(real['rank']),float(perm['rank']),rkb,
                 null_kind='同上,人内置换')
g.artifact_cannot_explain('#114 的实际贡献不能解释 Delta',
                          float(real.delta-G.loc['rating_corrected','delta']),float(real.delta),db)
g.require_resolvable_first('匹配后 corr(Delta,S) 可分辨【原始年龄】',abs(mcorr),cb,family='raw_age')
g.asserted('匹配把类别数差压下去了',mbal<0.1,f'残差 {mbal:.3f} sd')
d_sr=float(np.mean(dd[md])); n_sr=float(np.mean(dn[np.isfinite(dn)&KEEP]))
# 第二个量:在**扣掉时间表**的 Delta 上重问 corr(.,S),并按类别数匹配 —— 这才是 #129e 的对应版本
ms2=md&np.isfinite(S); j2=np.flatnonzero(ms2)
raw2=float(np.corrcoef(dd[ms2],S[ms2])[0,1])
cb2=float(np.std([np.corrcoef(dd[s_],S[s_])[0,1] for s_ in
                  (j2[rbd.integers(0,len(j2),len(j2))] for _ in range(200))]))
def matched2(seed):
    rg=np.random.default_rng(seed); m2=np.median(S[j2])
    hi=j2[S[j2]>m2]; lo=j2[S[j2]<=m2]
    c=(NC-NC[j2].mean())/NC[j2].std(); used=np.zeros(len(S),bool); P=[]
    for a in hi[rg.permutation(len(hi))]:
        dd_=np.abs(c[lo]-c[a]); dd_[used[lo]]=np.inf; k=int(np.argmin(dd_))
        if dd_[k]<0.25: used[lo[k]]=True; P.append((a,lo[k]))
    P=np.array(P); sel=np.r_[P[:,0],P[:,1]]
    return float(np.corrcoef(dd[sel],S[sel])[0,1]), abs(NC[P[:,0]].mean()-NC[P[:,1]].mean())/NC[j2].std()
mv2=[matched2(600+s) for s in range(5)]
mc2=float(np.mean([x[0] for x in mv2])); mb2=float(np.mean([x[1] for x in mv2]))
nl2=float(np.corrcoef(dn[np.isfinite(dn)&ms2],S[np.isfinite(dn)&ms2])[0,1])
print(f"  扣掉时间表后 corr(Delta,S) {raw2:+.4f}   匹配后 {mc2:+.4f}(保留 {100*mc2/raw2:.0f}%)"
      f"   置换零 {nl2:+.4f}   展布 {cb2:.4f} -> {abs(mc2)/cb2:.1f}x")
g.require_resolvable_first('匹配后 corr(Delta,S) 可分辨【扣掉时间表】',abs(mc2),cb2,family='deschedule')
g.negative_control('扣掉时间表后 corr(Delta,S) 对置换零',nl2,mc2,null_spread=cb2)
# 机制:argmin 由**左尾**决定,而题目均值对左尾是盲的。直接测,不叙述。
from scipy.stats import spearmanr as _sp
p10=np.array([np.nanpercentile(V[obs[:,j],j],10) for j in range(V.shape[1])])
p50=np.array([np.nanpercentile(V[obs[:,j],j],50) for j in range(V.shape[1])])
s10=_sp(rar,p10); s50=_sp(rar,p50); smn=_sp(rar,itm)
print(f"\n=== 机制:argmin 由左尾决定,而题目均值对左尾是盲的(题目 n={len(itm)})===")
print(f"  Spearman(稀有度, 题目起始年龄的 10 分位) = {s10.statistic:+.3f}  p={s10.pvalue:.4f}")
print(f"  Spearman(稀有度, 中位数)                = {s50.statistic:+.3f}  p={s50.pvalue:.4f}")
print(f"  Spearman(稀有度, 均值)                  = {smn.statistic:+.3f}  p={smn.pvalue:.4f}")
g.asserted('我的"左尾"机制被自己的检验否定了',abs(s10.statistic)<0.2,
           f'10 分位 rho={s10.statistic:+.3f}(p={s10.pvalue:.3f})—— 左尾是平的,不是它')
g.asserted('真正的机制:罕见类别在**中位数**上来得更晚,而**均值**对此完全盲',
           abs(s50.statistic)>0.3 and s50.pvalue<0.05 and abs(smn.statistic)<0.15,
           f'中位数 rho={s50.statistic:+.3f}(p={s50.pvalue:.4f}) vs 均值 rho={smn.statistic:+.3f}'
           f'(p={smn.pvalue:.3f})。#128 里"人群层几乎为零"是在**均值**上算的,汇总选错了')
g.require_resolvable_first('扣掉时间表后 Delta 可分辨',abs(d_sr-n_sr),db2,family='deschedule')
g.offset_control('扣掉时间表后 Delta 相对人内置换零',d_sr,n_sr,db2,
                 null_kind='人内置换起始年龄标签(在同一残差矩阵上)')
print(g)
D.to_csv(OUT/'grid.csv',index=False)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A14 R04 -- 66.852% 的时间表是用**均值**排的。换成中位数会怎样?

【CLOSURE,明确标注】本轮不开新战线,它保护本项目最老也最大的一个数字。

#130e 是本弧最便宜也最可迁移的东西:同一批 31 个题目,
    Spearman(稀有度, 起始年龄**均值**)   = +0.011  (p = 0.954)
    Spearman(稀有度, 起始年龄**中位数**) = +0.437  (p = 0.014)
同一组数据、同一个关系,换一个汇总统计量就从"完全没有"变成"显著"。

而 `#75` 的时间表(用它猜任意两个兴趣谁先谁后达 66.852 ± 0.191)与 `#63` 的
全人群上界 [60.5%, 66.5%],**都是从 `np.nanmean(V[tr],axis=0)` 排出来的**
(A03/R22 第 66 行)。所以这个问题不是学术的:如果中位数排出的顺序不同,那么
"人群共享一个发育顺序"这句话的**内容**就依赖于我选了哪个汇总量。

ESTIMAND        held-out 成对顺序准确率,作为**排序所用汇总统计量**的函数。
                以及:各汇总量排出的 31 个类别顺序之间的 Spearman。
IDENTIFICATION  训练/测试按人对半分,顺序只从训练半边学,准确率只在测试半边算。
SCOPE           报告 >=6 个类别起始年龄的人(与 A03/R22 同口径,便于直接比较)。
WORLDS          robust    各汇总量的准确率差 < 种子展布 -> #75 与汇总量的选择无关,
                          66.852% 是关于人的,不是关于我的代码的。
                fragile   有汇总量把准确率显著推高或推低 -> #75 的数字是"均值时间表"
                          的性质,标题必须加限定词,而更高的那个才是共享顺序的真实上界。
KILL            条件式:随机顺序必须在 50%,真顺序必须显著高于它,才读各汇总量之间的差。
POSITIVE CTRL   真顺序 vs 随机顺序。
NEGATIVE CTRL   把训练半边的顺序随机打乱。
NOISE FLOOR     8 seeds(与 A03/R22 同)。
MULTIPLICITY    6 个汇总量 x 8 seeds,整格发表。
IMPOSSIBLE      2 年分箱下大量并列,任何汇总量都受它限制;本轮只判**相对**差异。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_coverage

OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values; mask=~np.isnan(V)
prev=mask.mean(0); keep=np.flatnonzero(mask.sum(1)>=6)
print(f"eligible {len(keep):,}  categories {V.shape[1]}",flush=True)

SUM={'mean':      lambda A: np.nanmean(A,axis=0),
     'median':    lambda A: np.nanmedian(A,axis=0),
     'q25':       lambda A: np.nanpercentile(A,25,axis=0),
     'q40':       lambda A: np.nanpercentile(A,40,axis=0),
     'q60':       lambda A: np.nanpercentile(A,60,axis=0),
     'q75':       lambda A: np.nanpercentile(A,75,axis=0),
     'trimmed20': lambda A: np.array([np.nanmean(np.clip(A[:,j],
                       np.nanpercentile(A[:,j],10),np.nanpercentile(A[:,j],90)))
                       for j in range(A.shape[1])])}

def acc(order_vals,people,rng,mode='skip_ties'):
    """⚠ mode='skip_ties' 是 A03/R22 的原做法,它**跳过**排序量并列的对。
    分位数汇总会产生大量并列(中位数下 31 个类别里有 580 个有序对并列),于是那一臂
    是在一个**更容易的子集**上打分 —— 不同的输出上比准确率是无效的(#101b same_scale)。
    mode='half_credit' 在**同一批对**上打分,并列记 0.5,这才是可比的。"""
    right=0.;tot=0
    for i in people:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        for _ in range(min(10,len(j))):
            a,b=rng.choice(j,2,replace=False)
            if V[i,a]==V[i,b]: continue
            tie=order_vals[a]==order_vals[b]
            if tie and mode=='skip_ties': continue
            right+= 0.5 if tie else float((order_vals[a]<order_vals[b])==(V[i,a]<V[i,b]))
            tot+=1
    return 100*right/max(tot,1),tot

rows=[]; ORD={}
for seed in range(1,9):
    rng=np.random.default_rng(seed)
    p=rng.permutation(keep); tr,te=p[:len(p)//2],p[len(p)//2:]
    for nm,f in SUM.items():
        o=f(V[tr])
        a,n=acc(o,te,np.random.default_rng(seed*13),'skip_ties')       # 原做法
        a2,n2=acc(o,te,np.random.default_rng(seed*13),'half_credit')   # 同一批对
        rows.append(dict(seed=seed,summary=nm,acc=a,pairs=n,acc_fair=a2,pairs_fair=n2,
                         ties=int((o[:,None]==o[None,:]).sum()-len(o))))
        if seed==1: ORD[nm]=o
    o=SUM['mean'](V[tr]); po=rng.permutation(o)
    ar,_=acc(po,te,np.random.default_rng(seed*13),'skip_ties')
    ar2,_=acc(po,te,np.random.default_rng(seed*13),'half_credit')
    rows.append(dict(seed=seed,summary='random',acc=ar,pairs=0,acc_fair=ar2,pairs_fair=0,ties=0))
    print(f"  seed {seed}",flush=True)

D=pd.DataFrame(rows); D.to_csv(OUT/'grid.csv',index=False)
G=D.groupby('summary').acc.agg(['mean','std']); F=D.groupby('summary').acc_fair.agg(['mean','std'])
TI=D.groupby('summary').ties.mean(); NP=D.groupby('summary').pairs.mean(); NF=D.groupby('summary').pairs_fair.mean()
print("\n=== held-out 成对顺序准确率(8 seeds)===")
print(f"  {'汇总量':<11} {'跳过并列':>9} {'打分对数':>9} | {'同一批对':>9} {'sd':>6} {'vs mean':>9} {'类别并列':>8}")
base=G.loc['mean','mean']; basef=F.loc['mean','mean']
for nm in list(SUM)+['random']:
    print(f"  {nm:<11} {G.loc[nm,'mean']:>9.3f} {NP.loc[nm]:>9.0f} | {F.loc[nm,'mean']:>9.3f} "
          f"{F.loc[nm,'std']:>6.3f} {F.loc[nm,'mean']-basef:>+9.3f} {TI.loc[nm]:>8.0f}")
print(f"\n  ⚠ 跳过并列时 median 只在 {NP.loc['median']:.0f} 个对上打分,mean 在 {NP.loc['mean']:.0f} 个 —— "
      f"差 {100*(1-NP.loc['median']/NP.loc['mean']):.0f}%,是不同的输出")

from scipy.stats import spearmanr
print("\n=== 各汇总量排出的 31 类别顺序之间的 Spearman(seed 1)===")
ks=list(SUM)
print("  " + " ".join(f"{k:>10}" for k in ks))
for a in ks:
    print(f"  {a:<10}" + " ".join(f"{spearmanr(ORD[a],ORD[b]).statistic:>10.3f}" for b in ks))

sdm=float(F.loc['mean','std']); sdr=float(F.loc['random','std'])
best=max(SUM,key=lambda k: F.loc[k,'mean']); worst=min(SUM,key=lambda k: F.loc[k,'mean'])
spread=float(F.loc[best,'mean']-F.loc[worst,'mean'])
g=Gate('66.852% 是关于人的,还是关于我选了哪个汇总量的')
g.same_scale('两臂在同一批对上打分',float(NF.loc['mean']),float(NF.loc['median']),'打分对数')
g.asserted('随机顺序在 50%(按它**自己的**展布判,不按选定的阈值)',
           abs(F.loc['random','mean']-50)<2*sdr,
           f"{F.loc['random','mean']:.2f}% ± {sdr:.2f} -> 距 50 有 {abs(F.loc['random','mean']-50)/sdr:.2f}x")
g.asserted('真顺序显著高于随机',F.loc['mean','mean']-F.loc['random','mean']>5,
           f"{F.loc['mean','mean']:.3f}% vs {F.loc['random','mean']:.2f}%")
g.require_resolvable_first('同一批对上,汇总量之间的差是否可分辨',spread,sdm)
g.offset_control('同一批对上,最好汇总量 vs 均值',float(F.loc[best,'mean']),basef,sdm,
                 null_kind='同一设计下的均值时间表(不是零假设,是现行做法的基准)')
g.asserted('跳过并列会制造多大的假增益',True,
           f"median 跳过并列 {G.loc['median','mean']:.3f}% vs 同一批对 {F.loc['median','mean']:.3f}% "
           f"-> 假增益 {G.loc['median','mean']-F.loc['median','mean']:+.3f} 个百分点")
print(g)
print(f"\n  同一批对上:最好 {best} {F.loc[best,'mean']:.3f}%   最差 {worst} {F.loc[worst,'mean']:.3f}%   "
      f"极差 {spread:.3f} 个百分点 = {spread/sdm:.1f}x 种子 sd")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

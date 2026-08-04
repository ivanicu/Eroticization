import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A61 R279 -- C(递归)的正面检验,第二次 —— 这次用一把够锋利的尺子

`#226` 第一次问 C 失败,两个具名原因:量「内容」的尺子装反了一半(`#227`),
分层设计的分辨率不够(展布 ±0.10 压过 0.15 的量),而且**卡钳匹配根本没成**(中位 49/62/69)。
现在三条都变了:`#228a` 给了一个**零 = +0.0043、效应 = +0.4290** 的量(信噪比 ~100×),
`#233` 给了计数的正确读法,而这一轮把匹配做成**分箱抽样到同一分布**并在门里断言。

⚠ **换了因变量,而这是有理由的,不是换赛道**:`#226` 分层用的是位置分 `S`。
`#232b` 之后,**羞耻**才是两条罕见路径唯一同向的共同落点 —— 若赋值真的回流重塑表征,
它应当在**赋值最强的地方**(羞耻最高的人)留下最多的结构。

WORLDS          ① **C 递归**:赋值回流重塑表征 -> **羞耻越高的人,宽度剖面越有结构**
                ② **A / B**:表征的结构由内容/坐标决定,与这个人赋了多少值无关 -> 三层相同
ESTIMAND        按羞耻分三层,各层各算一次**人×块残差剖面的两半复现**(`#228a` 同款)。
KILL            **若复现随羞耻单调上升且最高层明显高于最低层(差 > 2× 展布)-> C 得到第一个正面证据;
                若三层相同 -> C 在这份数据上没有支持。**
⚠ 最强混杂         **勾选数** —— 勾得多的人格子测得更准,复现天然更高;而羞耻与勾选数相关。
                控制:三层**按勾选数分箱抽样到同一分布**,并在门里**断言三层中位相等**
                (`#226a` 的卡钳没成,这一次必须成)。
POSITIVE CTRL   种入「高羞耻层剖面更强」的结构,**强度扫描**,
                增量对着**实测基线**定价(`#228d`,不是对着零)。
NEGATIVE CTRL   两半各自独立块内跨人置换(`#228c`:必须在**劈开之后**)。
IMPOSSIBLE      横断面。「赋值回流塑造了表征」与「表征本来就更有结构的人更容易羞耻」
                在这份数据里**不可分**。能判的只有这个交互在不在。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl))
NB=len(MB); cov=np.zeros(NN); K=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1; K[ppl]+=M.sum(1)
ok=cov>=8; K=np.where(ok,K,np.nan)
SHN=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
sh=pd.to_numeric(d[SHN],errors='coerce').values.astype(float)
base=np.flatnonzero(ok&np.isfinite(sh)&np.isfinite(K))
print(f"块 {NB};有羞耻分且块覆盖 ≥8 的人 n = {len(base):,};"
      f"corr(羞耻, 勾选数) = {np.corrcoef(sh[base],K[base])[0,1]:+.4f}")

def profiles(rows, seed, plant=None, perm=False):
    rg=np.random.default_rng(seed); m=np.zeros(NN,bool); m[rows]=True
    A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        Mm=M if plant is None else np.clip(M+plant[ppl,b][:,None],0,1)
        o=rg.permutation(Mm.shape[1]); k=Mm.shape[1]//2
        ha=Mm[:,o[:k]].mean(1); hb=Mm[:,o[k:2*k]].mean(1)
        if perm: ha=ha[rg.permutation(len(ha))]; hb=hb[rg.permutation(len(hb))]
        A[b,ppl]=ha; B[b,ppl]=hb
    def prof(X):
        Fm=np.isfinite(X)&m[None,:]; Z=np.where(Fm,X,0.0); tot=Z.sum(0); ct=Fm.sum(0)
        Rr=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-Fm[b]>=6,(tot-Z[b])/np.maximum(ct-Fm[b],1),np.nan)
            Rr[b]=np.where(Fm[b],X[b]-lo,np.nan); Rr[b]=Rr[b]-np.nanmean(Rr[b])
        return Rr
    Ra,Rb=prof(A),prof(B); g=np.isfinite(Ra)&np.isfinite(Rb)&m[None,:]
    if g.sum()<3000: return np.nan,int(g.sum())
    return float(np.corrcoef(Ra[g],Rb[g])[0,1]), int(g.sum())

q1,q2=np.percentile(sh[base],[33.3,66.7])
TT={'低羞耻':base[sh[base]<=q1],'中羞耻':base[(sh[base]>q1)&(sh[base]<=q2)],'高羞耻':base[sh[base]>q2]}
print(f"三层 n = {[len(v) for v in TT.values()]};勾选数中位 = {[int(np.median(K[v])) for v in TT.values()]}")

# ⚠ 匹配:按勾选数十分位分箱,抽到三层共同的最小箱容量
rng=np.random.default_rng(20260804)
edges=np.percentile(K[base],np.arange(0,101,10))
MT={k:[] for k in TT}
for lo,hi in zip(edges[:-1],edges[1:]):
    pools={k:v[(K[v]>=lo)&(K[v]<hi if hi<edges[-1] else K[v]<=hi)] for k,v in TT.items()}
    n=min(len(p) for p in pools.values())
    for k in TT:
        if n>0: MT[k].append(rng.choice(pools[k],n,replace=False))
MT={k:np.concatenate(v) for k,v in MT.items()}
med=[int(np.median(K[v])) for v in MT.values()]
print(f"匹配后 n = {[len(v) for v in MT.values()]};勾选数中位 = **{med}**"
      f"(均值 {[round(float(np.mean(K[v])),1) for v in MT.values()]})")

def curve(T_, seeds=6, plant=None, perm=False):
    out={}
    for nm,rows in T_.items():
        vs=[profiles(rows,1100+s,plant,perm)[0] for s in range(seeds)]
        out[nm]=(float(np.nanmean(vs)),float(np.nanstd(vs)))
    return out
raw=curve(MT); nul=curve(MT,seeds=3,perm=True)
print(f"\n{'层':<8}{'剖面两半复现':>14}{'展布':>9}{'置换零':>10}")
rows=[]
for nm in MT:
    print(f"{nm:<8}{raw[nm][0]:>14.4f}{raw[nm][1]:>9.4f}{nul[nm][0]:>10.4f}")
    rows.append(dict(tertile=nm,n=len(MT[nm]),median_picks=int(np.median(K[MT[nm]])),
                     rep=raw[nm][0],sd=raw[nm][1],null=nul[nm][0]))
lo_,hi_=rows[0],rows[2]; gsd=float(np.hypot(lo_['sd'],hi_['sd']))
mono=rows[0]['rep']<=rows[1]['rep']<=rows[2]['rep']
print(f"高−低 = **{hi_['rep']-lo_['rep']:+.4f}** vs 2×展布 {2*gsd:.4f};单调 = {mono}")

# 正对照:强度扫描,增量对着实测基线定价
SW=[]
for gp in (0.02,0.05,0.10,0.20):
    P=np.zeros((NN,NB)); w=np.random.default_rng(55).standard_normal((NN,NB))
    P[MT['高羞耻']]=gp*w[MT['高羞耻']]; P[MT['中羞耻']]=0.4*gp*w[MT['中羞耻']]
    c=curve(MT,seeds=2,plant=P)
    SW.append((gp,c['高羞耻'][0]-c['低羞耻'][0]))
print(f"正对照强度扫描 g -> 高−低:"+' · '.join(f"{a:.2f}->{b:+.4f}" for a,b in SW)
      +f"  [基线 {hi_['rep']-lo_['rep']:+.4f} · 基线复现 {lo_['rep']:.4f}]")

T=pd.DataFrame(rows); check_columns(T,'R279')
T.to_csv(pathlib.Path(__file__).parent/'results'/'shame_tertiles.csv',index=False)

g=Gate('C 的正面检验(第二次)')
g.asserted('⚠ 最强混杂真的被控制住了:三层勾选数中位相等(`#226a` 那次没成)',
           len(set(med))==1, f"匹配前 {[int(np.median(K[v])) for v in TT.values()]} -> 匹配后 {med}")
g.asserted('正对照:强度扫描必须单调上升,且在某个 g 上超过实测基线的 2×展布(#228d)',
           SW[-1][1]>SW[0][1] and max(x[1] for x in SW)>2*gsd,
           ' · '.join(f"g={a:.2f} {b:+.4f}" for a,b in SW)+f";2×展布 {2*gsd:.4f}")
g.negative_control('高羞耻层的置换零',abs(float(hi_['null'])),abs(float(hi_['rep'])),
                   null_spread=None,
                   null_kind='两半各自独立块内跨人置换 —— 劈开之后打乱(#228c)')
g.offset_control('★ 高羞耻层 vs 低羞耻层',float(hi_['rep']),float(lo_['rep']),gsd,
                 null_kind='同一条管道在低羞耻层上的复现 —— 不是零假设,是「若没有交互,高层该落在哪」')
g.asserted('★ 注册的 kill:复现随羞耻单调上升且高层明显高于低层 -> C 得到第一个正面证据',
           mono and (hi_['rep']-lo_['rep']>2*gsd),
           f"三层 {[round(r['rep'],4) for r in rows]};单调={mono};"
           f"高−低 {hi_['rep']-lo_['rep']:+.4f} vs 2×展布 {2*gsd:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

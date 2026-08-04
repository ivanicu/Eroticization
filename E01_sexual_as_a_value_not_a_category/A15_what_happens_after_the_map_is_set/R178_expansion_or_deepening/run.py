import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A15 R01 -- 版图 17 岁后不再扩张,那成年后加进来的东西是什么样的?

#132d:29-32 岁的人,68.4% 的性兴趣是 17 岁前获得的,最晚的那个平均在 22.6 岁。
#132f:这给 Ivan 的模型 C(价值回流并重塑表征)划了一条时间边界。

那条边界之后仍然有 ~32% 的兴趣在进来。**它们和早期那些不一样吗?** 两个图景在心理学上
完全不同,而且它们对模型 C 的含义相反:

  DEEPEN   晚来的东西**靠近已有的中心** —— 成年后发生的是把既有的版图挖深,
           新条目是旧条目的邻居。C 被定域为"重排权重",不新增表征。
  EXPAND   晚来的东西**远离已有的中心** —— 成年后仍在开新的地。
           C 仍可以是一个持续的表征重塑过程。

第三种可能必须同时排除:
  ARTEFACT 晚来的只是**评分更低**的那些 —— 而 #114 已证明人把最爱的记得更早
           (-0.2000 年/评分 sd)。这条通路能**单独**造出"晚来的评分低",
           所以评分这条证据线天生是脏的,而**坐标位置那条不是**:
           #114 说的是"多喜欢",不是"喜欢哪个"。

ESTIMAND        对每个人:晚获得类别与他**早期集合**的平均共现相似度,减去把早/晚标签
                在他自己的类别集内置换后的期望。以及同样构造下的评分差与稀有度差。
IDENTIFICATION  零 = **人内置换早/晚标签**,精确保留这个人的类别集与晚的个数,
                只摧毁"哪些是晚的"。所以对"他喜欢多少东西"、"他整体早熟不早熟"、
                "他答了几个块"全部免疫。
SCOPE           29-32 岁档、报告 >=8 个类别起始年龄、且早/晚两侧都非空的人。
                分界 17.5 岁(= release 的 '17-18yo' 分箱下沿)。
WORLDS          DEEPEN / EXPAND / ARTEFACT(见上)
KILL            条件式:正对照必须开火**且**人内置换零必须为零,才读阈值。
POSITIVE CTRL   种植:强制把与早期集合最相似的 k 个类别标为"晚",相似度差必须随种植
                比例单调上升。g=0 必须逐位复现真实臂。
NEGATIVE CTRL   人内置换早/晚标签,5 个种子。
CONFOUND        #114 回忆偏差。按 #129i 的教训,喂给 artifact_cannot_explain 的必须是
                **实际贡献**(剥离该通道前后的差),不是单位幅度的上界。
                第二个:晚获得的类别在人群里可能本来就更罕见 -> 稀有度作为第三条线报出。
                第三个:类别数与晚的个数 —— 人内置换零已经把它们固定住了。
NOISE FLOOR     200 次按人自助 + 5 个置换种子。
MULTIPLICITY    3 条证据线 x 4 个种植水平 x 5 seeds x {含/不含评分校正},整格发表。
IMPOSSIBLE      因果与队列。29-32 岁档是一个横断面切片;"成年后加进来的"是**自报的**
                获得时间,而 #114/#119 已证明这个时间被系统性扭曲。本轮只判
                **早/晚两组在内容上是否不同**,不判它们何时真的发生。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

CUT=17.5
OLD=(age==30.5)                      # 29-32 档
E=(obs&(V<=CUT)); L=(obs&(V>CUT))
who=np.flatnonzero(KEEP&OLD&(E.sum(1)>=3)&(L.sum(1)>=1))
print(f"29-32 档、早/晚两侧都非空、>=8 个类别:{len(who):,} 人  "
      f"早 {E[who].sum(1).mean():.1f} 晚 {L[who].sum(1).mean():.1f}",flush=True)
base=np.flatnonzero(KEEP&OLD)
no_late=int((L[base].sum(1)==0).sum()); few_early=int(((L[base].sum(1)>0)&(E[base].sum(1)<3)).sum())
print(f"⚠ 显式记录被跳过的 {len(base)-len(who):,}/{len(base):,}(#118c):"
      f"**17 岁后一个新兴趣都没有 {no_late:,} 人 = {no_late/len(base):.1%}**,"
      f"早期不足 3 个 {few_early:,} 人 = {few_early/len(base):.1%}",flush=True)
check_coverage(len(who),int((KEEP&OLD).sum()),"早/晚都非空(35% 无晚期条目,已显式记录)",tol=0.40)

# ---- 类别-类别相似度:全人群共现(phi 系数),对角置零
Ob=obs.astype(float)
p=Ob.mean(0); C=(Ob.T@Ob)/len(Ob)
den=np.sqrt(np.outer(p*(1-p),p*(1-p))); den[den<1e-9]=1e-9
SIM=(C-np.outer(p,p))/den; np.fill_diagonal(SIM,0.)
print(f"共现相似度矩阵 {SIM.shape}  范围 {SIM.min():+.3f}..{SIM.max():+.3f}",flush=True)

# ⚠ 相似度可能只是稀有度的影子:罕见类别与**任何**类别的 phi 上界都低。
#   在配对层把稀有度回归掉,得到一个"去稀有度的连通性"。
iu=np.triu_indices(len(rar),1)
X=np.c_[np.ones(len(iu[0])),rar[iu[0]]+rar[iu[1]],rar[iu[0]]*rar[iu[1]],
        np.abs(rar[iu[0]]-rar[iu[1]])]
yv=SIM[iu]; bcoef=np.linalg.lstsq(X,yv,rcond=None)[0]
res=yv-X@bcoef
check_residualized(res,rar[iu[0]]+rar[iu[1]],"配对相似度对稀有度")
SIMR=np.zeros_like(SIM); SIMR[iu]=res; SIMR=SIMR+SIMR.T
print(f"稀有度解释了配对相似度的 R2 = {1-res.var()/yv.var():.3f};去稀有度后范围 "
      f"{SIMR[iu].min():+.3f}..{SIMR[iu].max():+.3f}",flush=True)

RM=np.full_like(V,np.nan)
for j,ri in best.items(): RM[:,j]=R[:,ri]

def lines(latemask,i):
    """给定'哪些是晚的',返回三条证据线。"""
    m=obs[i]; jj=np.flatnonzero(m)
    lt=jj[latemask]; er=jj[~latemask]
    if len(lt)==0 or len(er)==0: return None
    sim=SIM[np.ix_(lt,er)].mean()                       # 晚 -> 早期集合的平均相似度
    simr=SIMR[np.ix_(lt,er)].mean()                     # 同上,但配对稀有度已回归掉
    rat=np.nanmean(RM[i,lt])-np.nanmean(RM[i,er])       # 评分差
    rr =rar[lt].mean()-rar[er].mean()                   # 稀有度差
    return sim,rat,rr,simr

def run(assign,tag,seedbase=0):
    S3=[]
    for i in who:
        m=obs[i]; k=int(L[i,m].sum())
        lm=assign(i,m,k)
        r=lines(lm,i)
        if r is not None: S3.append(r)
    A=np.array(S3,dtype=float)
    return dict(tag=tag,n=len(A),sim=float(np.nanmean(A[:,0])),
                rating=float(np.nanmean(A[:,1])),rarity=float(np.nanmean(A[:,2])),
                simr=float(np.nanmean(A[:,3]))),A

real_row,A_real=run(lambda i,m,k: L[i,m],'real')
rows=[real_row]
for s_ in range(5):
    rg=np.random.default_rng(9200+s_)
    r,_=run(lambda i,m,k,rg=rg: np.isin(np.arange(m.sum()),rg.choice(m.sum(),k,replace=False)),'perm')
    rows.append({**r,'seed':s_})
# 正对照:把与早期集合最相似的 k 个标为晚
for g_ in [0.0,0.34,0.67,1.0]:
    def assign(i,m,k,g_=g_):
        base=L[i,m]
        if g_==0.: return base
        jj=np.flatnonzero(m); er=jj[~base]
        if len(er)==0: return base
        sc=SIM[np.ix_(jj,er)].mean(1)
        top=np.argsort(-sc)[:k]
        out=np.zeros(m.sum(),bool); out[top]=True
        n_sw=int(round(g_*k))
        mix=base.copy()
        if n_sw>0:
            on=np.flatnonzero(base); off=np.flatnonzero(out&~base)
            n_sw=min(n_sw,len(on),len(off))
            mix[on[:n_sw]]=False; mix[off[:n_sw]]=True
        return mix
    r,_=run(assign,f'plant{g_}'); rows.append(r)
    print(f"  plant {g_}",flush=True)

# #114 校正臂:把评分通道从起始年龄里剥掉,重新定早/晚
zr=(RM-np.nanmean(RM))/np.nanstd(RM)
f2=obs&np.isfinite(zr)
Vc=np.where(f2,V-(-0.2000)*zr,np.where(obs,V,np.nan))
Lc=(obs&(Vc>CUT))
r,_=run(lambda i,m,k: Lc[i,m],'rating_corrected'); rows.append(r)

D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
G=D.groupby('tag')[['sim','simr','rating','rarity']].mean()
print("\n=== 三条证据线(晚 − 早)===")
print(G.round(4).to_string())

rb=np.random.default_rng(4444); n=len(A_real); idx=np.arange(n)
B=np.array([np.nanmean(A_real[rb.integers(0,n,n)],axis=0) for _ in range(200)])
sb,rb_,rrb,srb=B.std(0)
real=G.loc['real']; perm=G.loc['perm']
print(f"\n按人自助(200):sim {sb:.4f}  simr {srb:.4f}  rating {rb_:.4f}  rarity {rrb:.4f}")

g=Gate('成年后加进来的是深化还是扩张')
g.degenerate_matches_reference('g=0 逐位复现 real',float(G.loc['plant0.0','sim']),float(real.sim))
mono=[float(G.loc[f'plant{q}','sim']) for q in [0.0,0.34,0.67,1.0]]
g.asserted('种植的相似度随比例单调上升',all(mono[i]<mono[i+1] for i in range(len(mono)-1)),
           " < ".join(f"{v:+.4f}" for v in mono))
g.require_resolvable_first('相似度差可分辨',abs(real.sim-perm.sim),sb,family='sim')
g.offset_control('相似度:晚的更靠近早期集合?',float(real.sim),float(perm.sim),sb,
                 null_kind='人内置换早/晚标签(保留类别集与晚的个数)')
g.artifact_cannot_explain('#114 的实际贡献不能解释相似度',
                          float(real.sim-G.loc['rating_corrected','sim']),float(real.sim-perm.sim),sb)
g.require_resolvable_first('去稀有度后的连通性差可分辨',abs(real.simr-perm.simr),srb,family='simr')
g.offset_control('去稀有度后:晚的仍然更不连通?',float(real.simr),float(perm.simr),srb,
                 null_kind='人内置换早/晚标签(保留类别集与晚的个数)')
g.artifact_cannot_explain('#114 的实际贡献不能解释去稀有度连通性',
                          float(real.simr-G.loc['rating_corrected','simr']),
                          float(real.simr-perm.simr),srb)
g.same_scale('去稀有度前后同一批人',float(D[D.tag=='real'].n.iloc[0]),float(D[D.tag=='real'].n.iloc[0]),'人数')
g.require_resolvable_first('评分差可分辨',abs(real.rating-perm.rating),rb_,family='rating')
g.offset_control('评分:晚的评分更低?',float(real.rating),float(perm.rating),rb_,
                 null_kind='同上,人内置换早/晚标签')
g.artifact_cannot_explain('#114 的实际贡献不能解释评分差',
                          float(real.rating-G.loc['rating_corrected','rating']),
                          float(real.rating-perm.rating),rb_)
g.require_resolvable_first('稀有度差可分辨',abs(real.rarity-perm.rarity),rrb,family='rarity')
g.offset_control('稀有度:晚的更罕见?',float(real.rarity),float(perm.rarity),rrb,
                 null_kind='同上,人内置换早/晚标签')
print(g)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

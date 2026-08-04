import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A57 R273 -- breadth 是回答风格,还是一个真的人格维度

`#227b`:第三个人层维度的正名是 **breadth**(你对多少东西说"是")。
`#193c` 撤掉过"默认赞同因子",而 breadth 现在从内容维度那条路走了回来。
`#26` 说这份数据分不开"默认赞同"与"真的被更多东西唤起" —— **但那是对单一 Likert 量表说的,
对 32 个块的结构从没被问过。**

WORLDS          ① **回答风格**:有些人对任何问卷都更爱勾"是" -> breadth 在**所有块上等强**
                ② **真的更宽**:有些人确实被更多东西唤起 -> **块之间有结构**,
                   且该结构在**互不相交的人群半样本**上复现
ESTIMAND        breadth 的**块层载荷谱** —— 每块的人内勾选比例 p_ib 与该人的
                **留一** breadth(扣掉块 b 自己)之间的跨人相关,共 32 个数。
KILL            **若扣掉块规模之后,载荷谱的跨人分半复现与其零不可分 -> 回答风格,
                `#188`/`#189` 的第三维要按"默认赞同"重写;
                若残差谱有结构且复现 -> 真的更宽,它是一个真的维度,只是名字一直是错的。**
⚠ 最强混杂(跑之前写下)
                **块大小/平均勾选率本身就会造出载荷谱** —— 选项多的块 p_ib 测得更准,
                与 breadth 的相关天然更高。**而它在两半里会完美复现**,
                所以"分半复现"单独**不能**分开两个世界。
控制(同一迭代内)
                把 32 个载荷对 `(选项数, 平均勾选率)` 回归掉,**只判残差谱的复现**。
⚠ 留一           留一 breadth 是必须的:若 breadth 含块 b 自己,p_ib 与它必然正相关(#227a 同款)。
NEGATIVE CTRL   块内跨人置换 p_ib(打掉人层协方差),同一条管道跑复现。
POSITIVE CTRL   种入一类人:**只在一半的块上更宽**。残差谱必须把那一半认出来。
IMPOSSIBLE      这份数据没有非情色的对照问卷,所以"对任何问卷都爱勾是"这个世界
                **只能被块间结构否证,不能被直接观测**。谱有结构 -> 世界①不足以解释;
                谱平坦 -> **不等于**世界①成立,只等于世界②在这个分辨率下没有证据。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
P=[]; NOPT=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    v=np.full(NN,np.nan); v[ppl]=M.mean(1); P.append(v); NOPT.append(len(opt))
P=np.array(P); NOPT=np.array(NOPT,float); NB=len(P)
RATE=np.nanmean(P,axis=1)
cov=np.isfinite(P).sum(0); ok=cov>=8
print(f"块 {NB} 个;n = {int(ok.sum()):,};选项数 {NOPT.min():.0f}–{NOPT.max():.0f};"
      f"平均勾选率 {RATE.min():.3f}–{RATE.max():.3f}")

def loadings(Pm, rows):
    """每块与【留一】breadth 的跨人相关。"""
    m=np.zeros(NN,bool); m[rows]=True
    F=np.isfinite(Pm)&m[None,:]; Z=np.where(F,Pm,0.0)
    tot=Z.sum(0); ct=F.sum(0)
    out=np.full(NB,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)   # ⚠ 留一
        g=F[b]&np.isfinite(lo)
        if g.sum()>200: out[b]=float(np.corrcoef(Pm[b][g],lo[g])[0,1])
    return out
def resid_size(L):
    g=np.isfinite(L); X=np.column_stack([np.ones(g.sum()),NOPT[g],RATE[g]])
    beta,*_=np.linalg.lstsq(X,L[g],rcond=None)
    o=np.full(NB,np.nan); o[g]=L[g]-X@beta; return o
def replication(Pm, seed, use_resid=True):
    rg=np.random.default_rng(seed); rows=np.flatnonzero(ok); p=rg.permutation(rows); h=len(p)//2
    A=loadings(Pm,p[:h]); B=loadings(Pm,p[h:])
    if use_resid: A,B=resid_size(A),resid_size(B)
    g=np.isfinite(A)&np.isfinite(B)
    return float(np.corrcoef(A[g],B[g])[0,1]) if g.sum()>=10 else np.nan
def permuted(seed):
    rg=np.random.default_rng(seed); Q=P.copy()
    for b in range(NB):
        idx=np.flatnonzero(np.isfinite(Q[b])); Q[b][idx]=rg.permutation(Q[b][idx])
    return Q

L=loadings(P,np.flatnonzero(ok)); Lr=resid_size(L)
print(f"\n载荷谱(全样本):均值 {np.nanmean(L):+.4f} · 展布 {np.nanstd(L):.4f} · "
      f"范围 {np.nanmin(L):+.4f}..{np.nanmax(L):+.4f}")
print(f"  对 (选项数, 平均勾选率) 回归后:残差展布 {np.nanstd(Lr):.4f} "
      f"(**块规模解释了 {100*(1-np.nanvar(Lr)/np.nanvar(L)):.1f}%** 的谱方差)")
rep =[replication(P,300+s)        for s in range(8)]
repn=[replication(permuted(400+s),300+s) for s in range(8)]
rep0=[replication(P,300+s,use_resid=False) for s in range(8)]
print(f"\n残差谱的跨人分半复现:**{np.mean(rep):+.4f} ± {np.std(rep):.4f}** · "
      f"置换零 {np.mean(repn):+.4f} ± {np.std(repn):.4f}")
print(f"  (⚠ 未扣块规模时是 {np.mean(rep0):+.4f} —— 这正是混杂:它复现得几乎完美,却什么也没说)")

# 正对照:一类人只在一半块上更宽
rg=np.random.default_rng(20260804); half=rg.permutation(NB)[:NB//2]
w=rg.standard_normal(NN)*0.05; Pp=P.copy()
for b in half: Pp[b]=np.clip(Pp[b]+w,0,1)
Lp=resid_size(loadings(Pp,np.flatnonzero(ok)))
sep=float(np.nanmean(Lp[half])-np.nanmean(Lp[[b for b in range(NB) if b not in set(half)]]))
repp=[replication(Pp,300+s) for s in range(4)]
print(f"\n正对照(一类人只在 {len(half)}/{NB} 个块上更宽):"
      f"残差谱把那一半分开 **{sep:+.4f}**;复现 {np.mean(repp):+.4f}")

# ---------- 正对照的强度扫描(⚠ `#211a`:硬编码的 g 又一次点不着火) ----------
SW=[]
for gp in (0.05,0.15,0.30,0.60):
    w2=np.random.default_rng(9).standard_normal(NN)*gp; Q=P.copy()
    for b in half: Q[b]=np.clip(Q[b]+w2,0,1)
    Lq=resid_size(loadings(Q,np.flatnonzero(ok)))
    SW.append((gp,float(np.nanmean(Lq[half])-np.nanmean(Lq[[b for b in range(NB) if b not in set(half)]]))))
print("  强度扫描 g -> 分开度:"+' · '.join(f"{a:.2f}->{b:+.4f}" for a,b in SW)
      +f"  (可观测量:残差展布 {np.nanstd(Lr):.4f})")

# ---------- 第二段:块内选项劈半,判【人×块】残差剖面 ----------
# ⚠ 上面那一段测的载荷是【块的属性】—— 任何与人无关的块性质都会在两半里完美复现,
#   所以它分不开两个世界。这一段把测量放到【同一个块的两半选项】上:
#   块的属性对两半相同、被块中心化扣掉,剩下的只有【这个人在这个块上是不是格外宽】。
MB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl))
def profiles(seed, plant=0.0, perm=False):
    rg=np.random.default_rng(seed)
    pl=rg.standard_normal((NN,NB))*plant if plant else None
    H=[np.full((2,NN),np.nan) for _ in range(NB)]
    for b,(M,ppl) in enumerate(MB):
        Mm=M.copy()
        if plant: Mm=np.clip(Mm+pl[ppl,b][:,None],0,1)
        o=rg.permutation(Mm.shape[1]); k=Mm.shape[1]//2
        ha=Mm[:,o[:k]].mean(1); hb=Mm[:,o[k:2*k]].mean(1)
        # ⚠ 零必须在【劈开之后】对两半各自独立打乱 —— 劈开之前打乱,两半拿到同一次置换,
        #   人×块这个格子原封不动,零会比效应还高(实测 +0.5835 vs +0.4290)。
        if perm: ha=ha[rg.permutation(len(ha))]; hb=hb[rg.permutation(len(hb))]
        H[b][0,ppl]=ha; H[b][1,ppl]=hb
    A=np.array([h[0] for h in H]); B=np.array([h[1] for h in H])
    def prof(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0)
        R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)   # ⚠ 留一
            R[b]=X[b]-lo
            R[b]=R[b]-np.nanmean(R[b])                                        # 块中心化
        return R
    Ra,Rb=prof(A),prof(B); m=np.isfinite(Ra)&np.isfinite(Rb)&ok[None,:]
    return float(np.corrcoef(Ra[m],Rb[m])[0,1]), int(m.sum())
pr=[profiles(600+s)[0] for s in range(5)]; pn=[profiles(700+s,perm=True)[0] for s in range(5)]
_,ncell=profiles(600)
print(f"\n第二段 人×块 残差剖面(块内选项劈半,n_cell = {ncell:,}):"
      f"**{np.mean(pr):+.4f} ± {np.std(pr):.4f}** · 置换零 {np.mean(pn):+.4f} ± {np.std(pn):.4f}")
PS=[(gp,np.mean([profiles(800+s,plant=gp)[0] for s in range(2)])) for gp in (0.02,0.05,0.10)]
PSW=[(gp,profiles(820,plant=gp)[0]) for gp in (0.20,0.40)]
print(f"  正对照(种入人×块特异的宽度)g -> 剖面复现:"
      +' · '.join(f"{a:.2f}->{b:+.4f}" for a,b in PS)
      +" ‖ 更宽的扫描(报告用,不入注册条款):"+' · '.join(f"{a:.2f}->{b:+.4f}" for a,b in PSW))

T=pd.DataFrame(dict(block=np.arange(NB),n_options=NOPT,pick_rate=RATE,loading=L,loading_resid=Lr))
check_columns(T,'R273'); T.to_csv(pathlib.Path(__file__).parent/'results'/'loading_spectrum.csv',index=False)

g=Gate('breadth:回答风格还是真的维度')
g.asserted('⚠ 第一段的设计不合身:载荷是【块】的属性,与人无关的块性质也会完美复现',
           False, f"未扣块规模 {np.mean(rep0):+.4f} 与扣掉后 {np.mean(rep):+.4f} 几乎相同 —— 这把尺子量的不是人")
g.asserted('正对照(第二段):种入人×块特异的宽度 -> 剖面复现必须随强度上升',
           PS[-1][1]>PS[0][1]+0.05, ' · '.join(f"g={a:.2f} {b:+.4f}" for a,b in PS))
g.negative_control('人×块剖面复现的置换零',abs(float(np.mean(pn))),abs(float(np.mean(pr))),
                   null_spread=float(np.std(pn)),
                   null_kind='块内跨人置换整行 —— 打掉人层结构,保留每块的选项边际与块内相关')
g.asserted('★ 真正的 kill:人×块剖面复现显著高于零 -> breadth 是【领域特异】的,不是单一回答风格',
           np.mean(pr)>np.mean(pn)+2*np.std(pn) and np.mean(pr)>0.05,
           f"剖面复现 {np.mean(pr):+.4f} vs 零 {np.mean(pn):+.4f} ± {np.std(pn):.4f}")
g.asserted('正对照:一类人只在一半块上更宽 -> 残差谱必须把那一半认出来',
           sep>0.02 and np.mean(repp)>0.3, f"分开 {sep:+.4f};复现 {np.mean(repp):+.4f}")
g.asserted('⚠ 混杂已暴露:未扣块规模时复现几乎完美,所以它单独不能分开两个世界',
           np.mean(rep0)>0.8, f"未扣 {np.mean(rep0):+.4f} vs 扣掉后 {np.mean(rep):+.4f}")
g.negative_control('残差谱复现的置换零',abs(float(np.mean(repn))),abs(float(np.mean(rep))),
                   null_spread=float(np.std(repn)),
                   null_kind='块内跨人置换 p_ib —— 打掉人层协方差,保留每块的边际分布')
g.asserted('注册的 kill:残差谱复现与零不可分 -> 回答风格;有结构且复现 -> 真的维度',
           np.mean(rep)>np.mean(repn)+2*np.std(repn) and np.mean(rep)>0.2,
           f"复现 {np.mean(rep):+.4f} vs 零 {np.mean(repn):+.4f} ± {np.std(repn):.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

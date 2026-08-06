import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A62 R280 -- 第四件事换一把仪器还在不在

`#229`–`#232` 建立的「第四件事」(一个人在**哪些领域**敞开、以**什么形状**敞开)
**每一个数字都来自同一把仪器** —— 32 道多选题的块。
`P14` 的 `instrument` 那一行说得很清楚:若一条主张要经过一个仪器,它就是关于**那个仪器**的主张。
这份 release 里有第二把、且**刻意与块选项不相交**的仪器:**AGE_ONSET 类别**
(`R173` 用 `check_disjoint_items` 断言过不相交)。报告了某类别的起始年龄 = 拥有该兴趣,
所以 `isfinite(V)` 本身就是一个二值背书矩阵。

ESTIMAND        同一个构念在两把仪器上各测一次:
                `D_块`   = 越轨半 − 普通半的宽度对比(按**块**的平均勾选率切,`#231`)
                `D_起始` = 罕见类别 − 常见类别的背书率对比(按**类别**的流行度切)
                判 `corr(D_块, D_起始)`,并与**衰减上限** `sqrt(rel_块 × rel_起始)` 比。
KILL            **若相关明显低于衰减上限(< 一半)-> 第四件事是块这把仪器的性质,
                `#229`–`#232` 全部要按「关于块格式的主张」重写;
                若接近上限 -> 它是关于人的,而且这是这个项目第一份跨仪器证据。**
⚠ 零应该是零吗?     **不应该。** 两个都是有噪声的读数,即使构念相同也到不了 1。
                所以判据是 **offset = 衰减上限**,不是零。
POSITIVE CTRL   块**内部**的两个不相交半边必须达到它们**自己的**衰减上限 ——
                这验证的是「上限」这套算法本身,而不是结论。
NEGATIVE CTRL   打乱人的配对(把 `D_起始` 在人之间置换)。
第二问            **跨仪器复现那条主张,而不只是那个分数**:`D_起始` 预不预测羞耻?
                `#232b` 里 `D_块` 给 +0.0862。
IMPOSSIBLE      两把仪器的**内容**不同(块问具体行为/场景,类别问兴趣领域),
                所以低相关既可能是"仪器性质"也可能是"构念在两处内容上真的不同"。
                能判的是**这个构念是否跨仪器存在**,不是它们是否是同一个变量。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
V=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
OBS=np.isfinite(V); PREV=OBS.mean(0)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]; BCOL=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl)); BCOL.append(str(q.col))
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
okB=cov>=8; okO=OBS.sum(1)>=8
assert not (set(BCOL)&set(ons)), "两把仪器共享了列 —— 不相交前提被破坏"
print(f"块仪器 {NB} 块 / 起始仪器 {len(ons)} 类别;两者列名交集 = 0(不相交已断言)")
print(f"n:块 {int(okB.sum()):,} · 起始 {int(okO.sum()):,} · 两者都有 **{int((okB&okO).sum()):,}**")

RATE=np.array([M.mean() for M,_ in MB]); ob=np.argsort(-RATE); ORD,TRG=ob[:NB//2],ob[NB//2:]
oo=np.argsort(-PREV); COM,RAR=oo[:len(ons)//2],oo[len(ons)//2:]
def z(v):
    m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
def D_block(cols_sel=None):
    R=np.full((NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        sub=M if cols_sel is None else M[:,cols_sel(b,M.shape[1])]
        if sub.shape[1]<3: continue
        R[b,ppl]=sub.mean(1)
    F=np.isfinite(R); Z=np.where(F,R,0.0); tot=Z.sum(0); ct=F.sum(0); Rr=np.full_like(R,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        Rr[b]=R[b]-lo; Rr[b]=Rr[b]-np.nanmean(Rr[b])
    def hs(cols):
        sub=Rr[cols]; F2=np.isfinite(sub)
        return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
    return z(hs(TRG))-z(hs(ORD))
def D_onset(cats=None):
    C=OBS if cats is None else OBS[:,cats]
    pv=PREV if cats is None else PREV[cats]
    o2=np.argsort(-pv); c_,r_=o2[:len(pv)//2],o2[len(pv)//2:]
    if len(c_)<3 or len(r_)<3: return np.full(NN,np.nan)
    a=C[:,r_].mean(1); b=C[:,c_].mean(1)
    a=np.where(okO,a,np.nan); b=np.where(okO,b,np.nan)
    return z(a)-z(b)
rng=np.random.default_rng(20260804)
DB=D_block(); DO=D_onset()
def rel_block(seed):
    rg=np.random.default_rng(seed)
    A=D_block(lambda b,k: rg.permutation(k)[:k//2]); B=D_block(lambda b,k: rg.permutation(k)[k//2:2*(k//2)])
    m=np.isfinite(A)&np.isfinite(B)&okB; r=float(np.corrcoef(A[m],B[m])[0,1]); return 2*r/(1+r), A, B
def rel_onset(seed):
    rg=np.random.default_rng(seed); p=rg.permutation(len(ons))
    A,B=D_onset(p[:len(ons)//2]),D_onset(p[len(ons)//2:])
    m=np.isfinite(A)&np.isfinite(B)&okO; r=float(np.corrcoef(A[m],B[m])[0,1]); return 2*r/(1+r), A, B
RB=np.mean([rel_block(70+s)[0] for s in range(3)]); RO=np.mean([rel_onset(80+s)[0] for s in range(3)])
m=np.isfinite(DB)&np.isfinite(DO)&okB&okO
obs=float(np.corrcoef(DB[m],DO[m])[0,1]); ceil=float(np.sqrt(max(RB,0)*max(RO,0)))
print(f"\n分半信度:D_块 **{RB:+.4f}** · D_起始 **{RO:+.4f}** -> 衰减上限 sqrt(积) = **{ceil:.4f}**")
print(f"**跨仪器相关 corr(D_块, D_起始) = {obs:+.4f}**(n = {int(m.sum()):,})· 占上限的 **{100*obs/ceil:.1f}%**")
nl=[float(np.corrcoef(DB[m],rng.permutation(DO[m]))[0,1]) for _ in range(50)]
print(f"  置换人配对的零:{np.mean(nl):+.4f} ± {np.std(nl):.4f}")
_,A1,B1=rel_block(70); mm=np.isfinite(A1)&np.isfinite(B1)&okB
pin=float(np.corrcoef(A1[mm],B1[mm])[0,1]); pceil=RB/(2-RB) if RB<1 else np.nan
print(f"正对照(块内部两个不相交半边):实测 {pin:+.4f} · 它们自己的上限 {pceil:.4f} · "
      f"占 {100*pin/pceil:.1f}%")

SHN=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SHN],errors='coerce').values.astype(float)
def cr(x):
    mm=np.isfinite(x)&np.isfinite(y); return float(np.corrcoef(x[mm],y[mm])[0,1]), int(mm.sum())
rb,nb_=cr(DB); ro,no_=cr(DO)
bt=lambda x: float(np.std([np.corrcoef(x[i],y[i])[0,1] for i in
    (rng.choice(np.flatnonzero(np.isfinite(x)&np.isfinite(y)),int(np.sum(np.isfinite(x)&np.isfinite(y))),True)
     for _ in range(200))]))
print(f"\n第二问 —— 跨仪器复现那条主张:")
print(f"  D_块  ↔ 羞耻 {rb:+.4f} ± {bt(DB):.4f}(n = {nb_:,};`#232b` 报 +0.0862)")
print(f"  **D_起始 ↔ 羞耻 {ro:+.4f} ± {bt(DO):.4f}(n = {no_:,})**")

# ⚠ 常驻混杂:D_起始 是不是就是「报了多少个类别」
NCAT=np.where(okO,OBS.sum(1),np.nan).astype(float)
mc=np.isfinite(DO)&np.isfinite(NCAT)
print(f"  ⚠ corr(D_起始, 报告类别数) = {np.corrcoef(DO[mc],NCAT[mc])[0,1]:+.4f}")
mm2=np.isfinite(DO)&np.isfinite(y)&np.isfinite(NCAT)
Xc=np.column_stack([np.ones(mm2.sum()),(DO[mm2]-DO[mm2].mean())/DO[mm2].std(),
                    (NCAT[mm2]-NCAT[mm2].mean())/NCAT[mm2].std()])
bb=np.linalg.lstsq(Xc,(y[mm2]-y[mm2].mean())/y[mm2].std(),rcond=None)[0]
sdb=float(np.std([np.linalg.lstsq(np.column_stack([np.ones(len(ix))]+[Xc[ix,1],Xc[ix,2]]),
    ((y[mm2]-y[mm2].mean())/y[mm2].std())[ix],rcond=None)[0][1]
    for ix in (rng.choice(mm2.sum(),mm2.sum(),True) for _ in range(200))]))
print(f"  **控制类别数后 D_起始 ↔ 羞耻 beta = {bb[1]:+.4f} ± {sdb:.4f}**(类别数自己 {bb[2]:+.4f})")
print(f"  解衰减(按各自信度):D_块 {rb/np.sqrt(RB):+.4f} · D_起始 {ro/np.sqrt(RO):+.4f}")

T=pd.DataFrame([dict(instrument='块',rel=RB,r_shame=rb,n=nb_),
                dict(instrument='起始',rel=RO,r_shame=ro,n=no_)])
check_columns(T,'R280'); T.to_csv(pathlib.Path(__file__).parent/'results'/'cross_instrument.csv',index=False)

g=Gate('第四件事换一把仪器还在不在')
g.asserted('两把仪器的题目不相交(#173 同款断言)',True,f"块列 {NB} / 起始列 {len(ons)};交集 0")
g.asserted('正对照:块内部两个不相交半边必须达到它们自己的衰减上限 —— 验证「上限」这套算法',
           pin>0.6*pceil, f"实测 {pin:+.4f} vs 上限 {pceil:.4f}({100*pin/pceil:.1f}%)")
g.negative_control('置换人配对',abs(float(np.mean(nl))),abs(obs),null_spread=float(np.std(nl)),
                   null_kind='把 D_起始 在人之间置换 —— 保留两个分布,只打掉配对')
g.offset_control('★ 跨仪器相关 vs 衰减上限',abs(obs),ceil,float(np.std(nl)),
                 null_kind='sqrt(两把仪器信度之积)—— 不是零假设,是「若两把仪器测的是同一个构念,相关最多能到哪」')
g.asserted('★ 注册的 kill:跨仪器相关达到上限的一半以上 -> 第四件事是关于人的,不是块格式的',
           abs(obs)>0.5*ceil, f"{obs:+.4f} vs 上限 {ceil:.4f}({100*obs/ceil:.1f}%)")
g.asserted('第二问:D_起始 复现 D_块 的羞耻方向',
           ro>0 and rb>0, f"D_块 {rb:+.4f} · D_起始 {ro:+.4f}")
g.asserted('⚠ 常驻混杂:控制报告类别数后 D_起始 的羞耻效应必须存活(> 2× 展布)',
           abs(bb[1])>2*sdb, f"beta {bb[1]:+.4f} ± {sdb:.4f};corr(D_起始, 类别数) = "
           f"{np.corrcoef(DO[mc],NCAT[mc])[0,1]:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A90 R335 -- `#179` 那个 +0.1185 拆开是多少

`#289b`:位置分 ↔ 羞耻 **+0.1185** 是在**合并样本**上算的,
而 `#289a` 说 **S 的预测在两组间不同**(留一:去掉 S 非不变性完全消失)。
**所以那个数是两个不同数字的平均,而这个项目从没报过它们各是多少。**

ESTIMAND        两组各自的 `corr(位置分 S, 羞耻)` 与自助展布;
                两者的差,与 **offset = 随机劈同样大小的两组之间的差**比。
⚠ 零不该是零      任意两组之间本来就有抽样差别。
KILL            **若两组的羞耻相关可分辨地不同 -> `#179` 那条最老的主张必须按组报,公开页要改;
                若相当 -> S 的非不变性不在羞耻这一格上,而那也要写清楚(它在别的格上)。**
⚠ 同报           两组各自的 n 与**组内**位置分信度(`#285b`:半样本上单个量会不稳)。
POSITIVE CTRL   两端:① 一个**已知按性别不同**的合成结局必须被判出差别;
                ② `agreeableness`(与性别只相关 −0.0749)当作结局必须判不出。
NEGATIVE CTRL   跨人置换羞耻(**只在有限值内**)。
IMPOSSIBLE      `biomale` 是二值自报(`#209` 同款登记)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl))
NB=len(MB); cov=np.zeros(NN); pos=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
ok=cov>=8; S=np.where(ok,pos/np.maximum(cov,1),np.nan)
SEX=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
SH=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SH],errors='coerce').values.astype(float)
def S_half(cols):
    p2=np.zeros(NN); c2=np.zeros(NN)
    for b in cols:
        M,ppl=MB[b]; rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
        p2[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0); c2[ppl]+=1
    return np.where(c2>=len(cols)//2,p2/np.maximum(c2,1),np.nan)
rngB=np.random.default_rng(20260804)
def cr(a,b,rows):
    m=np.zeros(NN,bool); m[rows]=True; m&=np.isfinite(a)&np.isfinite(b)
    if m.sum()<300: return np.nan,np.nan,0
    r=float(np.corrcoef(a[m],b[m])[0,1])
    sd=float(np.std([np.corrcoef(a[i],b[i])[0,1] for i in
        (rngB.choice(np.flatnonzero(m),int(m.sum()),True) for _ in range(300))]))
    return r,sd,int(m.sum())
have=ok&np.isfinite(SEX)&np.isfinite(y)
g0=np.flatnonzero(have&(SEX==0)); g1=np.flatnonzero(have&(SEX==1))
r_all,sd_all,n_all=cr(S,y,np.flatnonzero(have))
r0,sd0,n0=cr(S,y,g0); r1,sd1,n1=cr(S,y,g1)
print(f"合并样本:**{r_all:+.4f} ± {sd_all:.4f}**(n={n_all:,};`#179` 报 +0.1185)")
print(f"  `biomale=0`:**{r0:+.4f} ± {sd0:.4f}**(n={n0:,})")
print(f"  `biomale=1`:**{r1:+.4f} ± {sd1:.4f}**(n={n1:,})")
print(f"  **差 = {r1-r0:+.4f}**")
def relS(rows):
    vs=[]
    for s in range(3):
        p=np.random.default_rng(120+s).permutation(NB); h=NB//2
        a,b=S_half(p[:h]),S_half(p[h:])
        m=np.zeros(NN,bool); m[rows]=True; m&=np.isfinite(a)&np.isfinite(b)
        if m.sum()>400:
            r=float(np.corrcoef(a[m],b[m])[0,1]); vs.append(2*abs(r)/(1+abs(r)))
    return float(np.nanmean(vs))
print(f"  ⚠ 组内位置分信度:组 0 **{relS(g0):+.4f}** · 组 1 **{relS(g1):+.4f}**")
rngS=np.random.default_rng(4242)
RND=[]
for t in range(60):
    p=rngS.permutation(np.flatnonzero(have))
    a,_,_=cr(S,y,p[:len(g0)]); b,_,_=cr(S,y,p[len(g0):len(g0)+len(g1)])
    if np.isfinite(a) and np.isfinite(b): RND.append(b-a)
off=float(np.std(RND))
print(f"\n**offset(随机劈同样大小两组之间的差)= 0 ± {off:.4f}**"
      f"  -> 观测差 {r1-r0:+.4f} 是它的 **{abs(r1-r0)/max(off,1e-9):.1f}×**")
# ⚠ 第一版的正对照①是 `0.3*SEX + 噪声` —— **那只造成均值差,不造成相关的组间差**。
#    要让 `corr(S, 结局)` 在两组间不同,结局必须依赖 **SEX × S 的交互**,不是 SEX 的主效应。
#    这与 `#285c` 同族:**正对照的失败,几乎总是因为我给它的那个东西不具备被检验的性质。**
n_=rngS.standard_normal(NN)
zS=np.where(np.isfinite(S),(S-np.nanmean(S))/np.nanstd(S),0.0)
FK=np.where(have,0.35*np.nan_to_num(SEX)*zS+n_,np.nan)     # 交互项驱动
FK_BAD=np.where(have,0.3*np.nan_to_num(SEX)+n_,np.nan)     # 旧的那个,保留作对照
INV=np.asarray(pd.to_numeric(d['agreeablenessvariable'],errors='coerce').values,dtype=float)
def gap(yy):
    a,_,_=cr(S,yy,g0); b,_,_=cr(S,yy,g1); return b-a
print(f"正对照两端:① **交互项驱动**的合成结局 -> 差 **{gap(FK):+.4f}**(必须 > 2×offset)· "
      f"② `agreeableness` -> 差 **{gap(INV):+.4f}**(必须 ≈0)")
print(f"   ⚠ 对照:旧版正对照(SEX 主效应驱动)-> 差 **{gap(FK_BAD):+.4f}** —— "
      f"**一个只造成均值差的结局,造不出相关的组间差**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2)); z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[gap(perm_finite(y,300+i)) for i in range(20)]
print(f"负对照(置换羞耻,只在有限值内):{np.mean(nul):+.4f} ± {np.std(nul):.4f}")
T=pd.DataFrame([dict(grp='合并',r=r_all,sd=sd_all,n=n_all),
                dict(grp='biomale=0',r=r0,sd=sd0,n=n0),
                dict(grp='biomale=1',r=r1,sd=sd1,n=n1)])
check_columns(T,'R335'); T.to_csv(pathlib.Path(__file__).parent/'results'/'shame_by_group.csv',index=False)

g=Gate('`#179` 那个 +0.1185 拆开是多少')
g.asserted('⚠ 第一版正对照①不合身(SEX 主效应只造均值差,不造相关的组间差),已改为交互项驱动',
           abs(gap(FK_BAD))<2*off, f"旧版 {gap(FK_BAD):+.4f} —— 它本来就测不出东西")
g.asserted('正对照两端:交互项驱动的合成结局必须被判出、`agreeableness` 必须判不出',
           abs(gap(FK))>2*off and abs(gap(INV))<2*off,
           f"① {gap(FK):+.4f} · ② {gap(INV):+.4f} · offset ±{off:.4f}")
g.negative_control('置换羞耻',abs(float(np.mean(nul))),abs(r1-r0),
                   null_spread=float(np.std(nul)),null_kind='跨人置换羞耻(只在有限值内)—— 只打掉配对')
g.has_error_bar('组 0 的 S↔羞耻',r0,sd0,'bootstrap_人层')
g.has_error_bar('组 1 的 S↔羞耻',r1,sd1,'bootstrap_人层')
g.offset_control('★ 两组之差 vs 随机劈两组之差',float(r1-r0),0.0,off,
                 null_kind='随机劈同样大小两组之间的差 —— 不是零假设,是「若 S↔羞耻不随性别变,'
                           '两组该差多少」')
g.asserted('★ 注册的 kill:两组的羞耻相关可分辨地不同 -> `#179` 必须按组报',
           abs(r1-r0)>2*off,
           f"组 0 {r0:+.4f} ± {sd0:.4f} · 组 1 {r1:+.4f} ± {sd1:.4f};"
           f"差 {r1-r0:+.4f} vs offset ±{off:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

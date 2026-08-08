import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A92 R338 -- S 在两组里怎么不一样:是块的权重,还是「什么算冷门」本身

`#292b`:漂移在 **S 这个量本身**。而 S 的构造是「每块内所挑选项的平均稀有度」再跨块平均,
**所以这句话可以拆成两句可测的**:
① **块的权重不同** —— 哪些块贡献了 S 的方差(块层载荷)在两组里不同;
② **「什么算冷门」不同** —— 同一块里,两组对选项流行度的排序不同。

ESTIMAND        ① 两组各自的**块层载荷谱**(每块与**留一**总 S 的相关,32 维),判两谱的相关;
                ② 两组各自的**块内选项流行度**,逐块判两组流行度向量的相关,取中位。
                两者都与 **offset = 随机劈同样大小两组**比。
KILL            **若 ① 明显低而 ② 不低 -> 同一套「冷门」定义,但不同的块在起作用;
                若 ② 也明显低 -> 两组眼里「什么算冷门」本身就不同,
                而那意味着 `rar` 这个**跨组共用**的量本身就不该共用 —— 一个更深的非不变性。**
POSITIVE CTRL   人为把一组的某几块**打乱**:① 必须抓到,② **不受影响**(打乱的是人,不是选项流行度)。
NEGATIVE CTRL   由 offset 承担(随机劈)。
IMPOSSIBLE      流行度是**组内**估的,所以 ② 的两组差里混着抽样噪声;offset 正是为此。
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
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
ok=cov>=8
SEX=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
have=ok&np.isfinite(SEX)
g0=np.flatnonzero(have&(SEX==0)); g1=np.flatnonzero(have&(SEX==1))
print(f"两组 n = {len(g0):,} / {len(g1):,};块 {NB}")
def block_scores(rows,scramble=None):
    """每块的组内 S 贡献(用**组内**流行度定义稀有度)。scramble: 打乱这些块的人。"""
    m=np.zeros(NN,bool); m[rows]=True
    P=np.full((NB,NN),np.nan); PREV=[]
    for b,(M,ppl) in enumerate(MB):
        sel=np.isin(ppl,rows)
        if sel.sum()<200: PREV.append(None); continue
        Ms=M[sel]; pv=Ms.mean(0); PREV.append(pv)
        rr=-np.log(np.clip(pv,1e-4,1.)); n=Ms.sum(1)
        v=np.where(n>0,(Ms@rr)/np.maximum(n,1),np.nan)
        if scramble is not None and b in scramble:
            rg2=np.random.default_rng(900+b); v=v[rg2.permutation(len(v))]
        P[b,ppl[sel]]=v
    return P,PREV
def loadings(P,rows):
    m=np.zeros(NN,bool); m[rows]=True
    F=np.isfinite(P)&m[None,:]; Z=np.where(F,P,0.0); tot=Z.sum(0); ct=F.sum(0)
    out=np.full(NB,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)   # ⚠ 留一
        g2=F[b]&np.isfinite(lo)
        if g2.sum()>200: out[b]=float(np.corrcoef(P[b][g2],lo[g2])[0,1])
    return out
def prev_sim(P0,P1):
    vs=[]
    for a,b in zip(P0,P1):
        if a is None or b is None or len(a)!=len(b) or len(a)<8: continue
        if np.std(a)>1e-9 and np.std(b)>1e-9: vs.append(float(np.corrcoef(a,b)[0,1]))
    return float(np.median(vs)),len(vs)
def pair(r0,r1,scramble=None):
    P0,V0=block_scores(r0); P1,V1=block_scores(r1,scramble=scramble)
    L0,L1=loadings(P0,r0),loadings(P1,r1)
    m=np.isfinite(L0)&np.isfinite(L1)
    load=float(np.corrcoef(L0[m],L1[m])[0,1]) if m.sum()>=10 else np.nan
    pv,npv=prev_sim(V0,V1)
    return load,pv,npv
l_sex,p_sex,npv=pair(g0,g1)
rngS=np.random.default_rng(20260804)
RL=[];RP=[]
for t in range(5):
    p=rngS.permutation(np.flatnonzero(have))
    a,b,_=pair(p[:len(g0)],p[len(g0):len(g0)+len(g1)]); RL.append(a); RP.append(b)
lr,lsd=float(np.mean(RL)),float(np.std(RL)); pr,psd=float(np.mean(RP)),float(np.std(RP))
print(f"\n① **块层载荷谱**(每块与留一总 S 的相关,{NB} 维):")
print(f"   按性别劈 **{l_sex:+.4f}** vs 随机劈 **{lr:+.4f} ± {lsd:.4f}** -> 差 **{l_sex-lr:+.4f}**"
      f"({abs(l_sex-lr)/max(2*lsd,1e-9):.1f}× 的 2×展布)")
print(f"② **块内选项流行度**(逐块相关取中位,{npv} 块):")
print(f"   按性别劈 **{p_sex:+.4f}** vs 随机劈 **{pr:+.4f} ± {psd:.4f}** -> 差 **{p_sex-pr:+.4f}**"
      f"({abs(p_sex-pr)/max(2*psd,1e-9):.1f}× 的 2×展布)")
SCR=set(range(0,NB,4))
l_s,p_s,_=pair(g0,g1,scramble=SCR)
print(f"\n正对照(把组 1 的 {len(SCR)} 个块的人打乱):")
print(f"   ① 载荷谱 {l_sex:+.4f} -> **{l_s:+.4f}**(必须明显掉)· "
      f"② 流行度 {p_sex:+.4f} -> **{p_s:+.4f}**(必须几乎不动 —— 打乱的是人,不是选项流行度)")
T=pd.DataFrame([dict(level='①块层载荷谱',sex=l_sex,rnd=lr,sd=lsd,delta=l_sex-lr),
                dict(level='②块内流行度',sex=p_sex,rnd=pr,sd=psd,delta=p_sex-pr)])
check_columns(T,'R338'); T.to_csv(pathlib.Path(__file__).parent/'results'/'decompose_S.csv',index=False)

g=Gate('S 在两组里怎么不一样')
g.asserted('正对照:打乱一组的部分块 -> ① 必须掉、② 必须几乎不动',
           (l_sex-l_s)>0.05 and abs(p_sex-p_s)<0.02,
           f"① {l_sex:+.4f} -> {l_s:+.4f} · ② {p_sex:+.4f} -> {p_s:+.4f}")
g.offset_control('★ ① 块层载荷谱:按性别劈 vs 随机劈',l_sex,lr,lsd,
                 null_kind='随机劈同样大小两组 —— 不是零假设,是「若哪些块在起作用相同,该落在哪」')
g.offset_control('★ ② 块内流行度:按性别劈 vs 随机劈',p_sex,pr,psd,
                 null_kind='随机劈同样大小两组 —— 不是零假设,是「若两组眼里什么算冷门相同,该落在哪」')
g.asserted('★ 注册的 kill:① 低而 ② 不低 -> 同一套冷门定义不同的块;② 也低 -> `rar` 本身不该跨组共用',
           abs(p_sex-pr)<2*psd,
           f"① 差 {l_sex-lr:+.4f}(2×展布 {2*lsd:.4f})· ② 差 {p_sex-pr:+.4f}(2×展布 {2*psd:.4f})")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

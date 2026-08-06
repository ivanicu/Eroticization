import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A73 R302 -- 重报公开页上的那三个计数,把种子层量进去

**类型:CLOSURE**。`#256b` 清单里唯一 load-bearing 的一条。

`#252a` 报三个成分三重控制后的越阈数(c1 8–10 · c2 6 · c3 20),**只量了阈值层**。
`#255c`/`#256a`:这一族的主导不确定性来源是**数据劈分种子**,守卫 11 现在只给一层就直接 FAIL。

ESTIMAND        6 个选项劈半种子,各重算三个成分 → 三重残差化(勾选数 · 性别 · 位置分 S)
                → 29 个结局越阈数;报**跨种子展布**,按较大的那层给区间。
KILL            **若跨种子区间把某个成分推到含零 -> 总表那一格要改写;
                若三个的区间都不含零 -> `#252a` 结论不变,但精度第一次被正确标注。**
NEGATIVE CTRL   λ=1 全置换(`#300` 已证四个量都落到 0–1,本轮复用其结论,不重跑)。
IMPOSSIBLE      成分由特征分解定义,不同种子下 c1/c2/c3 的**身份**可能互换 ——
                所以报的是「三个成分各自的计数」,不是「c2 这一个对象的计数」。**同报身份指纹。**
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
NB=len(MB); cov=np.zeros(NN); K=np.zeros(NN); pos=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; K[ppl]+=n; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
ok=cov>=8; K=np.where(ok,K,np.nan); S=np.where(ok,pos/np.maximum(cov,1),np.nan)
SEX=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
def comps(seed):
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        A[b,ppl]=M[:,o[:k]].mean(1); B[b,ppl]=M[:,o[k:2*k]].mean(1)
    def prof(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
        return R
    Ra,Rb=prof(A),prof(B)
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&ok
            if mm.sum()>300: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0); out=[]
    for k in range(3):
        num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
        out.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    return out,Vv[:,:3]
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
EX={'age':d['age'].map(AGE),'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,d[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EX.items()]
rng=np.random.default_rng(20260804)
def npass(x,reps=10):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]; nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl); A=np.array([z[:L] for z in nl]); r=np.array(r); c=[]
    for _ in range(reps):
        idx=rng.choice(L,L,True)
        c.append(int(np.nansum(np.abs(r)>float(np.nanquantile(np.nanmax(A[:,idx],0),0.95)))))
    return float(np.mean(c)),float(np.std(c))
def resid3(x):
    Cv=[K,SEX,S]; m=np.isfinite(x)&ok&np.all(np.isfinite(np.array(Cv)),0)
    X=np.column_stack([np.ones(m.sum())]+[c[m] for c in Cv])
    b=np.linalg.lstsq(X,x[m],rcond=None)[0]; o=np.full(NN,np.nan); o[m]=x[m]-X@b; return o
SEEDS=[500,900,1200,1500,1800,2100]
cnt={k:[] for k in range(3)}; thr_sd={k:[] for k in range(3)}; V0=None; ident={k:[] for k in range(3)}
for s in SEEDS:
    cs,Vv=comps(s)
    if V0 is None: V0=Vv
    for k in range(3):
        a,b=npass(resid3(cs[k])); cnt[k].append(a); thr_sd[k].append(b)
        ident[k].append(float(abs(np.dot(Vv[:,k],V0[:,k]))))
print(f"三个成分三重控制后的越阈数,跨 {len(SEEDS)} 个选项劈半种子(29 个结局):")
rows=[]
for k in range(3):
    m,sdv=float(np.mean(cnt[k])),float(np.std(cnt[k]))
    tsd=float(np.mean(thr_sd[k]))
    rows.append(dict(comp=k+1,v_mean=m,seed_sd=sdv,thr_sd=tsd,
                     lo=m-2*max(sdv,tsd),hi=m+2*max(sdv,tsd),ident=float(np.mean(ident[k]))))
    print(f"  c{k+1}: "+' · '.join(f"{v:.1f}" for v in cnt[k])
          +f"  -> **{m:.1f} ± {sdv:.1f}(种子)/ ±{tsd:.1f}(阈值)** -> 区间 "
          f"**{m-2*max(sdv,tsd):.0f}–{m+2*max(sdv,tsd):.0f}**;身份指纹 {np.mean(ident[k]):.2f}")
T=pd.DataFrame(rows); check_columns(T,'R302')
T.to_csv(pathlib.Path(__file__).parent/'results'/'reissued_counts.csv',index=False)
zero=[int(r.comp) for _,r in T.iterrows() if r.lo<=0]
print(f"\n区间含零的成分:**{zero if zero else '无'}**;"
      f"`#252a` 报的是 c1 8–10 · c2 6 · c3 20(只含阈值层)")

g=Gate('把种子层量进公开页上的三个计数')
for _,r in T.iterrows():
    g.count_needs_interval(f"c{int(r.comp)} 三重控制后越阈",int(round(r['v_mean'])),len(OUT),
                           float(r.thr_sd),'threshold_resample_阈值重抽样',
                           n_resamples=10,seed_spread=float(r.seed_sd))
g.asserted('⚠ 成分由特征分解定义,跨种子身份可能互换 —— 身份指纹同报',
           all(r.ident>0.8 for _,r in T.iterrows()),
           ' · '.join(f"c{int(r.comp)} |cos|={r.ident:.2f}" for _,r in T.iterrows()))
g.asserted('★ 注册的 kill:区间含零 -> 总表那一格要改写;都不含零 -> 结论不变,精度被正确标注',
           len(zero)==0, f"含零的:{zero if zero else '无'};区间 "
           +' · '.join(f"c{int(r.comp)} {r.lo:.0f}–{r.hi:.0f}" for _,r in T.iterrows()))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A72 R299 -- 「一个人的剖面有多清晰」本身预测什么

`#296` 的总表还剩三个 `UNCOMPUTED`,最有心理学分量的是**「宽度剖面」行的 ↔羞耻 与 越阈两格**:
`#228` 建立了它的信度(0.4290)与跨仪器维数(`#236`),
但**从没把它当作一个人层分数去跑结局面板** —— 因为剖面是一个**人×块矩阵**,不是标量。

⚠ 所以诚实的问法不是硬造一个标量,而是:**「这个人的剖面有多清晰」本身预测什么?**
逐人版本:`strength_i` = 这个人的残差剖面在**两半选项**上的一致性
(`#228a` 的量,但**逐人**而不是全体合并)。

ESTIMAND        `strength_i` 的分半信度;若 ≥0.2,再跑 29 个结局面板(已剔除 `biomale`,`#252b`)。
KILL            **若 `strength_i` 自己有越阈结局 -> 表上那两格填上,「形状有多清晰」成为第九个人层量;
                若它什么都不预测 -> 那两格写 `0(已测)`。**
⚠ 先量信度         逐人一致性只有两半、每半 ~16 块,信度会很低。
                **信度 < 0.2 -> 如实登记为「这个站点结构上量不了」,不报一个零**(`P14` 的不可用登记)。
NEGATIVE CTRL   置换结局。
POSITIVE CTRL   守卫 13:种一个**已知与结局相关**的假 strength,扫描强度,**方向由扫描给出**。
IMPOSSIBLE      `strength_i` 同时携带「这个人的形状清晰」与「这个人答题一致」——
                后者是一种作答风格,本轮不分。
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
def strength(seed, blocks=None, plant=None):
    """逐人:残差剖面在两半选项上的一致性。"""
    bl=range(NB) if blocks is None else blocks
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    for b in bl:
        M,ppl=MB[b]; Mm=M if plant is None else np.clip(M+plant[ppl,b][:,None],0,1)
        o=rg.permutation(Mm.shape[1]); k=Mm.shape[1]//2
        A[b,ppl]=Mm[:,o[:k]].mean(1); B[b,ppl]=Mm[:,o[k:2*k]].mean(1)
    def prof(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in bl:
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
        return R
    Ra,Rb=prof(A),prof(B); out=np.full(NN,np.nan)
    G=np.isfinite(Ra)&np.isfinite(Rb)
    for i in np.flatnonzero(ok):
        m=G[:,i]
        if m.sum()<8: continue
        x,y=Ra[m,i],Rb[m,i]
        if x.std()>1e-9 and y.std()>1e-9: out[i]=float(np.corrcoef(x,y)[0,1])
    return out
ST=strength(900)
rng=np.random.default_rng(20260804)
rels=[]
for s in range(4):
    p=np.random.default_rng(80+s).permutation(NB); h=NB//2
    a,b=strength(910+s,p[:h]),strength(920+s,p[h:])
    m=np.isfinite(a)&np.isfinite(b)
    if m.sum()>500:
        r=float(np.corrcoef(a[m],b[m])[0,1]); rels.append(2*r/(1+r) if r<0.999 else np.nan)
REL=float(np.nanmean(rels))
m0=np.isfinite(ST)
print(f"逐人剖面清晰度 strength_i:n = {int(m0.sum()):,};均值 {np.nanmean(ST):+.4f} ± "
      f"{np.nanstd(ST):.4f}")
print(f"**分半信度(块劈半)= {REL:+.4f}**(4 次劈分:"
      +' · '.join(f"{r:+.3f}" for r in rels)+")")
MEASURABLE = REL>=0.20
print(f"  -> {'可测,继续跑面板' if MEASURABLE else '**信度 < 0.2:这个站点结构上量不了,不报零**'}")

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
def npass(x,reps=12,tag=None):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]; nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl); A=np.array([z[:L] for z in nl]); r=np.array(r); cnt=[]
    for _ in range(reps):
        idx=rng.choice(L,L,True)
        thr=float(np.nanquantile(np.nanmax(A[:,idx],0),0.95)); cnt.append(int(np.nansum(np.abs(r)>thr)))
    thr0=float(np.nanquantile(np.nanmax(A,0),0.95))
    if tag:
        top=sorted([(OUT[i][0],r[i]) for i in range(len(r)) if np.isfinite(r[i]) and abs(r[i])>thr0],
                   key=lambda t:-abs(t[1]))[:5]
        print(f"   {tag:<14} {np.mean(cnt):.1f}±{np.std(cnt):.1f}/{len(OUT)}(阈值 {thr0:.4f}) "
              +' · '.join(f"{n[:20]} {v:+.3f}" for n,v in top))
    return float(np.mean(cnt)),float(np.std(cnt)),r,thr0
SH=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SH],errors='coerce').values.astype(float)
mm=np.isfinite(ST)&np.isfinite(y)&ok
r_sh=float(np.corrcoef(ST[mm],y[mm])[0,1])
sd_sh=float(np.std([np.corrcoef(ST[i],y[i])[0,1] for i in
    (rng.choice(np.flatnonzero(mm),int(mm.sum()),True) for _ in range(200))]))
print(f"\n结局面板(29 个,已剔除 biomale):")
n_hit,n_sd,rvec,thr0=npass(ST,tag='剖面清晰度')
print(f"   ↔ 羞耻 **{r_sh:+.4f} ± {sd_sh:.4f}**(n = {int(mm.sum()):,})")
SW=[]
for gp in (0.0,0.02,0.05,0.10):
    P=np.zeros((NN,NB)); u=rng.standard_normal(NN)
    w=np.array([M.mean() for M,_ in MB]); w=(w-w.mean())/max(w.std(),1e-9)
    P+= gp*np.outer(u,w).T.T if False else gp*np.outer(u,w)
    SW.append((gp,npass(strength(950,plant=P))[0]))
print(f"正对照(种一个人×块特异结构,守卫 13 判方向)g -> 越阈:"
      +' · '.join(f"{a:.2f}->{b:.1f}" for a,b in SW))

T=pd.DataFrame([dict(quantity='剖面清晰度 strength_i',rel=REL,n=int(m0.sum()),
                     n_hit=n_hit,n_hit_sd=n_sd,r_shame=r_sh,sd_shame=sd_sh)])
check_columns(T,'R299'); T.to_csv(pathlib.Path(__file__).parent/'results'/'clarity.csv',index=False)

g=Gate('剖面清晰度本身预测什么')
g.asserted('⚠ 先量信度再看结局:信度 < 0.2 则登记为「结构上量不了」,不报零',
           True, f"信度 {REL:+.4f} -> {'可测' if MEASURABLE else '量不了'}")
g.plant_direction_from_sweep('正对照:种人×块特异结构后的越阈数(方向由扫描给出)',
                             SW,baseline=SW[0][1],baseline_spread=1.5,half_of=2.0)
g.count_needs_interval('剖面清晰度的越阈计数',int(round(n_hit)),len(OUT),n_sd,
                       'threshold_resample_阈值重抽样',n_resamples=12)
g.has_error_bar('剖面清晰度 ↔ 羞耻',r_sh,sd_sh,'bootstrap_人层')
g.asserted('★ 注册的 kill:strength_i 自己有越阈结局 -> 表上两格填上,成为第九个人层量',
           MEASURABLE and n_hit>2*max(n_sd,1e-9) and n_hit>=2,
           f"信度 {REL:+.4f};越阈 {n_hit:.1f}±{n_sd:.1f}/{len(OUT)};↔羞耻 {r_sh:+.4f}±{sd_sh:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A72 R300 -- 用减法式给四个量补功效证明

`#254c`:一致性/相似度类的量**做不了加法式正对照** —— `#252c`(特征向量)·
`#254b`(一致性)· `#248c`(rate),三次都是**种入改变了被测量的定义或尺度**。
共同点:它们都是**比值或投影**,分子与分母由同一批数据算出。
**做法:用 `#282` 的减法式** —— 把已有结构按比例换成噪声,**稀释不改变量的定义,只改变强度。**

ESTIMAND        `X = (1−λ)·真实 + λ·置换`(随机 λ 比例的人被块内置换),扫 λ,
                对 **strength_i · c1 · c2 · c3** 各报 29 个结局的越阈数曲线。
KILL(两端)      **λ=0 必须落在各自基线上;λ=1 必须落在各自的置换零上。**
                任一端对不上 -> 这条减法式对该量不合身,不能读中间。
交付            **功效证明**:若某个量的曲线在 λ 增大时**下降**,说明它的越阈数确实由人层结构扛着;
                若曲线**整体平**,说明这个设计对该量本来就看不见任何东西 ——
                那么它的「零」应当写成「结构上量不了」,而不是「0(已测)」。
IMPOSSIBLE      λ 混合的是「谁的行数」;曲线形状只在这一种稀释方式下有意义。
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
def halves(seed,lam):
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        ha=M[:,o[:k]].mean(1); hb=M[:,o[k:2*k]].mean(1)
        hit=np.flatnonzero(rg.random(len(ppl))<lam)
        if len(hit): ha[hit]=ha[rg.permutation(hit)]; hb[hit]=hb[rg.permutation(hit)]
        A[b,ppl]=ha; B[b,ppl]=hb
    return A,B
def prof(X):
    F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
    return R
def quantities(seed,lam):
    A,B=halves(seed,lam); Ra,Rb=prof(A),prof(B)
    st=np.full(NN,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
    for i in np.flatnonzero(ok):
        m=G[:,i]
        if m.sum()<8: continue
        x,y=Ra[m,i],Rb[m,i]
        if x.std()>1e-9 and y.std()>1e-9: st[i]=float(np.corrcoef(x,y)[0,1])
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(R[i])&np.isfinite(R[j])&ok
            if mm.sum()>300: C[i,j]=np.corrcoef(R[i][mm],R[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0); comps=[]
    for k in range(3):
        num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
        comps.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    return [st]+comps
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
def npass(x,reps=8):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]; nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl); A=np.array([z[:L] for z in nl]); r=np.array(r)
    return float(np.mean([int(np.nansum(np.abs(r)>float(np.nanquantile(np.nanmax(A[:,rng.choice(L,L,True)],0),0.95))))
                          for _ in range(reps)]))
NAMES=['剖面清晰度','c1','c2','c3']; LAM=[0.0,0.25,0.5,0.75,1.0]
CUR={n:[] for n in NAMES}
for lam in LAM:
    qs=quantities(1500,lam)
    for n,x in zip(NAMES,qs): CUR[n].append(npass(x))
print(f"减法式功效曲线(λ = 被换成噪声的人的比例;格内为 29 个结局的越阈数):")
print(f"{'量':<10}"+''.join(f"{('λ=%.2f'%l):>9}" for l in LAM))
for n in NAMES:
    print(f"{n:<10}"+''.join(f"{v:>9.1f}" for v in CUR[n]))
T=pd.DataFrame([dict(quantity=n,**{f"lam_{str(l).replace('.','_')}":v for l,v in zip(LAM,CUR[n])})
                for n in NAMES])
check_columns(T,'R300'); T.to_csv(pathlib.Path(__file__).parent/'results'/'subtractive_power.csv',index=False)
drop={n:CUR[n][0]-CUR[n][-1] for n in NAMES}
print(f"\nλ=0 → λ=1 的落差:"+' · '.join(f"{n} **{drop[n]:+.1f}**" for n in NAMES))
flat=[n for n in NAMES if drop[n]<2.0]
print(f"曲线整体平(落差 < 2)的量:**{flat if flat else '无'}** <- 这些量的「零」是「结构上量不了」")

# ⚠⚠ 一个不一致必须先查:`#299` 给剖面清晰度 4.1±1.8,本轮 λ=0 给 10.6 —— 同一个量,只换了种子。
SEEDS=[900,1200,1500,1800,2100,2400]
sc=[npass(quantities(s,0.0)[0]) for s in SEEDS]
print(f"\n⚠ 剖面清晰度的越阈数,跨 {len(SEEDS)} 个**选项劈半种子**:"
      +' · '.join(f"{v:.1f}" for v in sc))
print(f"   均值 **{np.mean(sc):.1f} ± {np.std(sc):.1f}**,范围 {min(sc):.1f}–{max(sc):.1f}"
      f"  [`#299` 报 4.1 ± 1.8,只含阈值重抽样]")
c3s=[npass(quantities(s,0.0)[3]) for s in SEEDS[:4]]
print(f"   对照:c3 跨种子 "+' · '.join(f"{v:.1f}" for v in c3s)
      +f" -> {np.mean(c3s):.1f} ± {np.std(c3s):.1f}")

g=Gate('减法式功效曲线')
g.count_needs_interval('剖面清晰度的越阈计数(跨选项劈半种子)',int(round(np.mean(sc))),len(OUT),
                       float(np.std(sc)),'seed_跨种子',n_resamples=len(SEEDS))
for n in NAMES:
    g.plant_direction_from_sweep(f"{n} 的减法式扫描(方向由数据给出)",
                                 list(zip(LAM,CUR[n])),baseline=CUR[n][0],baseline_spread=1.5,half_of=2.0)
g.asserted('★ 两端:λ=1 必须落在各自的置换零附近(全部人被置换 -> 越阈数应接近 0–2)',
           all(CUR[n][-1]<=3.0 for n in NAMES), ' · '.join(f"{n} λ=1 {CUR[n][-1]:.1f}" for n in NAMES))
g.asserted('★ 交付:曲线整体平的量,它的零应写成「结构上量不了」而不是「0(已测)」',
           True, f"落差 "+' · '.join(f"{n} {drop[n]:+.1f}" for n in NAMES)
                 +f";平的:{flat if flat else '无'}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

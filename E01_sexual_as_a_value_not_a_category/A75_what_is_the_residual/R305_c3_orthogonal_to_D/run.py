import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A75 R305 -- c3 里扣掉 D 之后,还剩什么

`#259a`:c3 与 D 的剖面差 **14%** 是可分辨的,**但没说差在哪里**。

ESTIMAND        ① 逐个结局比两条剖面,列出 |r_c3| − |r_D| 最大的几个;
                ② 把 **D 从 c3 里回归掉**,用 `c3⊥D` 重跑 29 个结局面板(**两层区间**)。
KILL            **若 `c3⊥D` 仍有越阈结局(区间不含零)-> 那 14% 是独立的东西,值得单独命名;
                若残差什么都不剩 -> c3 只是「更灵敏的 D」,总表该合并成一行。**
POSITIVE CTRL   把 D 从「D 自己的带噪声复制」里回归掉,**残差必须什么都不剩** ——
                证明这个减法有效。
NEGATIVE CTRL   置换结局。
⚠ 口径           **全仪器**构造(不劈块),并同报全仪器的 `corr(c3, D)` ——
                `#259b`:半块与全仪器不是同一个口径。
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
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
def build(seed):
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
    def z(v):
        m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
    def hs(R,cols):
        sub=R[cols]; F2=np.isfinite(sub)
        return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
    D=(z(hs(Ra,TRG))-z(hs(Ra,ORD))+z(hs(Rb,TRG))-z(hs(Rb,ORD)))/2
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&ok
            if mm.sum()>300: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(Vv[:,2][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,2])[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan),D
C3,D=build(500)
m0=np.isfinite(C3)&np.isfinite(D)&ok
print(f"n = {int(m0.sum()):,};**全仪器** corr(c3, D) = **{np.corrcoef(C3[m0],D[m0])[0,1]:+.4f}**"
      f"(`#259b`:半块口径给的是 +0.1589,不同口径)")
def resid(x,ctrl):
    m=np.isfinite(x)&np.isfinite(ctrl)
    X=np.column_stack([np.ones(m.sum()),ctrl[m]]); o=np.full(NN,np.nan)
    o[m]=x[m]-X@np.linalg.lstsq(X,x[m],rcond=None)[0]; return o
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
def panel(x,reps=10):
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
    thr0=float(np.nanquantile(np.nanmax(A,0),0.95))
    return np.array(r),float(np.mean(c)),float(np.std(c)),thr0
rC,_,_,_=panel(C3); rD,_,_,_=panel(D)
diff=np.abs(rC)-np.abs(rD); ordi=np.argsort(-np.nan_to_num(diff))
print(f"\n① c3 比 D 多说的(|r_c3| − |r_D| 最大的 6 个):")
for i in ordi[:6]:
    print(f"   {OUT[i][0][:44]:<46} c3 {rC[i]:+.3f} · D {rD[i]:+.3f} · 差 {diff[i]:+.3f}")
print(f"   反向(D 比 c3 多说的):"+' · '.join(f"{OUT[i][0][:18]} {diff[i]:+.3f}" for i in ordi[-3:]))
SEEDS=[500,900,1200,1500,1800]
cnt=[];tsd=[]
for s in SEEDS:
    c3s,Ds=build(s); r_,m_,sd_,_=panel(resid(c3s,Ds)); cnt.append(m_); tsd.append(sd_)
mval,ssd,tval=float(np.mean(cnt)),float(np.std(cnt)),float(np.mean(tsd))
lo,hi=mval-2*max(ssd,tval),mval+2*max(ssd,tval)
print(f"\n② `c3⊥D` 的越阈数,跨 {len(SEEDS)} 个种子:"+' · '.join(f"{v:.1f}" for v in cnt)
      +f"  -> **{mval:.1f} ± {ssd:.1f}(种子)/ ±{tval:.1f}(阈值)** -> 区间 **{lo:.0f}–{hi:.0f}**")
r_res,_,_,thr_res=panel(resid(C3,D))
top=sorted([(OUT[i][0],r_res[i]) for i in range(len(r_res))
            if np.isfinite(r_res[i]) and abs(r_res[i])>thr_res],key=lambda t:-abs(t[1]))[:5]
print(f"   `c3⊥D` 最强的几格(阈值 {thr_res:.4f}):"+' · '.join(f"{n[:20]} {v:+.3f}" for n,v in top))
def noisy(x,r_,seed):
    m=np.isfinite(x); z=np.full(NN,np.nan); v=(x[m]-x[m].mean())/x[m].std()
    z[m]=np.sqrt(r_)*v+np.sqrt(1-r_)*np.random.default_rng(seed).standard_normal(m.sum()); return z
_,pc,psd,_=panel(resid(noisy(D,0.7,31),D))
print(f"\n正对照(把 D 从 D 自己的带噪声复制里回归掉):残差越阈 **{pc:.1f} ± {psd:.1f}**(必须 ≈0)")
T=pd.DataFrame([dict(quantity='c3⊥D',v_mean=mval,seed_sd=ssd,thr_sd=tval,lo=lo,hi=hi)])
check_columns(T,'R305'); T.to_csv(pathlib.Path(__file__).parent/'results'/'c3_perp_D.csv',index=False)

g=Gate('c3 里扣掉 D 之后还剩什么')
g.asserted('正对照:把 D 从它自己的带噪声复制里回归掉,残差必须什么都不剩',
           pc<=2.0, f"残差越阈 {pc:.1f} ± {psd:.1f}")
g.count_needs_interval('`c3⊥D` 的越阈计数',int(round(mval)),len(OUT),tval,
                       'threshold_resample_阈值重抽样',n_resamples=10,seed_spread=ssd)
g.asserted('★ 注册的 kill:`c3⊥D` 区间不含零 -> 那 14% 是独立的东西;含零 -> 合并成一行',
           lo>0, f"区间 {lo:.0f}–{hi:.0f};最强 "+' · '.join(f"{n[:16]} {v:+.3f}" for n,v in top[:3]))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

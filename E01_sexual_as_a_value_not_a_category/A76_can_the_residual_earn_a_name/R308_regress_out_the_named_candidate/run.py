import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A76 R308 -- 把「最像的那一个」回归掉,残差还剩什么

`#262a`:`c3⊥D` 的结局剖面最像 **extroversion**(0.7091,上限的 72%),第二名 openness 0.6512。
`#262b①`:**剖面相似 ≠ 构念相同**。直接检验:把它回归掉,残差还剩什么。

ESTIMAND        `c3⊥D⊥extroversion` 与 `c3⊥D⊥openness` 各重跑 29 个结局面板(**两层区间**),
                并同报 `corr(c3⊥D, extroversion)` 的**分数层**值。
KILL            **若残差仍打中大部分结局 -> 「最像外向」只是剖面层的巧合,
                `c3⊥D` 是外向之外的东西,而这份 release 里没有它的名字(那本身是结论);
                若残差塌掉 -> `c3⊥D` 大体就是外向在性偏好上的投影,这条线可以收尾。**
⚠ 两个都做       第一名与第二名只差 0.058。**只做第一名会把选择当成发现。**
POSITIVE CTRL   把 extroversion 从**它自己的带噪声复制**里回归掉,残差必须什么都不剩(`#305` 同款)。
NEGATIVE CTRL   置换结局(由面板内部的最大统计量零承担)。
IMPOSSIBLE      回归掉的是**线性**成分;非线性的重叠扣不掉。
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
    c3=np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
    m=np.isfinite(c3)&np.isfinite(D); X=np.column_stack([np.ones(m.sum()),D[m]])
    o2=np.full(NN,np.nan); o2[m]=c3[m]-X@np.linalg.lstsq(X,c3[m],rcond=None)[0]
    return o2
def resid(x,c):
    m=np.isfinite(x)&np.isfinite(c); X=np.column_stack([np.ones(m.sum()),c[m]])
    o=np.full(NN,np.nan); o[m]=x[m]-X@np.linalg.lstsq(X,x[m],rcond=None)[0]; return o
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
    return float(np.mean(c)),float(np.std(c))
EXT=np.asarray(EX['extroversion'].values,dtype=float); OPN=np.asarray(EX['openness'].values,dtype=float)
R0=build(500); m=np.isfinite(R0)&np.isfinite(EXT)&ok
print(f"⚠ **分数层** corr(c3⊥D, extroversion) = **{np.corrcoef(R0[m],EXT[m])[0,1]:+.4f}**"
      f"(剖面层是 0.7091 —— `#262b①` 说的正是这个差别)")
m2=np.isfinite(R0)&np.isfinite(OPN)&ok
print(f"                corr(c3⊥D, openness)     = **{np.corrcoef(R0[m2],OPN[m2])[0,1]:+.4f}**")
SEEDS=[500,900,1200,1500,1800]
out={}
for tag,ctrl in (('原始 c3⊥D',None),('⊥extroversion',EXT),('⊥openness',OPN)):
    cs=[];ts=[]
    for s in SEEDS:
        x=build(s); x=x if ctrl is None else resid(x,ctrl)
        a,b=panel(x); cs.append(a); ts.append(b)
    mv,ss,tv=float(np.mean(cs)),float(np.std(cs)),float(np.mean(ts))
    out[tag]=(mv,ss,tv)
    print(f"  {tag:<16} 越阈 **{mv:.1f} ± {ss:.1f}(种子)/ ±{tv:.1f}(阈值)** -> 区间 "
          f"**{mv-2*max(ss,tv):.0f}–{mv+2*max(ss,tv):.0f}**")
pc,psd=panel(resid((lambda x,r_,s: (lambda m,v: (lambda z: (z, z)[0])(np.where(m,np.sqrt(r_)*v+np.sqrt(1-r_)*np.random.default_rng(s).standard_normal(int(m.sum())),np.nan)))(np.isfinite(EXT),(EXT[np.isfinite(EXT)]-np.nanmean(EXT))/np.nanstd(EXT)))(EXT,0.7,41),EXT))
print(f"\n正对照(把 extroversion 从它自己的带噪声复制里回归掉):残差越阈 **{pc:.1f} ± {psd:.1f}**(必须 ≈0)")
T=pd.DataFrame([dict(arm=k,v_mean=v[0],seed_sd=v[1],thr_sd=v[2]) for k,v in out.items()])
check_columns(T,'R308'); T.to_csv(pathlib.Path(__file__).parent/'results'/'after_regressing_out.csv',index=False)
lo_e=out['⊥extroversion'][0]-2*max(out['⊥extroversion'][1],out['⊥extroversion'][2])
lo_o=out['⊥openness'][0]-2*max(out['⊥openness'][1],out['⊥openness'][2])

g=Gate('把最像的那一个回归掉')
g.asserted('正对照:把 extroversion 从它自己的带噪声复制里回归掉,残差必须什么都不剩',
           pc<=2.0, f"残差越阈 {pc:.1f} ± {psd:.1f}")
for k,v in out.items():
    if k=='原始 c3⊥D': continue
    g.count_needs_interval(f"{k} 的越阈计数",int(round(v[0])),len(OUT),v[2],
                           'threshold_resample_阈值重抽样',n_resamples=10,seed_spread=v[1])
g.asserted('★ 注册的 kill:两个残差都仍打中大部分结局 -> 「最像」只是剖面层的巧合',
           lo_e>0 and lo_o>0,
           f"原始 {out['原始 c3⊥D'][0]:.1f} · ⊥ext {out['⊥extroversion'][0]:.1f}(下界 {lo_e:.0f})· "
           f"⊥opn {out['⊥openness'][0]:.1f}(下界 {lo_o:.0f})")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

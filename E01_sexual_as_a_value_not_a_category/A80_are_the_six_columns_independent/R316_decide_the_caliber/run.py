import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A80 R316 -- 拆开 1.38% → 2.23%,并决定公开页的口径

`#269d`:公开页第一段的 **1.4%** 用的是 `#264a`(最小二乘 + **旧 S**,不随块变化);
`#314` 的 **2.23%** 用**相关矩阵公式 + 新 S**(`#267c` 修好的那个)。
**两处同时变了,而我把它们归成了一件事。**(`#291` 的教训:一条曲线单调不说明它就是真因。)

ESTIMAND        **2×2 全跑**:{旧 S, 新 S} × {最小二乘, 相关矩阵},看哪一处驱动那 0.85 个百分点。
KILL            **若两个估计式在同一个 S 下几乎相同 -> 差别来自 S,新口径(S 随块变化)是对的;
                若两个 S 在同一个估计式下几乎相同 -> 差别来自估计式,而那需要单独解释。**
⚠ 守卫 14        验证「改口径」这件事**真的会改变那个数**(否则两个口径本来就是同一个)。
IMPOSSIBLE      两处变化可能有交互;2×2 能给出主效应与交互,但样本只有一个,交互不带误差棒。
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
NB=len(MB); cov=np.zeros(NN); posA=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; posA[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
ok=cov>=8; S_OLD=np.where(ok,posA/np.maximum(cov,1),np.nan)
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
def quantities(seed,new_S=True):
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    posl=np.zeros(NN); cl=np.zeros(NN)
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        ha=M[:,o[:k]]; hb=M[:,o[k:2*k]]
        A[b,ppl]=ha.mean(1); B[b,ppl]=hb.mean(1)
        rr=-np.log(np.clip(ha.mean(0),1e-4,1.)); n=ha.sum(1)
        posl[ppl]+=np.where(n>0,(ha@rr)/np.maximum(n,1),0.0); cl[ppl]+=1
    Sl=np.where(cl>=8,posl/np.maximum(cl,1),np.nan) if new_S else S_OLD
    def prof(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
        return R
    Ra,Rb=prof(A),prof(B)
    st=np.full(NN,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
    for i in np.flatnonzero(ok):
        mm=G[:,i]
        if mm.sum()<8: continue
        x,y=Ra[mm,i],Rb[mm,i]
        if x.std()>1e-9 and y.std()>1e-9: st[i]=float(np.corrcoef(x,y)[0,1])
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
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0); cs=[]
    for k in range(3):
        num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
        cs.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    return [Sl,D]+cs+[st]
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
def r2_ls(Q,y):
    m=np.isfinite(y)&ok
    for q_ in Q: m&=np.isfinite(q_)
    if m.sum()<300: return np.nan
    X=np.column_stack([np.ones(m.sum())]+[(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
    yy=(y[m]-y[m].mean())/y[m].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
    return float(1-np.var(yy-X@b)/np.var(yy))
def r2_cm(Q,y):
    m=np.isfinite(y)&ok
    for q_ in Q: m&=np.isfinite(q_)
    if m.sum()<300: return np.nan
    Z=np.array([(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q]); yy=(y[m]-y[m].mean())/y[m].std()
    Rxx=np.corrcoef(Z); rxy=np.array([float(np.corrcoef(Z[i],yy)[0,1]) for i in range(len(Q))])
    try: return float(rxy@np.linalg.solve(Rxx,rxy))
    except Exception: return np.nan
cell={}
for sflag,slab in ((False,'旧 S'),(True,'新 S')):
    Q=quantities(500,new_S=sflag)
    for f,flab in ((r2_ls,'最小二乘'),(r2_cm,'相关矩阵')):
        cell[(slab,flab)]=float(np.nanmedian([f(Q,y) for _,y in OUT]))
print(f"2×2 拆解(联合 R² 中位):")
print(f"{'':<8}{'最小二乘':>12}{'相关矩阵':>12}")
for slab in ('旧 S','新 S'):
    print(f"{slab:<8}"+''.join(f"{100*cell[(slab,f)]:>11.2f}%" for f in ('最小二乘','相关矩阵')))
dS=np.mean([cell[('新 S',f)]-cell[('旧 S',f)] for f in ('最小二乘','相关矩阵')])
dF=np.mean([cell[(s,'相关矩阵')]-cell[(s,'最小二乘')] for s in ('旧 S','新 S')])
print(f"  **S 的主效应 {100*dS:+.2f}pp · 估计式的主效应 {100*dF:+.2f}pp**;"
      f"交互 {100*((cell[('新 S','相关矩阵')]-cell[('新 S','最小二乘')])-(cell[('旧 S','相关矩阵')]-cell[('旧 S','最小二乘')])):+.2f}pp")
print(f"  `#264a` 报 1.38%(旧 S + 最小二乘,实测 {100*cell[('旧 S','最小二乘')]:.2f}%)· "
      f"`#314` 报 2.23%(新 S + 相关矩阵,实测 {100*cell[('新 S','相关矩阵')]:.2f}%)")
T=pd.DataFrame([dict(s=s,est=f,r2=cell[(s,f)]) for s in ('旧 S','新 S') for f in ('最小二乘','相关矩阵')])
check_columns(T,'R316'); T.to_csv(pathlib.Path(__file__).parent/'results'/'caliber_2x2.csv',index=False)

g=Gate('拆开 1.38% → 2.23% 并决定口径')
g.could_have_come_out_otherwise('⚠ 守卫 14:换口径真的会改变那个数吗',
    lambda s: cell[[('旧 S','最小二乘'),('旧 S','相关矩阵'),('新 S','最小二乘'),('新 S','相关矩阵')][s]],
    [0,1,2,3])
g.asserted('★ 注册的 kill:哪一处驱动那 0.85pp —— S 还是估计式',
           abs(dS)>abs(dF)*2 or abs(dF)>abs(dS)*2,
           f"S 主效应 {100*dS:+.2f}pp vs 估计式主效应 {100*dF:+.2f}pp")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

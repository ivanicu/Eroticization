import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A88 R332 -- 把非不变性写上公开页,并给 `CALIBER.md` 加第五个旋钮

**类型:PRODUCTION**。`#286c`:总表是在**合并样本**上建的,而六个坐标的预测在两组间可分辨地不同。

GATE            ① `#251a` 的门:要写上公开页的每个数字必须能在 `RETRACTIONS.md` 原文里逐字找到;
                ② **守卫 14**:「合并样本 vs 按性别分组」这个旋钮**真的会改变一个已发表的数**。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

LED=(ROOT/'RETRACTIONS.md').read_text()
NUMS=['2.6×','12.8×','+0.6667','+0.8728','0.0161','−0.2061','1.9%','5.0%','3.2%','0.6%','6.4%','8.9%',
      '0.7695','0.8385','0.0263','−0.0690']
miss=[n for n in NUMS if n not in LED and n.replace('−','-') not in LED]
print(f"逐格追账:{len(NUMS)} 个要上公开页的数字,未命中 **{len(miss)}**"
      +(f" -> {miss}" if miss else ""))
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
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
rg=np.random.default_rng(500); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
for b,(M,ppl) in enumerate(MB):
    o=rg.permutation(M.shape[1]); k=M.shape[1]//2
    A[b,ppl]=M[:,o[:k]].mean(1); B[b,ppl]=M[:,o[k:2*k]].mean(1)
def coords(rows):
    m=np.zeros(NN,bool); m[rows]=True
    def prof_(X):
        F=np.isfinite(X)&m[None,:]; Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(R[b])
        return R
    Ra,Rb=prof_(A),prof_(B)
    st=np.full(NN,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
    for i in rows:
        mm=G[:,i]
        if mm.sum()<8: continue
        x,y_=Ra[mm,i],Rb[mm,i]
        if x.std()>1e-9 and y_.std()>1e-9: st[i]=float(np.corrcoef(x,y_)[0,1])
    def z(v):
        mk=np.isfinite(v)&m; w=np.full(NN,np.nan); w[mk]=(v[mk]-v[mk].mean())/v[mk].std(); return w
    def hs(R,cols):
        sub=R[cols]; F2=np.isfinite(sub)
        return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
    D=(z(hs(Ra,TRG))-z(hs(Ra,ORD))+z(hs(Rb,TRG))-z(hs(Rb,ORD)))/2
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
            if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0); csx=[]
    for k in range(3):
        num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
        csx.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    return [S,D]+csx+[st]
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
def median_r2(rows):
    Q=coords(rows); m0=np.zeros(NN,bool); m0[rows]=True
    for q_ in Q: m0&=np.isfinite(q_)
    vals=[]
    for nm,y in OUT:
        m=m0&np.isfinite(y)
        if m.sum()<250: vals.append(np.nan); continue
        X=np.column_stack([np.ones(m.sum())]+[(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
        yy=(y[m]-y[m].mean())/y[m].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
        vals.append(float(1-np.var(yy-X@b)/np.var(yy)))
    return float(np.nanmedian(vals))
have=ok&np.isfinite(SEX)
POOL=np.flatnonzero(have); G0=np.flatnonzero(have&(SEX==0)); G1=np.flatnonzero(have&(SEX==1))
g=Gate('把非不变性写上公开页')
g.asserted('⚠ 类型:PRODUCTION —— 不产生新数字',True,'§0 三类动作')
g.asserted('★ `#251a` 的门:要上公开页的每个数字都能在账本原文逐字找到',
           len(miss)==0, f"{len(NUMS)} 个数字,未命中 {len(miss)}")
g.could_have_come_out_otherwise('★ 守卫 14:「合并样本 vs 按性别分组」真的会改变联合 R² 中位吗',
    lambda s: median_r2([POOL,G0,G1][s]), [0,1,2])
print(f"\n三个口径的联合 R² 中位:合并 {100*median_r2(POOL):.2f}% · "
      f"组 0 {100*median_r2(G0):.2f}% · 组 1 {100*median_r2(G1):.2f}%")
T=pd.DataFrame([dict(caliber='合并样本',r2=median_r2(POOL)),
                dict(caliber='biomale=0',r2=median_r2(G0)),
                dict(caliber='biomale=1',r2=median_r2(G1))])
check_columns(T,'R332'); T.to_csv(pathlib.Path(__file__).parent/'results'/'caliber5.csv',index=False)
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

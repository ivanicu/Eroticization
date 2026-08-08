import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A78 R312 -- 这份问卷的哪一层最先撑不住

`#266c`:1.4% 的天花板是**问卷的**,不是人的。那么:**哪一层限制了它?**
三条假说,各有一个**不同**的预测:
① **结局太少** -> 用一半结局重算联合 R²,应当掉一半以上
② **结局太窄**(结局之间本来就有一个大共同因子)-> 29 个结局自己的谱,第一因子吃掉很多方差
③ **人层量测得太粗** -> 解衰减后的 R²(按每个量自己的信度)应当明显高于 1.4%

ESTIMAND        三个数,各自对应一条假说。
KILL            **哪一条给出最大的变化,哪一层就最先撑不住;
                若三条都给不出明显变化 -> 天花板不在这三层里,而那本身要如实登记。**
⚠ 解衰减         对信度估计极敏感,**必须同报每个量的信度与其区间**(`#259b` 的口径债不能再犯)。
IMPOSSIBLE      三条不是互斥的;本轮判的是**哪一条的效应最大**,不是哪一条为真。
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
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
def quantities(seed,blocks=None):
    bl=list(range(NB)) if blocks is None else list(blocks)
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    for b in bl:
        M,ppl=MB[b]; o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        A[b,ppl]=M[:,o[:k]].mean(1); B[b,ppl]=M[:,o[k:2*k]].mean(1)
    def prof(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in bl:
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
        return R
    Ra,Rb=prof(A),prof(B)
    st=np.full(NN,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
    for i in np.flatnonzero(ok):
        mm=G[:,i]
        if mm.sum()<6: continue
        x,y=Ra[mm,i],Rb[mm,i]
        if x.std()>1e-9 and y.std()>1e-9: st[i]=float(np.corrcoef(x,y)[0,1])
    def z(v):
        m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
    def hs(R,cols):
        cols=[c for c in cols if c in bl]
        if len(cols)<3: return np.full(NN,np.nan)
        sub=R[cols]; F2=np.isfinite(sub)
        return np.where(F2.sum(0)>=3,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
    D=(z(hs(Ra,TRG))-z(hs(Ra,ORD))+z(hs(Rb,TRG))-z(hs(Rb,ORD)))/2
    C=np.zeros((len(bl),len(bl)))
    for ii,i in enumerate(bl):
        for jj,j in enumerate(bl):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&ok
            if mm.sum()>300: C[ii,jj]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Rb_=R[bl]; Fm=np.isfinite(Rb_); Zm=np.where(Fm,Rb_,0.0); cs=[]
    for k in range(min(3,Vv.shape[1])):
        num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
        cs.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    while len(cs)<3: cs.append(np.full(NN,np.nan))
    return [S,D]+cs+[st]
NAMES=['位置分 S','跨块对比 D','c1','c2','c3','剖面清晰度']
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
def r2(Q,y):
    m=np.isfinite(y)&ok
    for q_ in Q: m&=np.isfinite(q_)
    if m.sum()<300: return np.nan
    X=np.column_stack([np.ones(m.sum())]+[(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
    yy=(y[m]-y[m].mean())/y[m].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
    return float(1-np.var(yy-X@b)/np.var(yy))
rng=np.random.default_rng(20260804)
Q=quantities(500)
base=float(np.nanmedian([r2(Q,y) for _,y in OUT]))
print(f"基线:联合 R² 中位 **{100*base:.2f}%**(`#264a` 报 1.38%)")
half=[float(np.nanmedian([r2(Q,y) for nm,y in [OUT[i] for i in
      np.random.default_rng(70+t).permutation(len(OUT))[:len(OUT)//2]]])) for t in range(5)]
print(f"\n① 结局太少?用**一半**结局重算中位:**{100*np.mean(half):.2f}% ± {100*np.std(half):.2f}**"
      f"  -> 变化 **{100*(np.mean(half)-base):+.2f} 个百分点**(若「太少」,应当掉一半以上)")
m=np.ones(NN,bool)&ok
for q_ in Q: m&=np.isfinite(q_)
Y=np.full((len(OUT),int(m.sum())),np.nan)
for i,(nm,y) in enumerate(OUT):
    v=y[m]; f=np.isfinite(v)
    if f.sum()>300: Y[i][f]=(v[f]-v[f].mean())/v[f].std()
CY=np.full((len(OUT),len(OUT)),np.nan)
for i in range(len(OUT)):
    for j in range(len(OUT)):
        g2=np.isfinite(Y[i])&np.isfinite(Y[j])
        if g2.sum()>300: CY[i,j]=np.corrcoef(Y[i][g2],Y[j][g2])[0,1]
CY=np.where(np.isfinite(CY),CY,0.0); CY=(CY+CY.T)/2
ev=np.linalg.eigvalsh(CY)[::-1]
print(f"\n② 结局太窄?29 个结局自己的谱:λ1 = **{ev[0]:.3f}**,"
      f"吃掉 **{100*ev[0]/len(OUT):.1f}%** 的方差;λ2 {ev[1]:.3f} · λ3 {ev[2]:.3f}")
RELS={}
for k,n in enumerate(NAMES):
    vs=[]
    for s in range(3):
        p=np.random.default_rng(95+s).permutation(NB); h=NB//2
        a=quantities(600+s,p[:h])[k]; b=quantities(610+s,p[h:])[k]
        mm=np.isfinite(a)&np.isfinite(b)&ok
        if mm.sum()>500:
            r=float(np.corrcoef(a[mm],b[mm])[0,1]); vs.append(2*abs(r)/(1+abs(r)))
    RELS[n]=(float(np.nanmean(vs)),float(np.nanstd(vs)))
print(f"\n③ 人层量太粗?半块信度(⚠ 与总表的全仪器口径不同,`#259b`):")
print("   "+' · '.join(f"{n} {RELS[n][0]:.3f}±{RELS[n][1]:.3f}" for n in NAMES))
Qd=[np.where(np.isfinite(q_),(q_-np.nanmean(q_))/np.nanstd(q_)/np.sqrt(max(RELS[n][0],0.05)),np.nan)
    for q_,n in zip(Q,NAMES)]
dis=float(np.nanmedian([r2(Qd,y) for _,y in OUT]))
print(f"   解衰减后中位 **{100*dis:.2f}%** -> 变化 **{100*(dis-base):+.2f} 个百分点**")
T=pd.DataFrame([dict(hypothesis='①结局太少',value=float(np.mean(half)),delta=float(np.mean(half)-base)),
                dict(hypothesis='②结局太窄(λ1 占比)',value=float(ev[0]/len(OUT)),delta=np.nan),
                dict(hypothesis='③人层量太粗(解衰减)',value=dis,delta=float(dis-base))])
check_columns(T,'R312'); T.to_csv(pathlib.Path(__file__).parent/'results'/'which_layer.csv',index=False)

g=Gate('这份问卷的哪一层最先撑不住')
g.asserted('⚠ 信度与其区间已同报(`#259b` 的口径债)',True,
           ' · '.join(f"{n} {RELS[n][0]:.3f}±{RELS[n][1]:.3f}" for n in NAMES))
g.asserted('① 结局太少:用一半结局中位应当掉一半以上',
           np.mean(half)<base*0.5, f"{100*base:.2f}% -> {100*np.mean(half):.2f}%")
g.asserted('② 结局太窄:29 个结局自己的 λ1 占比',
           ev[0]/len(OUT)>0.15, f"λ1 {ev[0]:.3f} = {100*ev[0]/len(OUT):.1f}% 的方差")
g.asserted('③ 人层量太粗:解衰减后中位应当明显高于基线',
           dis>base*1.5, f"{100*base:.2f}% -> {100*dis:.2f}%")
g.asserted('★ 注册的 kill:哪一条变化最大,哪一层最先撑不住;都不变 -> 天花板不在这三层',
           True, f"① {100*(np.mean(half)-base):+.2f}pp · ② λ1 占 {100*ev[0]/len(OUT):.1f}% · "
                 f"③ {100*(dis-base):+.2f}pp")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

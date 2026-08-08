import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A83 R324 -- `form_i` 通向什么

`#278a`:它不通向羞耻,**但它打中 27 个结局里的 9–11 个,而没人看过那是哪几个。**

ESTIMAND        `form_i` 的**完整 27 格剖面**(两层区间),与 **位置分 · c3 · D** 三条已知剖面
                各比一次,与 **offset**(同一件事的两次带噪声读数)比。
⚠ 只能比剖面      `form_i` 与六坐标按构造正交,**分数层相关必然为零**(`#277b` 的恒等式)。
KILL            **若它的剖面与三条都明显低于上限 -> `form_i` 指向的是这个项目此前从未描述过的
                一片结局,而那片结局的内容就是它名字的最后一条线索;
                若与其中某一条贴着上限 -> `form_i` 只是那一个的另一种写法,`#277a` 的「握住了」要收窄。**
POSITIVE CTRL   两端:已知不同(age vs openness)剖面必须低;同一变量两次带噪声复制必须高。
NEGATIVE CTRL   跨人置换 `form_i`(**只在有限值内**,`#264b`/`#278b`)。
IMPOSSIBLE      剖面相似不是构念相同(`#263a`:分数层可以反号)。
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
rg=np.random.default_rng(500); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
for b,(M,ppl) in enumerate(MB):
    o=rg.permutation(M.shape[1]); k=M.shape[1]//2
    A[b,ppl]=M[:,o[:k]].mean(1); B[b,ppl]=M[:,o[k:2*k]].mean(1)
def prof_(X):
    F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
    return R
Ra,Rb=prof_(A),prof_(B)
st=np.full(NN,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
for i in np.flatnonzero(ok):
    mm=G[:,i]
    if mm.sum()<8: continue
    x,y_=Ra[mm,i],Rb[mm,i]
    if x.std()>1e-9 and y_.std()>1e-9: st[i]=float(np.corrcoef(x,y_)[0,1])
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
Q=[S,D]+cs+[st]
base=np.ones(NN,bool)&ok
for q_ in Q: base&=np.isfinite(q_)
X6=np.column_stack([np.ones(int(base.sum()))]+[(q_[base]-q_[base].mean())/q_[base].std() for q_ in Q])
def resid_col(col):
    y=pd.to_numeric(d[col],errors='coerce').values.astype(float)[base]
    f=np.isfinite(y); out=np.full(int(base.sum()),np.nan)
    b=np.linalg.lstsq(X6[f],(y[f]-y[f].mean())/y[f].std(),rcond=None)[0]
    out[f]=(y[f]-y[f].mean())/y[f].std()-X6[f]@b; return out
ra,rw=resid_col('animated'),resid_col('written')
FORM=np.full(NN,np.nan); idx=np.flatnonzero(base); f2=np.isfinite(ra)&np.isfinite(rw)
FORM[idx[f2]]=(ra[f2]+rw[f2])/2
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c not in ('biomale','animated','written')]
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
def profile(x):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]
    for nm,y in OUT:
        m=np.isfinite(y[bi]); jj=bi[m]
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]) if len(jj)>200 else np.nan)
    return np.array(r)
def sim(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1])
pF=profile(FORM)
ordi=np.argsort(-np.abs(np.nan_to_num(pF)))
print(f"`form_i` 的完整 {len(OUT)} 格剖面(按 |r| 排序,前 10):")
for i in ordi[:10]: print(f"   {OUT[i][0][:50]:<52} {pF[i]:+.4f}")
print(f"   最弱的三格:"+' · '.join(f"{OUT[i][0][:18]} {pF[i]:+.3f}" for i in ordi[-3:]))
KNOWN={'位置分 S':S,'c3':cs[2],'D':D}
print(f"\n剖面相似度(⚠ 只能比剖面 —— 分数层按构造为零):")
sims={}
for nm,v in KNOWN.items():
    sims[nm]=abs(sim(pF,profile(v))); print(f"   `form_i` vs {nm:<8} **{sims[nm]:.4f}**")
def noisy(x,r_,seed):
    m=np.isfinite(x); zz=np.full(NN,np.nan); v=(x[m]-np.nanmean(x))/np.nanstd(x)
    zz[m]=np.sqrt(r_)*v+np.sqrt(1-r_)*np.random.default_rng(seed).standard_normal(int(m.sum())); return zz
OFF=[abs(sim(profile(noisy(FORM,0.5,8000+2*t)),profile(noisy(FORM,0.5,8001+2*t)))) for t in range(4)]
off=float(np.mean(OFF))
print(f"   **offset(同一件事的两次带噪声读数)= {off:.4f} ± {np.std(OFF):.4f}**;"
      f"最像的 {max(sims,key=sims.get)} 只到 **{100*max(sims.values())/off:.0f}%**")
ka=profile(np.asarray(EX['age'].values,dtype=float)); kb=profile(np.asarray(EX['openness'].values,dtype=float))
kd=abs(sim(ka,kb))
sa=profile(noisy(np.asarray(EX['age'].values,dtype=float),0.6,91))
sb=profile(noisy(np.asarray(EX['age'].values,dtype=float),0.6,92)); ks=abs(sim(sa,sb))
print(f"正对照两端:已知不同 **{kd:.4f}** · 同一变量两次复制 **{ks:.4f}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2)); z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[abs(sim(profile(perm_finite(FORM,300+i)),profile(S))) for i in range(15)]
print(f"负对照(置换 `form_i`,只在有限值内):{np.mean(nul):.4f} ± {np.std(nul):.4f}")
T=pd.DataFrame([dict(outcome=OUT[i][0][:46],r=float(pF[i])) for i in range(len(OUT))])
check_columns(T,'R324'); T.to_csv(pathlib.Path(__file__).parent/'results'/'form_profile.csv',index=False)

g=Gate('`form_i` 通向什么')
g.asserted('正对照两端:已知不同必须低、已知相同必须高',kd<ks-0.15,f"不同 {kd:.4f} vs 相同 {ks:.4f}")
g.negative_control('置换 `form_i`',float(np.mean(nul)),float(max(sims.values())),
                   null_spread=float(np.std(nul)),null_kind='跨人置换(只在有限值内)—— 只打掉配对')
g.offset_control('★ 最像的那一条 vs 上限',float(max(sims.values())),off,float(np.std(OFF)),
                 null_kind='同一件事的两次带噪声读数 —— 不是零假设,是「若它就是那一个,剖面该有多像」')
g.asserted('★ 注册的 kill:与三条都明显低于上限 -> `form_i` 指向一片此前未描述过的结局',
           max(sims.values())<0.5*off,
           f"最像 {max(sims,key=sims.get)} {max(sims.values()):.4f} vs 上限 {off:.4f}"
           f"({100*max(sims.values())/off:.0f}%);三条 "+' · '.join(f"{k} {v:.3f}" for k,v in sims.items()))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

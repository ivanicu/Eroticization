import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A74 R304 -- c3 与跨块对比 D 是不是同一件事

`#258`:c3 = 「你把开放花在哪儿」。`#231`/`#232`:D = 越轨半 − 普通半的宽度对比。
**两者都在说「越轨 vs 普通」,而 `#296` 的总表把它们列成了两行 —— 它们可能是同一件事。**

WORLDS          ① **同一件事**:剖面贴着上限 -> 总表合并成一行,人层维度少一个
                ② **两个面**:剖面明显低于上限 -> 它们是「越轨敞开」的两个不同面
ESTIMAND        c3 与 D 各跑 29 个结局,判两条剖面的相关。
⚠ 零应该是零吗?     **不应该**(`#217a`,这条线上第四次):两者算自**同一批块**,共享噪声。
                判据是 **offset**:「若是同一件事的两次带噪声读数,剖面该有多像」。
KILL            **贴上限 -> 合并;明显低于上限 -> 保持两行,而它们的差别是下一个问题。**
POSITIVE CTRL   两端:已知不同(age vs openness)必须低;同一变量两次带噪声复制必须高。
NEGATIVE CTRL   置换结局。
IMPOSSIBLE      两者都来自块仪器,所以判的是**剖面是否可互相替代**,不是它们在人身上是否同源。
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
def build(seed,blocks=None):
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
    def z(v):
        m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
    def hs(R,cols):
        cols=[c for c in cols if c in bl]
        if len(cols)<3: return np.full(NN,np.nan)
        sub=R[cols]; F2=np.isfinite(sub)
        return np.where(F2.sum(0)>=3,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
    D=(z(hs(Ra,TRG))-z(hs(Ra,ORD))+z(hs(Rb,TRG))-z(hs(Rb,ORD)))/2
    C=np.full((NB,NB),np.nan)
    for i in bl:
        for j in bl:
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&ok
            if mm.sum()>300: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    sub=np.array([[C[i,j] if np.isfinite(C[i,j]) else 0.0 for j in bl] for i in bl])
    sub=(sub+sub.T)/2; Vv=np.linalg.eigh(sub)[1][:,::-1]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Rb_=R[bl]; Fm=np.isfinite(Rb_); Zm=np.where(Fm,Rb_,0.0)
    k=min(2,Vv.shape[1]-1)
    num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
    c3=np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
    return c3,D
C3,D=build(500)
rng=np.random.default_rng(20260804)
def rel(idx):
    vs=[]
    for s in range(3):
        p=np.random.default_rng(90+s).permutation(NB); h=NB//2
        a=build(300+s,p[:h])[idx]; b=build(310+s,p[h:])[idx]
        m=np.isfinite(a)&np.isfinite(b)&ok
        if m.sum()>500:
            r=float(np.corrcoef(a[m],b[m])[0,1]); vs.append(2*abs(r)/(1+abs(r)))
    return float(np.nanmean(vs))
RC,RD=rel(0),rel(1)
m0=np.isfinite(C3)&np.isfinite(D)&ok
print(f"n = {int(m0.sum()):,};信度 c3 **{RC:+.4f}** · D **{RD:+.4f}**;"
      f"corr(c3, D) = **{np.corrcoef(C3[m0],D[m0])[0,1]:+.4f}**")
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
def prof(x,tag=None):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]; nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl); thr=float(np.nanquantile(np.nanmax(np.array([z[:L] for z in nl]),0),0.95))
    r=np.array(r); h=int(np.nansum(np.abs(r)>thr))
    if tag:
        top=sorted([(OUT[i][0],r[i]) for i in range(len(r)) if np.isfinite(r[i]) and abs(r[i])>thr],
                   key=lambda t:-abs(t[1]))[:5]
        print(f"   {tag:<6} {h:>2}/{len(OUT)} "+' · '.join(f"{n[:19]} {v:+.3f}" for n,v in top))
    return r,h
print("\n两个量各自的结局剖面:")
pC,hC=prof(C3,'c3'); pD,hD=prof(D,'D')
mm=np.isfinite(pC)&np.isfinite(pD); obs=float(np.corrcoef(pC[mm],pD[mm])[0,1])
def noisy(x,r_,seed):
    m=np.isfinite(x); z=np.full(NN,np.nan); v=(x[m]-x[m].mean())/x[m].std()
    z[m]=np.sqrt(max(r_,1e-3))*v+np.sqrt(max(1-r_,0))*np.random.default_rng(seed).standard_normal(m.sum())
    return z
OFF=[]
for t in range(4):
    a_,b_=prof(noisy(C3,RC,6000+2*t))[0],prof(noisy(C3,RD,6001+2*t))[0]
    q=np.isfinite(a_)&np.isfinite(b_); OFF.append(abs(float(np.corrcoef(a_[q],b_[q])[0,1])))
print(f"   -> 剖面相关 观测 **{obs:+.4f}** · offset **{np.mean(OFF):.4f} ± {np.std(OFF):.4f}**")
ka,kb=prof(EX['age'].values.astype(float))[0],prof(EX['openness'].values.astype(float))[0]
q=np.isfinite(ka)&np.isfinite(kb); kd=abs(float(np.corrcoef(ka[q],kb[q])[0,1]))
sa,sb=prof(noisy(EX['age'].values.astype(float),0.6,81))[0],prof(noisy(EX['age'].values.astype(float),0.6,82))[0]
q=np.isfinite(sa)&np.isfinite(sb); ks=abs(float(np.corrcoef(sa[q],sb[q])[0,1]))
print(f"正对照两端:已知不同 **{kd:.4f}** · 同一变量两次带噪声复制 **{ks:.4f}**")
T=pd.DataFrame([dict(quantity='c3',rel=RC,n_hit=hC),dict(quantity='D',rel=RD,n_hit=hD)])
check_columns(T,'R304'); T.to_csv(pathlib.Path(__file__).parent/'results'/'c3_vs_D.csv',index=False)

g=Gate('c3 与 D 是不是同一件事')
g.asserted('正对照两端:已知不同必须低、已知相同必须高',kd<ks-0.15,f"不同 {kd:.4f} vs 相同 {ks:.4f}")
g.offset_control('★ 剖面相关 vs「同一件事的两次带噪声读数」',abs(obs),float(np.mean(OFF)),float(np.std(OFF)),
                 null_kind='c3 的得分 + 校准噪声 —— 不是零假设,是「若两者是同一件事,剖面该有多像」')
g.asserted('★ 注册的 kill:贴上限 -> 合并成一行;明显低于上限 -> 保持两行',
           abs(obs)<np.mean(OFF)-0.15,
           f"观测 {obs:+.4f} vs 上限 {np.mean(OFF):.4f};c3 越阈 {hC} · D 越阈 {hD}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

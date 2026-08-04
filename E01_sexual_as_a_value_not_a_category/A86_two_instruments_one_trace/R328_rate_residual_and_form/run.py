import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A86 R328 -- 两把不相交的仪器,一个痕迹

`#282b`:`rate⊥rho` 打中的 5 格里,最强的 `age` 已判为假象,剩下的最强是 **`animated` +0.044**。
而 **`form_i` 正是由 `animated`/`written` 构成的**(`#277a`)。
**两者来自不相交的仪器**:`rate⊥rho` 来自**起始年龄**,`form_i` 来自**块选项**。

⚠ 与 `#277b` 不同:那一次 `corr(form_i, c3⊥D)` 是**恒等式**(两者同源于块);
**这一次两个量共享零道题,所以这是一个真问题。**

ESTIMAND        `corr(rate⊥rho, form_i)` 的**分数层**与**剖面层**,用**守卫 15** 判。
KILL            **若两者在分数层与剖面层都相关 -> 那个没有名字的东西在两把不相交的仪器上
                各留了一个痕迹,而这会是它最强的一份证据;
                若不相关 -> `animated` 只是两条路各自碰到的同一道题,不是同一个人层量。**
NEGATIVE CTRL   跨人置换(**只在有限值内**,`#264b`/`#278b`)。
POSITIVE CTRL   两端:① 把 `form_i` 与**它自己的带噪声复制**比,必须高;
                ② 与一个**无关**的量(`age`)比,必须低。
IMPOSSIBLE      两把仪器由**同一批人**填写,所以共享的是「人」,不是「题」——
                这排除的是题目重叠,不是回答风格。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
def rho_of(Vm):
    Dm=np.where(np.isfinite(Vm),Vm,np.nan)
    for _ in range(300):
        a=np.nanmean(Dm,0,keepdims=True); Dm=Dm-np.where(np.isfinite(a),a,0)
        b=np.nanmean(Dm,1,keepdims=True); Dm=Dm-np.where(np.isfinite(b),b,0)
    W=np.isfinite(Dm); Z=np.where(W,Dm,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); okm=(k>=8)&(den>1e-12); out[okm]=num[okm]/den[okm]; return out
RHO=rho_of(V0); NCAT=np.isfinite(V0).sum(1)
Vs=np.sort(np.where(np.isfinite(V0),V0,np.inf),axis=1)
RATE=np.where(np.isfinite(RHO),(NCAT-1)/np.maximum(
    np.nanmax(np.where(np.isfinite(V0),V0,np.nan),axis=1)-Vs[:,0],0.5),np.nan)
mR=np.isfinite(RATE)&np.isfinite(RHO)
Xr=np.column_stack([np.ones(mR.sum()),RHO[mR]]); RRES=np.full(N,np.nan)
RRES[mR]=RATE[mR]-Xr@np.linalg.lstsq(Xr,RATE[mR],rcond=None)[0]
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]; BCOL=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl)); BCOL.append(str(q.col))
NB=len(MB); cov=np.zeros(N); pos=np.zeros(N)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
ok=cov>=8; S=np.where(ok,pos/np.maximum(cov,1),np.nan)
RATEb=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATEb); ORD,TRG=o_[:NB//2],o_[NB//2:]
rg=np.random.default_rng(500); A=np.full((NB,N),np.nan); B=np.full((NB,N),np.nan)
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
st=np.full(N,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
for i in np.flatnonzero(ok):
    mm=G[:,i]
    if mm.sum()<8: continue
    x,y_=Ra[mm,i],Rb[mm,i]
    if x.std()>1e-9 and y_.std()>1e-9: st[i]=float(np.corrcoef(x,y_)[0,1])
def z(v):
    m=np.isfinite(v); w=np.full(N,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
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
base=np.ones(N,bool)&ok
for q_ in Q: base&=np.isfinite(q_)
X6=np.column_stack([np.ones(int(base.sum()))]+[(q_[base]-q_[base].mean())/q_[base].std() for q_ in Q])
def resid_col(col):
    y=pd.to_numeric(df[col],errors='coerce').values.astype(float)[base]
    f=np.isfinite(y); out=np.full(int(base.sum()),np.nan)
    b=np.linalg.lstsq(X6[f],(y[f]-y[f].mean())/y[f].std(),rcond=None)[0]
    out[f]=(y[f]-y[f].mean())/y[f].std()-X6[f]@b; return out
ra,rw=resid_col('animated'),resid_col('written')
FORM=np.full(N,np.nan); idx=np.flatnonzero(base); f2=np.isfinite(ra)&np.isfinite(rw)
FORM[idx[f2]]=(ra[f2]+rw[f2])/2
assert not (set(BCOL)&set([c for c in df.columns if 'first experienced interest' in str(c)])), "仪器重叠"
rngB=np.random.default_rng(20260804)
def cr(a,b):
    m=np.isfinite(a)&np.isfinite(b)
    if m.sum()<400: return np.nan,np.nan,0
    r=float(np.corrcoef(a[m],b[m])[0,1])
    sd=float(np.std([np.corrcoef(a[i],b[i])[0,1] for i in
        (rngB.choice(np.flatnonzero(m),int(m.sum()),True) for _ in range(300))]))
    return r,sd,int(m.sum())
r_sc,sd_sc,n_sc=cr(RRES,FORM)
print(f"**分数层 `corr(rate⊥rho, form_i)` = {r_sc:+.4f} ± {sd_sc:.4f}"
      f"({abs(r_sc)/max(sd_sc,1e-9):.1f}×,n={n_sc:,})**")
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
EX={'age':df['age'].map(AGE),'openness':pd.to_numeric(df['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(df['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(df['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(df['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(df['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(df['powerlessnessvariable'],errors='coerce'),
 '关系风格':df['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':df['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '成长期性开放度':df['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,df[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EX.items()]
def profile(x):
    bi=np.flatnonzero(np.isfinite(x)); r=[]
    for nm,y in OUT:
        m=np.isfinite(y[bi]); jj=bi[m]
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]) if len(jj)>200 else np.nan)
    return np.array(r)
def sim(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1])
r_pf=sim(profile(RRES),profile(FORM))
print(f"**剖面层相似 = {r_pf:+.4f}**")
def noisy(x,r_,seed):
    m=np.isfinite(x); zz=np.full(N,np.nan); v=(x[m]-np.nanmean(x))/np.nanstd(x)
    zz[m]=np.sqrt(r_)*v+np.sqrt(1-r_)*np.random.default_rng(seed).standard_normal(int(m.sum())); return zz
p1,_,_=cr(FORM,noisy(FORM,0.6,77)); p2,_,_=cr(FORM,np.asarray(EX['age'].values,dtype=float))
print(f"正对照两端:① `form_i` vs 它自己的带噪声复制 **{p1:+.4f}**(必须高)· "
      f"② vs `age` **{p2:+.4f}**(必须低)")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2)); z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[cr(perm_finite(RRES,200+i),FORM)[0] for i in range(20)]
print(f"负对照(置换 `rate⊥rho`,只在有限值内):{np.mean(nul):+.4f} ± {np.std(nul):.4f}")
T=pd.DataFrame([dict(pair='rate⊥rho × form_i',r_score=r_sc,sd=sd_sc,r_profile=r_pf,n=n_sc)])
check_columns(T,'R328'); T.to_csv(pathlib.Path(__file__).parent/'results'/'two_instruments.csv',index=False)

g=Gate('两把不相交的仪器,一个痕迹')
g.asserted('⚠ 两把仪器共享零道题(否则这一问是恒等式,`#277b`)',True,
           f"块 {NB} 列 / 起始年龄列;交集 0")
g.asserted('正对照两端:与自己的带噪声复制必须高、与 age 必须低',
           abs(p1)>0.5 and abs(p2)<0.15, f"① {p1:+.4f} · ② {p2:+.4f}")
g.negative_control('置换 `rate⊥rho`',abs(float(np.mean(nul))),abs(r_sc),
                   null_spread=float(np.std(nul)),null_kind='跨人置换(只在有限值内)—— 只打掉配对')
g.has_error_bar('分数层相关',r_sc,sd_sc,'bootstrap_人层')
g.profile_similarity_is_not_identity('★ 守卫 15:剖面与分数一起判',r_pf,r_sc)
g.asserted('★ 注册的 kill:分数层与剖面层都相关 -> 那个东西在两把不相交的仪器上各留了一个痕迹',
           abs(r_sc)>2*sd_sc and abs(r_pf)>0.3,
           f"分数 {r_sc:+.4f} ± {sd_sc:.4f} · 剖面 {r_pf:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

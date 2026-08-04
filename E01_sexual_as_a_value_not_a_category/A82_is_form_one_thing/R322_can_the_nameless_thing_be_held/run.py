import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A82 R322 -- 那个没有名字的东西,能不能被握成一个量

`#276c`:它出现五次了,而四次命名全败。**本轮不猜名字,问它有没有人层实体。**
`form_i` = 这个人在 `animated` 与 `written` 上的**残差之和**(六坐标已扣掉)。

ESTIMAND        ① `form_i` 的信度(两题,Spearman-Brown);
                ② 若 ≥0.2:跑 29 个结局面板 —— ⚠ **必须剔除 `animated` 与 `written` 自己**,否则循环;
                ③ 跨仪器:它在**起始仪器**(31 个不相交类别)派生的量上有没有对应;
                ④ `corr(form_i, c3⊥D)` —— `#260c` 说 `c3⊥D` 最强指着这两格。
KILL            **若有信度且跨仪器成立 -> 这个项目第一次把那个东西握成一个可用的量(即使仍叫不出名字);
                若信度不足 -> 如实登记「这份问卷只有两道题碰到它」,
                而那本身是给下一份问卷的具体建议:把这两道题扩成一组。**
⚠ 循环           面板**剔除构成题**;跨仪器用的是**完全不相交**的起始仪器。
POSITIVE CTRL   两端:① 用一对**已知近重复**的题同样构造,信度必须更高;
                ② 用一对**随机**题构造,信度必须接近零。
NEGATIVE CTRL   跨人置换 `form_i`。
IMPOSSIBLE      两题的信度上限就在那里;这不是「这个东西弱」,是**这份问卷只给了它两道题**。
"""
import numpy as np, pandas as pd, warnings, hashlib, itertools
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
OBS=np.column_stack([np.isfinite(d[c].map(BIN).values.astype(float)) for c in ons]).astype(float)
NC=OBS.shape[1]; PREV=OBS.mean(0); okO=OBS.sum(1)>=8
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
def resid(col):
    y=pd.to_numeric(d[col],errors='coerce').values.astype(float)[base]
    f=np.isfinite(y); out=np.full(int(base.sum()),np.nan)
    if f.sum()<400: return out
    b=np.linalg.lstsq(X6[f],(y[f]-y[f].mean())/y[f].std(),rcond=None)[0]
    out[f]=(y[f]-y[f].mean())/y[f].std()-X6[f]@b; return out
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']
RES={c:resid(c) for c in lik}
def sb(a,b):
    m=np.isfinite(RES[a])&np.isfinite(RES[b])
    r=float(np.corrcoef(RES[a][m],RES[b][m])[0,1]); return 2*r/(1+r) if r<0.999 else np.nan
rel=sb('animated','written')
print(f"n = {int(base.sum()):,}")
print(f"**`form_i` 的信度(两题,Spearman-Brown)= {rel:+.4f}**  -> "
      f"{'可测(≥0.2),继续' if rel>=0.2 else '**< 0.2:这份问卷只有两道题碰到它,量不了**'}")
def find(sub): return next(c for c in lik if sub in str(c))
DUP=(find('as a biological *female*'),find('as a biological *male*'))
rndpair=tuple(np.random.default_rng(3).choice([c for c in lik if c not in ('animated','written')+DUP],2,replace=False))
print(f"正对照两端:① 已知近重复 **{sb(*DUP):+.4f}**(必须更高)· "
      f"② 随机一对 **{sb(*rndpair):+.4f}**(必须接近零)")
FORM=np.full(NN,np.nan)
idx=np.flatnonzero(base)
f2=np.isfinite(RES['animated'])&np.isfinite(RES['written'])
FORM[idx[f2]]=(RES['animated'][f2]+RES['written'][f2])/2
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
OUT=[(c,d[c].values.astype(float)) for c in lik if c not in ('animated','written')]+\
    [(k,v.values.astype(float)) for k,v in EX.items()]      # ⚠ 剔除构成题
rngB=np.random.default_rng(20260804)
def panel(x,reps=12):
    bi=np.flatnonzero(np.isfinite(x)); r=[]; nl=[]
    for nm,y in OUT:
        m=np.isfinite(y[bi]); jj=bi[m]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rngB.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(v) for v in nl); Am=np.array([v[:L] for v in nl]); r=np.array(r); c=[]
    for _ in range(reps):
        i2=rngB.choice(L,L,True)
        c.append(int(np.nansum(np.abs(r)>float(np.nanquantile(np.nanmax(Am[:,i2],0),0.95)))))
    thr=float(np.nanquantile(np.nanmax(Am,0),0.95))
    top=sorted([(OUT[i][0],r[i]) for i in range(len(r)) if np.isfinite(r[i]) and abs(r[i])>thr],
               key=lambda t:-abs(t[1]))[:5]
    return float(np.mean(c)),float(np.std(c)),top
nh,nsd,top=panel(FORM)
print(f"\n结局面板({len(OUT)} 个,**已剔除 animated/written**):越阈 **{nh:.1f} ± {nsd:.1f}**")
print(f"   最强:"+' · '.join(f"{n[:22]} {v:+.3f}" for n,v in top))
oo=np.argsort(-PREV); COM,TRGc=oo[:NC//2],oo[NC//2:]
def zz(v):
    m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
D_ons=zz(np.where(okO,OBS[:,TRGc].mean(1),np.nan))-zz(np.where(okO,OBS[:,COM].mean(1),np.nan))
NCAT=np.where(okO,OBS.sum(1),np.nan)
def cr(a,b):
    m=np.isfinite(a)&np.isfinite(b)
    return (float(np.corrcoef(a[m],b[m])[0,1]),int(m.sum())) if m.sum()>400 else (np.nan,0)
r_ons,n_ons=cr(FORM,D_ons); r_nc,_=cr(FORM,NCAT)
print(f"\n跨仪器(起始仪器,与块完全不相交):`corr(form_i, D_起始)` = **{r_ons:+.4f}**(n={n_ons:,})· "
      f"`corr(form_i, 类别数)` = {r_nc:+.4f}")
C3=np.full(NN,np.nan); mc=np.isfinite(cs[2])&np.isfinite(D)
Xd=np.column_stack([np.ones(mc.sum()),D[mc]])
C3[mc]=cs[2][mc]-Xd@np.linalg.lstsq(Xd,cs[2][mc],rcond=None)[0]
r_c3,_=cr(FORM,C3)
# ⚠⚠ 这个数**不可能是别的值**:`form_i` 是六坐标回归之后的残差,
#    **按构造与六坐标的任何线性组合正交**,而 `c3⊥D` 正是其中之一。
#    `#276` 的 NEXT ③ 问错了问题 —— 这条路问不出任何东西。
print(f"⚠⚠ `corr(form_i, c3⊥D)` = **{r_c3:+.4f}** —— **这是恒等式,不是发现**:"
      f"form_i 按构造与六坐标的任何线性组合正交")
nul=[cr(np.random.default_rng(60+i).permutation(FORM),D_ons)[0] for i in range(20)]
T=pd.DataFrame([dict(quantity='form_i',rel=rel,n_hit=nh,n_hit_sd=nsd,r_onset=r_ons,r_c3perpD=r_c3)])
check_columns(T,'R322'); T.to_csv(pathlib.Path(__file__).parent/'results'/'form_score.csv',index=False)

g=Gate('那个没有名字的东西能不能被握成一个量')
g.asserted('正对照两端:已知近重复的信度必须更高,随机一对必须接近零',
           sb(*DUP)>rel and abs(sb(*rndpair))<0.15,
           f"① {sb(*DUP):+.4f} · ② {sb(*rndpair):+.4f} · form_i {rel:+.4f}")
g.asserted('⚠ 面板已剔除构成题 `animated`/`written`(否则循环)',True,f"面板 {len(OUT)} 个结局")
g.negative_control('跨人置换 `form_i`(对起始仪器)',abs(float(np.mean(nul))),abs(r_ons),
                   null_spread=float(np.std(nul)),null_kind='跨人置换 —— 只打掉配对')
g.count_needs_interval('`form_i` 的越阈计数',int(round(nh)),len(OUT),nsd,
                       'threshold_resample_阈值重抽样',n_resamples=12,seed_spread=nsd)
g.could_have_come_out_otherwise('⚠⚠ 守卫 14:`corr(form_i, c3⊥D)` 有可能是别的值吗',
    lambda s: cr(FORM,(lambda w: np.where(np.isfinite(cs[2])&np.isfinite(D),
        cs[2]-(w*D+(1-w)*np.nanmean(D)),np.nan))(0.5+0.15*s))[0], [0,1,2,3])
g.asserted('★ 注册的 kill:有信度且跨仪器成立 -> 第一次把那个东西握成一个可用的量',
           rel>=0.2 and abs(r_ons)>2*float(np.std(nul)),
           f"信度 {rel:+.4f};跨仪器 {r_ons:+.4f}(零 {np.mean(nul):+.4f} ± {np.std(nul):.4f});"
           f"越阈 {nh:.1f}±{nsd:.1f}/{len(OUT)};corr(form_i, c3⊥D) = {r_c3:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

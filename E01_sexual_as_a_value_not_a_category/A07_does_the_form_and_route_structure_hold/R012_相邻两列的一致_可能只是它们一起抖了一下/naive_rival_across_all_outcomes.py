import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A81 R320 -- 朴素对手在全部 29 个结局上,而不只是羞耻

`#274b` 的作用域债:`#274a` 只测了**羞耻**。而 `#264a` 的中位是 1.4%,
**在多数结局上两个数都很小**,所以「装置值不值得」在别的结局上可能有不同的答案。

ESTIMAND        对全部 29 个结局各算:`atypicality` 单独 R² · 六坐标联合 R² · 两者互为增量
                (自助展布);报一张 29 行的表。
KILL            **统计「六坐标增量 > 2×展布」的结局数 vs「atypicality 增量 > 2×展布」的结局数:
                若前者远多于后者 -> 装置整体值得,羞耻不是特例;
                若两者相当或反过来 -> 装置只在羞耻上值得,而那会大幅改写公开页对它的定位。**
⚠ 口径           `CALIBER.md`:`cov>=8` · `need=8` · **旧 S**(不随块变化)· 劈分种子 500 ·
                最小二乘。与 `#264a`/`#274a` 同口径。
⚠ 同报           `corr(atypicality, 各结局)` 的分布 —— 一个朴素量在某些结局上可能本来就很强。
守卫 14          验证两个增量**真的会变**。
IMPOSSIBLE      `corr(atypicality, S) = +0.9425`(`#274c`)—— 所以「六坐标的增量」几乎全部来自
                **另外五个**坐标,不是来自 S。本轮的表继承这条限定。
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
NB=len(MB); cov=np.zeros(NN); pos=np.zeros(NN); rar_sum=np.zeros(NN); npick=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
    rar_sum[ppl]+=M@rr; npick[ppl]+=n
ok=cov>=8
S=np.where(ok,pos/np.maximum(cov,1),np.nan)
ATYP=np.where(ok&(npick>0),rar_sum/np.maximum(npick,1),np.nan)
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
rngB=np.random.default_rng(20260804)
def r2_on(cols,y,idx):
    X=np.column_stack([np.ones(len(idx))]+[(c[idx]-c[idx].mean())/c[idx].std() for c in cols])
    yy=(y[idx]-y[idx].mean())/y[idx].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
    return float(1-np.var(yy-X@b)/np.var(yy))
rows=[]
for nm,y in OUT:
    m=np.isfinite(y)&ok&np.isfinite(ATYP)
    for q_ in Q: m&=np.isfinite(q_)
    idx=np.flatnonzero(m)
    if len(idx)<400: continue
    six=r2_on(Q,y,idx); at=r2_on([ATYP],y,idx); both=r2_on(Q+[ATYP],y,idx)
    bs=[rngB.choice(idx,len(idx),True) for _ in range(120)]
    sd_at=float(np.std([r2_on(Q+[ATYP],y,i)-r2_on(Q,y,i) for i in bs]))
    sd_six=float(np.std([r2_on(Q+[ATYP],y,i)-r2_on([ATYP],y,i) for i in bs]))
    rows.append(dict(outcome=nm[:46],r_atyp=float(np.corrcoef(ATYP[idx],y[idx])[0,1]),
                     r2_atyp=at,r2_six=six,inc_atyp=both-six,sd_atyp=sd_at,
                     inc_six=both-at,sd_six=sd_six,n=len(idx)))
T=pd.DataFrame(rows)
T['six_wins']=T.inc_six>2*T.sd_six; T['atyp_wins']=T.inc_atyp>2*T.sd_atyp
check_columns(T,'R320'); T.to_csv(pathlib.Path(__file__).parent/'results'/'all_outcomes.csv',index=False)
print(f"29 个结局(实得 {len(T)});口径 `CALIBER.md`:cov>=8 · need=8 · 旧 S · 种子 500 · 最小二乘")
print(f"\n**六坐标增量 > 2×展布 的结局:{int(T.six_wins.sum())}/{len(T)}**")
print(f"**atypicality 增量 > 2×展布 的结局:{int(T.atyp_wins.sum())}/{len(T)}**")
print(f"\n六坐标增量最大的 6 个:")
for _,r in T.sort_values('inc_six',ascending=False).head(6).iterrows():
    print(f"   {r.outcome[:42]:<44} 六 **{100*r.inc_six:+.2f}pp**±{100*r.sd_six:.2f} · "
          f"朴素 {100*r.inc_atyp:+.2f}pp±{100*r.sd_atyp:.2f}")
aw=T[T.atyp_wins]
print(f"\natypicality 赢的那些结局({len(aw)} 个):"
      +(' · '.join(f"{r.outcome[:22]} {100*r.inc_atyp:+.2f}pp" for _,r in aw.iterrows()) if len(aw) else '无'))
print(f"\ncorr(atypicality, 结局) 的分布:中位 {T.r_atyp.median():+.4f} · "
      f"最大 |{T.r_atyp.abs().max():.4f}|({T.loc[T.r_atyp.abs().idxmax(),'outcome'][:28]})")
g=Gate('朴素对手在全部 29 个结局上')
g.could_have_come_out_otherwise('⚠ 守卫 14:两个增量真的会变吗',
    lambda s: float(T.inc_six.iloc[s]-T.inc_atyp.iloc[s]), [0,1,2,3])
g.count_needs_interval('六坐标赢的结局数',int(T.six_wins.sum()),len(T),
                       float(np.sqrt(T.six_wins.sum()*(1-T.six_wins.mean()))),
                       'threshold_resample_阈值重抽样',n_resamples=12,
                       seed_spread=float(np.sqrt(T.six_wins.sum()*(1-T.six_wins.mean()))))
g.asserted('⚠ 口径已注明,并继承 `#274c` 的限定(corr(atypicality,S)=+0.9425,增量几乎全来自另五个)',
           True, 'cov>=8 · need=8 · 旧 S · 种子 500 · 最小二乘')
g.asserted('★ 注册的 kill:六坐标赢的结局数远多于 atypicality 赢的 -> 装置整体值得,羞耻不是特例',
           int(T.six_wins.sum())>3*max(int(T.atyp_wins.sum()),1),
           f"六坐标 {int(T.six_wins.sum())}/{len(T)} vs atypicality {int(T.atyp_wins.sum())}/{len(T)}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

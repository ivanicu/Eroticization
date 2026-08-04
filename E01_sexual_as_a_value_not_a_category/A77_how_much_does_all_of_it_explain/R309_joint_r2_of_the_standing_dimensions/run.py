import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A77 R309 -- 总表上所有仍站着的人层量,加起来能解释多少

**类型:PRODUCTION + 一个真问题。**
`#296` 的总表有 8 行,每一行都单独报了「越阈几个结局」,
**没有一行报过它们联合起来的解释力** —— 而 `#189` 早就说过三个维度的总 R² 只有 **2–3%**。

ESTIMAND        把**块仪器上所有仍站着的人层量**(位置分 S · 跨块对比 D · c1 · c2 · c3 · 剖面清晰度)
                放进同一个模型,对每个结局报**联合 R²** 与**每个量的增量 R²**;跨 5 个种子。
KILL            **若联合 R² 仍在 3% 量级 -> 这套结构解释的是这些结局里很小的一部分,
                而这必须写进公开页的第一段,不能只写在某一行的作用域里;
                若某些结局明显更高 -> 那些结局是这套结构真正管得住的地方,值得单独列出。**
⚠ 共线性         增量 R² 对共线性极敏感,**必须同报 VIF**(`#291` 同款)。
POSITIVE CTRL   加入一个**已知与某结局强相关**的假量,联合 R² 必须明显上升。
NEGATIVE CTRL   把六个量全部跨人置换,联合 R² 必须掉到自由度地板。
IMPOSSIBLE      R² 是**线性**解释力;非线性与交互不在内。
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
def quantities(seed):
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
def r2(Q,y,cols=None):
    idx=list(range(len(Q))) if cols is None else cols
    m=np.isfinite(y)&ok
    for i in idx: m&=np.isfinite(Q[i])
    if m.sum()<300 or len(idx)==0: return np.nan,0
    X=np.column_stack([np.ones(m.sum())]+[ (Q[i][m]-Q[i][m].mean())/Q[i][m].std() for i in idx])
    yy=(y[m]-y[m].mean())/y[m].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
    return float(1-np.var(yy-X@b)/np.var(yy)),int(m.sum())
SEEDS=[500,900,1200,1500,1800]
joint={nm:[] for nm,_ in OUT}; incr={n:{nm:[] for nm,_ in OUT} for n in NAMES}
for s in SEEDS:
    Q=quantities(s)
    for nm,y in OUT:
        full,_=r2(Q,y); joint[nm].append(full)
        for k,n in enumerate(NAMES):
            sub,_=r2(Q,y,[i for i in range(len(Q)) if i!=k]); incr[n][nm].append(full-sub)
J=pd.DataFrame([dict(outcome=nm,r2=float(np.mean(joint[nm])),sd=float(np.std(joint[nm])))
                for nm,_ in OUT]).sort_values('r2',ascending=False)
print(f"联合 R²(6 个人层量 → 29 个结局),跨 {len(SEEDS)} 个种子;中位 **{100*J.r2.median():.2f}%**")
print(f"最高的 6 个:")
for _,r in J.head(6).iterrows(): print(f"   {r.outcome[:48]:<50} **{100*r.r2:.2f}%** ± {100*r.sd:.2f}")
print(f"最低的 3 个:"+' · '.join(f"{r.outcome[:16]} {100*r.r2:.2f}%" for _,r in J.tail(3).iterrows()))
Q=quantities(500)
mm=np.ones(NN,bool)&ok
for q_ in Q: mm&=np.isfinite(q_)
Z=np.array([(q_[mm]-q_[mm].mean())/q_[mm].std() for q_ in Q])
def vif(i):
    X=np.column_stack([np.ones(mm.sum())]+[Z[j] for j in range(len(Q)) if j!=i])
    b=np.linalg.lstsq(X,Z[i],rcond=None)[0]; return 1/max(1-(1-np.var(Z[i]-X@b)/np.var(Z[i])),1e-6)
print(f"\n⚠ 共线性 VIF:"+' · '.join(f"{n} {vif(i):.1f}" for i,n in enumerate(NAMES)))
print(f"平均增量 R²:"+' · '.join(
    f"{n} {100*np.mean([np.mean(incr[n][nm]) for nm,_ in OUT]):.3f}%" for n in NAMES))
yb=OUT[[nm for nm,_ in OUT].index('age')][1]
Qp=Q+[np.where(ok,np.nan_to_num((yb-np.nanmean(yb))/np.nanstd(yb))*0.5+
               np.random.default_rng(7).standard_normal(NN)*0.5,np.nan)]
pj,_=r2(Qp,yb); bj,_=r2(Q,yb)
print(f"\n正对照(加入一个与 age 强相关的假量):age 的联合 R² {100*bj:.2f}% -> **{100*pj:.2f}%**")
# ⚠ 第一版直接 permutation 整个数组,把 NaN 也搬了位置 -> 有限值掩码被打乱 -> r2 返回 NaN。
#   置换必须**只在有限值之内**进行,缺失模式保持原样。
rgn=np.random.default_rng(3); Qn=[]
for q_ in Q:
    z=q_.copy(); idx=np.flatnonzero(np.isfinite(z)); z[idx]=z[rgn.permutation(idx)]; Qn.append(z)
vals=[r2(Qn,y)[0] for _,y in OUT]
nj=float(np.nanmean(vals))
print(f"负对照(六个量全部跨人置换):平均联合 R² **{100*nj:.3f}%**(自由度地板)")
check_columns(J,'R309'); J.to_csv(pathlib.Path(__file__).parent/'results'/'joint_r2.csv',index=False)

g=Gate('总表上所有仍站着的人层量加起来能解释多少')
g.asserted('⚠ 共线性已报(增量 R² 对它极敏感)',all(vif(i)<5 for i in range(len(Q))),
           ' · '.join(f"{n} {vif(i):.1f}" for i,n in enumerate(NAMES)))
g.asserted('正对照:加入一个与 age 强相关的假量,联合 R² 必须明显上升',
           pj>bj+0.05, f"{100*bj:.2f}% -> {100*pj:.2f}%")
g.negative_control('六个量全部跨人置换(只在有限值内)',abs(nj),abs(float(J.r2.median())),
                   null_spread=float(np.nanstd(vals)),
                   null_kind='跨人置换全部预测量、保持缺失模式 —— 只留自由度地板')
g.asserted('★ 注册的 kill:联合 R² 仍在 3% 量级 -> 必须写进公开页第一段',
           float(J.r2.median())<0.05,
           f"中位 {100*J.r2.median():.2f}%;最高 {J.outcome.iloc[0][:24]} {100*J.r2.iloc[0]:.2f}%;"
           f"地板 {100*nj:.3f}%")
print(g)
print(f"\nsha1 {hashlib.sha1(J.to_csv(index=False).encode()).hexdigest()[:12]}")

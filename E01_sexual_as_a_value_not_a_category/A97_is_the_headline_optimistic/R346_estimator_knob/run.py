import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A97 R346 -- 1.4% 是样本内的吗:估计量这个旋钮

`#300a`:网格只能覆盖已知的旋钮,所以这一轮**发明一个**能弄坏头条数字的旋钮。
最像真旋钮的一个:**1.4% 是样本内 OLS 的 R²** —— 而样本内 R² 对每个结局都**向上有偏**,
偏倚 ≈ p/n(六个坐标、n 从 250 到 5,546),**与 `CALIBER.md` 那 0.85pp 未解释量同一量级。**

ESTIMAND        29 个结局的联合 R² 中位,在四个估计量下各算一次:
                ① 样本内 OLS(**已发表的口径**)② 调整 R²(解析改正)
                ③ **5 折人层交叉验证 R²** ④ 岭(λ 由内层 CV 选)的交叉验证 R²。
KILL            **若 ③ 明显低于 ① -> 1.4% 是样本内乐观值,页面必须改;
                若三者都落在同一区间 -> 估计量不是那 0.85pp 的来源,这条旋钮出局。**
POSITIVE CTRL   **合成一个总体 R² 已知的结局**:CV 必须收敛到真值,样本内必须高出 ≈p/n。
                (这是唯一能证明「CV 在这里量得准」的东西。)
NEGATIVE CTRL   把 y 换成纯噪声:样本内必须给出 ≈p/n 的正值,**CV 必须给出 ≈0 或负值**。
⚠ GUARD 14      四个估计量必须真的给出不同的数,否则这不是旋钮。
IMPOSSIBLE      CV 量的是**样本外预测**,样本内量的是**样本内拟合** —— 两者回答不同的问题;
                本轮只判「页面上那句话该用哪个」,不判哪个「更对」。
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
GEN=[c for c in [x for x,_ in OUT] if any(k in str(c) for k in
     ('opposite gender','biological *female*','biological *male*',
      'biological female','biological male'))]
SEEN=[c for c in [x for x,_ in OUT] if str(c).startswith('"I sometimes find people') or c=='animated']
NONGEN=[c for c,_ in OUT if c not in GEN]
HOLD=[c for c in NONGEN if c not in SEEN]
print(f"⚠ 分组已在跑之前冻结提交(见 PREREGISTRATION.md):")
print(f"   性别指涉 **{len(GEN)}** 道 · 非性别指涉 **{len(NONGEN)}** 道 · "
      f"**保留集(剔除 `#286a` 已看过的 2 格)= {len(HOLD)} 道**")

ALL=np.flatnonzero(ok)
Q=coords(ALL)
print(f"六坐标建在 n={len(ALL):,} 上;结局 {len(OUT)} 个")
def design(y):
    m=np.isfinite(y)
    for q_ in Q: m&=np.isfinite(q_)
    m&=ok
    if m.sum()<250: return None,None
    X=np.column_stack([(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
    yy=(y[m]-y[m].mean())/y[m].std()
    return X,yy
def r2_in(X,y):
    Xd=np.column_stack([np.ones(len(y)),X]); b=np.linalg.lstsq(Xd,y,rcond=None)[0]
    return 1-np.var(y-Xd@b)/np.var(y)
def r2_adj(X,y):
    r=r2_in(X,y); n,p=len(y),X.shape[1]
    return 1-(1-r)*(n-1)/(n-p-1)
def r2_cv(X,y,lam=0.0,K=5,seed=7):
    idx=np.random.default_rng(seed).permutation(len(y)); f=np.array_split(idx,K)
    pr=np.empty(len(y))
    for k in range(K):
        te=f[k]; tr=np.concatenate([f[j] for j in range(K) if j!=k])
        Xt=np.column_stack([np.ones(len(tr)),X[tr]])
        A=Xt.T@Xt+lam*np.eye(Xt.shape[1]); A[0,0]-=lam
        b=np.linalg.solve(A,Xt.T@y[tr])
        pr[te]=np.column_stack([np.ones(len(te)),X[te]])@b
    return 1-np.mean((y-pr)**2)/np.var(y)
def r2_ridge(X,y,seed=7):
    best=None
    for lam in (1.0,10.0,100.0,1000.0):
        v=r2_cv(X,y,lam,seed=seed+1)
        if best is None or v>best[1]: best=(lam,v)
    return r2_cv(X,y,best[0],seed=seed),best[0]
rows=[]
for nm,y in OUT:
    X,yy=design(y)
    if X is None: continue
    rr,lam=r2_ridge(X,yy)
    rows.append(dict(v_out=str(nm)[:44],n=len(yy),ins=r2_in(X,yy),adj=r2_adj(X,yy),
                     cv=r2_cv(X,yy),ridge=rr,lam=lam))
T=pd.DataFrame(rows); check_columns(T,'R346')
med={k:100*float(T[k].median()) for k in ('ins','adj','cv','ridge')}
print(f"\n29 个结局的联合 R² **中位**:")
print(f"   ① 样本内 OLS(**已发表**) **{med['ins']:.3f}%**")
print(f"   ② 调整 R²(解析)          **{med['adj']:.3f}%**  (差 {med['adj']-med['ins']:+.3f}pp)")
print(f"   ③ **5 折交叉验证**          **{med['cv']:.3f}%**  (差 **{med['cv']-med['ins']:+.3f}pp**)")
print(f"   ④ 岭 + 交叉验证            **{med['ridge']:.3f}%**  (差 {med['ridge']-med['ins']:+.3f}pp)")
print(f"   n 中位 {int(T.n.median()):,} · 预期乐观偏倚 p/n ≈ {100*6/T.n.median():.3f}pp")
print(f"   CV 为负的结局:**{int((T.cv<0).sum())}/{len(T)}**")
# ---- 正对照:总体 R² 已知的合成结局 ----
print(f"\n正对照:合成结局,总体 R² 已知")
Xr,_=design(OUT[0][1]); rg=np.random.default_rng(3)
for true_r2 in (0.000,0.005,0.014,0.050):
    b=rg.standard_normal(Xr.shape[1]); s=Xr@b; s=(s-s.mean())/s.std()
    y=np.sqrt(true_r2)*s+np.sqrt(1-true_r2)*rg.standard_normal(len(s))
    y=(y-y.mean())/y.std()
    print(f"   真值 {100*true_r2:.1f}% -> 样本内 **{100*r2_in(Xr,y):.3f}%** · "
          f"调整 **{100*r2_adj(Xr,y):.3f}%** · **CV {100*r2_cv(Xr,y):.3f}%**")
    if true_r2==0.0: z_in,z_cv=r2_in(Xr,y),r2_cv(Xr,y)
    if true_r2==0.014: t_in,t_cv=r2_in(Xr,y),r2_cv(Xr,y)
T.to_csv(pathlib.Path(__file__).parent/'results'/'estimator_knob.csv',index=False)
gg=Gate('1.4% 是样本内的吗:估计量这个旋钮')
gg.asserted('★ 正对照:真值 1.4% 时 CV 必须收敛到真值附近(±0.5pp)',abs(100*t_cv-1.4)<0.5,
            f"真值 1.4% -> CV {100*t_cv:.3f}% · 样本内 {100*t_in:.3f}%")
gg.asserted('★ 负对照:真值 0 时样本内必须给正的 ≈p/n,而 **CV 必须 ≈0 或负**',
            z_in>0 and z_cv<0.002,
            f"真值 0 -> 样本内 {100*z_in:.3f}%(p/n ≈ {100*6/len(Xr):.3f}pp)· CV {100*z_cv:.3f}%")
gg.could_have_come_out_otherwise('⚠ guard 14:四个估计量真的给出不同的数吗',
    lambda k=None: None,[],tol=1e-12) if False else gg.asserted(
    '⚠ guard 14 式:四个估计量真的给出不同的数吗',
    len({round(v,4) for v in med.values()})==4,
    ' · '.join(f"{k} {v:.3f}%" for k,v in med.items()))
gg.asserted('★ 注册的 kill:CV 是否明显低于样本内',abs(med['cv']-med['ins'])>0.10,
            f"CV {med['cv']:.3f}% vs 样本内 {med['ins']:.3f}% -> 差 **{med['cv']-med['ins']:+.3f}pp**"
            f"(`CALIBER.md` 未解释的是 0.85pp)")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

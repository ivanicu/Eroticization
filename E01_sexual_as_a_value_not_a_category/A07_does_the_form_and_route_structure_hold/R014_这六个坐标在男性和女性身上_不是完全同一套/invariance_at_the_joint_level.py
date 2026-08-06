import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A88 R331 -- 结局层的不变性,改用联合 R² 曲线

`#285b`:结局层不可用 —— 随机劈的展布 **0.4247**,因为**每个坐标各自的 29 格剖面本身就不稳**
(c1 甚至反号)。
**改比六个坐标合起来的联合 R² 曲线** —— 联合量稳得多(`#264a` 的联合 R² 跨种子只 ±0.4pp)。

ESTIMAND        每组各算 29 个结局的**联合 R²**,得两条 29 维曲线,判它们的相关;
                offset = **随机劈同样大小**。
KILL            **若联合 R² 曲线在两组间也可分辨地不同 -> 结构非不变是稳的,总表必须按组拆开;
                若不可分辨 -> 非不变只出现在载荷层,那是一个关于**构造**的差别,
                不是关于**预测**的差别 —— 两者的实践含义完全不同。**
POSITIVE CTRL   ① 一个**已知按性别不同**的合成量;
                ② **一个真实存在但已知不随性别变的量**(先报它与 `biomale` 的相关以证明这一点)
                   —— **不能再用纯噪声**(`#285c`)。
IMPOSSIBLE      联合 R² 抹掉了「哪个坐标在起作用」;它只判「合起来预测得一样不一样」。
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
def curve(rows,extra=None):
    Q=coords(rows)+([extra] if extra is not None else [])
    m0=np.zeros(NN,bool); m0[rows]=True
    for q_ in Q: m0&=np.isfinite(q_)
    out=[]
    for nm,y in OUT:
        m=m0&np.isfinite(y)
        if m.sum()<250: out.append(np.nan); continue
        X=np.column_stack([np.ones(m.sum())]+[(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
        yy=(y[m]-y[m].mean())/y[m].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
        out.append(float(1-np.var(yy-X@b)/np.var(yy)))
    return np.array(out)
def sim(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>=10 else np.nan
have=ok&np.isfinite(SEX)
g0=np.flatnonzero(have&(SEX==0)); g1=np.flatnonzero(have&(SEX==1))
c0,c1=curve(g0),curve(g1); r_sex=sim(c0,c1)
rngS=np.random.default_rng(20260804)
RND=[]
for t in range(4):
    p=rngS.permutation(np.flatnonzero(have))
    RND.append(sim(curve(p[:len(g0)]),curve(p[len(g0):len(g0)+len(g1)])))
r_rnd,sd_rnd=float(np.mean(RND)),float(np.std(RND))
print(f"两组 n = {len(g0):,} / {len(g1):,}")
print(f"\n**联合 R² 曲线的相关:按性别劈 {r_sex:+.4f} vs 随机劈 {r_rnd:+.4f} ± {sd_rnd:.4f}"
      f"  -> 差 {r_sex-r_rnd:+.4f}(2×展布 {2*sd_rnd:.4f})**")
print(f"   曲线本身:性别组 0 中位 {100*np.nanmedian(c0):.2f}% · 组 1 中位 {100*np.nanmedian(c1):.2f}%")
dif=np.abs(c0-c1); oi=np.argsort(-np.nan_to_num(dif))
print(f"   两组差最大的三格:"+' · '.join(
    f"{OUT[i][0][:22]} {100*c0[i]:.1f}% vs {100*c1[i]:.1f}%" for i in oi[:3]))
for nm in ('conscientiousness','agreeableness','openness'):
    v=np.asarray(EX[nm].values,dtype=float); m=np.isfinite(v)&np.isfinite(SEX)
    print(f"   ⚠ 正对照②候选 `{nm}` 与 biomale 的相关 = {np.corrcoef(v[m],SEX[m])[0,1]:+.4f}")
INV=np.asarray(EX['agreeableness'].values,dtype=float)
FAKE=np.where(have,np.nan_to_num(SEX)*2.0+rngS.standard_normal(NN)*0.5,np.nan)
def pc(extra):
    a=curve(g0,extra); b=curve(g1,extra); s_sex=sim(a,b)
    p=rngS.permutation(np.flatnonzero(have))
    s_rnd=sim(curve(p[:len(g0)],extra),curve(p[len(g0):len(g0)+len(g1)],extra))
    return s_sex-s_rnd
p1=pc(FAKE); p2=pc(INV)
print(f"\n正对照两端:① 加入已知按性别不同的量 -> 差 **{p1:+.4f}**(必须明显 <0)")
print(f"           ② 加入真实存在但几乎不随性别变的 `agreeableness` -> 差 **{p2:+.4f}**(必须 ≈0)")
T=pd.DataFrame([dict(outcome=OUT[i][0][:44],r2_g0=float(c0[i]),r2_g1=float(c1[i])) for i in range(len(OUT))])
check_columns(T,'R331'); T.to_csv(pathlib.Path(__file__).parent/'results'/'joint_curves.csv',index=False)

g=Gate('结局层的不变性(联合 R² 曲线)')
g.asserted('正对照两端:已知按性别不同的量必须被判出;真实但不随性别变的量必须判不出',
           p1<-0.05 and abs(p2)<0.10, f"① {p1:+.4f} · ② {p2:+.4f}(⚠ ② 用的是真实变量,不是纯噪声,`#285c`)")
g.offset_control('★ 联合 R² 曲线:按性别劈 vs 随机劈',r_sex,r_rnd,sd_rnd,
                 null_kind='随机劈同样大小两组 —— 不是零假设,是「若预测不变,按性别劈该落在哪」')
g.asserted('★ 注册的 kill:曲线可分辨地不同 -> 非不变是稳的,总表按组拆;不可分辨 -> 非不变只在构造层',
           abs(r_sex-r_rnd)<2*sd_rnd,
           f"性别劈 {r_sex:+.4f} vs 随机劈 {r_rnd:+.4f} ± {sd_rnd:.4f};差 {r_sex-r_rnd:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

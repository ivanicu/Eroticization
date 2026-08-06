import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A91 R337 -- 非不变性是题目层的,还是这个量本身的

`#291c③`:支配/臣服没被 `#333` 的**措辞**分组挡住,而它们**在内容上**与性别高度相关。
**所以该问:按「内容上与性别相关」重新分组,非不变性还剩多少?**

ESTIMAND        用一个**与我无关的外部判据**给 29 个结局排序 —— **每道题自己与 `biomale` 的 |r|**
                (数据给的,不是我读题读出来的);按它从高到低**逐步剔除**,
                重画「按性别劈 vs 随机劈」的差,得到一条**剔除曲线**。
KILL            **若剔掉与性别最相关的少数几道之后差就塌了 -> 非不变性是「几道题与性别相关」的直接后果;
                若曲线要剔到很深才塌 -> S 的漂移不是题目层的,是这个量本身的。**
⚠ 无需污染声明     排序**由数据给出**,不由我挑。
⚠ 每一步          都要**重跑随机劈的 offset**(结局数变了,展布会变)。
POSITIVE CTRL   沿用 `#290b` 的**交互项驱动**版本。
IMPOSSIBLE      「与性别相关」用的是**边际相关**;一道题可以边际无关而交互有关(支配那一格就可能是)。
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
have=ok&np.isfinite(SEX)
rank=[]
for nm,y in OUT:
    m=have&np.isfinite(y)
    rank.append((nm,abs(float(np.corrcoef(y[m],SEX[m])[0,1])) if m.sum()>300 else 0.0))
rank.sort(key=lambda t:-t[1])
print(f"与 `biomale` 边际相关最高的 6 道(**由数据排序,不是我挑的**):")
for nm,v in rank[:6]: print(f"   {nm[:50]:<52} |r| = {v:.4f}")
print(f"   ⚠ 支配那一格的 |r| = "
      +f"{[v for n,v in rank if 'dominant' in str(n)][0]:.4f}"
      +f"(排第 {[i for i,(n,_) in enumerate(rank) if 'dominant' in str(n)][0]+1})")
g0=np.flatnonzero(have&(SEX==0)); g1=np.flatnonzero(have&(SEX==1))
rngS=np.random.default_rng(20260804)
def curve(rows,only,extra=None):
    Q=coords(rows)+([extra] if extra is not None else [])
    m0=np.zeros(NN,bool); m0[rows]=True
    for q_ in Q: m0&=np.isfinite(q_)
    out=[]
    for nm,y in OUT:
        if nm not in only: continue
        m=m0&np.isfinite(y)
        if m.sum()<250: out.append(np.nan); continue
        X=np.column_stack([np.ones(m.sum())]+[(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
        yy=(y[m]-y[m].mean())/y[m].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
        out.append(float(1-np.var(yy-X@b)/np.var(yy)))
    return np.array(out)
def sim(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>=10 else np.nan
def gap(only,extra=None):
    rs=sim(curve(g0,only,extra),curve(g1,only,extra)); rr=[]
    for _ in range(4):
        p=rngS.permutation(np.flatnonzero(have))
        rr.append(sim(curve(p[:len(g0)],only,extra),curve(p[len(g0):len(g0)+len(g1)],only,extra)))
    return rs-float(np.mean(rr)),float(np.std(rr))
KS=[0,2,4,6,8,10,12]
rows=[]
print(f"\n剔除曲线(逐步剔掉与性别最相关的前 k 道;每步重跑 offset):")
for k in KS:
    only=set(n for n,_ in rank[k:])
    if len(only)<11: break
    dd,sd=gap(only)
    rows.append(dict(k=k,n_left=len(only),delta=dd,sd=sd,ratio=abs(dd)/max(2*sd,1e-9)))
    print(f"   剔掉 {k:>2} 道(剩 {len(only)}):差 **{dd:+.4f}** ± {sd:.4f} "
          f"({abs(dd)/max(2*sd,1e-9):.1f}× 的 2×展布)")
T=pd.DataFrame(rows); check_columns(T,'R337')
T.to_csv(pathlib.Path(__file__).parent/'results'/'removal_curve.csv',index=False)
zS=np.where(np.isfinite(S),(S-np.nanmean(S))/np.nanstd(S),0.0)
FK=np.where(have,0.35*np.nan_to_num(SEX)*zS+rngS.standard_normal(NN),np.nan)
pc,_=gap(set(n for n,_ in rank[12:]),extra=FK)
base_deep=[r for r in rows if r['k']==12]
print(f"\n正对照(在剔掉 12 道之后加入交互项驱动的量):差 **{pc:+.4f}** vs 该步基线 "
      f"{base_deep[0]['delta']:+.4f} -> 变动 **{pc-base_deep[0]['delta']:+.4f}**(必须明显 <0)")
first_dead=next((r['k'] for r in rows if r['ratio']<1.0),None)
print(f"**曲线第一次掉到不可分辨(<1× 的 2×展布)是在剔掉 {first_dead} 道时**"
      if first_dead is not None else "**曲线全程都可分辨**")
g=Gate('非不变性是题目层的还是这个量本身的')
g.asserted('⚠ 排序由数据给出(每道题与 `biomale` 的边际 |r|),不由我挑 -> 无需污染声明',
           True, f"前三:"+' · '.join(f"{n[:16]} {v:.3f}" for n,v in rank[:3]))
g.asserted('正对照:在最深一步加入交互项驱动的量,差必须明显变负',
           (pc-base_deep[0]['delta'])<-0.05, f"{base_deep[0]['delta']:+.4f} -> {pc:+.4f}")
g.asserted('★ 注册的 kill:剔掉少数几道后差就塌 -> 题目层;要剔很深才塌 -> S 这个量本身',
           first_dead is not None and first_dead<=4,
           f"第一次不可分辨在 k={first_dead};曲线 "
           +' · '.join(f"k={r['k']} {r['delta']:+.3f}({r['ratio']:.1f}×)" for r in rows))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

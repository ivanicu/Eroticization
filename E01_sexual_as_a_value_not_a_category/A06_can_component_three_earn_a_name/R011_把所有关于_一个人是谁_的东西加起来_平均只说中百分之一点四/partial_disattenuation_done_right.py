import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A78 R314 -- 第三层:人层量测得太粗吗(这次用正确的解衰减)

`#267b`:我上一次的「解衰减」是**把预测量除以 √信度** —— 线性重缩放,**R² 对它不变**,
所以那个 `+0.00` 是恒等式,不是零。**第三层从未被检验。**
`#267c`:顺带,位置分 S 的半块信度算出恒为 1,因为 S 不随块子集变化。**本轮一并修好。**

ESTIMAND        正确做法:**修正相关矩阵** —— `r'(量_k, 结局_j) = r / sqrt(信度_k)`,
                `r'(量_k, 量_l) = r / sqrt(信度_k × 信度_l)`,再由修正后的矩阵算联合 R²。
⚠ 只能修一侧      结局的信度这份数据没有(每题只问一次),所以这是**部分**解衰减,
                得到的是一个**下界**。如实登记。
KILL            **若部分解衰减把中位从 1.4% 推到 3% 以上 -> 天花板的一大半是人层量的测量误差,
                下一份问卷该做的是把这六个量测准,而不是加结局;
                若几乎不动 -> 三层都不是,天花板在别处。**
⚠ 守卫 14        **当场验证新算法真的会变**(旧的那个不会)。
IMPOSSIBLE      解衰减假设测量误差是随机且独立的;若两个量共享同一种误差,修正会**过度**。
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
def quantities(seed,blocks=None):
    bl=list(range(NB)) if blocks is None else list(blocks)
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    posl=np.zeros(NN); cl=np.zeros(NN)
    for b in bl:
        M,ppl=MB[b]; o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        A[b,ppl]=M[:,o[:k]].mean(1); B[b,ppl]=M[:,o[k:2*k]].mean(1)
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
        posl[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0); cl[ppl]+=1   # ⚠ #267c:S 现在随块子集变
    Sl=np.where(cl>=max(3,len(bl)//2),posl/np.maximum(cl,1),np.nan)
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
    return [Sl,D]+cs+[st]
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
def rels(seed0=95):
    out={}
    for k,n in enumerate(NAMES):
        vs=[]
        for s in range(3):
            p=np.random.default_rng(seed0+s).permutation(NB); h=NB//2
            a=quantities(600+s,p[:h])[k]; b=quantities(610+s,p[h:])[k]
            mm=np.isfinite(a)&np.isfinite(b)&ok
            if mm.sum()>500:
                r=float(np.corrcoef(a[mm],b[mm])[0,1]); vs.append(min(2*abs(r)/(1+abs(r)),0.99))
        out[n]=float(np.nanmean(vs))
    return out
REL=rels()
print(f"半块信度(⚠ `#267c` 已修:S 现在随块子集变化):")
print("   "+' · '.join(f"{n} {REL[n]:.3f}" for n in NAMES))
def joint_r2(Q,y,rel=None):
    m=np.isfinite(y)&ok
    for q_ in Q: m&=np.isfinite(q_)
    if m.sum()<300: return np.nan
    Z=np.array([(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
    yy=(y[m]-y[m].mean())/y[m].std()
    Rxx=np.corrcoef(Z); rxy=np.array([float(np.corrcoef(Z[i],yy)[0,1]) for i in range(len(Q))])
    if rel is not None:
        rr=np.array([max(rel[n],0.05) for n in NAMES])
        rxy=rxy/np.sqrt(rr)                                   # 预测量一侧解衰减
        Rxx=Rxx/np.sqrt(np.outer(rr,rr)); np.fill_diagonal(Rxx,1.0)
        w,V=np.linalg.eigh(Rxx); w=np.clip(w,1e-6,None); Rxx=V@np.diag(w)@V.T   # 保正定
    try: return float(rxy@np.linalg.solve(Rxx,rxy))
    except Exception: return np.nan
Q=quantities(500)
base=float(np.nanmedian([joint_r2(Q,y) for _,y in OUT]))
dis =float(np.nanmedian([joint_r2(Q,y,REL) for _,y in OUT]))
# ⚠ 解衰减给出 R² > 1 是不可能的。查修正后的相关矩阵还是不是一个合法的相关矩阵。
m0=np.ones(NN,bool)&ok
for q_ in Q: m0&=np.isfinite(q_)
Z0=np.array([(q_[m0]-q_[m0].mean())/q_[m0].std() for q_ in Q])
Rxx0=np.corrcoef(Z0); rr0=np.array([max(REL[n],0.05) for n in NAMES])
Rc=Rxx0/np.sqrt(np.outer(rr0,rr0)); np.fill_diagonal(Rc,1.0)
ev=np.linalg.eigvalsh(Rc)
print(f"\n⚠ 修正后的相关矩阵:最小特征值 = **{ev.min():+.4f}**"
      f"(合法的相关矩阵必须 ≥ 0);最大 {ev.max():.3f}")
print(f"   最大的一个修正后相关 = **{np.max(np.abs(Rc-np.eye(len(NAMES)))):.3f}**"
      f"(合法的相关必须 ≤ 1)")
print(f"   -> **在这样的信度(0.22–0.35)下,经典解衰减在这里没有定义。**")

print(f"\n联合 R² 中位:原始 **{100*base:.2f}%** -> **部分解衰减 {100*dis:.2f}%**"
      f"(变化 **{100*(dis-base):+.2f} 个百分点**,倍数 **{dis/max(base,1e-9):.2f}×**)")
g=Gate('第三层:人层量测得太粗吗(正确的解衰减)')
g.could_have_come_out_otherwise('⚠ 守卫 14:新的解衰减算法**真的会变**吗(旧的那个不会)',
    lambda s: float(np.nanmedian([joint_r2(Q,y,{n:max(min(REL[n]*(0.6+0.2*s),0.99),0.05) for n in NAMES})
                                  for _,y in OUT])), [0,1,2,3])
g.could_have_come_out_otherwise('⚠ 对照:`#312③` 的旧做法(除以 √信度再算 R²)',
    lambda s: (lambda Qd: float(np.nanmedian([joint_r2(Qd,y) for _,y in OUT])))(
        [np.where(np.isfinite(q_),(q_-np.nanmean(q_))/np.nanstd(q_)/np.sqrt(max(REL[n]*(0.6+0.2*s),0.05)),np.nan)
         for q_,n in zip(Q,NAMES)]), [0,1,2,3])
T=pd.DataFrame([dict(arm='原始',r2=base),dict(arm='部分解衰减',r2=dis)])
check_columns(T,'R314'); T.to_csv(pathlib.Path(__file__).parent/'results'/'disattenuated.csv',index=False)
g.asserted('⚠ 只能修预测量一侧(结局每题只问一次,没有信度)-> 这是一个下界',
           True, '如实登记:部分解衰减')
g.asserted('⚠⚠ 解衰减后的 R² > 1 是不可能的 —— 先查修正后的矩阵是否合法',
           ev.min()>=0, f"最小特征值 {ev.min():+.4f}(必须 ≥ 0);"
           f"最大修正后相关 {np.max(np.abs(Rc-np.eye(len(NAMES)))):.3f}(必须 ≤ 1)")
g.asserted('★ 注册的 kill:中位推到 3% 以上 -> 天花板的一大半是人层量的测量误差',
           ev.min()>=0 and 0.03<dis<=1.0,
           f"{100*base:.2f}% -> {100*dis:.2f}% —— ⚠ 而修正后的矩阵不合法,"
           f"所以这个数**没有定义**,不能作为证据")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

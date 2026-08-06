import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A78 R311 -- 六个坐标之后,剩下的 98.6% 是噪声还是结构

`#264a`:六个人层量联合解释中位 1.4%。**但它没说剩下的是噪声还是结构。**
这直接决定这个项目该不该继续找第七个维度。

ESTIMAND        把六个量从**每一个结局**里回归掉 -> 29 条残差;
                对 **29×29 残差相关矩阵**取特征谱。
⚠ 零不该是零      29 个结局本身高度相关(大多是情色 Likert 题),所以判据是
                「**扣掉六个坐标之后,还剩多少人层共同结构**」——
                **offset = 把六个坐标换成六个随机量之后的同一条谱。**
KILL            **若残差谱的头几个特征值明显超过 offset -> 存在这六个坐标没抓到的人层因子,
                「找第七个维度」有依据,而且本轮给出它的大小;
                若残差谱贴着 offset -> 这六个坐标已经吃掉了这份数据里可被线性人层结构解释的部分,
                1.4% 就是这个仪器的天花板 —— 那是一条比任何新维度都重要的结论。**
POSITIVE CTRL   守卫 13:种入一个已知的第七维(强度扫描),残差谱必须把它认出来,
                **方向由扫描给出,我不写**。
IMPOSSIBLE      「线性人层结构」——非线性与交互不在内(`#264` 同款登记)。
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
NO=len(OUT)
def spectrum(Q,plant=None,g=0.0):
    m=np.ones(NN,bool)&ok
    for q_ in Q: m&=np.isfinite(q_)
    X=np.column_stack([np.ones(m.sum())]+[(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
    Rz=np.full((NO,int(m.sum())),np.nan)
    for i,(nm,y) in enumerate(OUT):
        yy=y.copy()
        if plant is not None and g: yy=yy+g*np.nan_to_num(plant)*(1 if i%2==0 else -1)
        v=yy[m]; f=np.isfinite(v)
        if f.sum()<300: continue
        b=np.linalg.lstsq(X[f],(v[f]-v[f].mean())/v[f].std(),rcond=None)[0]
        r=np.full(int(m.sum()),np.nan); r[f]=(v[f]-v[f].mean())/v[f].std()-X[f]@b
        Rz[i]=r
    C=np.full((NO,NO),np.nan)
    for i in range(NO):
        for j in range(NO):
            g2=np.isfinite(Rz[i])&np.isfinite(Rz[j])
            if g2.sum()>300: C[i,j]=np.corrcoef(Rz[i][g2],Rz[j][g2])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2
    return np.linalg.eigvalsh(C)[::-1]
rng=np.random.default_rng(20260804)
Q=quantities(500); obs=spectrum(Q)
OFFS=[]
for t in range(4):
    Qr=[np.where(ok,np.random.default_rng(400+t*7+i).standard_normal(NN),np.nan) for i in range(6)]
    OFFS.append(spectrum(Qr))
off=np.mean(OFFS,0); offsd=np.std(OFFS,0)
print(f"残差谱前 6(29 个结局,六个坐标已回归掉):")
print(f"{'':>4}{'观测':>10}{'offset(六个随机量)':>22}{'超出':>10}")
for k in range(6):
    print(f"  λ{k+1}{obs[k]:>10.4f}{off[k]:>16.4f} ± {offsd[k]:.4f}{obs[k]-off[k]:>10.4f}")
excess=float(np.sum(np.clip(obs[:6]-off[:6],0,None)))
print(f"  前 6 个特征值的总超出 = **{excess:.4f}**;λ1 超出 **{obs[0]-off[0]:+.4f}** vs 2×展布 {2*offsd[0]:.4f}")
plant=np.where(ok,rng.standard_normal(NN),np.nan)
SW=[]
for gp in (0.0,0.05,0.15,0.40):
    sp=spectrum(Q,plant=plant,g=gp); SW.append((gp,float(sp[0]-off[0])))
print(f"正对照(种入一个已知第七维,守卫 13 判方向)g -> λ1 超出:"
      +' · '.join(f"{a:.2f}->{b:+.4f}" for a,b in SW))
T=pd.DataFrame([dict(k=i+1,obs=float(obs[i]),off=float(off[i]),off_sd=float(offsd[i]),
                     excess=float(obs[i]-off[i])) for i in range(8)])
check_columns(T,'R311'); T.to_csv(pathlib.Path(__file__).parent/'results'/'residual_spectrum.csv',index=False)

g=Gate('六个坐标之后剩下的是噪声还是结构')
g.plant_direction_from_sweep('正对照:种入已知第七维后 λ1 的超出(方向由扫描给出)',
                             SW,baseline=SW[0][1],baseline_spread=float(offsd[0]),half_of=2*float(offsd[0]))
g.offset_control('★ λ1 观测 vs offset(六个随机量)',float(obs[0]),float(off[0]),float(offsd[0]),
                 null_kind='把六个坐标换成六个随机量之后的同一条谱 —— 不是零假设,'
                           '是「若这六个坐标什么都没吃掉,λ1 该在哪」')
g.asserted('★ 注册的 kill:头几个特征值明显超过 offset -> 存在第七个维度;贴着 offset -> 1.4% 是天花板',
           (obs[0]-off[0])>2*offsd[0],
           f"λ1 {obs[0]:.4f} vs offset {off[0]:.4f} ± {offsd[0]:.4f};前 6 总超出 {excess:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

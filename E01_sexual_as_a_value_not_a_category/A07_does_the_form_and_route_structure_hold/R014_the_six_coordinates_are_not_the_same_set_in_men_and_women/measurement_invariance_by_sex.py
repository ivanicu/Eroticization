import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A88 R330 -- 这套结构对不同的人是不是同一套

这个项目**控制**过性别(`#252a`),**但从没问过这六个坐标对不同性别是不是同一个东西**。
这不是控制问题,是**测量不变性**问题,而它改变整套结构的地位:
**若结构因人群而异,那它不是「人的坐标」,是「这个样本里占多数那群人的坐标」。**

ESTIMAND        按 `biomale` 分两组,各自**独立**重建六个坐标,判:
                ① **载荷层**:两组的块×块残差相关矩阵的前 k 维**子空间复现**;
                ② **结局层**:两组各自的 29 格剖面相关。
⚠ 零不该是零      任意两个半样本的载荷/剖面本来就高度相关。
                **判据是「按性别劈」vs「随机劈同样大小」** —— offset,不是零。
KILL            **若按性别劈与随机劈无差别 -> 结构不变,这是这个项目第一次证明它测的是「人」;
                若明显更低 -> 六个坐标在两个人群里不是同一个东西,总表必须按组拆开报。**
POSITIVE CTRL   两端:① 一个**已知按性别不同**的合成量必须被判出差别;
                ② 一个**已知不随性别变**的量必须判不出。
IMPOSSIBLE      `biomale` 是二值自报;它与其它未测的人口学变量混在一起(`#209` 同款登记)。
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
def prof_(X,rows):
    m=np.zeros(NN,bool); m[rows]=True
    F=np.isfinite(X)&m[None,:]; Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(R[b])
    return R
def coords(rows):
    Ra,Rb=prof_(A,rows),prof_(B,rows)
    m=np.zeros(NN,bool); m[rows]=True
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
    return [S,D]+csx+[st],Vv
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
def profile(x,rows):
    m=np.zeros(NN,bool); m[rows]=True; bi=np.flatnonzero(np.isfinite(x)&m); r=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]) if len(jj)>150 else np.nan)
    return np.array(r)
def sim(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>=10 else np.nan
def compare(rows1,rows2):
    Q1,V1=coords(rows1); Q2,V2=coords(rows2)
    sub=float(np.linalg.norm((V1[:,:3]@V1[:,:3].T)@V2[:,:3],'fro')**2/3)
    pf=[sim(profile(Q1[k],rows1),profile(Q2[k],rows2)) for k in range(6)]
    return sub,float(np.nanmean(pf)),pf
have=ok&np.isfinite(SEX)
g0=np.flatnonzero(have&(SEX==0)); g1=np.flatnonzero(have&(SEX==1))
print(f"两组 n = {len(g0):,} / {len(g1):,}(biomale=0 / 1)")
sub_sex,pf_sex,pfs=compare(g0,g1)
rngS=np.random.default_rng(20260804)
RND=[]
for t in range(4):
    p=rngS.permutation(np.flatnonzero(have)); h=len(g0)
    RND.append(compare(p[:h],p[h:h+len(g1)]))
sub_rnd=float(np.mean([r[0] for r in RND])); sub_sd=float(np.std([r[0] for r in RND]))
pf_rnd=float(np.mean([r[1] for r in RND])); pf_sd=float(np.std([r[1] for r in RND]))
print(f"\n① 载荷层(前 3 维子空间复现):**按性别劈 {sub_sex:.4f}** vs "
      f"**随机劈 {sub_rnd:.4f} ± {sub_sd:.4f}**  -> 差 **{sub_sex-sub_rnd:+.4f}**")
print(f"② 结局层(六个坐标各自的 29 格剖面相关,取均值):**按性别劈 {pf_sex:+.4f}** vs "
      f"**随机劈 {pf_rnd:+.4f} ± {pf_sd:.4f}**  -> 差 **{pf_sex-pf_rnd:+.4f}**")
print(f"   六个坐标分别:"+' · '.join(f"{v:+.3f}" for v in pfs))
FAKE=np.where(have,np.nan_to_num(SEX)*2.0+rngS.standard_normal(NN)*0.5,np.nan)
FLAT=np.where(have,rngS.standard_normal(NN),np.nan)
def one(x,r1,r2): return sim(profile(x,r1),profile(x,r2))
p_sex_fake=one(FAKE,g0,g1); p_rnd_fake=np.mean([one(FAKE,rngS.permutation(np.flatnonzero(have))[:len(g0)],
                                                   rngS.permutation(np.flatnonzero(have))[:len(g1)]) for _ in range(3)])
p_sex_flat=one(FLAT,g0,g1); p_rnd_flat=np.mean([one(FLAT,rngS.permutation(np.flatnonzero(have))[:len(g0)],
                                                    rngS.permutation(np.flatnonzero(have))[:len(g1)]) for _ in range(3)])
print(f"\n正对照两端:① 已知按性别不同的合成量:性别劈 {p_sex_fake:+.4f} vs 随机劈 {p_rnd_fake:+.4f} "
      f"-> 差 **{p_sex_fake-p_rnd_fake:+.4f}**(必须明显 <0)")
print(f"           ② 已知不随性别变的量:性别劈 {p_sex_flat:+.4f} vs 随机劈 {p_rnd_flat:+.4f} "
      f"-> 差 **{p_sex_flat-p_rnd_flat:+.4f}**(必须 ≈0)")
T=pd.DataFrame([dict(level='载荷层',sex=sub_sex,rnd=sub_rnd,sd=sub_sd,delta=sub_sex-sub_rnd),
                dict(level='结局层',sex=pf_sex,rnd=pf_rnd,sd=pf_sd,delta=pf_sex-pf_rnd)])
check_columns(T,'R330'); T.to_csv(pathlib.Path(__file__).parent/'results'/'invariance.csv',index=False)

g=Gate('这套结构对不同的人是不是同一套')
g.asserted('正对照两端:已知按性别不同的量必须被判出差别、已知不变的量必须判不出',
           (p_sex_fake-p_rnd_fake)<-0.10 and abs(p_sex_flat-p_rnd_flat)<0.15,
           f"① {p_sex_fake-p_rnd_fake:+.4f} · ② {p_sex_flat-p_rnd_flat:+.4f}")
g.offset_control('★ 载荷层:按性别劈 vs 随机劈',sub_sex,sub_rnd,sub_sd,
                 null_kind='随机劈同样大小两组 —— 不是零假设,是「若结构不变,按性别劈该落在哪」')
g.offset_control('★ 结局层:按性别劈 vs 随机劈',pf_sex,pf_rnd,pf_sd,
                 null_kind='随机劈同样大小两组 —— 同上')
g.asserted('★ 注册的 kill:按性别劈与随机劈无差别 -> 结构不变;明显更低 -> 总表必须按组拆开',
           abs(sub_sex-sub_rnd)<2*sub_sd and abs(pf_sex-pf_rnd)<2*pf_sd,
           f"载荷层 {sub_sex-sub_rnd:+.4f}(2×展布 {2*sub_sd:.4f})· "
           f"结局层 {pf_sex-pf_rnd:+.4f}(2×展布 {2*pf_sd:.4f})")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

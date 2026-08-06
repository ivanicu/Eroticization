import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A76 R307 -- 不猜名字,直接问 `c3⊥D` 的剖面与哪一个已知外部量最像

三次命名失败的共同点:**我先看载荷,想出名字,再去找题目验证它**(`#201` `#202` `#306`)。
第四次换一条路:**候选集是全部 10 个非情色字段,不是我挑的** —— 所以本轮不需要污染声明。

ESTIMAND        每个非情色字段各算它自己的 29 维结局剖面(**剔除自身那一格**,否则 r=1 自我膨胀),
                判 `c3⊥D` 的剖面与哪一个相关最高,与 **offset** 比。
KILL            **若某一个明显最像(超过 offset 的一半,且领先第二名 > 2×展布)
                -> 那个外部量是名字的第一份候选,而它是数据给的不是我猜的;
                若全部低于 offset 的一半 -> `c3⊥D` 在这份 release 里没有对应的外部量,
                而那本身是一条结论:它需要一个这份数据没测的变量。**
POSITIVE CTRL   把某个外部量的**带噪声复制**当作待判量 -> 必须指认它自己(两端:另一个必须不被指认)。
NEGATIVE CTRL   跨人置换 `c3⊥D`。
IMPOSSIBLE      「最像」是剖面层的相似,不是同一构念的证明(`#259a`:分数层可以几乎不相关)。
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
def build(seed=500):
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
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(Vv[:,2][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,2])[:,None]).sum(0)
    c3=np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
    m=np.isfinite(c3)&np.isfinite(D); X=np.column_stack([np.ones(m.sum()),D[m]])
    o2=np.full(NN,np.nan); o2[m]=c3[m]-X@np.linalg.lstsq(X,c3[m],rcond=None)[0]
    return o2
RES=build()
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
ATT={'Significantly less attractive':-3,'Moderately less attractive':-2,'Slightly less attractive':-1,
     'About average attractiveness':0,'Slightly more attractive':1,'Moderately more attractive':2,
     'Significantly more attractive':3}
EX={'age':d['age'].map(AGE),'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '自评吸引力':d['Compared to other people of your same gender and age range, you are (yh6d44s)'].map(ATT),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,d[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EX.items()]
rng=np.random.default_rng(20260804)
def profile(x,skip=None):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]
    for nm,y in OUT:
        if nm==skip: r.append(np.nan); continue
        mm=np.isfinite(y[bi]); jj=bi[mm]
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]) if len(jj)>200 else np.nan)
    return np.array(r)
def sim(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>=10 else np.nan
pR=profile(RES)
rows=[]
for nm in EX:
    y=np.asarray(EX[nm].values,dtype=float)
    rows.append((nm,abs(sim(profile(y,skip=nm),profile(RES,skip=nm)))))
rows.sort(key=lambda t:-t[1])
print(f"`c3⊥D` 的剖面与 10 个非情色字段的剖面相似度(已剔除自身那一格):")
for nm,v in rows: print(f"   {nm:<18} **{v:.4f}**")
def noisy(x,r_,seed):
    m=np.isfinite(x); z=np.full(NN,np.nan); v=(x[m]-np.nanmean(x))/np.nanstd(x)
    z[m]=np.sqrt(r_)*v+np.sqrt(1-r_)*np.random.default_rng(seed).standard_normal(m.sum()); return z
OFF=[]
for t in range(4):
    a=noisy(RES,0.6,7000+2*t); b=noisy(RES,0.6,7001+2*t)
    OFF.append(abs(sim(profile(a),profile(b))))
off=float(np.mean(OFF))
print(f"\noffset(同一件事的两次带噪声读数)= **{off:.4f} ± {np.std(OFF):.4f}**;"
      f"最像的 {rows[0][0]} 只到 **{100*rows[0][1]/off:.0f}%**")
lead=rows[0][1]-rows[1][1]
sdl=float(np.std([abs(sim(profile(noisy(np.asarray(EX[rows[0][0]].values,dtype=float),0.95,900+i),skip=rows[0][0]),
                          profile(RES,skip=rows[0][0]))) for i in range(8)]))
print(f"第一名 {rows[0][0]} {rows[0][1]:.4f} · 第二名 {rows[1][0]} {rows[1][1]:.4f} · "
      f"领先 **{lead:+.4f}** vs 2×展布 {2*sdl:.4f}")
q=np.asarray(EX['neuroticism'].values,dtype=float)
pc=[(nm,abs(sim(profile(np.asarray(EX[nm].values,dtype=float),skip=nm),
                profile(noisy(q,0.8,55),skip=nm)))) for nm in EX]
pc.sort(key=lambda t:-t[1])
print(f"\n正对照(把 neuroticism 的带噪声复制当作待判量):指认 **{pc[0][0]}**"
      f"({'✅' if pc[0][0]=='neuroticism' else '❌'})· 第二名 {pc[1][0]} {pc[1][1]:.3f}")
nul=[abs(sim(profile(rng.permutation(RES)),profile(np.asarray(EX[rows[0][0]].values,dtype=float)))) for _ in range(20)]
print(f"置换 `c3⊥D` 的零:{np.mean(nul):.4f} ± {np.std(nul):.4f}")
T=pd.DataFrame([dict(field=nm,similarity=v) for nm,v in rows])
check_columns(T,'R307'); T.to_csv(pathlib.Path(__file__).parent/'results'/'name_candidates.csv',index=False)

g=Gate('让数据自己指名字')
g.asserted('正对照:把某个外部量的带噪声复制当作待判量,必须指认它自己',
           pc[0][0]=='neuroticism', f"指认 {pc[0][0]};第二名 {pc[1][0]} {pc[1][1]:.3f}")
g.negative_control('置换 `c3⊥D`',float(np.mean(nul)),float(rows[0][1]),
                   null_spread=float(np.std(nul)),null_kind='跨人置换 —— 只打掉配对')
g.offset_control('★ 最像的那一个 vs 上限',float(rows[0][1]),off,float(np.std(OFF)),
                 null_kind='同一件事的两次带噪声读数 —— 不是零假设,是「若它就是那个外部量,该有多像」')
g.asserted('★ 注册的 kill:某一个超过 offset 的一半且领先第二名 > 2×展布 -> 名字的第一份候选',
           rows[0][1]>0.5*off and lead>2*sdl,
           f"最像 {rows[0][0]} {rows[0][1]:.4f}(上限的 {100*rows[0][1]/off:.0f}%)· "
           f"领先 {lead:+.4f} vs 2×展布 {2*sdl:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

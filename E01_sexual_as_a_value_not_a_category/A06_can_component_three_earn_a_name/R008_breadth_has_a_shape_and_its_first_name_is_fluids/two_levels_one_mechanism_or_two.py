import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A59 R277 -- 罕见的两个层次:一个机制的两次投影,还是两个机制

`#231d`:罕见分裂成了两个层次 —— **块内**你挑的选项有多冷门(位置分 `S`,`#179`),
和**你根本在哪些块上敞开**(越轨半 − 普通半的宽度对比 `D`,`#231`)。
两者只相关约 +0.2,却都通向羞耻。**但没说它们是不是同一个东西的两次投影。**

WORLDS          ① **一个机制,两个投影** —— 真正起作用的是"你偏离常模多远",
                   `S` 与 `D` 只是它的两个读数 ->
                   两条结局剖面**接近"同一件事"的上限**,且共同因子能**吸收**两者对羞耻的效应
                ② **两个机制** —— `S` 是"口味的稀有度",`D` 是"越界的意愿" ->
                   剖面明显低于上限,且至少一个**残差效应存活**
ESTIMAND        (a) `S` 与 `D` 各跑 31 个结局的剖面相关;
                (b) 取两者的第一主成分作为共同因子,判各自对羞耻的**残差**效应。
KILL            **若剖面接近上限 且 两个残差效应都塌进展布 -> 一个机制;
                若剖面明显低于上限 且 至少一个残差效应存活 -> 两个机制。**
⚠ 零应该是零吗?        **不应该**(`#230d` 同款):两者共享同一批人、同一批块、同一条管道。
                判据是 **offset**:「若是同一件事的两次带噪声读数,剖面该有多像」,
                噪声按**各自的分半信度**校准。
NEGATIVE CTRL   置换结局。
POSITIVE CTRL   两端:已知不同(age vs openness)剖面必须低;
                同一变量的两次带噪声复制必须高。
IMPOSSIBLE      「一个机制」与「两个高度相关的机制」在横断面上不可分 ——
                能判的是**这份数据里它们是否可互相替代**,不是它们在人身上是否同源。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
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

def build(seed):
    """一次选项劈半 -> (S_A,S_B, D_A,D_B)。两个层次都用同一次劈分,保证信度可比。"""
    rg=np.random.default_rng(seed); S=[np.zeros(NN),np.zeros(NN)]; ct=[np.zeros(NN),np.zeros(NN)]
    Rh=[np.full((NB,NN),np.nan),np.full((NB,NN),np.nan)]
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        for side,cols in ((0,o[:k]),(1,o[k:2*k])):
            sub=M[:,cols]; rr=-np.log(np.clip(sub.mean(0),1e-4,1.))
            n=sub.sum(1); S[side][ppl]+=np.where(n>0,(sub@rr)/np.maximum(n,1),0.0)
            ct[side][ppl]+=(n>0); Rh[side][b,ppl]=sub.mean(1)
    out=[]
    for side in (0,1):
        out.append(np.where(ct[side]>=6,S[side]/np.maximum(ct[side],1),np.nan))
    Ds=[]
    for side in (0,1):
        X=Rh[side]; F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); c2=F.sum(0)
        R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(c2-F[b]>=6,(tot-Z[b])/np.maximum(c2-F[b],1),np.nan)
            R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
        def hs(cols):
            sub=R[cols]; F2=np.isfinite(sub)
            return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
        t_,o2=hs(TRG),hs(ORD)
        def z(v):
            m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
        Ds.append(z(t_)-z(o2))
    return out[0],out[1],Ds[0],Ds[1]
SA,SB,DA,DB=build(500)
def rel(a,b):
    m=np.isfinite(a)&np.isfinite(b)&ok; r=float(np.corrcoef(a[m],b[m])[0,1]); return 2*r/(1+r)
RS,RD=rel(SA,SB),rel(DA,DB)
S=(SA+SB)/2; D=(DA+DB)/2
m0=np.isfinite(S)&np.isfinite(D)&ok
print(f"块 {NB};n = {int(m0.sum()):,}")
print(f"位置分 S 分半信度 **{RS:+.4f}** · 跨块对比 D **{RD:+.4f}** · corr(S,D) = **{np.corrcoef(S[m0],D[m0])[0,1]:+.4f}**")

lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
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
rng=np.random.default_rng(913)
def prof(x,tag=None):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]; nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl); thr=float(np.nanquantile(np.nanmax(np.array([z[:L] for z in nl]),0),0.95))
    r=np.array(r); h=int(np.nansum(np.abs(r)>thr))
    if tag:
        top=sorted([(OUT[i][0],r[i]) for i in range(len(r)) if np.isfinite(r[i]) and abs(r[i])>thr],
                   key=lambda t:-abs(t[1]))[:4]
        print(f"   {tag:<12} {h:>2}/{len(OUT)}(阈值 {thr:.4f}) "+' · '.join(f"{n[:19]} {v:+.3f}" for n,v in top))
    return r,h
def hits(x):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]; nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl); thr=float(np.nanquantile(np.nanmax(np.array([z[:L] for z in nl]),0),0.95))
    r=np.array(r); return {OUT[i][0] for i in range(len(r)) if np.isfinite(r[i]) and abs(r[i])>thr}, r, thr
def pc(a,b):
    ma=np.isfinite(a)&np.isfinite(b); m2=np.isfinite(a)&np.isfinite(b)
    return float(np.corrcoef(a[m2],b[m2])[0,1])
print("\n两个层次各自的结局剖面:")
(pS,hS),(pD,hD)=prof(S,'位置分 S'),prof(D,'跨块对比 D')
mm=np.isfinite(pS)&np.isfinite(pD); obs=float(np.corrcoef(pS[mm],pD[mm])[0,1])
def noisy(x,r_,seed):
    m=np.isfinite(x); z=np.full(NN,np.nan); v=(x[m]-x[m].mean())/x[m].std()
    z[m]=np.sqrt(max(r_,1e-3))*v+np.sqrt(max(1-r_,0))*np.random.default_rng(seed).standard_normal(m.sum())
    return z
OFF=[]
for t in range(4):
    a_,b_=prof(noisy(S,RS,3000+2*t))[0],prof(noisy(S,RD,3001+2*t))[0]
    q=np.isfinite(a_)&np.isfinite(b_); OFF.append(abs(float(np.corrcoef(a_[q],b_[q])[0,1])))
print(f"   -> 剖面相关 观测 **{obs:+.4f}** · offset(同一件事的两次带噪声读数)**{np.mean(OFF):.4f} ± {np.std(OFF):.4f}**")

# 共同因子:两者的第一主成分,判各自对羞耻的残差效应
SHN=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SHN],errors='coerce').values.astype(float)
m=np.isfinite(S)&np.isfinite(D)&np.isfinite(y)&ok
zs=(S[m]-S[m].mean())/S[m].std(); zd=(D[m]-D[m].mean())/D[m].std(); zy=(y[m]-y[m].mean())/y[m].std()
Zm=np.column_stack([zs,zd]); g1=Zm@np.linalg.eigh(np.cov(Zm,rowvar=False))[1][:,-1]
g1=(g1-g1.mean())/g1.std()
def beta(cols):
    X=np.column_stack([np.ones(m.sum())]+cols); return np.linalg.lstsq(X,zy,rcond=None)[0]
b_alone_s=float(np.corrcoef(zs,zy)[0,1]); b_alone_d=float(np.corrcoef(zd,zy)[0,1])
bg=beta([g1]); bs=beta([g1,zs]); bd=beta([g1,zd]); bb=beta([zs,zd])
boot=lambda cols,i: np.std([np.linalg.lstsq(np.column_stack([np.ones(len(ix))]+[c[ix] for c in cols]),zy[ix],rcond=None)[0][i]
                            for ix in (np.random.default_rng(9000+k).choice(m.sum(),m.sum(),True) for k in range(150))])
sd_s=float(boot([g1,zs],2)); sd_d=float(boot([g1,zd],2))
print(f"\n共同因子(S 与 D 的第一主成分)对羞耻 beta = {bg[1]:+.4f}")
print(f"  加回 S 后:S 的残差 beta **{bs[2]:+.4f}** ± {sd_s:.4f}(单独时 {b_alone_s:+.4f})")
print(f"  加回 D 后:D 的残差 beta **{bd[2]:+.4f}** ± {sd_d:.4f}(单独时 {b_alone_d:+.4f})")
print(f"  两者同时:S {bb[1]:+.4f} · D {bb[2]:+.4f}(n = {m.sum():,})")

HS,rS,tS=hits(S); HD,rD,tD=hits(D)
both=sorted(HS&HD)
print(f"\n两条剖面的交集({len(HS)} ∩ {len(HD)} = **{len(both)}** 个结局):")
for nm in sorted(both,key=lambda n:-abs(rS[[o[0] for o in OUT].index(n)])):
    i=[o[0] for o in OUT].index(nm)
    print(f"   {nm[:52]:<52} S {rS[i]:+.4f} · D {rD[i]:+.4f}")

kd_a,kd_b=prof(EX['age'].values.astype(float))[0],prof(EX['openness'].values.astype(float))[0]
q=np.isfinite(kd_a)&np.isfinite(kd_b); kd=abs(float(np.corrcoef(kd_a[q],kd_b[q])[0,1]))
ks_a,ks_b=prof(noisy(EX['age'].values.astype(float),0.6,41))[0],prof(noisy(EX['age'].values.astype(float),0.6,42))[0]
q=np.isfinite(ks_a)&np.isfinite(ks_b); ks=abs(float(np.corrcoef(ks_a[q],ks_b[q])[0,1]))
print(f"\n正对照两端:已知不同(age vs openness)**{kd:.4f}** · 同一变量两次带噪声复制 **{ks:.4f}**")

T=pd.DataFrame([dict(level='S_块内位置',rel=RS,n_hit=hS,r_shame_alone=b_alone_s,beta_resid=bs[2],sd=sd_s),
                dict(level='D_跨块对比',rel=RD,n_hit=hD,r_shame_alone=b_alone_d,beta_resid=bd[2],sd=sd_d)])
check_columns(T,'R277'); T.to_csv(pathlib.Path(__file__).parent/'results'/'two_levels.csv',index=False)

g=Gate('罕见的两个层次:一个机制还是两个')
g.asserted('正对照两端:已知不同必须低、已知相同必须高',kd<ks-0.15,f"已知不同 {kd:.4f} vs 已知相同 {ks:.4f}")
g.offset_control('剖面相关 vs「同一件事的两次带噪声读数」',abs(obs),float(np.mean(OFF)),float(np.std(OFF)),
                 null_kind='S 的得分 + 校准噪声 —— 不是零假设,是「若两个层次是同一件事,剖面该有多像」')
g.asserted('★ 注册的 kill:剖面明显低于上限 且 至少一个残差效应存活 -> 两个机制',
           abs(obs)<np.mean(OFF)-0.15 and (abs(bs[2])>2*sd_s or abs(bd[2])>2*sd_d),
           f"剖面 {obs:+.4f} vs 上限 {np.mean(OFF):.4f};残差 S {bs[2]:+.4f}±{sd_s:.4f} · D {bd[2]:+.4f}±{sd_d:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

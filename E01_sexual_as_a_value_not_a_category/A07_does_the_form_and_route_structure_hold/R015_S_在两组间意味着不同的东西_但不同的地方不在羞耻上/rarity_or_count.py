import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A93 R340 -- 漂的是「有多冷门」,还是「选了几个」

`#292b` 排除了**题**;`#294a` 排除了**稀有度的定义**。
S = 「每块内所挑选项的**平均**稀有度」再跨块平均 —— 剩下没查的是**「平均」这一步的分母**:
**两组在每块里选了几个。**

ESTIMAND        把 S 的位置分别换成 ① S(基线)② `Nq` = 跨块平均的**块内选项数** ③ `S ⊥ Nq`
                ④ S 与 `Nq` **都放进去**(七坐标),各测 `#286a` 的联合 R² 曲线相关
                (按性别劈 vs 随机劈,**同口径**)。
KILL            **若 ③ `S⊥Nq` 的差塌掉而 ② `Nq` 的差大 -> 非不变性根本不是关于「冷门」的,
                是关于「选了多少」的 —— `#286`–`#294` 整条线的心理学读法要换
                (选择广度 ≠ 稀有度偏好);
                若 ③ 保住 -> 真的是稀有度,`Nq` 只是同行的另一个量。**
POSITIVE CTRL   `SEX × S` 交互项强度扫描(`#290b`),过 guard 13。
NEGATIVE CTRL   由 offset 承担(随机劈,同口径)。
⚠ CONTROL       **`Nq` 与 S 在合并样本上本就相关** —— 先报这个相关,再报 partial,
                否则 ② 和 ① 的差别读不出来是不是同一个东西的两次测量。
IMPOSSIBLE      四臂共用同一批随机劈种子,所以臂间差是配对的;但每臂自己的展布仍只有 4 次劈。
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

# ---- 把 S 拆成两半 ----
cvn=np.zeros(NN); nsum=np.zeros(NN)
for M,ppl in MB:
    cvn[ppl]+=1; nsum[ppl]+=M.sum(1)
NQ=np.where(ok,nsum/np.maximum(cvn,1),np.nan)            # 跨块平均的块内选项数
mfin=np.isfinite(S)&np.isfinite(NQ)
r_SN=float(np.corrcoef(S[mfin],NQ[mfin])[0,1])
print(f"⚠ 先报这个:`corr(S, Nq)` = **{r_SN:+.4f}**(n={mfin.sum():,}) —— "
      f"两者{'不是' if abs(r_SN)<0.30 else '是'}同一个东西的两次测量\n")
def resid(a,b,rows):
    m=np.zeros(NN,bool); m[rows]=True; m&=np.isfinite(a)&np.isfinite(b)
    out=np.full(NN,np.nan)
    x=b[m]; x=(x-x.mean())/x.std(); y=a[m]
    out[m]=y-np.polyval(np.polyfit(x,y,1),x)
    return out
_C={}
def coordsX(rows,arm):
    k=rows.tobytes()
    if k not in _C: _C[k]=coords(rows)
    base=_C[k]; rest=list(base[1:])
    if arm=='S':      return [base[0]]+rest
    if arm=='Nq':     return [NQ]+rest
    if arm=='S⊥Nq':   return [resid(base[0],NQ,rows)]+rest
    if arm=='S+Nq':   return [base[0],NQ]+rest
def curve(rows,arm,outs=None):
    Q=coordsX(rows,arm)
    m0=np.zeros(NN,bool); m0[rows]=True
    for q_ in Q: m0&=np.isfinite(q_)
    out=[]
    for nm,y in (outs if outs is not None else OUT):
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
rngS=np.random.default_rng(20260804); RND=[]
for _ in range(4):
    p=rngS.permutation(np.flatnonzero(have)); RND.append((p[:len(g0)].copy(),p[len(g0):len(g0)+len(g1)].copy()))
print(f"两组 n = {len(g0):,} / {len(g1):,};四臂共用同一批随机劈种子\n")
RES={}
for arm,lab in (('S','① S(基线 `#286a`)'),('Nq','② `Nq` 块内选项数'),
                ('S⊥Nq','③ **`S ⊥ Nq`**'),('S+Nq','④ S 与 `Nq` 都放')):
    s=sim(curve(g0,arm),curve(g1,arm))
    rr=[sim(curve(a,arm),curve(b,arm)) for a,b in RND]
    rm,rd=float(np.mean(rr)),float(np.std(rr)); RES[arm]=(s,rm,rd)
    print(f"  {lab:<22} 性别劈 **{s:+.4f}** vs 随机劈 **{rm:+.4f} ± {rd:.4f}** -> "
          f"差 **{s-rm:+.4f}**({abs(s-rm)/max(2*rd,1e-9):.1f}× 的 2×展布)")
b0=RES['S'][0]-RES['S'][1]; b3=RES['S⊥Nq'][0]-RES['S⊥Nq'][1]; b2=RES['Nq'][0]-RES['Nq'][1]
print(f"\n③ 相对基线保留 **{100*b3/b0:.1f}%** · ② 相对基线 **{100*b2/b0:.1f}%**")
rgP=np.random.default_rng(77)
def plant_outs(g):
    Sz=np.where(np.isfinite(S),(S-np.nanmean(S))/np.nanstd(S),np.nan)
    return [(f'plant{k}',(0.15+0.05*k)*Sz*(1.0+g*SEX)
             +0.30*np.random.default_rng(1000+k).standard_normal(NN)) for k in range(12)]
SW=[]
for g in (0.0,0.6,1.2,2.0):
    po=plant_outs(g)
    s_=sim(curve(g0,'S',po),curve(g1,'S',po))
    r_=float(np.mean([sim(curve(a,'S',po),curve(b,'S',po)) for a,b in RND[:2]]))
    SW.append((g,s_-r_)); print(f"  正对照 g={g:.1f}: **{s_-r_:+.4f}**")
T=pd.DataFrame([dict(arm=k,sex=v[0],rnd=v[1],sd=v[2],delta=v[0]-v[1]) for k,v in RES.items()])
check_columns(T,'R340'); T.to_csv(pathlib.Path(__file__).parent/'results'/'rarity_or_count.csv',index=False)
gg=Gate('漂的是「有多冷门」,还是「选了几个」')
gg.plant_direction_from_sweep('正对照:`SEX × S` 交互项强度扫描',SW,SW[0][1],baseline_spread=2*RES['S'][2])
gg.asserted(f'⚠ 控制:`corr(S, Nq)` = {r_SN:+.4f} —— 两者不是同一个量的两次测量',abs(r_SN)<0.60,
            f"|{r_SN:+.4f}| < 0.60 才能把 ② 和 ① 读成两件事")
for arm,lab in (('S','① S 基线'),('Nq','② `Nq`'),('S⊥Nq','③ **`S ⊥ Nq`**')):
    s,rm,rd=RES[arm]
    gg.offset_control(f'{lab}:按性别劈 vs 随机劈',s,rm,rd,
        null_kind='随机劈同样大小两组,同一批种子 —— 不是零假设,是「若结构相同该落在哪」')
gg.component_difference_is_not_mechanism('★ 注册的 kill:去掉 `Nq` 之后 S 还漂不漂',
    b0,b3,2*RES['S⊥Nq'][2],'`Nq` 块内选项数')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

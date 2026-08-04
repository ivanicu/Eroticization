import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A92 R339 -- 换成组内的尺子,非不变性还剩多少

`#293b`:`rar = −log(选项流行度)` 在**合并样本**上算一次给每个人用,而它在两组里不是同一个向量
(块内流行度按性别劈 +0.9039 vs 随机劈 +0.9973)。**那就换成组内的尺子,看非不变性还剩多少。**

⚠ **只动这一个旋钮。** 另外五个坐标(D · c1 · c2 · c3 · 清晰度)本来就是**组内**算的
(`prof_` 的残差、特征向量、清晰度都带 `m` 掩码),**唯独 S 用的 `rar` 是合并的**。
所以「合并尺子 vs 组内尺子」是一个干净的单变量对照,而且两臂在**同一次运行、同一批种子**下测
(`#316` 的教训:跨运行比 = 两个变化混在一起)。

ESTIMAND        `#286a` 的联合 R² 曲线相关(按性别劈 vs 随机劈),在**两种 S** 下各测一次;
                ⚠ **offset 也必须用同样的口径** —— 随机劈的两组也各自算自己的 `rar`。
KILL            **若换成组内尺子后差塌到不可分辨 -> 非不变性就是尺子问题,而修法是现成的:
                每个人群用自己的 `rar`;
                若仍可分辨 -> 尺子只是其中一部分,还有别的东西在漂。**
POSITIVE CTRL   `SEX × S` **交互项**驱动的合成结局族,强度 g 扫描(`#290b`:主效应造不出组间差),
                过 guard 13:g=0 必须落回随机劈水平,且随 g 单调下降。
NEGATIVE CTRL   由 offset 承担(随机劈,同口径)。
IMPOSSIBLE      组内 `rar` 用组内样本估,所以它自带估计噪声;offset 用同口径正是为此。
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

# ---- 唯一被动的旋钮:S 的尺子 ----
def S_within(rows):
    """用 **组内** 流行度定义稀有度,重建 S。其余五个坐标本来就是组内算的。"""
    m=np.zeros(NN,bool); m[rows]=True
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        sel=m[ppl]
        if sel.sum()<200: continue
        Ms=M[sel]
        rr=-np.log(np.clip(Ms.mean(0),1e-4,1.))            # ⚠ 组内流行度
        n=Ms.sum(1); v=np.where(n>0,(Ms@rr)/np.maximum(n,1),0.0)
        cv[ppl[sel]]+=1; ps[ppl[sel]]+=v
    return np.where(cv>=8,ps/np.maximum(cv,1),np.nan)

_C={}
def coords6(rows,kind):
    k=(id(rows),rows.tobytes()[:64],len(rows))
    if k not in _C: _C[k]=coords(rows)
    base=_C[k]
    return ([base[0]] if kind=='pooled' else [S_within(rows)])+list(base[1:])

def curve(rows,kind,outs=None):
    Q=coords6(rows,kind)
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
POOL=np.flatnonzero(have)
rngS=np.random.default_rng(20260804)
RND=[]
for _ in range(4):
    p=rngS.permutation(POOL); RND.append((p[:len(g0)].copy(),p[len(g0):len(g0)+len(g1)].copy()))
print(f"两组 n = {len(g0):,} / {len(g1):,};结局 {len(OUT)} 道;随机劈 {len(RND)} 次")
print(f"⚠ 唯一被动的旋钮是 S 的 `rar`;D·c1·c2·c3·清晰度本来就是组内算的。\n")

RES={}
for kind,lab in (('pooled','合并尺子(`#286a` 的口径)'),('within','**组内尺子**')):
    s=sim(curve(g0,kind),curve(g1,kind))
    rr=[sim(curve(a,kind),curve(b,kind)) for a,b in RND]
    rm,rd=float(np.mean(rr)),float(np.std(rr))
    RES[kind]=(s,rm,rd)
    print(f"{lab:<26} 性别劈 **{s:+.4f}** vs 随机劈 **{rm:+.4f} ± {rd:.4f}** -> "
          f"差 **{s-rm:+.4f}**({abs(s-rm)/max(2*rd,1e-9):.1f}× 的 2×展布)")
sp,rp,dp=RES['pooled']; sw,rw,dw=RES['within']
print(f"\n换尺子后,差从 **{sp-rp:+.4f}** 变成 **{sw-rw:+.4f}** —— "
      f"保留了 **{100*(sw-rw)/(sp-rp):.1f}%**")

# ---- 正对照:SEX × S 交互项驱动的合成结局族,强度 g 扫描(#290b) ----
rgP=np.random.default_rng(77)
def plant_outs(g):
    Sz=np.where(np.isfinite(S),(S-np.nanmean(S))/np.nanstd(S),np.nan)
    o=[]
    for k in range(12):
        a=0.15+0.05*k
        y=a*Sz*(1.0+g*SEX)+0.30*np.random.default_rng(1000+k).standard_normal(NN)
        o.append((f'plant{k}',y))
    return o
SW=[]
for g in (0.0,0.6,1.2,2.0):
    po=plant_outs(g)
    s_=sim(curve(g0,'pooled',po),curve(g1,'pooled',po))
    r_=float(np.mean([sim(curve(a,'pooled',po),curve(b,'pooled',po)) for a,b in RND[:2]]))
    SW.append((g,s_-r_)); print(f"  正对照 g={g:.1f}: 性别劈−随机劈 = **{s_-r_:+.4f}**")

T=pd.DataFrame([dict(ruler='合并',sex=sp,rnd=rp,sd=dp,delta=sp-rp),
                dict(ruler='组内',sex=sw,rnd=rw,sd=dw,delta=sw-rw)])
check_columns(T,'R339'); T.to_csv(pathlib.Path(__file__).parent/'results'/'within_group_ruler.csv',index=False)

g_=Gate('换成组内的尺子,非不变性还剩多少')
g_.plant_direction_from_sweep('正对照:`SEX × S` 交互项强度扫描(#290b:主效应造不出组间差)',
    SW, SW[0][1], baseline_spread=2*dp)
g_.offset_control('合并尺子:按性别劈 vs 随机劈',sp,rp,dp,
    null_kind='随机劈同样大小两组,**同用合并尺子** —— 不是零假设,是「若结构相同该落在哪」')
g_.offset_control('★ **组内尺子**:按性别劈 vs 随机劈',sw,rw,dw,
    null_kind='随机劈同样大小两组,**两组各自算自己的 `rar`** —— 同口径,否则是两个变化混在一起')
g_.asserted('★ 注册的 kill:换成组内尺子后差是否塌到不可分辨',
    abs(sw-rw)<2*dw, f"组内 差 {sw-rw:+.4f} vs 2×展布 {2*dw:.4f};合并 差 {sp-rp:+.4f}")
print(g_)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

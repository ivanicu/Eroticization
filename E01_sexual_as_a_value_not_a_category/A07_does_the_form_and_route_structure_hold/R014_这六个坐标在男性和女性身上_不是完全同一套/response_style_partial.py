import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A94 R341 -- 差在预测量上,还是在结局自己身上

`#295c`:题、尺子、分母三个零件都换过了,差一直在 —— 而那正是 frontier §3 的**盆地形状**。
**所以这一轮是我不希望它成立的那一个。**

联合 R² 曲线是**跨 29 个结局**比的。一个组间的曲线差,可以来自**预测量**(我一直在查的),
**也可以来自结局自己**:若某组在结局 k 上全距更窄、或答题风格不同,
它在那格的 R² 就会因为**与 S 无关的原因**变化。**这个混淆我从 `#286` 起一次都没控制过。**

ESTIMAND        结局侧参考曲线 = 只用**答题风格**(跨所有 Likert 题的**均值 · 标准差 · 作答数**,
                三个纯测量侧变量,都不在结局里)算的 29 维 R² 曲线;
                把它**逐组偏出**六坐标曲线后重测「按性别劈 vs 随机劈」。
KILL            **若偏出后差塌掉 -> `#286`–`#295` 整条线测的是结局侧的测量差,不是结构差,
                这十轮的心理学读法整个换掉;
                若差还在 -> 它真的在预测量上,而这十轮才算站住。**
POSITIVE CTRL   人为给一组的 6 个结局加噪(信度下降),**偏出必须把这个人为差收回去**。
NEGATIVE CTRL   由 offset 承担(随机劈,**同样做偏出**)。
⚠ CONTROL       **必须同时报「偏出后随机劈的展布有没有一起缩」** ——
                偏出会同时压缩两边,不然就是拿更小的差比更小的展布。
IMPOSSIBLE      答题风格是测量侧的**代理**,不是信度本身;单题 Likert 没有重复测量,信度不可直接估。
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

# ---- 结局侧:答题风格(纯测量侧,都不在 OUT 里) ----
LK=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in lik])
STY=[('风格·均值',np.nanmean(LK,1)),('风格·标准差',np.nanstd(LK,1)),
     ('风格·作答数',np.isfinite(LK).sum(1).astype(float))]
print(f"答题风格 {len(STY)} 个变量,建在 {LK.shape[1]} 道 Likert 上;都不在 {len(OUT)} 个结局里")
_C={}
def coords6(rows):
    k=rows.tobytes()
    if k not in _C: _C[k]=coords(rows)
    return _C[k]
def _r2(rows,Q,outs):
    m0=np.zeros(NN,bool); m0[rows]=True
    for q_ in Q: m0&=np.isfinite(q_)
    out=[]
    for nm,y in outs:
        m=m0&np.isfinite(y)
        if m.sum()<250: out.append(np.nan); continue
        X=np.column_stack([np.ones(m.sum())]+[(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
        yy=(y[m]-y[m].mean())/y[m].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
        out.append(float(1-np.var(yy-X@b)/np.var(yy)))
    return np.array(out)
def curve(rows,outs=None):  return _r2(rows,coords6(rows),outs if outs is not None else OUT)
def refcurve(rows,outs=None): return _r2(rows,[v for _,v in STY],outs if outs is not None else OUT)
def sim(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>=10 else np.nan
def sim_p(a,b,ra,rb):
    """⚠ 逐组把结局侧参考曲线偏出,再比残差曲线。"""
    m=np.isfinite(a)&np.isfinite(b)&np.isfinite(ra)&np.isfinite(rb)
    if m.sum()<10: return np.nan
    def rs(y,x):
        x=(x[m]-x[m].mean())/x[m].std(); return y[m]-np.polyval(np.polyfit(x,y[m],1),x)
    A_,B_=rs(a,ra),rs(b,rb)
    return float(np.corrcoef(A_,B_)[0,1])
have=ok&np.isfinite(SEX)
g0=np.flatnonzero(have&(SEX==0)); g1=np.flatnonzero(have&(SEX==1))
rngS=np.random.default_rng(20260804); RND=[]
for _ in range(4):
    p=rngS.permutation(np.flatnonzero(have)); RND.append((p[:len(g0)].copy(),p[len(g0):len(g0)+len(g1)].copy()))
print(f"两组 n = {len(g0):,} / {len(g1):,}\n")
c0,c1_=curve(g0),curve(g1); r0,r1=refcurve(g0),refcurve(g1)
print(f"⚠ 先报结局侧本身:参考曲线的组间相关 **{sim(r0,r1):+.4f}** "
      f"—— 答题风格在两组间{'也不一样' if sim(r0,r1)<0.85 else '基本一样'}")
print(f"⚠ 六坐标曲线与参考曲线的相关:组0 **{sim(c0,r0):+.4f}** · 组1 **{sim(c1_,r1):+.4f}**")
def arm(fn,tag):
    s=fn(g0,g1); rr=[fn(a,b) for a,b in RND]
    rm,rd=float(np.mean(rr)),float(np.std(rr))
    print(f"  {tag:<20} 性别劈 **{s:+.4f}** vs 随机劈 **{rm:+.4f} ± {rd:.4f}** -> "
          f"差 **{s-rm:+.4f}**({abs(s-rm)/max(2*rd,1e-9):.1f}× 的 2×展布)")
    return s,rm,rd
raw=arm(lambda a,b: sim(curve(a),curve(b)),'① 原样(`#286a`)')
par=arm(lambda a,b: sim_p(curve(a),curve(b),refcurve(a),refcurve(b)),'② **偏出结局侧后**')
print(f"\n偏出后差从 **{raw[0]-raw[1]:+.4f}** 变成 **{par[0]-par[1]:+.4f}** —— "
      f"保留 **{100*(par[0]-par[1])/(raw[0]-raw[1]):.1f}%**")
print(f"⚠ 展布有没有一起缩:随机劈展布 **{raw[2]:.4f}** -> **{par[2]:.4f}** "
      f"({100*par[2]/raw[2]:.0f}%)")
# ---- 正对照:给组1 的 6 个结局加噪,偏出必须把这个人为差收回去 ----
rgN=np.random.default_rng(4242); HIT=set(range(0,len(OUT),5))
def noisy(rows,lvl):
    m=np.zeros(NN,bool); m[rows]=True
    o=[]
    for k,(nm,y) in enumerate(OUT):
        yy=y.copy()
        if k in HIT and lvl>0:
            z=lvl*np.nanstd(y)*np.random.default_rng(700+k).standard_normal(NN)
            yy=np.where(m,y+z,y)
        o.append((nm,yy))
    return o
print(f"\n正对照:给**组1**的 {len(HIT)} 个结局加噪(信度下降),偏出必须把人为差收回去")
for lvl in (0.0,1.0,2.0):
    oo0,oo1=noisy(g1,0.0),noisy(g1,lvl)   # 组0 不动,组1 加噪
    a,b=curve(g0,oo0),curve(g1,oo1); ra,rb=refcurve(g0,oo0),refcurve(g1,oo1)
    print(f"  lvl={lvl:.1f}: 原样 **{sim(a,b):+.4f}** · 偏出后 **{sim_p(a,b,ra,rb):+.4f}**")
    if lvl==0.0: base_raw,base_par=sim(a,b),sim_p(a,b,ra,rb)
    if lvl==2.0: hi_raw,hi_par=sim(a,b),sim_p(a,b,ra,rb)
rec=abs(hi_par-base_par)/max(abs(hi_raw-base_raw),1e-9)   # ⚠ 两个量都是负的;上一版用 max(x,1e-9) 夹了一个负数,把比值毁成 4.3e8 而守卫仍 PASS
print(f"  收回比例:原样掉 **{hi_raw-base_raw:+.4f}**,偏出后掉 **{hi_par-base_par:+.4f}** "
      f"—— 收回 **{100*(1-rec):+.1f}%**(>0 才叫收回)")
T=pd.DataFrame([dict(arm='原样',sex=raw[0],rnd=raw[1],sd=raw[2],delta=raw[0]-raw[1]),
                dict(arm='偏出结局侧',sex=par[0],rnd=par[1],sd=par[2],delta=par[0]-par[1])])
check_columns(T,'R341'); T.to_csv(pathlib.Path(__file__).parent/'results'/'outcome_side.csv',index=False)
gg=Gate('差在预测量上,还是在结局自己身上')
gg.bounded_statistic_out_of_range('⚠ guard 17:收回比例先过自己的定义域',1-rec,-2,2,'收回比例')
gg.asserted('★ 正对照:给组1 的 6 个结局加噪,偏出必须把人为差收回去',(1-rec)>0.40,
            f"原样掉 {hi_raw-base_raw:+.4f} -> 偏出后掉 {hi_par-base_par:+.4f}(收回 {100*(1-rec):+.1f}%)—— 不到 +40% 则这个代理**测不动信度**,整臂 UNVERIFIED")
gg.asserted('⚠ 展布没有被偏出压掉(否则是拿更小的差比更小的展布)',par[2]>0.4*raw[2],
            f"随机劈展布 {raw[2]:.4f} -> {par[2]:.4f}({100*par[2]/raw[2]:.0f}%)")
gg.offset_control('① 原样:按性别劈 vs 随机劈',raw[0],raw[1],raw[2],
    null_kind='随机劈同样大小两组 —— 不是零假设,是「若结构相同该落在哪」')
gg.offset_control('★ ② **偏出结局侧后**:按性别劈 vs 随机劈',par[0],par[1],par[2],
    null_kind='随机劈同样大小两组,**同样做偏出** —— 同口径,否则是两个变化混在一起')
gg.component_difference_is_not_mechanism('★ 注册的 kill:偏出结局侧后差是否塌掉',
    raw[0]-raw[1],par[0]-par[1],2*par[2],'结局侧答题风格')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A97 R347 -- 嵌套 CV:连六个坐标一起留出去

`#301d`:1.087% 是**上界**,因为六个坐标是在**全体人**身上估的,只有回归系数被折叠了。
而**坐标才是这个方法的主要产物** —— 它们的过拟合才是真问题。

ESTIMAND        5 折人层嵌套 CV:每折**只用训练集**重估全部六个坐标
                (块内留一剖面 · 越轨/普通对比 · 三个特征向量 · 清晰度),
                测试集的人**投影**到训练集估出来的结构上,读测试集 R²。
KILL            **若明显低于 1.087% -> 页面上那个数还要再降一档,并且要说清「坐标本身也是估出来的」;
                若几乎不动 -> 六个坐标很稳,那本身是一条关于它们的好消息。**
POSITIVE CTRL   合成一个**总体 R² 已知**的结局,嵌套 CV 必须仍收敛到真值 ——
                否则读到的是**投影的偏差**,不是过拟合(投影是新代码,这是它唯一的验证)。
NEGATIVE CTRL   纯噪声结局:嵌套 CV 必须给 ≈0 或负。
⚠ 符号          特征向量符号不必跨折对齐 —— 回归会吸收它。这一条是**代码简化的理由**,写下来
                是因为 `R210:73` 那次符号坑正是在**求和**而不是回归里出的。
IMPOSSIBLE      嵌套 CV 量的是「这套流程用在没见过的人身上」;它**不**告诉你坐标的真值有多准。
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

def fit_apply(tr,te):
    """在 tr 上估全部六个坐标,把 te 的人**投影**上去。返回 (Q_tr, Q_te)。"""
    mt=np.zeros(NN,bool); mt[tr]=True
    # ① S:稀有度用**训练集**流行度
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        s=mt[ppl]
        if s.sum()<200: continue
        rr=-np.log(np.clip(M[s].mean(0),1e-4,1.)); n=M.sum(1)
        v=np.where(n>0,(M@rr)/np.maximum(n,1),np.nan)      # 训练集的尺子,量所有人
        g=np.isfinite(v); cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    S_=np.where(cv>=8,ps/np.maximum(cv,1),np.nan)
    # ② 剖面:块的留一均值也只用训练集的人
    def prof_(X):
        # ⚠ 留一块均值是**人内**量(这个人自己的其它块),**不是从别人身上估的** ——
        # 所以它不能加训练掩码。第一版加了,于是测试集的剖面被整片抹成 NaN
        # (诊断的告密者:坐标的有限数**恰好等于** |训练集|)。
        # 真正从别人身上估的只有**中心化常数 mu**,那一个才用训练集。
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0)
        R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo,np.nan)
            mu=np.nanmean(np.where(F[b]&mt,R[b],np.nan))   # ← 只有这一个来自训练集
            R[b]=R[b]-mu
        return R
    Ra,Rb=prof_(A),prof_(B)
    st=np.full(NN,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
    for i in np.concatenate([tr,te]):
        mm=G[:,i]
        if mm.sum()<8: continue
        x,y_=Ra[mm,i],Rb[mm,i]
        if x.std()>1e-9 and y_.std()>1e-9: st[i]=float(np.corrcoef(x,y_)[0,1])
    def z(v):
        k=np.isfinite(v)&mt                                 # 标准化常数来自训练集
        w=np.full(NN,np.nan); g=np.isfinite(v)
        w[g]=(v[g]-v[k].mean())/v[k].std(); return w
    def hs(R,cols):
        sub=R[cols]; F2=np.isfinite(sub)
        return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
    D_=(z(hs(Ra,TRG))-z(hs(Ra,ORD))+z(hs(Rb,TRG))-z(hs(Rb,ORD)))/2
    C=np.full((NB,NB),np.nan)                               # 特征向量只用训练集的人估
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&mt
            if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0); cs=[]
    for k in range(3):                                      # 投影:所有人,训练集的特征向量
        num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
        cs.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    return [S_,D_]+cs+[st]

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

rgF=np.random.default_rng(20260804)
ALL=np.flatnonzero(ok); idx=rgF.permutation(ALL); FOLD=np.array_split(idx,5)
print(f"n={len(ALL):,},5 折 · 每折重估六个坐标")
QT=[]
for k in range(5):
    te=FOLD[k]; tr=np.concatenate([FOLD[j] for j in range(5) if j!=k])
    QT.append((tr,te,fit_apply(tr,te)))
    print(f"  折 {k+1}/5 完成(训练 {len(tr):,} / 测试 {len(te):,})")
def nested_r2(y):
    num=0.0; den=0.0; tot=[]
    for tr,te,Q in QT:
        mt=np.isfinite(y); mte=mt.copy()
        for q_ in Q: mt&=np.isfinite(q_); 
        m_tr=np.zeros(NN,bool); m_tr[tr]=True; m_te=np.zeros(NN,bool); m_te[te]=True
        a=mt&m_tr; b=mt&m_te
        if a.sum()<200 or b.sum()<50: continue
        Xa=np.column_stack([np.ones(a.sum())]+[(q_[a]-q_[a].mean())/max(q_[a].std(),1e-12) for q_ in Q])
        ya=y[a]; mu,sd=ya.mean(),ya.std()
        bcoef=np.linalg.lstsq(Xa,(ya-mu)/sd,rcond=None)[0]
        Xb=np.column_stack([np.ones(b.sum())]+[(q_[b]-q_[a].mean())/max(q_[a].std(),1e-12) for q_ in Q])
        pred=Xb@bcoef; yb=(y[b]-mu)/sd
        num+=float(((yb-pred)**2).sum()); den+=float(((yb-yb.mean())**2).sum()); tot.append(b.sum())
    return (1-num/den) if den>0 else np.nan, int(np.sum(tot))
tr0,te0,Q0d=QT[0]
FIN=[int(np.isfinite(q_).sum()) for q_ in Q0d]
print("  诊断:每个坐标的有限数 ->", FIN, f"(|训练集| {len(tr0):,} / 全体 {len(ALL):,})")
_y=OUT[0][1].astype(float); _m=np.isfinite(_y)
for q_ in Q0d: _m&=np.isfinite(q_)
print(f"  诊断:结局0 与六坐标同时有限的人 = {int(_m.sum()):,}")
rows=[]
for nm,y in OUT:
    r,n_=nested_r2(y.astype(float))
    if np.isfinite(r): rows.append(dict(v_out=str(nm)[:44],n=n_,nested=r))
T=pd.DataFrame(rows); check_columns(T,'R347')
med=100*float(T.nested.median())
print(f"\n**嵌套 CV 的联合 R² 中位 = {med:.3f}%**  (n 中位 {int(T.n.median()):,};"
      f"结局 {len(T)} 个;为负的 **{int((T.nested<0).sum())}/{len(T)}**)")
print(f"   对比:样本内 **1.390%** · 固定坐标的 5 折 CV **1.087%** · **嵌套 {med:.3f}%**")
Q0=QT[0][2]
print(f"\n正对照/负对照:合成结局,总体 R² 已知")
mm=np.ones(NN,bool)
for q_ in Q0: mm&=np.isfinite(q_)
mm&=ok; base=np.column_stack([q_[mm] for q_ in Q0]); rg=np.random.default_rng(3)
CTL={}
for true_r2 in (0.000,0.014,0.050):
    bb=rg.standard_normal(base.shape[1]); s=base@bb; s=(s-s.mean())/s.std()
    yv=np.full(NN,np.nan); yv[mm]=np.sqrt(true_r2)*s+np.sqrt(1-true_r2)*rg.standard_normal(len(s))
    r,_=nested_r2(yv); CTL[true_r2]=100*r
    print(f"   真值 {100*true_r2:.1f}% -> **嵌套 CV {100*r:.3f}%**")
T.to_csv(pathlib.Path(__file__).parent/'results'/'nested_cv.csv',index=False)
gg=Gate('嵌套 CV:连六个坐标一起留出去')
gg.apply_reached_the_test_set('⚠ guard 18:投影这一步真的跑到测试集了吗',FIN,len(tr0),len(ALL),
                              labels=['S','D','c1','c2','c3','清晰度'])
gg.asserted('★ 正对照:真值 1.4% 时嵌套 CV 必须仍收敛到真值附近(±0.6pp)',abs(CTL[0.014]-1.4)<0.6,
            f"真值 1.4% -> 嵌套 {CTL[0.014]:.3f}% —— 这是**投影这段新代码**唯一的验证")
gg.asserted('★ 负对照:真值 0 时嵌套 CV 必须 ≈0 或负',CTL[0.0]<0.2,f"真值 0 -> {CTL[0.0]:.3f}%")
gg.asserted('★ 注册的 kill:嵌套 CV 是否明显低于固定坐标的 1.087%',abs(med-1.087)>0.10,
            f"嵌套 **{med:.3f}%** vs 固定坐标 **1.087%** vs 样本内 **1.390%**")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

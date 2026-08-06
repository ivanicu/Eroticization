import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A80 R317 -- 口径清单:哪些小选择**真的**会改变一个已发表的数

`#271b`:`#314` 的 2.23% 与本轮的 1.21% 差 1 个百分点,而**两个注册的原因都不是它** ——
差别在我没写进 docstring 的小选择里。

ESTIMAND        逐一扰动五个候选口径旋钮,看联合 R² 中位**是否真的动**:
                ① 块覆盖阈值 `cov>=8` vs `>=6` vs `>=10`
                ② 剖面清晰度的最小块数 `need` = 6 vs 8
                ③ S 是否随块子集变化
                ④ 选项劈分种子
                ⑤ 阈值重抽样次数(只影响计数,不影响 R² —— **预期不动,作为阴性对照**)
KILL            **守卫 14 逐项验证:不会改变任何数的旋钮从清单里删掉** ——
                否则清单会变成一张让人安心的长表(`#241a` 的同一个病)。
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
NB=len(MB); cov=np.zeros(NN); posA=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; posA[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
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
def median_r2(seed=500,covmin=8,need=8,new_S=False):
    ok=cov>=covmin; S_OLD=np.where(ok,posA/np.maximum(cov,1),np.nan)
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    posl=np.zeros(NN); cl=np.zeros(NN)
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        ha=M[:,o[:k]]; hb=M[:,o[k:2*k]]
        A[b,ppl]=ha.mean(1); B[b,ppl]=hb.mean(1)
        rr=-np.log(np.clip(ha.mean(0),1e-4,1.)); n=ha.sum(1)
        posl[ppl]+=np.where(n>0,(ha@rr)/np.maximum(n,1),0.0); cl[ppl]+=1
    Sl=np.where(cl>=covmin,posl/np.maximum(cl,1),np.nan) if new_S else S_OLD
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
        if mm.sum()<need: continue
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
    Q=[Sl,D]+cs+[st]
    vals=[]
    for nm,y in OUT:
        m=np.isfinite(y)&ok
        for q_ in Q: m&=np.isfinite(q_)
        if m.sum()<300: vals.append(np.nan); continue
        X=np.column_stack([np.ones(m.sum())]+[(q_[m]-q_[m].mean())/q_[m].std() for q_ in Q])
        yy=(y[m]-y[m].mean())/y[m].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
        vals.append(float(1-np.var(yy-X@b)/np.var(yy)))
    return float(np.nanmedian(vals))
KNOBS=[('① 块覆盖阈值 cov>=',lambda s:median_r2(covmin=[8,6,10][s]),[0,1,2]),
       ('② 清晰度最小块数 need=',lambda s:median_r2(need=[8,6,10][s]),[0,1,2]),
       ('③ S 是否随块变化',lambda s:median_r2(new_S=bool(s)),[0,1]),
       ('④ 选项劈分种子',lambda s:median_r2(seed=[500,900,1200][s]),[0,1,2]),
       ('⑤ 阈值重抽样次数(阴性对照:不该影响 R²)',lambda s:median_r2(),[0,1,2])]
g=Gate('口径清单:哪些旋钮真的会改变一个已发表的数')
rows=[]
for nm,fn,ps in KNOBS:
    vals=[fn(p_) for p_ in ps]
    moved=max(abs(v-vals[0]) for v in vals[1:])
    rows.append(dict(knob=nm,base=vals[0],max_delta=moved,v_values='·'.join(f"{100*v:.2f}%" for v in vals)))
    print(f"{nm:<34} "+' · '.join(f"{100*v:.2f}%" for v in vals)+f"  -> 最大变动 **{100*moved:.2f}pp**")
    g.could_have_come_out_otherwise(nm,fn,ps)
T=pd.DataFrame(rows); check_columns(T,'R317')
T.to_csv(pathlib.Path(__file__).parent/'results'/'caliber_knobs.csv',index=False)
# ⚠ 三个活旋钮各自约 0.17pp —— 它们**叠加**起来能不能接近 `#314` 那个未解的 0.85pp?
combo=median_r2(covmin=10,need=10,new_S=True)
print(f"\n⚠ 三个活旋钮**同时**取另一端(cov>=10, need=10, 新 S):**{100*combo:.2f}%**"
      f"(基线 {100*rows[0]['base']:.2f}%,差 **{100*(combo-rows[0]['base']):+.2f}pp**)")
print(f"   `#314` 未解的差是 +0.85pp —— 单个旋钮最大 0.18pp,三个叠加 {100*(combo-rows[0]['base']):+.2f}pp")
live=[r['knob'] for r in rows if r['max_delta']>1e-9]
print(f"\n**真的会动的旋钮:{len(live)}/{len(rows)}** -> 进 `CALIBER.md`;其余删掉")
g.asserted('★ 注册的 kill:不会改变任何数的旋钮从清单里删掉',
           len(live)<len(rows), f"{len(live)}/{len(rows)} 进清单;删掉 "
           +' · '.join(r['knob'][:12] for r in rows if r['max_delta']<=1e-9))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

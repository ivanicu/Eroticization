import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A81 R319 -- 一个最朴素的「你离常模有多远」,打不打得赢这六个坐标

这个项目有**六个坐标、14 个守卫、一张口径清单**,而它对**羞耻**的联合解释是 **2.77%**(`#264a`)。
**从没有人问过一个更朴素的对手。**

ESTIMAND        `atypicality` = 一个人**所有**块上所选选项的**平均稀有度**
                —— 不残差化、不特征分解、不劈半、不留一。**最朴素的那一个。**
                判:它单独对羞耻的相关;六坐标的联合 R²;两者互为增量。
KILL            **若 `atypicality` 单独就接近或超过六坐标的联合解释 -> 这整套装置对羞耻是过度设计的,
                必须写进公开页;
                若六坐标明显更强且 `atypicality` 的增量接近零 -> 这套装置抓到了朴素量抓不到的东西。**
⚠ 必须同报        `atypicality` 与位置分 S **高度相关**(S 就是它的留一残差化版本),
                所以**用增量 R² 下结论,不用单独 R²**。
守卫 14          验证两个增量**真的会变**。
IMPOSSIBLE      只对**羞耻**这一个结局判;别的结局不在本轮范围。
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
NB=len(MB); cov=np.zeros(NN); pos=np.zeros(NN); rar_sum=np.zeros(NN); npick=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
    rar_sum[ppl]+=M@rr; npick[ppl]+=n                       # ⚠ atypicality:最朴素的那一个
ok=cov>=8
S=np.where(ok,pos/np.maximum(cov,1),np.nan)
ATYP=np.where(ok&(npick>0),rar_sum/np.maximum(npick,1),np.nan)
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
def quantities(seed=500):
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
SH=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SH],errors='coerce').values.astype(float)
Q=quantities()
m0=np.isfinite(y)&ok&np.isfinite(ATYP)
for q_ in Q: m0&=np.isfinite(q_)
print(f"n = {int(m0.sum()):,}")
print(f"⚠ 必须先报:corr(atypicality, 位置分 S) = "
      f"**{np.corrcoef(ATYP[m0],S[m0])[0,1]:+.4f}**(S 就是它的留一残差化版本)")
def r2(cols):
    if not cols: return 0.0
    X=np.column_stack([np.ones(m0.sum())]+[(c[m0]-c[m0].mean())/c[m0].std() for c in cols])
    yy=(y[m0]-y[m0].mean())/y[m0].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
    return float(1-np.var(yy-X@b)/np.var(yy))
r_at=float(np.corrcoef(ATYP[m0],y[m0])[0,1])
six=r2(Q); at=r2([ATYP]); both=r2(Q+[ATYP])
print(f"\n对羞耻:")
print(f"  `atypicality` 单独相关 **{r_at:+.4f}** · 单独 R² **{100*at:.2f}%**")
print(f"  六坐标联合 R² **{100*six:.2f}%**(`#264a` 报 2.77%)")
print(f"  两者合起来 **{100*both:.2f}%**")
print(f"  **`atypicality` 的增量 {100*(both-six):+.2f}pp · 六坐标的增量 {100*(both-at):+.2f}pp**")
rng=np.random.default_rng(20260804)
bt=lambda f: float(np.std([f(i) for i in (rng.choice(np.flatnonzero(m0),int(m0.sum()),True) for _ in range(150))]))
def r2_on(cols,idx):
    X=np.column_stack([np.ones(len(idx))]+[(c[idx]-c[idx].mean())/c[idx].std() for c in cols])
    yy=(y[idx]-y[idx].mean())/y[idx].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
    return float(1-np.var(yy-X@b)/np.var(yy))
sd_at=bt(lambda i: r2_on(Q+[ATYP],i)-r2_on(Q,i))
sd_six=bt(lambda i: r2_on(Q+[ATYP],i)-r2_on([ATYP],i))
print(f"  自助展布:atypicality 增量 ±{100*sd_at:.2f}pp · 六坐标增量 ±{100*sd_six:.2f}pp")
T=pd.DataFrame([dict(model='atypicality 单独',r2=at),dict(model='六坐标',r2=six),
                dict(model='两者',r2=both)])
check_columns(T,'R319'); T.to_csv(pathlib.Path(__file__).parent/'results'/'naive_rival.csv',index=False)

g=Gate('朴素对手打不打得赢六个坐标')
g.could_have_come_out_otherwise('⚠ 守卫 14:两个增量真的会变吗',
    lambda s: r2_on(Q+[ATYP],rng.choice(np.flatnonzero(m0),int(m0.sum()),True))
              -r2_on(Q,rng.choice(np.flatnonzero(m0),int(m0.sum()),True)), [0,1,2,3])
g.has_error_bar('atypicality 的增量 R²',float(both-six),sd_at,'bootstrap_人层')
g.has_error_bar('六坐标的增量 R²',float(both-at),sd_six,'bootstrap_人层')
g.asserted('★ 注册的 kill:atypicality 单独接近或超过六坐标 -> 整套装置对羞耻是过度设计的',
           at>=six*0.9, f"atypicality {100*at:.2f}% vs 六坐标 {100*six:.2f}%")
g.asserted('★ 另一支:六坐标明显更强且 atypicality 增量 ≈0 -> 装置抓到了朴素量抓不到的',
           six>at*1.5 and abs(both-six)<2*sd_at,
           f"六坐标 {100*six:.2f}% vs atypicality {100*at:.2f}%;"
           f"atypicality 增量 {100*(both-six):+.2f} ± {100*sd_at:.2f}pp")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

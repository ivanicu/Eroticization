import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A80 R315 -- 这六列是不是六个独立的测量

`#269b`:六个坐标全部由**同一次选项劈半**构造 —— 它们**共享测量误差**,
而这正是解衰减失效(修正后相关 1.794)的原因。
`#269b` 还指出:**分数层近乎正交(VIF ≤ 1.6)不等于测量误差正交**。

ESTIMAND        两次**独立**的选项劈半 A / B,各算一遍六个量:
                - 同名跨劈分 `corr(x_k^A, x_k^B)` = **信度**
                - 异名**跨**劈分 `corr(x_k^A, x_l^B)` = 不含共享误差的**真分数相关**
                - 异名**同**劈分 `corr(x_k^A, x_l^A)` = 真分数相关 **+ 共享误差**
                **共享误差痕迹 = 同劈分异名 − 跨劈分异名**,报一张 6×6 的表。
KILL            **若某些格子的痕迹明显大于零(> 2×跨种子展布)-> 这六列不是六个独立的测量,
                总表读法必须加这一句;若全部贴零 -> 六列独立,`#269b` 的推论要收窄。**
POSITIVE CTRL   守卫 14:这张表**必须随劈分种子变化**(否则它是构造出来的)。
NEGATIVE CTRL   把 B 劈分的人**跨人置换** -> 跨劈分相关必须塌到零,痕迹 = 同劈分相关本身。
IMPOSSIBLE      「共享误差」在这里 = 同一次抽签造成的共同扰动;
                它与「两个量本来就相关」在**同劈分**里不可分,**只有跨劈分能分**。
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
NAMES=['位置分 S','跨块对比 D','c1','c2','c3','剖面清晰度']
def quantities(seed):
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    posl=np.zeros(NN); cl=np.zeros(NN)
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        ha=M[:,o[:k]]; hb=M[:,o[k:2*k]]
        A[b,ppl]=ha.mean(1); B[b,ppl]=hb.mean(1)
        rr=-np.log(np.clip(ha.mean(0),1e-4,1.)); n=ha.sum(1)
        posl[ppl]+=np.where(n>0,(ha@rr)/np.maximum(n,1),0.0); cl[ppl]+=1
    Sl=np.where(cl>=8,posl/np.maximum(cl,1),np.nan)
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
    return [Sl,D]+cs+[st]
def cor(a,b,perm=False,seed=0):
    m=np.isfinite(a)&np.isfinite(b)&ok
    if m.sum()<500: return np.nan
    x,y=a[m],b[m]
    if perm: y=np.random.default_rng(seed).permutation(y)
    return float(np.corrcoef(x,y)[0,1])
def trace_matrix(sA,sB,perm=False):
    QA,QB=quantities(sA),quantities(sB); n=len(NAMES)
    same=np.array([[cor(QA[i],QA[j]) for j in range(n)] for i in range(n)])
    cross=np.array([[cor(QA[i],QB[j],perm,seed=sA+j) for j in range(n)] for i in range(n)])
    return same,cross,np.array([cor(QA[i],QB[i],perm,seed=sA+i) for i in range(n)])
S1,X1,rel1=trace_matrix(500,900)
print(f"信度(同名跨劈分):"+' · '.join(f"{NAMES[i]} {rel1[i]:+.3f}" for i in range(6)))
tr=S1-X1
print(f"\n共享误差痕迹 = 同劈分异名 − 跨劈分异名(对角无意义,置空):")
print(f"{'':<12}"+''.join(f"{n[:6]:>9}" for n in NAMES))
for i,n in enumerate(NAMES):
    print(f"{n[:11]:<12}"+''.join('       —' if i==j else f"{tr[i,j]:>+9.3f}" for j in range(6)))
off=np.array([tr[i,j] for i in range(6) for j in range(6) if i!=j])
print(f"\n非对角痕迹:均值 **{off.mean():+.4f}** · 最大 |{np.abs(off).max():.4f}| · "
      f"绝对值 > 0.05 的格子 **{int((np.abs(off)>0.05).sum())}/30**")
TR=[]
for s in range(4):
    a,b,_=trace_matrix(1000+7*s,1100+7*s)
    TR.append(np.array([(a-b)[i,j] for i in range(6) for j in range(6) if i!=j]))
sd=np.std(TR,0)
big=int(np.sum(np.abs(off)>2*sd))
print(f"跨 4 对劈分的展布:中位 ±{np.median(sd):.4f};"
      f"**明显大于零(>2×展布)的格子 {big}/30**")
Sp,Xp,_=trace_matrix(500,900,perm=True)
print(f"负对照(把 B 劈分的人跨人置换):跨劈分异名相关均值 "
      f"{np.mean([Xp[i,j] for i in range(6) for j in range(6) if i!=j]):+.4f}(必须 ≈0)")
T=pd.DataFrame(tr,index=NAMES,columns=NAMES).reset_index().rename(columns={'index':'v_row'})
check_columns(T,'R315'); T.to_csv(pathlib.Path(__file__).parent/'results'/'shared_error.csv',index=False)

g=Gate('这六列是不是六个独立的测量')
g.could_have_come_out_otherwise('⚠ 守卫 14:这张痕迹表必须随劈分种子变化',
    lambda s: float(np.mean(np.abs(TR[s]))) if s<len(TR) else float('nan'), list(range(4)))
g.negative_control('把 B 劈分的人跨人置换',
    abs(float(np.mean([Xp[i,j] for i in range(6) for j in range(6) if i!=j]))),
    abs(float(np.mean([X1[i,j] for i in range(6) for j in range(6) if i!=j]))),
    null_spread=float(np.std([Xp[i,j] for i in range(6) for j in range(6) if i!=j])),
    null_kind='跨人置换 B 劈分 —— 打掉配对,保留两边的分布')
g.asserted('★ 注册的 kill:有格子明显大于零 -> 这六列不是六个独立的测量',
           big>0, f"{big}/30 个格子 > 2×展布;最大 |{np.abs(off).max():.4f}|;"
                  f"均值 {off.mean():+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

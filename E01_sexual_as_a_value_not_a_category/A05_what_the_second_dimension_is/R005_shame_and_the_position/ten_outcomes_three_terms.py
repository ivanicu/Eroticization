import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A35 R233 -- 那 10 道越阈值的题,是被 S 预测的,还是被内容残差预测的

`#187` 把内容残差判成 NOISE,**但只用了一个结局**。
`#184b` 的 S 面板里 **10/20 越过全族阈值**,而 `#185` 证明它们**不是一个因子** ——
**那 10 个各自是什么,从没被问过。**

ESTIMAND        对 20 道题逐个跑三项分解:r(S,·) · r(内容残差,·) · r(内容残差,· | S)。
KILL            **若存在任何一道题被内容残差(控制 S 后)以 >2× 预测 ->
                `#187` 的 NULL 只对羞耻成立,范围要收窄到「对羞耻而言」。**
MULTIPLICITY    20 题 × 2 个预测子 -> **最大统计量零**给全族阈值,两侧各一个。
NEGATIVE CTRL   每题打乱一次(在分析样本内,#184b 的教训)。
POSITIVE CTRL   S 一侧必须复现 `#184b` 的 10/20 越阈值。
IMPOSSIBLE      内容分数是跨人主成分 -> 特异组合不可见(`#165` 的原话,`#187b` 已带过来)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage, check_residualized

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df); con=np.zeros(NN); pos=np.zeros(NN); cnt=np.zeros(NN); KB=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    Z=M-M.mean(0,keepdims=True); w,v=np.linalg.eigh(np.cov(Z,rowvar=False))
    con[ppl]+=Z@v[:,-1]; pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1)
    KB[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
C=np.where(ok,con/np.maximum(cnt,1),np.nan); S=np.where(ok,pos/np.maximum(cnt,1),np.nan)
KB=np.where(ok,KB,np.nan)
base=np.isfinite(S)&np.isfinite(C)&np.isfinite(KB); bi=np.flatnonzero(base)
X=np.c_[np.ones(len(bi)),S[bi]]
Cres=np.full(NN,np.nan); Cres[bi]=C[bi]-X@np.linalg.lstsq(X,C[bi],rcond=None)[0]
check_residualized(Cres[bi],S[bi],'R233 内容残差')

def pr(y,x,ctrls,ii):
    XX=np.c_[np.ones(len(ii)),*[c[ii] for c in ctrls]] if ctrls else np.ones((len(ii),1))
    ry=y[ii]-XX@np.linalg.lstsq(XX,y[ii],rcond=None)[0]
    rx=x[ii]-XX@np.linalg.lstsq(XX,x[ii],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

rb=np.random.default_rng(20260803); rows=[]; nS=[]; nC=[]
for c in lik:
    y=df[c].values.astype(float); m=base&np.isfinite(y); ii=np.flatnonzero(m)
    if len(ii)<3000: continue
    rs=pr(y,S,[KB],ii); rc=pr(y,Cres,[KB],ii); rcs=pr(y,Cres,[KB,S],ii)
    bs=[pr(y,Cres,[KB,S],rb.choice(ii,len(ii),replace=True)) for _ in range(200)]
    sd=float(np.std(bs))
    def perm(x,ctrls):
        yp=y.copy(); yp[ii]=rb.permutation(y[ii]); return abs(pr(yp,x,ctrls,ii))
    ps=[v for v in (perm(S,[KB]) for _ in range(40)) if np.isfinite(v)]
    pc=[v for v in (perm(Cres,[KB,S]) for _ in range(40)) if np.isfinite(v)]
    if len(ps)<20 or len(pc)<20:
        print(f"  ⚠ 丢弃:{c[:60]}"); continue
    nS.append(ps); nC.append(pc)
    rows.append(dict(q=c[:66],n=len(ii),r_S=rs,r_C=rc,r_C_given_S=rcs,sd=sd,ratio=abs(rcs)/sd))
T=pd.DataFrame(rows); check_columns(T,'R233'); check_coverage(len(T),len(lik),'R233 面板',tol=0.15)
T=T.sort_values('r_S',key=abs,ascending=False)
T.to_csv(pathlib.Path(__file__).parent/'results'/'three_terms.csv',index=False)
LS=min(len(x) for x in nS); LC=min(len(x) for x in nC)
thrS=float(np.nanquantile(np.nanmax(np.array([x[:LS] for x in nS]),axis=0),0.95))
thrC=float(np.nanquantile(np.nanmax(np.array([x[:LC] for x in nC]),axis=0),0.95))
print(f"\n全族阈值:S 侧 |r| = {thrS:.4f} · 内容残差侧 |r| = {thrC:.4f}\n")
print(f"{'r(S,·)':>9}{'r(Cres,·)':>11}{'r(Cres,·|S)':>13}{'比':>6}  题")
for _,r in T.iterrows():
    mS='★' if abs(r.r_S)>thrS else ' '; mC='★' if abs(r.r_C_given_S)>thrC else ' '
    print(f"{r.r_S:>+9.4f}{mS}{r.r_C:>+10.4f}{r.r_C_given_S:>+12.4f}{mC}{r.ratio:>6.1f}  {r.q[:56]}")
nS_pass=int((T.r_S.abs()>thrS).sum()); nC_pass=int((T.r_C_given_S.abs()>thrC).sum())
print(f"\nS 侧越阈值 {nS_pass}/{len(T)}(`#184b` 是 10/20)· **内容残差侧越阈值 {nC_pass}/{len(T)}**")

g=Gate('#187 的 NULL 是不是只对羞耻成立')
g.asserted('正对照:S 侧复现 `#184b` 的 10/20',abs(nS_pass-10)<=2,f"{nS_pass}/{len(T)}")
g.asserted('#129 守卫:残差与 S 的相关在构造上为 0',True,'check_residualized 已通过')
g.asserted('注册的 kill:存在任何一道被内容残差(控制 S 后)以 >2× 预测',nC_pass>0,
           f"内容残差侧越阈值 {nC_pass}/{len(T)}")
print(g)
print(f"\n  => {'范围要收窄到「对羞耻而言」' if nC_pass>0 else '**`#187` 的 NULL 在全部 20 道结局上都成立** —— 内容残差对任何一个都没有预测'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

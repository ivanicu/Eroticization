import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A83 R323 -- `form_i` 和羞耻

这条线上每一个人层量最后都被拿去问羞耻(`#179` 位置分 +0.1185 · `#230b` c3 +0.1286 ·
`#231` 两半反号 · `#235b` 跨仪器 +0.1384),**而 `#277a` 刚握住的 `form_i` 还没被问过。**

ESTIMAND        `corr(form_i, 羞耻)`;`form_i` 在六坐标之上的**增量 R²**(自助展布);
                以及它与 `#231` 的**越轨/普通两半**的关系。
⚠ 不循环          羞耻**不是** `form_i` 的构成题(构成题是 `animated`/`written`)。
⚠ 不是恒等式       `form_i` 与六坐标按构造正交,但它的**增量**不因此为零 ——
                它等于**羞耻残差上的那一份**,是一个真的可变的量(守卫 14 验证)。
KILL            **若 `form_i` 有独立于六坐标的羞耻增量 -> 这个无名的东西是通向羞耻的**第三条**路,
                `#232b`「两条路只在羞耻处相交」要扩写成三条;
                若增量为零 -> 它是一个真实但**与羞耻无关**的维度,
                而那会是这条线上**第一个不通向羞耻**的人层量,本身也值得写。**
POSITIVE CTRL   守卫 14;外加把一个**已知与羞耻相关**的量(位置分 S)当作 `form_i` 放进同一条管道,
                增量必须明显 > 0(证明这条管道测得出增量)。
NEGATIVE CTRL   跨人置换 `form_i`。
IMPOSSIBLE      横断面自报;`form_i` 只由两道题构成,信度 0.3382,所以任何零都是**衰减过的**。
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
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
rg=np.random.default_rng(500); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
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
    x,y_=Ra[mm,i],Rb[mm,i]
    if x.std()>1e-9 and y_.std()>1e-9: st[i]=float(np.corrcoef(x,y_)[0,1])
def z(v):
    m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
def hs(R,cols):
    sub=R[cols]; F2=np.isfinite(sub)
    return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
ORD_S=hs(Ra,ORD); TRG_S=hs(Ra,TRG)                       # #231 的两半
D=(z(TRG_S)-z(ORD_S)+z(hs(Rb,TRG))-z(hs(Rb,ORD)))/2
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
Q=[S,D]+cs+[st]
base=np.ones(NN,bool)&ok
for q_ in Q: base&=np.isfinite(q_)
X6=np.column_stack([np.ones(int(base.sum()))]+[(q_[base]-q_[base].mean())/q_[base].std() for q_ in Q])
def resid_col(col):
    y=pd.to_numeric(d[col],errors='coerce').values.astype(float)[base]
    f=np.isfinite(y); out=np.full(int(base.sum()),np.nan)
    b=np.linalg.lstsq(X6[f],(y[f]-y[f].mean())/y[f].std(),rcond=None)[0]
    out[f]=(y[f]-y[f].mean())/y[f].std()-X6[f]@b; return out
ra,rw=resid_col('animated'),resid_col('written')
FORM=np.full(NN,np.nan); idx=np.flatnonzero(base); f2=np.isfinite(ra)&np.isfinite(rw)
FORM[idx[f2]]=(ra[f2]+rw[f2])/2
SH=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SH],errors='coerce').values.astype(float)
m0=np.isfinite(FORM)&np.isfinite(y)&base
print(f"n = {int(m0.sum()):,}")
r_sh=float(np.corrcoef(FORM[m0],y[m0])[0,1])
rngB=np.random.default_rng(20260804)
sd_r=float(np.std([np.corrcoef(FORM[i],y[i])[0,1] for i in
    (rngB.choice(np.flatnonzero(m0),int(m0.sum()),True) for _ in range(300))]))
print(f"**`corr(form_i, 羞耻)` = {r_sh:+.4f} ± {sd_r:.4f}({abs(r_sh)/max(sd_r,1e-9):.1f}×)**")
print(f"   对照:位置分 +0.1185(`#179`)· c3 +0.1286(`#230b`)· D_起始 +0.1384(`#235b`)")
def r2_on(cols,idx_):
    X=np.column_stack([np.ones(len(idx_))]+[(c[idx_]-c[idx_].mean())/c[idx_].std() for c in cols])
    yy=(y[idx_]-y[idx_].mean())/y[idx_].std(); b=np.linalg.lstsq(X,yy,rcond=None)[0]
    return float(1-np.var(yy-X@b)/np.var(yy))
ii=np.flatnonzero(m0)
six=r2_on(Q,ii); both=r2_on(Q+[FORM],ii); inc=both-six
bs=[rngB.choice(ii,len(ii),True) for _ in range(200)]
sd_inc=float(np.std([r2_on(Q+[FORM],i)-r2_on(Q,i) for i in bs]))
print(f"\n六坐标 R² **{100*six:.2f}%** -> 加上 `form_i` **{100*both:.2f}%**;"
      f"**增量 {100*inc:+.2f} ± {100*sd_inc:.2f}pp({abs(inc)/max(sd_inc,1e-12):.1f}×)**")
# ⚠ 我原本写的正对照是「把位置分 S 放进同一条管道」—— **S 已在六坐标里**,
#    设计矩阵精确奇异,SVD 直接不收敛。**一个正对照如果用的是模型里已有的量,它连跑都跑不起来。**
#    改用两个**不在**六坐标里的量。
AT=np.full(NN,np.nan); rs=np.zeros(NN); np_=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); rs[ppl]+=M@rr; np_[ppl]+=M.sum(1)
AT=np.where(ok&(np_>0),rs/np.maximum(np_,1),np.nan)
# ⚠ `AT` 在 `ii` 上有缺失 -> 设计矩阵含 NaN -> SVD 不收敛。正对照要用**它自己的索引集**。
ii_at=np.flatnonzero(m0&np.isfinite(AT))
inc_at=r2_on(Q+[AT],ii_at)-r2_on(Q,ii_at)
print(f"   正对照(改用 atypicality,不在六坐标里但 `#274a` 已知增量 ≈0):"
      f"{100*inc_at:+.2f}pp(n={len(ii_at):,})")
# 另一端:一个**已知与羞耻强相关**且不在六坐标里的量 —— 用羞耻自己的带噪声复制
_rgp=np.random.default_rng(11); _m=np.isfinite(y)&m0
FAKE=np.full(NN,np.nan)
FAKE[np.flatnonzero(_m)]=0.7*((y[_m]-y[_m].mean())/y[_m].std())+0.7*_rgp.standard_normal(int(_m.sum()))
ii_fk=np.flatnonzero(m0&np.isfinite(FAKE))
inc_fk=r2_on(Q+[FAKE],ii_fk)-r2_on(Q,ii_fk)
print(f"   正对照(一个与羞耻强相关的假量,必须给出明显增量):**{100*inc_fk:+.2f}pp**")
# ⚠⚠ 又是 `#264b` 那个错:直接 permutation 整个数组会把 NaN 也搬位置,缺失模式被打乱。
#    **我在 `#264b` 亲手写下这条教训,而这里重犯了** —— 置换只在有限值之内进行。
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2)); z2[j]=z2[np.random.default_rng(seed).permutation(j)]
    return z2
nul=[r2_on(Q+[perm_finite(FORM,70+i)],ii)-six for i in range(30)]
print(f"负对照(跨人置换 `form_i`):增量 {100*np.mean(nul):+.3f} ± {100*np.std(nul):.3f}pp")
def cr(a,b):
    m=np.isfinite(a)&np.isfinite(b)&ok
    return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>400 else np.nan
print(f"\n与 `#231` 的两半:`corr(form_i, 普通半)` = **{cr(FORM,ORD_S):+.4f}** · "
      f"`corr(form_i, 越轨半)` = **{cr(FORM,TRG_S):+.4f}**")
T=pd.DataFrame([dict(quantity='form_i',r_shame=r_sh,sd_r=sd_r,inc_r2=inc,sd_inc=sd_inc,
                     r_ord=cr(FORM,ORD_S),r_trg=cr(FORM,TRG_S),n=int(m0.sum()))])
check_columns(T,'R323'); T.to_csv(pathlib.Path(__file__).parent/'results'/'form_shame.csv',index=False)

g=Gate('`form_i` 和羞耻')
g.could_have_come_out_otherwise('⚠ 守卫 14:`form_i` 的羞耻增量真的会变吗',
    lambda s: r2_on(Q+[FORM],rngB.choice(ii,len(ii),True))-r2_on(Q,rngB.choice(ii,len(ii),True)),
    [0,1,2,3])
g.has_error_bar('`corr(form_i, 羞耻)`',r_sh,sd_r,'bootstrap_人层')
g.has_error_bar('`form_i` 的增量 R²',inc,sd_inc,'bootstrap_人层')
g.negative_control('跨人置换 `form_i`',abs(float(np.mean(nul))),abs(inc),
                   null_spread=float(np.std(nul)),null_kind='跨人置换 —— 只打掉配对')
g.asserted('正对照两端:与羞耻强相关的假量必须给出明显增量;atypicality 必须接近零(`#274a`)',
           inc_fk>0.05 and abs(inc_at)<0.01, f"假量 {100*inc_fk:+.2f}pp · atypicality {100*inc_at:+.2f}pp")
g.asserted('★ 注册的 kill:`form_i` 有独立于六坐标的羞耻增量 -> 通向羞耻的第三条路',
           inc>2*sd_inc, f"增量 {100*inc:+.2f} ± {100*sd_inc:.2f}pp;"
           f"corr {r_sh:+.4f} ± {sd_r:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

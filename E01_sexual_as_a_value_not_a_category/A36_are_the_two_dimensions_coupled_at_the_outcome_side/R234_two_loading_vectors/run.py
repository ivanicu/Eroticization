import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A36 R234 -- 两个维度在人层正交,在**结局层**耦不耦合

`#188` 立起「两个正交维度,各自的外部相关物不同」——
但**正交性只在构造上成立**:`corr(S, Cres) = 0` 是回归的定义,不是发现(`#129`)。
没被检验的是:**两条外部相关向量,是不是也彼此独立。**

    ORTHOGONAL  corr(a, b) 落在零里 -> 「两个正交维度」在结局层也成立
    COUPLED     corr(a, b) 显著非零 -> 人层正交,**结局层耦合**,措辞要改

其中 a_j = r(S, y_j | 勾选数),b_j = r(Cres, y_j | 勾选数),j 跑遍 20 道 Likert。

ESTIMAND        corr(a, b),**单位是结局(n=20),不是人**(`#20` 的教训:估计量的 n 是它自己的单位数)。
KILL            条件式:先要**结局协方差对照开火**(见下),再判 corr(a,b)。
STRONGEST CONFOUND(跑之前写下)
                **20 道结局本身彼此相关。** 即使 S ⟂ Cres,两条向量也会被**同一个结局协方差**
                抹平而显得耦合。**这不是发现,是结局侧的几何。**
CONTROL(同一迭代内)
                用**两个随机正交的人层向量**过同一条管道 —— 结局协方差原封不动,只毁掉人层信号。
                若随机对也给出同样大的 corr(a,b),那么真实值不可读。
NEGATIVE CTRL   打乱配对(把 b 的结局标签打乱)-> 保留两条边际分布,毁掉配对。
NOISE FLOOR     人层 bootstrap 300 次,重算两条向量再算 corr。
IMPOSSIBLE      n=20 个结局。**任何 |corr| < 0.44 在 n=20 上都不显著**(粗略),
                所以本设计只能判**大的**耦合。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_residualized

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
X0=np.c_[np.ones(len(bi)),S[bi]]
Cres=np.full(NN,np.nan); Cres[bi]=C[bi]-X0@np.linalg.lstsq(X0,C[bi],rcond=None)[0]
check_residualized(Cres[bi],S[bi],'R234 内容残差')
Y=df[lik].values.astype(float)
print(f"结局 {len(lik)} 道;S/Cres 有效 {len(bi):,} 人")

def vecs(x1,x2,ii):
    a=[];b=[]
    for j in range(Y.shape[1]):
        y=Y[:,j]; m=np.isfinite(y[ii])
        jj=ii[m]
        XX=np.c_[np.ones(len(jj)),KB[jj]]
        ry=y[jj]-XX@np.linalg.lstsq(XX,y[jj],rcond=None)[0]
        r1=x1[jj]-XX@np.linalg.lstsq(XX,x1[jj],rcond=None)[0]
        r2=x2[jj]-XX@np.linalg.lstsq(XX,x2[jj],rcond=None)[0]
        a.append(np.corrcoef(ry,r1)[0,1]); b.append(np.corrcoef(ry,r2)[0,1])
    return np.array(a),np.array(b)

a,b=vecs(S,Cres,bi)
rho=float(np.corrcoef(a,b)[0,1])
print(f"\ncorr(a, b) 跨 {len(a)} 道结局 = **{rho:+.4f}**")
print(f"  a(S 侧)范围 {a.min():+.4f}..{a.max():+.4f}   b(Cres 侧)范围 {b.min():+.4f}..{b.max():+.4f}")

rb=np.random.default_rng(20260803)
# 负对照:打乱配对
pair_null=[float(np.corrcoef(a,rb.permutation(b))[0,1]) for _ in range(2000)]
# 混杂对照:两个随机正交的人层向量,结局协方差不变
ctrl=[]
for _ in range(60):
    u1=rb.standard_normal(NN); u2=rb.standard_normal(NN)
    u2b=u2[bi]-np.c_[np.ones(len(bi)),u1[bi]]@np.linalg.lstsq(np.c_[np.ones(len(bi)),u1[bi]],u2[bi],rcond=None)[0]
    u2f=np.full(NN,np.nan); u2f[bi]=u2b
    aa,bb=vecs(u1,u2f,bi); ctrl.append(float(np.corrcoef(aa,bb)[0,1]))
# 人层 bootstrap
boot=[]
for _ in range(120):
    s_=rb.choice(bi,len(bi),replace=True)
    aa,bb=vecs(S,Cres,s_); boot.append(float(np.corrcoef(aa,bb)[0,1]))
print(f"\n负对照(打乱配对,2000 次)   {np.mean(pair_null):+.4f} ± {np.std(pair_null):.4f}")
print(f"混杂对照(随机正交对,60 次) {np.mean(ctrl):+.4f} ± {np.std(ctrl):.4f}   <- 结局协方差原封不动")
print(f"人层 bootstrap(120 次)     {np.mean(boot):+.4f} ± {np.std(boot):.4f}")
loo=[float(np.corrcoef(np.delete(a,k),np.delete(b,k))[0,1]) for k in range(len(a))]
print(f"留一(20 道)               {min(loo):+.4f} .. {max(loo):+.4f}")

T=pd.DataFrame(dict(q=[c[:60] for c in lik],a_S=a,b_Cres=b))
check_columns(T,'R234'); T.to_csv(pathlib.Path(__file__).parent/'results'/'vectors.csv',index=False)
g=Gate('两个维度在结局层耦不耦合')
g.asserted('#129 守卫:人层正交是构造出来的,已确认',True,'check_residualized 通过')
g.offset_control('corr(a,b) vs 随机正交对',rho,float(np.mean(ctrl)),float(np.std(boot)),
                 null_kind='两个随机正交人层向量过同一条管道 —— **结局协方差原封不动**,'
                           '只毁掉人层信号;这不是零假设,是结局侧几何的基线')
g.negative_control('打乱配对',float(abs(np.mean(pair_null))),rho,null_spread=float(np.std(pair_null)))
g.resolvable('corr(a,b) 本身',rho,float(np.std(boot)))
g.no_sign_crossing('留一 20 次同号',loo)
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- ⚠ 跑完才想到的机械解释,必须当场查 ---------------------------------------
# b_j = r(Cres, y_j) = [c_j − ρ_CS·a_j] / √(1−ρ_CS²),其中 c_j = r(C, y_j)、ρ_CS = corr(C,S)。
# **只要 ρ_CS ≠ 0,残差化就在构造上给 b 灌进一个对 a 的负依赖。**
# 我上面的"随机正交对"对照抓不到它 —— 那一对本来 ρ=0,所以没有可灌的东西。
# 正确的判据是**原始**内容向量 c 与 a 的相关:c 没有被 S 残差化过。
print("\n---- 机械解释:残差化本身会不会造出这个负相关 ----")
rho_CS=float(np.corrcoef(C[bi],S[bi])[0,1])
c_raw=[]
for j in range(Y.shape[1]):
    y=Y[:,j]; m=np.isfinite(y[bi]); jj=bi[m]
    XX=np.c_[np.ones(len(jj)),KB[jj]]
    ry=y[jj]-XX@np.linalg.lstsq(XX,y[jj],rcond=None)[0]
    rc=C[jj]-XX@np.linalg.lstsq(XX,C[jj],rcond=None)[0]
    c_raw.append(np.corrcoef(ry,rc)[0,1])
c_raw=np.array(c_raw)
rho_raw=float(np.corrcoef(a,c_raw)[0,1])
pred=-rho_CS/np.sqrt(max(1-rho_CS**2,1e-9))         # 残差化灌进去的斜率
print(f"  corr(C, S) 人层 = {rho_CS:+.4f}")
print(f"  corr(a, c_raw)  = {rho_raw:+.4f}   <- **c 没有被残差化过,这一格才是判据**")
print(f"  corr(a, b)      = {rho:+.4f}")
print(f"  残差化灌进去的斜率 −ρ_CS/√(1−ρ_CS²) = {pred:+.4f}")
bootr=[]
for _ in range(120):
    s_=rb.choice(bi,len(bi),replace=True)
    aa,_=vecs(S,Cres,s_)
    cc=[]
    for j in range(Y.shape[1]):
        y=Y[:,j]; m=np.isfinite(y[s_]); jj=s_[m]
        XX=np.c_[np.ones(len(jj)),KB[jj]]
        ry=y[jj]-XX@np.linalg.lstsq(XX,y[jj],rcond=None)[0]
        rc=C[jj]-XX@np.linalg.lstsq(XX,C[jj],rcond=None)[0]
        cc.append(np.corrcoef(ry,rc)[0,1])
    bootr.append(float(np.corrcoef(aa,np.array(cc))[0,1]))
sd_raw=float(np.std(bootr))
print(f"  corr(a, c_raw) 的人层 bootstrap sd = {sd_raw:.4f}  -> {abs(rho_raw)/sd_raw:.1f}×")

g2=Gate('这个负耦合是真的,还是残差化造的')
g2.asserted('可判前提:corr(C,S) 确实非零,机械解释确实在场',abs(rho_CS)>0.05,
            f"corr(C,S) = {rho_CS:+.4f} -> 残差化会灌进斜率 {pred:+.4f}")
g2.resolvable('corr(a, c_raw)(未残差化)',rho_raw,sd_raw)
g2.asserted('判据:未残差化的版本仍为负',rho_raw<0,f"{rho_raw:+.4f}")
g2.offset_control('corr(a,c_raw) vs 残差化灌进的斜率',rho_raw,pred,sd_raw,
                  null_kind='残差化在构造上给 b 灌进的对 a 的依赖 —— 不是零假设,是纯代数的量')
print(g2)
print(f"\n  => {'COUPLED(真的)' if (rho_raw<0 and abs(rho_raw)>2*sd_raw) else 'ARTIFACT / UNVERIFIED'}")

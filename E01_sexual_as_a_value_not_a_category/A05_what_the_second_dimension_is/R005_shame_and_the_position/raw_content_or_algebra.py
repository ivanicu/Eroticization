import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A36 R235 -- `#188` 的 6/20,是内容侧的发现,还是 `#189b` 那个代数

`#189b`:`b_j = [c_j − ρ_CS·a_j]/√(1−ρ_CS²)`,而 `ρ_CS = +0.4323`。
所以只要 `a_j` 大,`b_j` 就会被推向负 —— **一个 `c_j` 本身为零的结局也能有一个大的 `b_j`。**
`#188` 的 6/20 全部是 `b_j`,**没有一个是 `c_j`。**

    GENUINE   `c_j`(未残差化)本身可分辨、且与 `b_j` 同号 -> 内容侧真的预测它
    ALGEBRA   `c_j` 不可分辨而 `b_j` 大 -> 那个格子是 `−ρ_CS·a_j` 项造的

ESTIMAND        对 20 道结局同时给 `a_j = r(S,·)`、`c_j = r(C,·)`、`b_j = r(Cres,·)`,
                **各自的全族阈值**(最大统计量零),逐格判 GENUINE / ALGEBRA。
KILL            **若 `#188` 的 6 道里有 ≥4 道判 ALGEBRA -> `#188` 的 6/20 要重判。**
POSITIVE CTRL   `animated`:`#233` 已给 `c = −0.1259`,比 `b = −0.1063` **还大** ->
                它必须判 GENUINE。判不到就说明这一轮的 c 侧管道坏了。
NEGATIVE CTRL   构造一个 `c_j ≡ 0` 而 `a_j` 大的合成结局(= S + 噪声,与 C 无关的部分),
                它必须判 ALGEBRA。
IMPOSSIBLE      `C` 是跨人主成分 -> 特异组合不可见(`#165` 原话,一路带下来)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_residualized, check_coverage

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
check_residualized(Cres[bi],S[bi],'R235 内容残差')
rho_CS=float(np.corrcoef(C[bi],S[bi])[0,1])
print(f"corr(C,S) = {rho_CS:+.4f};n = {len(bi):,}")

# 负对照结局:与 C 无关、只由 S 驱动的合成题
rb=np.random.default_rng(20260803)
synth=np.full(NN,np.nan); synth[bi]=S[bi]+rb.standard_normal(len(bi))*S[bi].std()*1.5
OUT=[(c,df[c].values.astype(float)) for c in lik]+[('【负对照】纯 S 驱动的合成题',synth)]

def rr_(y,x,ii):
    m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
    XX=np.c_[np.ones(len(jj)),KB[jj]]
    ry=y[jj]-XX@np.linalg.lstsq(XX,y[jj],rcond=None)[0]
    rx=x[jj]-XX@np.linalg.lstsq(XX,x[jj],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1]), len(jj)

rows=[]; nulls={'a':[],'c':[],'b':[]}
for name,y in OUT:
    a,n1=rr_(y,S,bi); c,_=rr_(y,C,bi); b,_=rr_(y,Cres,bi)
    sd=float(np.std([rr_(y,C,rb.choice(bi,len(bi),replace=True))[0] for _ in range(200)]))
    for k,x in (('a',S),('c',C),('b',Cres)):
        ps=[]
        for _ in range(40):
            yp=y.copy(); yp[bi]=rb.permutation(y[bi])
            v,_=rr_(yp,x,bi)
            if np.isfinite(v): ps.append(abs(v))
        if len(ps)>=20: nulls[k].append(ps)
    rows.append(dict(q=name[:60],n=n1,a_S=a,c_C=c,b_Cres=b,sd_c=sd))
T=pd.DataFrame(rows); check_columns(T,'R235')
check_coverage(len(T),len(OUT),'R235 面板',tol=0.10)
thr={}
for k in nulls:
    L=min(len(x) for x in nulls[k])
    thr[k]=float(np.nanquantile(np.nanmax(np.array([x[:L] for x in nulls[k]]),axis=0),0.95))
print(f"\n全族阈值:a(S) {thr['a']:.4f} · **c(C 原始) {thr['c']:.4f}** · b(Cres) {thr['b']:.4f}\n")

six=['animated','biomale','written','submissive','allrollidentity','biological *female*']
def is_six(q): return any(k.lower() in q.lower() for k in six)
T['verdict']=np.where(T.b_Cres.abs()>thr['b'],
                      np.where((T.c_C.abs()>thr['c'])&(np.sign(T.c_C)==np.sign(T.b_Cres)),'GENUINE','ALGEBRA'),'—')
T.to_csv(pathlib.Path(__file__).parent/'results'/'abc.csv',index=False)
print(f"{'a(S)':>9}{'c(C)':>10}{'b(Cres)':>10}  {'判':<9} 题")
for _,r in T.sort_values('b_Cres',key=abs,ascending=False).iterrows():
    if abs(r.b_Cres)>thr['b'] or '负对照' in r.q or 'ashamed' in r.q:
        print(f"{r.a_S:>+9.4f}{r.c_C:>+10.4f}{r.b_Cres:>+10.4f}  {r.verdict:<9} {r.q[:52]}")
SIX=T[(T.b_Cres.abs()>thr['b'])&(~T.q.str.contains('负对照'))]
n_alg=int((SIX.verdict=='ALGEBRA').sum())
print(f"\n越 b 阈值的真实结局 {len(SIX)} 道,其中 ALGEBRA {n_alg} 道")

ani=T[T.q.str.contains('animated')].iloc[0]; neg=T[T.q.str.contains('负对照')].iloc[0]
g=Gate('#188 的 6/20 是发现还是代数')
g.asserted('正对照:animated 的 c 比 b 还大 -> 必须判 GENUINE',ani.verdict=='GENUINE',
           f"c={ani.c_C:+.4f} b={ani.b_Cres:+.4f} -> {ani.verdict}")
g.asserted('负对照:纯 S 驱动的合成题必须判 ALGEBRA(或不越阈值)',
           neg.verdict in ('ALGEBRA','—'),
           f"a={neg.a_S:+.4f} c={neg.c_C:+.4f} b={neg.b_Cres:+.4f} -> {neg.verdict}")
g.asserted('注册的 kill:6 道里 >=4 道判 ALGEBRA -> `#188` 的 6/20 要重判',n_alg>=4,
           f"ALGEBRA {n_alg}/{len(SIX)}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

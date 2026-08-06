import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A38 R238 -- 勾选数自己预测什么

`#192b`:「位置贴羞耻、内容贴媒介」的干净分工**只在给定勾选数之后才出现**。
而**勾选数从来只是一个控制项** —— `#104` 证过位置倾向去掉它后保留 67%,
**但没人问过它自己预测什么。**

ESTIMAND        把勾选数 K 当**主变量**跑 20 道 Likert 面板(最大统计量零给全族阈值);
                再把 `r(S, 羞耻)` 分解成 K 能解释的部分与不能的部分。
KILL            **若 K 自己就强预测羞耻(越全族阈值,且 |r| 达到 S 的一半以上)->
                「位置↔羞耻」要重新分解,它可能有一大半是"勾得多的人更容易有一样让自己羞耻的东西"。**
NEGATIVE CTRL   每题在分析样本内打乱(`#184b` 的教训)。
POSITIVE CTRL   把 S 当主变量跑同一条面板,必须复现 `#184b` 的 羞耻 排名 1。
DECOMPOSITION   `r(S,羞耻)` raw · 给定 K · 而 `r(K,羞耻)` 给定 S —— 三个数一起报。
IMPOSSIBLE      K 与 S 在人层相关(`#100` 已报 +0.608,而零是 +0.719)——
                两者的"各自贡献"在共线下不可分,只能报三个数,不得声称分开了(`#182b` 的教训)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
SHAME=next(c for c in lik if 'ashamed' in c)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df); pos=np.zeros(NN); cnt=np.zeros(NN); K=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); K[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
S=np.where(ok,pos/np.maximum(cnt,1),np.nan); K=np.where(ok,K,np.nan)
base=np.isfinite(S)&np.isfinite(K); bi=np.flatnonzero(base)
rho_SK=float(np.corrcoef(S[bi],K[bi])[0,1])
print(f"n = {len(bi):,};corr(S, K) = {rho_SK:+.4f}   (`#100` 报的是 +0.608,零 +0.719)")

def mr(y,x,ii,ctrls=()):
    m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
    XX=np.c_[np.ones(len(jj)),*[c[jj] for c in ctrls]] if ctrls else np.ones((len(jj),1))
    ry=y[jj]-XX@np.linalg.lstsq(XX,y[jj],rcond=None)[0]
    rx=x[jj]-XX@np.linalg.lstsq(XX,x[jj],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1]), len(jj)

rb=np.random.default_rng(20260803); rows=[]; nK=[]; nS=[]
for c in lik:
    y=df[c].values.astype(float)
    rK,n=mr(y,K,bi); rS,_=mr(y,S,bi)
    sd=float(np.std([mr(y,K,rb.choice(bi,len(bi),replace=True))[0] for _ in range(200)]))
    for store,x in ((nK,K),(nS,S)):
        ps=[]
        for _ in range(40):
            yp=y.copy(); yp[bi]=rb.permutation(y[bi]); v,_=mr(yp,x,bi)
            if np.isfinite(v): ps.append(abs(v))
        if len(ps)>=20: store.append(ps)
    rows.append(dict(q=c[:62],n=n,r_K=rK,r_S=rS,sd=sd,ratio=abs(rK)/sd))
T=pd.DataFrame(rows); check_columns(T,'R238'); check_coverage(len(T),len(lik),'R238 面板',tol=0.10)
thr=lambda nl:(lambda L: float(np.nanquantile(np.nanmax(np.array([x[:L] for x in nl]),axis=0),0.95)))(min(len(x) for x in nl))
thrK,thrS=thr(nK),thr(nS)
T=T.sort_values('r_K',key=abs,ascending=False)
T.to_csv(pathlib.Path(__file__).parent/'results'/'breadth_panel.csv',index=False)
print(f"\n全族阈值:K 侧 |r| = {thrK:.4f} · S 侧 |r| = {thrS:.4f}\n")
print(f"{'r(K,·)':>9}{'r(S,·)':>10}{'比':>7}  题")
for _,r in T.head(10).iterrows():
    print(f"{r.r_K:>+9.4f}{'★' if abs(r.r_K)>thrK else ' '}{r.r_S:>+9.4f}{'★' if abs(r.r_S)>thrS else ' '}{r.ratio:>6.1f}  {r.q[:56]}")
nK_pass=int((T.r_K.abs()>thrK).sum()); nS_pass=int((T.r_S.abs()>thrS).sum())
sh_K=float(T[T.q.str.contains('ashamed')].r_K.iloc[0]); sh_S=float(T[T.q.str.contains('ashamed')].r_S.iloc[0])
rk_S=int(T.reset_index().index[T.reset_index().q.str.contains('ashamed')][0])+1
print(f"\nK 侧越阈值 {nK_pass}/{len(T)} · S 侧 {nS_pass}/{len(T)}")
print(f"羞耻:r(K,·) = {sh_K:+.4f}(K 侧排名 {rk_S})· r(S,·) = {sh_S:+.4f}")

# 三项分解(共线,不得声称分开)
y=df[SHAME].values.astype(float)
a_raw,_=mr(y,S,bi); a_gK,_=mr(y,S,bi,(K,))
k_raw,_=mr(y,K,bi); k_gS,_=mr(y,K,bi,(S,))
print(f"\n三项分解(共线 corr(S,K) = {rho_SK:+.4f},**不得声称分开了**):")
print(f"  r(S,羞耻) raw {a_raw:+.4f}  ->  给定 K {a_gK:+.4f}")
print(f"  r(K,羞耻) raw {k_raw:+.4f}  ->  给定 S {k_gS:+.4f}")

g=Gate('勾选数自己预测什么')
g.asserted('正对照:S 侧羞耻仍是最强之一(复现 `#184b`)',abs(sh_S)>thrS,f"{sh_S:+.4f} vs 阈值 {thrS:.4f}")
g.asserted('可判前提:多数题在 K 侧为零(否则 K 对什么都出结果)',nK_pass<=len(T)//2,
           f"K 侧越阈值 {nK_pass}/{len(T)}")
g.threshold_outside_noise('r(K, 羞耻) vs 全族阈值',abs(sh_K),thrK,float(T[T.q.str.contains('ashamed')].sd.iloc[0]))
g.asserted('注册的 kill:K 自己越阈值且达到 S 的一半以上',
           (abs(sh_K)>thrK) and (abs(sh_K)>=0.5*abs(sh_S)),
           f"|r(K,羞耻)| = {abs(sh_K):.4f} · 阈值 {thrK:.4f} · S 的一半 = {0.5*abs(sh_S):.4f}")
g.asserted('共线性已量化,三项分解不得读成"各自贡献"',True,
           f"corr(S,K) = {rho_SK:+.4f};S 给定 K {a_gK:+.4f} · K 给定 S {k_gS:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A33 R230 -- S 预测的是羞耻,还是一个「广义性态度」因子

`#184b`:S↔羞耻 排名 1/20 并越过全族阈值,**但 10/20 题都越过**。
所以「羞耻是最大的那一个」成立,而「S 预测的是羞耻」还没被检验。

ESTIMAND        对 S 面板里越阈值的那些题做主成分;判 **S↔羞耻 在控制 PC1 之后还剩多少**。
KILL            条件式:先要 **PC1 确实是一个共同因子**(解释率 > 30%,且羞耻在它上面有载荷);
                再判:**控制 PC1 后 S↔羞耻 掉到 2× 以下 -> `#179` 的措辞要从「羞耻」改成那个因子。**
NEGATIVE CTRL   用**没越阈值**的那些题做同样的 PC1,控制它应当**几乎不动** S↔羞耻。
POSITIVE CTRL   控制羞耻本身,S↔羞耻 必须归零(管道自检)。
IMPOSSIBLE      PC1 由这些题**自己**定义,而羞耻是其中之一 -> 控制 PC1 会**部分地控制羞耻自己**。
                所以"掉多少"有一个**结构性下界**;本轮同时报「把羞耻排除在 PC1 之外」的版本。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
P=pd.read_csv(pathlib.Path(__file__).resolve().parents[1]/'R229_discriminant_panel'/'results'/'panel_S.csv')
thr=0.0393                                    # `#184b` 的全族阈值,钉住
top=[q for q in P[P.r.abs()>thr].q]; rest=[q for q in P[P.r.abs()<=thr].q]
print(f"越阈值 {len(top)} 题 · 未越 {len(rest)} 题")

lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
def match(short): return next((c for c in lik if c[:70]==short), None)
TOP=[match(q) for q in top]; TOP=[c for c in TOP if c]
REST=[match(q) for q in rest]; REST=[c for c in REST if c]
SHAME=next(c for c in lik if 'ashamed' in c)
print(f"匹配回原列名:越阈值 {len(TOP)} · 未越 {len(REST)}")

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df); pos=np.zeros(NN); cnt=np.zeros(NN); KB=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); KB[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8; S=np.where(ok,pos/np.maximum(cnt,1),np.nan); KB=np.where(ok,KB,np.nan)

def pc1(cols, mask):
    Y=df[cols].values.astype(float)[mask]
    Y=np.where(np.isfinite(Y),Y,np.nanmean(Y,axis=0))
    Z=(Y-Y.mean(0))/np.maximum(Y.std(0),1e-9)
    w,v=np.linalg.eigh(np.cov(Z,rowvar=False))
    ev=w[::-1]/w.sum(); return Z@v[:,-1], float(ev[0]), v[:,-1]

def pr(y,x,ctrls,idx):
    X=np.c_[np.ones(len(idx)),*[c[idx] for c in ctrls]] if ctrls else np.ones((len(idx),1))
    ry=y[idx]-X@np.linalg.lstsq(X,y[idx],rcond=None)[0]
    rx=x[idx]-X@np.linalg.lstsq(X,x[idx],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

sh=df[SHAME].values.astype(float)
m=np.isfinite(sh)&np.isfinite(S)&np.isfinite(KB); idx=np.flatnonzero(m)
F_top,ev_top,load_top=pc1(TOP,m); F_rest,ev_rest,_=pc1(REST,m)
F_ex,ev_ex,_=pc1([c for c in TOP if c!=SHAME],m)          # 把羞耻排除在 PC1 之外
Ft=np.full(NN,np.nan); Ft[idx]=F_top
Fr=np.full(NN,np.nan); Fr[idx]=F_rest
Fe=np.full(NN,np.nan); Fe[idx]=F_ex
li=TOP.index(SHAME); print(f"\nPC1(越阈值 {len(TOP)} 题)解释率 {100*ev_top:.1f}%,羞耻载荷 {load_top[li]:+.3f}")
print(f"PC1(未越阈值 {len(REST)} 题)解释率 {100*ev_rest:.1f}%")
print(f"PC1(越阈值但**排除羞耻**)解释率 {100*ev_ex:.1f}%")

rb=np.random.default_rng(20260803); rows=[]
for name,ctrls in [('raw',[KB]),('+PC1(含羞耻)',[KB,Ft]),('+PC1(排除羞耻)',[KB,Fe]),
                   ('+PC1(未越阈值)',[KB,Fr]),('+羞耻自身',[KB,sh])]:
    r=pr(sh,S,ctrls,idx)
    bs=[pr(sh,S,ctrls,rb.choice(idx,len(idx),replace=True)) for _ in range(300)]
    sd=float(np.std(bs)); rows.append(dict(model=name,r=r,sd=sd,ratio=abs(r)/sd if sd>0 else np.nan))
T=pd.DataFrame(rows); check_columns(T,'R230'); T.to_csv(pathlib.Path(__file__).parent/'results'/'factor.csv',index=False)
print(f"\n{'模型':<18}{'r(S, 羞耻)':>12}{'sd':>9}{'比':>7}")
for _,r in T.iterrows(): print(f"{r.model:<18}{r.r:>+12.4f}{r.sd:>9.4f}{r.ratio:>7.1f}")

raw=float(T[T.model=='raw'].r.iloc[0])
ex =float(T[T.model=='+PC1(排除羞耻)'].r.iloc[0]); ex_sd=float(T[T.model=='+PC1(排除羞耻)'].sd.iloc[0])
g=Gate('S 预测的是羞耻还是一个广义因子')
g.asserted('可判前提:PC1 确实是一个共同因子(解释率 >30%)',ev_ex>0.30,f"{100*ev_ex:.1f}%")
g.asserted('正对照:控制羞耻自身,S↔羞耻 必须归零',
           abs(float(T[T.model=='+羞耻自身'].r.iloc[0]))<0.02,
           f"{float(T[T.model=='+羞耻自身'].r.iloc[0]):+.6f}")
g.asserted('负对照:控制**未越阈值**题的 PC1 应几乎不动',
           abs(float(T[T.model=='+PC1(未越阈值)'].r.iloc[0])-raw)<0.02,
           f"{float(T[T.model=='+PC1(未越阈值)'].r.iloc[0]):+.4f} vs raw {raw:+.4f}")
g.resolvable('控制排除羞耻的 PC1 之后,S↔羞耻',ex,ex_sd)
g.asserted('注册的 kill:控制 PC1 后掉到 2× 以下 -> 措辞要改',abs(ex)/ex_sd<2,
           f"{abs(ex)/ex_sd:.1f}×;保留 {100*ex/raw:.0f}%")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

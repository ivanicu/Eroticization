import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A33 R229 -- 判别效度:羞耻在 20 道 Likert 里排第几

`#183` 的 NEXT 注册的是:找 2–3 个**非性的负性自我评价**题跑同一条管道,
若它们与羞耻的 +0.027 无法区分,「性的羞耻」要降级成「负性自我评价」。

⚠ **这个检验做不了。** 这份 release 的 −3..+3 Likert 只有 **20 道**,
逐条读过:**没有一道是非性的负性自我评价**(`biomale` / `highenergy` 是一般特质,
其余全是性的)。**不拿代理冒充** —— 注册的检验判 **STRUCTURALLY IMPOSSIBLE**。

**替代设计(明说是替代)**:把同一条管道跑遍**全部 20 道**,给出羞耻的**排名**与整条谱。
它答不了"是不是一般性自我否定",但能答:**这条关联是不是对任意 Likert 题都出结果。**

ESTIMAND        r(截距, 每一道 Likert),同一批人、同一条管道、同一个控制(类别数)。
KILL            条件式:先要**多数题目为零**(否则管道对什么都出结果 -> 整轮不可读);
                再判:**羞耻若不在前三,那么它不特殊。**
NEGATIVE CTRL   每题各自打乱一次。
MULTIPLICITY    20 题 -> 用**最大统计量零**给全族阈值(而不是逐题 2×)。
IMPOSSIBLE      见上:非性的负性自我评价题**不存在于本数据**。
                所以「性的羞耻 vs 一般性自我否定」这个分岔,**在这份 release 上永远判不了**。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
print(f"−3..+3 Likert 题 {len(lik)} 道")
O=pd.read_csv('data/derived/onset.csv')
onset=[c for c in O.columns if re.search(r'How old were you when you first',c)]
A_=O[onset].apply(pd.to_numeric,errors='coerce').values
A_=np.where((A_>=2)&(A_<=60),A_,np.nan)
assert np.isfinite(A_).sum()>10000
have=np.isfinite(A_); rar=-np.log(np.clip(have.mean(0),1e-4,1.)); K=have.sum(1).astype(float)
n_=A_.shape[0]; I0=np.full(n_,np.nan)
for i in range(n_):
    idx=np.flatnonzero(np.isfinite(A_[i]))
    if len(idx)<6: continue
    x=A_[i,idx].astype(float)
    if x.std()<1e-9: continue
    I0[i]=np.linalg.lstsq(np.c_[np.ones(len(x)),x],rar[idx],rcond=None)[0][0]

def pr(y,x,idx):
    X=np.c_[np.ones(len(idx)),K[idx]]
    ry=y[idx]-X@np.linalg.lstsq(X,y[idx],rcond=None)[0]
    rx=x[idx]-X@np.linalg.lstsq(X,x[idx],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

rb=np.random.default_rng(20260803); rows=[]; nulls=[]
for c in lik:
    y=df[c].values.astype(float)
    m=np.isfinite(y)&np.isfinite(I0)&np.isfinite(K); idx=np.flatnonzero(m)
    if len(idx)<3000: continue
    r=pr(y,I0,idx)
    bs=[pr(y,I0,rb.choice(idx,len(idx),replace=True)) for _ in range(200)]
    # ⚠ 打乱后某些题可能给出 nan(子样本里取值退化)-> 最大统计量必须 nan-safe,
    #   否则一个 nan 让整个全族阈值变成 nan,而下游只会静静地不标 ★。
    # ⚠ #184b:第一版打乱的是**整条数组**(含 NaN)—— 打乱后 NaN 落进分析样本,
    #   相关变 nan,于是 20 道题里 17 道被 `continue` 静默丢掉,面板只剩 3 道。
    #   **打乱必须在分析样本内做。** 静默截断吃掉了 85% 的面板(#118c 又一次)。
    def perm_r():
        yp=y.copy(); yp[idx]=rb.permutation(y[idx])
        return abs(pr(yp,I0,idx))
    per=[v for v in (perm_r() for _ in range(50)) if np.isfinite(v)]
    if len(per)<25:
        print(f"  ⚠ 丢弃(打乱后有效次数 {len(per)}/50):{c[:60]}"); continue
    nulls.append(per)
    rows.append(dict(q=c[:70],n=len(idx),r=r,sd=float(np.std(bs)),ratio=abs(r)/float(np.std(bs))))
# ⚠ #184a:第一版列名叫 `item` —— pandas Series 自带 `.item()`,`r.q[:74]` 取到的是方法。
#   与 `#156` 的 `G.shift`、`#166c` 的 `T.cov` **同一类,第四次**。而 `check_columns`
#   就是为这一类写的,我**第三次**没有调用它。现在调用。
T=pd.DataFrame(rows); check_columns(T,'R229 面板'); T=T.sort_values('r',key=abs,ascending=False)
T.to_csv(pathlib.Path(__file__).parent/'results'/'panel.csv',index=False)

# 最大统计量零:每一轮 permutation 取 20 题里的最大 |r|
L=min(len(x) for x in nulls); Nl=np.array([x[:L] for x in nulls])
maxstat=np.nanmax(Nl,axis=0)
thr=float(np.nanquantile(maxstat,0.95))
assert np.isfinite(thr), '全族阈值是 nan —— 不得继续'
print(f"\n最大统计量零(20 题,50 次):全族 95% 阈值 |r| = {thr:.4f}\n")
print(f"{'排名':<5}{'|r|':>8}{'比':>7}  题目")
for k,(_,r) in enumerate(T.iterrows(),1):
    mark='★' if abs(r.r)>thr else ' '
    print(f"{k:<5}{r.r:>+8.4f}{r.ratio:>7.1f} {mark} {r.q[:74]}")
sh_row=T[T.q.str.contains('ashamed')].iloc[0]
rank=int(T.reset_index().index[T.reset_index().q.str.contains('ashamed')][0])+1
n_pass=int((T.r.abs()>thr).sum())
from lib.gates import check_coverage
check_coverage(len(T),len(lik),'R229 面板覆盖',tol=0.15)
print(f"\n羞耻排名 {rank}/{len(T)};超过全族阈值的题 {n_pass}/{len(T)}")

g=Gate('羞耻在这条谱里特不特殊')
g.asserted('注册的检验 STRUCTURALLY IMPOSSIBLE —— 明说,不用代理冒充',True,
           '本 release 的 20 道 −3..+3 Likert 里没有一道是**非性的负性自我评价**;'
           '「性的羞耻 vs 一般性自我否定」在这份数据上永远判不了')
g.asserted('可判前提:多数题为零(否则管道对什么都出结果)',n_pass<=len(T)//2,
           f"超阈值 {n_pass}/{len(T)}")
g.threshold_outside_noise('羞耻 vs 全族阈值',float(abs(sh_row.r)),thr,float(sh_row.sd))
g.asserted('注册的 kill:羞耻若不在前三,它就不特殊',rank<=3,f"排名 {rank}/{len(T)}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 同一条面板,换成 `#179` 的块层位置分 S ------------------------------------
# 上面的面板测的是**截距**(`#180`/`#182`/`#183` 那条线),不是 `#179` 那个 9× 的 S。
# **一条声明的多重性必须在它自己的量上判**,不能借另一条的结果。
print("\n" + "="*70 + "\n同一条面板,换成 #179 的块层位置分 S\n" + "="*70)
import zlib
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df); pos=np.zeros(NN); cntb=np.zeros(NN); KB=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); KB[ppl]+=M.sum(1); cntb[ppl]+=1
okb=cntb>=8
S=np.where(okb,pos/np.maximum(cntb,1),np.nan); KB=np.where(okb,KB,np.nan)
print(f"S 有效 {np.isfinite(S).sum():,} 人")

def prS(y,x,idx):
    X=np.c_[np.ones(len(idx)),KB[idx]]
    ry=y[idx]-X@np.linalg.lstsq(X,y[idx],rcond=None)[0]
    rx=x[idx]-X@np.linalg.lstsq(X,x[idx],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

rows2=[]; nulls2=[]
for c in lik:
    y=df[c].values.astype(float)
    m=np.isfinite(y)&np.isfinite(S)&np.isfinite(KB); idx=np.flatnonzero(m)
    if len(idx)<3000: continue
    r=prS(y,S,idx)
    bs=[prS(y,S,rb.choice(idx,len(idx),replace=True)) for _ in range(200)]
    def perm2():
        yp=y.copy(); yp[idx]=rb.permutation(y[idx]); return abs(prS(yp,S,idx))
    per=[v for v in (perm2() for _ in range(50)) if np.isfinite(v)]
    if len(per)<25: continue
    nulls2.append(per); rows2.append(dict(q=c[:70],n=len(idx),r=r,sd=float(np.std(bs)),
                                          ratio=abs(r)/float(np.std(bs))))
T2=pd.DataFrame(rows2); check_columns(T2,'R229 S 面板'); T2=T2.sort_values('r',key=abs,ascending=False)
T2.to_csv(pathlib.Path(__file__).parent/'results'/'panel_S.csv',index=False)
L2=min(len(x) for x in nulls2); thr2=float(np.nanquantile(np.nanmax(np.array([x[:L2] for x in nulls2]),axis=0),0.95))
check_coverage(len(T2),len(lik),'R229 S 面板覆盖',tol=0.15)
print(f"\n全族 95% 阈值 |r| = {thr2:.4f}\n")
for k,(_,r) in enumerate(T2.head(8).iterrows(),1):
    print(f"{k:<4}{r.r:>+8.4f}{r.ratio:>7.1f} {'★' if abs(r.r)>thr2 else ' '} {r.q[:70]}")
sh2=T2[T2.q.str.contains('ashamed')].iloc[0]
rank2=int(T2.reset_index().index[T2.reset_index().q.str.contains('ashamed')][0])+1
n2=int((T2.r.abs()>thr2).sum())
print(f"\nS 面板:羞耻 {sh2.r:+.4f}({sh2.ratio:.1f}×),排名 {rank2}/{len(T2)};超阈值 {n2}/{len(T2)}")

g2=Gate('S ↔ 羞耻 在多重性下站不站得住')
g2.asserted('可判前提:S 面板覆盖了大部分题目',len(T2)>=len(lik)*0.85,f"{len(T2)}/{len(lik)}")
g2.threshold_outside_noise('S↔羞耻 vs 全族阈值',float(abs(sh2.r)),thr2,float(sh2.sd))
g2.asserted('S↔羞耻 越过全族阈值',abs(sh2.r)>thr2,f"{abs(sh2.r):.4f} vs {thr2:.4f}")
g2.asserted('截距↔羞耻 越过全族阈值',abs(sh_row.r)>thr,f"{abs(sh_row.r):.4f} vs {thr:.4f}")
print(g2)

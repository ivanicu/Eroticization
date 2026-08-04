import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A35 R232 -- `#165` 那个去位置后的内容残差,是信号还是噪声

`#165`:内容分半信度 +0.26,减去泄漏底噪 +0.076 = **+0.187,1.6× -> UNVERIFIED**,
而 `#165` 自己写着"本设计答不了"。
**`#179`–`#186` 造出了一个当时没有的工具:一个与块零 item 重叠的外部结局(羞耻)。**
一个**只是噪声**的残差,不可能预测一个外部变量。

    SIGNAL  内容残差预测羞耻,且与 S 的预测**不重叠** -> `#165` 的 UNVERIFIED 可以收
    NOISE   内容残差对羞耻**毫无预测** -> 那个残差是噪声,从 UNVERIFIED 降为 NULL
    SHADOW  内容残差预测羞耻,但控制 S 后归零 -> 它只是 S 的影子,`#165` 的说法不变

ESTIMAND        r(内容残差, 羞耻),raw 与控制 S 之后;并判 S 与内容残差**各自的增量**。
KILL            条件式:先要 **S↔羞耻 在同一批人上复现 `#179`**(否则管道变了,整轮不可读);
                再判:**内容残差 raw < 2× -> NOISE;raw >2× 而控制 S 后 <2× -> SHADOW;
                两者都 >2× -> SIGNAL。**
NEGATIVE CTRL   跨人打乱羞耻。
POSITIVE CTRL   把 S 本身当"内容残差"塞进同一条管道,必须强测到(证明管道能测到东西)。
IMPOSSIBLE      内容分数由跨人主成分定义(`#165` 的 IMPOSSIBLE 原话)——
                完全特异的组合对它不可见,所以 NOISE 判定只覆盖**共享的**内容维度。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_residualized

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
sh=df['"I am ashamed or embarrassed about at least some of what arouses me" (7cw1ziu)'].values.astype(float)
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
    Z=M-M.mean(0,keepdims=True)
    w,v=np.linalg.eigh(np.cov(Z,rowvar=False)); pc=v[:,-1]
    con[ppl]+=Z@pc; pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1)
    KB[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
C=np.where(ok,con/np.maximum(cnt,1),np.nan); S=np.where(ok,pos/np.maximum(cnt,1),np.nan)
KB=np.where(ok,KB,np.nan)
print(f"S 与 C 有效 {int(np.isfinite(S).sum()):,} 人")

m=np.isfinite(sh)&np.isfinite(S)&np.isfinite(C)&np.isfinite(KB); idx=np.flatnonzero(m)
# `#165` 的"去位置后的内容":把 C 对 S 回归掉
X=np.c_[np.ones(len(idx)),S[idx]]
Cres=np.full(NN,np.nan); Cres[idx]=C[idx]-X@np.linalg.lstsq(X,C[idx],rcond=None)[0]
check_residualized(Cres[idx],S[idx],'R232 内容残差')      # #129:残差与被回归掉的协变量相关必须是 0

def pr(y,x,ctrls,ii):
    XX=np.c_[np.ones(len(ii)),*[c[ii] for c in ctrls]] if ctrls else np.ones((len(ii),1))
    ry=y[ii]-XX@np.linalg.lstsq(XX,y[ii],rcond=None)[0]
    rx=x[ii]-XX@np.linalg.lstsq(XX,x[ii],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

rb=np.random.default_rng(20260803); rows=[]
for nm,x,ctrls in [('S(复现 #179)',S,[KB]),('内容残差 raw',Cres,[KB]),
                   ('内容残差 +S',Cres,[KB,S]),('正对照:把 S 当残差塞进来',S,[KB])]:
    r=pr(sh,x,ctrls,idx)
    bs=[pr(sh,x,ctrls,rb.choice(idx,len(idx),replace=True)) for _ in range(400)]
    sd=float(np.std(bs)); rows.append(dict(term=nm,r=r,sd=sd,ratio=abs(r)/sd))
T=pd.DataFrame(rows); check_columns(T,'R232'); T.to_csv(pathlib.Path(__file__).parent/'results'/'terms.csv',index=False)
print(f"\n{'项':<24}{'r(·, 羞耻)':>12}{'sd':>9}{'比':>7}")
for _,r in T.iterrows(): print(f"{r.term:<24}{r.r:>+12.4f}{r.sd:>9.4f}{r.ratio:>7.1f}")

shp=sh.copy(); shp[idx]=rb.permutation(sh[idx])
r_null=pr(shp,Cres,[KB],idx)
print(f"\n跨人打乱羞耻:内容残差 r = {r_null:+.4f}")

r_s=float(T[T.term.str.startswith('S(')].r.iloc[0]); rt_s=float(T[T.term.str.startswith('S(')].ratio.iloc[0])
r_c=float(T[T.term=='内容残差 raw'].r.iloc[0]); rt_c=float(T[T.term=='内容残差 raw'].ratio.iloc[0])
r_cs=float(T[T.term=='内容残差 +S'].r.iloc[0]); rt_cs=float(T[T.term=='内容残差 +S'].ratio.iloc[0])
g=Gate('#165 的内容残差:信号还是噪声')
g.asserted('可判前提:S↔羞耻 在同一批人上复现 `#179` 的 +0.1155',abs(r_s-0.1155)<0.02,
           f"{r_s:+.4f} vs #179 的 +0.1155")
g.negative_control('跨人打乱羞耻(内容残差)',abs(r_null),max(abs(r_c),1e-6))
g.asserted('#129 守卫:残差与 S 的相关在构造上为 0',True,'check_residualized 已通过')
g.require_resolvable_first('内容残差 raw',r_c,float(T[T.term=='内容残差 raw'].sd.iloc[0]),family='content')
verdict=('SIGNAL' if (rt_c>2 and rt_cs>2) else ('SHADOW' if rt_c>2 else 'NOISE'))
g.asserted(f'注册的判定:{verdict}',True,
           f"raw {rt_c:.1f}× · 控制 S 后 {rt_cs:.1f}× -> **{verdict}**")
print(g)
print(f"\n  => **{verdict}**" + ("  —— `#165` 的那个残差对一个外部结局毫无预测,应从 UNVERIFIED 降为 NULL"
      if verdict=='NOISE' else ""))
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

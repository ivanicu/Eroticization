import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A39 R240 -- 这条线测到的是「羞耻」,还是一个「自我—身份」维度

`#194` 的 NEXT:`allrollidentity` 是这 20 道里**唯一一个关于"身份"而非"欲望"**的量,
它在内容侧越阈值(−0.0423)、去性别后保留 89%,**而从没被追过。**
若它与羞耻**共享同一个法则网络画像**,那么 `#179` 的构念名要重估。

ESTIMAND        **法则网络画像**:对每一道题,算它与**其余 18 道**的相关向量;
                判 `corr(画像_羞耻, 画像_allrollidentity)`,n = 18。
KILL            **若两者画像相关 > 0.60 -> 它们在这份数据里测的是同一个东西,
                `#179` 的构念名要从"羞耻"重估。**
POSITIVE CTRL   `animated` 与 `written`(同为媒介题)的画像相关必须**高** ——
                低了说明这个画像方法测不出"同一个东西",整轮不可读。
NEGATIVE CTRL   `biomale` 与 `highenergy`(一个人口学、一个一般特质)画像相关必须**低**。
NOISE FLOOR     人层 bootstrap 300 次重算整条画像。
IMPOSSIBLE      n = 18 个共同题。|corr| < 0.47 在 n=18 上不显著,所以只能判**大的**相似。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
SHAME=next(c for c in lik if 'ashamed' in c)
IDENT=next(c for c in lik if c=='allrollidentity')
ANI=next(c for c in lik if c=='animated'); WRI=next(c for c in lik if c=='written')
BIO=next(c for c in lik if c=='biomale'); HIG=next(c for c in lik if c=='highenergy')
Y=df[lik].values.astype(float)
ok=np.isfinite(Y).all(axis=1); bi=np.flatnonzero(ok)
print(f"20 道全非缺失的人 {len(bi):,}")

def profile(col_a, col_b, ii):
    """a 与 b 各自对**其余 18 道**(排除 a、b)的相关向量。"""
    others=[c for c in lik if c not in (col_a,col_b)]
    A=df[col_a].values.astype(float)[ii]; B=df[col_b].values.astype(float)[ii]
    pa=[];pb=[]
    for c in others:
        z=df[c].values.astype(float)[ii]
        pa.append(np.corrcoef(A,z)[0,1]); pb.append(np.corrcoef(B,z)[0,1])
    return np.array(pa),np.array(pb),others

def sim(col_a,col_b,ii):
    pa,pb,_=profile(col_a,col_b,ii)
    return float(np.corrcoef(pa,pb)[0,1]), len(pa)

rb=np.random.default_rng(20260803)
PAIRS=[('羞耻 vs allrollidentity',SHAME,IDENT),
       ('【正对照】animated vs written',ANI,WRI),
       ('【负对照】biomale vs highenergy',BIO,HIG),
       ('羞耻 vs animated',SHAME,ANI),
       ('羞耻 vs biomale',SHAME,BIO)]
rows=[]
for name,a,b in PAIRS:
    s,n=sim(a,b,bi)
    bs=[sim(a,b,rb.choice(bi,len(bi),replace=True))[0] for _ in range(300)]
    sd=float(np.std(bs))
    rows.append(dict(pair=name,n_items=n,similarity=s,sd=sd,ratio=abs(s)/sd))
T=pd.DataFrame(rows); check_columns(T,'R240'); T.to_csv(pathlib.Path(__file__).parent/'results'/'profiles.csv',index=False)
print(f"\n{'画像相似度':>10}{'sd':>9}{'比':>7}  对(n = 18 道共同题)")
for _,r in T.iterrows(): print(f"{r.similarity:>+10.4f}{r.sd:>9.4f}{r.ratio:>7.1f}  {r.pair}")

# 逐项:羞耻与 allrollidentity 在哪些题上像 / 不像
pa,pb,others=profile(SHAME,IDENT,bi)
d=pd.DataFrame(dict(q=[c[:54] for c in others],shame=pa,ident=pb,gap=pa-pb)).sort_values('gap',key=abs,ascending=False)
d.to_csv(pathlib.Path(__file__).parent/'results'/'itemwise.csv',index=False)
print(f"\n差最大的 5 道:")
for _,r in d.head(5).iterrows(): print(f"  羞耻 {r.shame:+.3f} · 身份 {r.ident:+.3f} · 差 {r.gap:+.3f}  {r.q[:52]}")

main=T[T.pair.str.startswith('羞耻 vs allroll')].iloc[0]
pos =T[T.pair.str.contains('正对照')].iloc[0]
neg =T[T.pair.str.contains('负对照')].iloc[0]
g=Gate('羞耻与身份是不是同一个东西')
g.asserted('正对照:animated 与 written 的画像必须高度相似',pos.similarity>0.60,
           f"{pos.similarity:+.4f}")
g.asserted('负对照:biomale 与 highenergy 的画像必须不相似',abs(neg.similarity)<0.60,
           f"{neg.similarity:+.4f}")
g.resolvable('羞耻 vs 身份 的画像相似度',float(main.similarity),float(main.sd))
g.threshold_outside_noise('相似度 vs 注册阈值 0.60',float(main.similarity),0.60,float(main.sd))
g.asserted('注册的 kill:画像相关 > 0.60 -> 构念名要重估',main.similarity>0.60,
           f"{main.similarity:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

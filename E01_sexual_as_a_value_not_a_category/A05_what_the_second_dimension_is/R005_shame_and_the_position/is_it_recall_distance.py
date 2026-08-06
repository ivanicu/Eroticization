import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A33 R228 -- 「起点比轨迹重要」还剩一个对手:年龄

`#182`:截距(早期冷门度)预测羞耻 +0.0270(2.9×,n=10,567),两次独立复现。
**但截距是"年龄轴原点上的冷门度"**,而一个人**现在**多大、他回溯了多远,
会同时影响截距与羞耻。`#162` 已证回忆偏差是**人群规律**,
**但没人问过它是否也搬运了这条关联。**

ESTIMAND        r(截距, 羞耻) 在 ① 控制当前年龄之后 ② 按年龄段分层之内。
                年龄只有 5 个段(14-17 / 18-20 / 21-24 / 25-28 / 29-32),用段中点。
KILL            **跑之前注册**(`#182` 的 NEXT 原话):
                「若关联在最年长的一层消失,那么它可能是回溯距离的产物,而不是"起点"的产物。」
                ⚠ 补一句它没写全的:回溯距离若是机制,关联应当**随年龄单调变化**;
                只在最老一层消失而中间层不动,是**功效**问题不是机制问题 -> 那时判 UNVERIFIED。
NEGATIVE CTRL   段内打乱羞耻。
NOISE FLOOR     每层各自 bootstrap 500。
MULTIPLICITY    5 层 × {raw, +年龄} + 全样本 2 格 = 12 格,整格发表。
IMPOSSIBLE      年龄段只有 5 个且最老一段是 29-32 —— **本数据里没有中年以上的人**,
                所以"回溯 30 年"这个区间根本不在样本里。这一条限制结论的外推。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
sh=df['"I am ashamed or embarrassed about at least some of what arouses me" (7cw1ziu)'].values.astype(float)
BAND={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
age=df['age'].map(BAND).values.astype(float)
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
    X=np.c_[np.ones(len(x)),x]
    I0[i]=np.linalg.lstsq(X,rar[idx],rcond=None)[0][0]
print(f"截距有效 {np.isfinite(I0).sum():,};年龄有效 {np.isfinite(age).sum():,}")
print(f"corr(截距, 当前年龄) = {np.corrcoef(I0[np.isfinite(I0)&np.isfinite(age)],age[np.isfinite(I0)&np.isfinite(age)])[0,1]:+.4f}")

def pr(y,x,ctrls,idx):
    X=np.c_[np.ones(len(idx)),*[c[idx] for c in ctrls]] if ctrls else np.ones((len(idx),1))
    ry=y[idx]-X@np.linalg.lstsq(X,y[idx],rcond=None)[0]
    rx=x[idx]-X@np.linalg.lstsq(X,x[idx],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

rb=np.random.default_rng(20260803); rows=[]
base=np.isfinite(sh)&np.isfinite(I0)&np.isfinite(K)&np.isfinite(age)
for name,mask,ctrls in ([('全样本 raw',base,[K]),('全样本 +年龄',base,[K,age])]+
                        [(f'{b} 段内',base&(age==v),[K]) for b,v in BAND.items()]):
    idx=np.flatnonzero(mask)
    if len(idx)<300: continue
    r=pr(sh,I0,ctrls,idx)
    bs=[pr(sh,I0,ctrls,rb.choice(idx,len(idx),replace=True)) for _ in range(500)]
    sd=float(np.std(bs))
    null=pr(rb.permutation(sh),I0,ctrls,idx)
    rows.append(dict(stratum=name,n=len(idx),r=r,sd=sd,ratio=abs(r)/sd,null=null))
T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'strata.csv',index=False)
print(f"\n{'层':<14}{'n':>7}{'r(截距, 羞耻)':>14}{'sd':>9}{'|r|/sd':>8}{'打乱零':>9}")
for _,r in T.iterrows():
    print(f"{r.stratum:<14}{int(r.n):>7,}{r.r:>+14.4f}{r.sd:>9.4f}{r.ratio:>8.1f}{r.null:>+9.4f}")

st=T[T.stratum.str.contains('段内')]
old=st.iloc[-1]; young=st.iloc[0]
mono=np.corrcoef([BAND[s.split(' ')[0]] for s in st.stratum],st.r)[0,1]
print(f"\n  按年龄的单调性:corr(年龄段中点, r) = {mono:+.4f}")

g=Gate('这条关联是不是回溯距离的产物')
g.asserted('可判前提:全样本 raw 仍可分辨(否则无从谈分层)',
           float(T[T.stratum=='全样本 raw'].ratio.iloc[0])>2,
           f"{float(T[T.stratum=='全样本 raw'].r.iloc[0]):+.4f}, "
           f"{float(T[T.stratum=='全样本 raw'].ratio.iloc[0]):.1f}×")
g.negative_control('段内打乱羞耻(全样本)',float(abs(T[T.stratum=='全样本 raw'].null.iloc[0])),
                   float(T[T.stratum=='全样本 raw'].r.iloc[0]))
g.asserted('控制当前年龄后不塌',
           float(T[T.stratum=='全样本 +年龄'].ratio.iloc[0])>2,
           f"{float(T[T.stratum=='全样本 +年龄'].r.iloc[0]):+.4f}, "
           f"{float(T[T.stratum=='全样本 +年龄'].ratio.iloc[0]):.1f}×")
g.asserted('注册的 kill:关联在最年长一层消失',float(old.ratio)<2,
           f"{old.stratum} r={old.r:+.4f} ({old.ratio:.1f}×)")
g.asserted('⚠ 补注册:回溯距离若是机制,r 应随年龄**单调**变化',abs(mono)>0.7,
           f"corr(年龄, r) = {mono:+.4f} —— 不单调则只在一层消失是**功效**问题,不是机制问题")
g.no_sign_crossing('五层同号',[float(x) for x in st.r])
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

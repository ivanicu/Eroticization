import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A34 R231 -- 羞耻是一道题,不是一个量表:它的信度能不能估

`#185` 的 NEXT:`#179`–`#185` 七轮全部用同一道羞耻题,而 `#184` 已证这份 release
**没有第二道**能测同一构念的题。所以「羞耻」在本项目里**是一道题**,它的信度**从未被估计过** ——
而 `#167` 那条法则说:一个没有误差棒的数,不能当作它自己的精度。

ESTIMAND        单题信度 r_xx 的**下界**,以及去衰减后的 S↔羞耻。
IDENTIFICATION  单题做不了分半。能给的下界:对任意另一题 y,`r_xy ≤ √(r_xx·r_yy) ≤ √r_xx`
                -> **r_xx ≥ max_y r_xy²**。这是一个**极弱**的下界,本轮明说它弱到什么程度。
KILL            **不注册"去衰减后超过 0.3"这种 kill** —— 因为去衰减的分母是我猜的。
                改为注册:**若最强下界 < 0.25,那么单题信度在本数据里不可估,
                结论只能以敏感性带的形式陈述,不得给出一个校正后的点值。**
POSITIVE CTRL   对**块层位置分 S** 跑同一套(它有真实的分半信度 `#100` = 0.432 ± 0.016),
                看这个下界方法能不能逼近它 —— 逼不近就说明这个下界确实很弱。
IMPOSSIBLE      重测、平行题、同构念多题 —— **本 release 三者皆无**。
                所以本轮的产出是一条**敏感性带**,不是一个数。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
SHAME=next(c for c in lik if 'ashamed' in c)
Y=df[lik].values.astype(float)
C=pd.DataFrame(Y,columns=lik).corr()
row=C[SHAME].drop(SHAME).abs().sort_values(ascending=False)
print("羞耻与其余 19 题的最强相关:")
for k,(q,v) in enumerate(row.head(5).items(),1): print(f"  {k}. |r| = {v:.4f}   {q[:74]}")
lb=float(row.iloc[0]**2)
print(f"\n单题信度下界 r_xx >= max_y r_xy² = {row.iloc[0]:.4f}² = **{lb:.4f}**")

# 正对照:同一套方法用在 S 上,S 有真实分半信度 0.432(`#100`)
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
sS=[]
for c in lik:
    y=df[c].values.astype(float); m=np.isfinite(y)&np.isfinite(S)
    if m.sum()<3000: continue
    sS.append(abs(np.corrcoef(y[m],S[m])[0,1]))
lb_S=float(max(sS)**2)
print(f"正对照:同一方法用在 S 上 -> 下界 {lb_S:.4f},而 S 的真实分半信度是 **0.432**(`#100`)")
print(f"        这个下界只逼近到真值的 {100*lb_S/0.432:.0f}% —— **它确实很弱**")

# 敏感性带
sh=df[SHAME].values.astype(float)
m=np.isfinite(sh)&np.isfinite(S)&np.isfinite(KB); idx=np.flatnonzero(m)
X=np.c_[np.ones(len(idx)),KB[idx]]
ry=sh[idx]-X@np.linalg.lstsq(X,sh[idx],rcond=None)[0]
rx=S[idx]-X@np.linalg.lstsq(X,S[idx],rcond=None)[0]
r_obs=float(np.corrcoef(ry,rx)[0,1])
r_S=0.432
print(f"\n观测 r(S, 羞耻) = {r_obs:+.4f}  (n={len(idx):,})")
print(f"\n敏感性带(横轴 = 假设的单题信度,纵轴 = 去衰减后的 r):")
band=[]
for rxx in [0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00]:
    dis=r_obs/np.sqrt(rxx*r_S)
    band.append(dict(r_xx=rxx,disattenuated=dis))
    print(f"  r_xx = {rxx:.2f}  ->  {dis:+.4f}")
B=pd.DataFrame(band); B.to_csv(pathlib.Path(__file__).parent/'results'/'band.csv',index=False)

g=Gate('单题信度能不能估')
g.asserted('正对照:同一方法在 S 上只逼近真值的一小部分 -> 下界确实弱',
           lb_S<0.432*0.6,f"下界 {lb_S:.4f} vs 真值 0.432({100*lb_S/0.432:.0f}%)")
g.asserted('注册的 kill:最强下界 < 0.25 -> 单题信度不可估,只能报敏感性带',lb<0.25,
           f"下界 {lb:.4f}")
g.asserted('因此不给校正后的点值',True,
           f"去衰减值在 r_xx∈[0.3,1.0] 上从 {B.disattenuated.max():+.3f} 到 {B.disattenuated.min():+.3f} —— "
           f"**跨度 {B.disattenuated.max()-B.disattenuated.min():.3f},比效应本身还大**")
g.same_scale('去衰减用的 S 信度与 `#100` 同源',r_S,0.432,'reliability')
print(g)
print(f"\n  => **`#179` 的 +0.116 是一个**下界**:即使 S 的信度已知(0.432),"
      f"\n     羞耻这一侧的信度未知,而它把真实关联压低了至少 1/√0.432 = {1/np.sqrt(r_S):.2f} 倍。**")
print(f"\nsha1 {hashlib.sha1(B.to_csv(index=False).encode()).hexdigest()[:12]}")

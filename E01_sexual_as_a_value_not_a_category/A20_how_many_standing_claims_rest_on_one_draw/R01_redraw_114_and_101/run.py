import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A20 R01 -- 有多少现存声明的零,是一次抽样?

#152e:**当一个比较的基准是一次随机实现时,判它的分母是那个基准的实现展布,
不是被比较量的抽样展布。** `#148` 因为这个把 1.4× 报成了 20.0×。

扫 139 个持久化结果文件:53 个含零臂,**17 个的零臂没有 seed 列**。分诊后,
两个支撑**现存声明**的零确认是单次抽样:

  #114(A12/R10,「人把最爱的性兴趣记得更早,约九个月」,19.8× SE)
        `rp = np.random.default_rng(9)` —— **一次**同类别内打乱。
        README 写的「同类别内打乱零 = 效应的 0%」出自那一次。
  #101/#102(A11/R14,「唯一挂得住的外部锚是性别 +0.093」)
        `curveball(M, default_rng(8100))` —— **一次**保边际实现。

**重抽它们的零 20 次,看那"0%"是不是运气。**

ESTIMAND        零的实现分布:20 次独立重抽的均值与展布,以及"零占效应的比例"的区间。
IDENTIFICATION  与原轮**逐字相同**的零构造,只把种子换成 20 个。
SCOPE           原轮的口径。
WORLDS          LUCKY  零的实现展布很大,原来的一次抽样落在低端 -> 那条 negative_control
                       的判定是运气,声明要加限定
                SOLID  实现展布很小,任何一次抽样都会给出同样的判定 -> 一次抽样在这里够用,
                       而这本身是值得记录的(不是所有单抽都有问题)
KILL            条件式:效应臂必须复现原轮的值(逐位或在自助展布内),才读零的重抽。
POSITIVE CTRL   效应臂复现(见上)。
NEGATIVE CTRL   —— 本轮**就是**在造零对照,不再套一层。
NOISE FLOOR     20 次重抽。
MULTIPLICITY    2 条声明 x 20 抽,整格发表。
IMPOSSIBLE      重抽只检验**实现方差**;若零的**构造**本身有偏(而不是吵),重抽看不见。
"""
import numpy as np, pandas as pd, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').values
P=(np.nan_to_num(R)>0).astype(float); breadth=P.sum(1)
def norm(s): return re.sub(r'[^a-z]',' ',s.lower())
best={}
for j,c in enumerate(ons):
    m=re.search(r'interest in ([a-z /-]+)',norm(c))
    if not m: continue
    ws=set(w for w in m.group(1).split() if len(w)>4)
    if not ws: continue
    sc=[(len(ws&set(norm(rc).split())),i) for i,rc in enumerate(rate)]
    s,i=max(sc)
    if s>=1: best[j]=i
prec=np.nanmean(V,axis=1)
rows=[]
for j,ri in best.items():
    col=V[:,j]; y=R[:,ri]
    m=np.isfinite(col)&np.isfinite(y)&(np.isfinite(V).sum(1)>=6)&np.isfinite(prec)
    idx=np.flatnonzero(m)
    if len(idx)<400: continue
    s=col[idx].sum(); n=len(idx); loo=(s-col[idx])/(n-1)
    dev=(col[idx]-loo)-(prec[idx]-np.nanmean(prec[idx])); dev=dev-dev.mean()
    rows.append(pd.DataFrame(dict(cat=j,person=idx,dev=dev,rating=y[idx],
                                  breadth=breadth[idx],prec=prec[idx])))
L=pd.concat(rows,ignore_index=True)
print(f"#114 长表 {len(L):,} 行  {L.person.nunique():,} 人  {L.cat.nunique()} 类别",flush=True)

def slope(Lx):
    x=(Lx.rating.values-Lx.rating.mean())/Lx.rating.std()
    Z=np.c_[np.ones(len(Lx)),x,(Lx.breadth.values-Lx.breadth.mean())/Lx.breadth.std(),
            (Lx.prec.values-Lx.prec.mean())/Lx.prec.std()]
    return float(np.linalg.lstsq(Z,Lx.dev.values,rcond=None)[0][1])
b_real=slope(L)
jj=L.person.unique(); gi={p:np.flatnonzero(L.person.values==p) for p in jj}
rbb=np.random.default_rng(0)
se_real=float(np.std([slope(L.iloc[np.concatenate([gi[p] for p in rbb.choice(jj,len(jj))])])
                      for _ in range(200)]))
print(f"  效应臂 {b_real:+.4f} ± {se_real:.4f}({abs(b_real)/se_real:.1f}×;#114 报 −0.2000,19.8×)",flush=True)

draws=[]
for s in range(20):
    rp=np.random.default_rng(9 if s==0 else 1000+s)     # s=0 复现原轮的那一次
    Ls=L.copy()
    Ls['rating']=Ls.groupby('cat').rating.transform(lambda v: rp.permutation(v.values))
    draws.append(slope(Ls))
draws=np.array(draws)
print(f"\n#114 的零,重抽 20 次:")
print(f"  原轮那一次(seed 9)= {draws[0]:+.5f}   占效应 {100*abs(draws[0]/b_real):.1f}%")
print(f"  20 抽:均值 {draws.mean():+.5f}  sd {draws.std():.5f}  "
      f"范围 [{draws.min():+.5f}, {draws.max():+.5f}]")
print(f"  占效应的比例范围 {100*abs(draws).min()/abs(b_real):.1f}% .. "
      f"{100*abs(draws).max()/abs(b_real):.1f}%")

T=pd.DataFrame(dict(draw=range(20),null=draws))
T.to_csv(pathlib.Path(__file__).parent/'results'/'redraw114.csv',index=False)
g=Gate('#114 的"零 = 效应的 0%"是不是一次抽样的运气')
g.asserted('效应臂复现原轮(正对照)',abs(abs(b_real)-0.2000)<0.03,
           f"{b_real:+.4f} vs #114 报的 −0.2000")
g.asserted('零的实现展布相对效应有多大',True,
           f"sd(零) {draws.std():.5f} = 效应的 {100*draws.std()/abs(b_real):.1f}%;"
           f"最大的一抽占效应 {100*abs(draws).max()/abs(b_real):.1f}%")
g.negative_control('零(20 抽里最不利的一抽)对效应',float(abs(draws).max()),float(abs(b_real)),
                   null_spread=float(draws.std()))
g.equivalent_within('零的分布是否被界在效应的 10% 内',float(np.abs(draws).mean()),
                    float(draws.std()),0.10*abs(b_real))
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

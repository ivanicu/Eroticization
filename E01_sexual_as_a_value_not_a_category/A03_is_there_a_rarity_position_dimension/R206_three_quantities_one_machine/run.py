import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A21 R03 -- 「何时」的人层成分有多小:三个量放进同一台机器

#160c:「什么」人侧信度 0.57–0.67,「何时」只有 0.055–0.108,差 5.8–10.4 倍。
#160 的 NEXT 问:这和 `#63` 的上界是不是同一件事的两种说法。

而账本里已经有第三个可比的数:**`Entry 8`** —— 「时间表贴合度是一种个人特质」被杀掉,
**跨类别分半 Spearman-Brown = 0.214,"mostly measurement noise"**。

**三个量,一台机器**(同一批人、同一种劈分、同一个 k、同一个 SB):
  S  你**喜欢**多罕见的东西                        —— #100 / #160
  z  你的罕见兴趣比时间表**早多少**                —— #159 / #160
  a  共享时间表**贴合你**的程度(成对顺序准确率)  —— Entry 8,SB 0.214

ESTIMAND        三个量各自的分半信度,作为每半单位数 k 的函数,同一台机器。
IDENTIFICATION  三个量都在"人内把单位劈成不相交两半 -> 各算一次 -> 跨人相关 -> SB"上跑。
SCOPE           同时满足三个量单位数要求的人。
WORLDS          SAME   a 落在 z 附近(~0.1)-> 「何时」的人层成分小,两条独立的路同一个数
                SPLIT  a 明显高于 z -> "贴合时间表的程度"与"偏离方向"是两件不同的事,
                       而 Entry 8 的 0.214 是对的
KILL            条件式:三个量的人特异种植都必须把信度推上去,才读三条曲线的相对位置。
POSITIVE CTRL   三个种植。
NEGATIVE CTRL   人内置换各自的单位标签。
NOISE FLOOR     5 个劈分种子。
MULTIPLICITY    k ∈ {4,5,6,7} x 3 个量 x {真实, 置换, 种植},整格发表。
IMPOSSIBLE      三个量的单位信息量不同(块 vs 类别 vs 类别对),匹配 k 拉不平这一条。
                与 `#160d` 同一条范围,继续携带。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A21_who_owns_the_personal_20_percent'
          /'R205_what_versus_when_same_design'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('KS=[4,5,6,7]')[0])   # 跨轮依赖显式声明(P16)

pop=np.array([np.nanmean(V[obs[:,j],j]) for j in range(V.shape[1])])   # #75 的共享时间表

def adherence(cols,i,perm,rg):
    """共享时间表在这个人的类别对上的成对顺序准确率。"""
    j=np.intersect1d(np.flatnonzero(obs[i]),cols)
    if len(j)<3: return np.nan
    y=V[i,j].copy()
    if perm: y=y[rg.permutation(len(y))]
    p=pop[j]; right=0; tot=0
    for a in range(len(j)):
        for b in range(a+1,len(j)):
            if y[a]==y[b] or p[a]==p[b]: continue
            right+=((p[a]<p[b])==(y[a]<y[b])); tot+=1
    return right/tot if tot>=3 else np.nan

def half_rel3(kind,k,seed,perm=False,plant=0.):
    rg=np.random.default_rng(seed); tie=np.random.default_rng(20260803)
    u=rg.standard_normal(len(V)); A=[];B=[]
    for i in np.flatnonzero(KEEP&okA):
        if kind=='S':
            av=np.flatnonzero(okB[i])
            if len(av)<2*k: continue
            p_=rg.permutation(av); h1,h2=p_[:k],p_[k:2*k]
            x1=np.nanmean(BMc[i,h1])+plant*u[i]; x2=np.nanmean(BMc[i,h2])+plant*u[i]
            if perm:
                x1=np.nanmean(BMc[rg.integers(len(BMc)),h1]); x2=np.nanmean(BMc[rg.integers(len(BMc)),h2])
        else:
            av=np.flatnonzero(obs[i])
            if len(av)<2*k: continue
            p_=rg.permutation(av); h1,h2=p_[:k],p_[k:2*k]
            f=z_on if kind=='z' else adherence
            x1=f(h1,i,perm,rg,tie) if kind=='z' else f(h1,i,perm,rg)
            x2=f(h2,i,perm,rg,tie) if kind=='z' else f(h2,i,perm,rg)
            if plant: x1=(x1 or np.nan)+plant*u[i]; x2=(x2 or np.nan)+plant*u[i]
        if np.isfinite(x1) and np.isfinite(x2): A.append(x1); B.append(x2)
    A=np.array(A); B=np.array(B)
    if len(A)<300: return np.nan,np.nan,len(A)
    r=float(np.corrcoef(A,B)[0,1])
    return (2*r/(1+r) if r>-0.99 else np.nan), r, len(A)

KS=[4,5,6,7]; rows=[]
print(f"\n{'k':<4}{'量':<4}{'n':>7}{'分半 r':>10}{'SB 信度':>10}{'置换零':>10}{'种植':>9}")
for k in KS:
    for kind in ['S','z','a']:
        sd_=zlib.crc32(f'{kind}{k}x'.encode())%9973
        sb,r,n=half_rel3(kind,k,sd_)
        sbn,_,_=half_rel3(kind,k,sd_+1,perm=True)
        pl=1.0 if kind!='a' else 0.15
        sbp,_,_=half_rel3(kind,k,sd_+2,plant=pl)
        rows.append(dict(k=k,kind=kind,n=n,r=r,sb=sb,sb_null=sbn,sb_plant=sbp))
        print(f"{k:<4}{kind:<4}{n:>7,}{r:>+10.4f}{sb:>+10.4f}{sbn:>+10.4f}{sbp:>+9.4f}",flush=True)

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'three.csv',index=False)
P=T.pivot(index='k',columns='kind',values='sb')
print(f"\n  三条曲线(SB 信度):")
print("   " + P.round(4).to_string())
print(f"\n  a(时间表贴合度)均值 {P['a'].mean():+.4f}   z 均值 {P['z'].mean():+.4f}   "
      f"S 均值 {P['S'].mean():+.4f}")
print(f"  Entry 8 报的 a = 0.214;#160 报的 z = 0.055–0.108")

g=Gate('「何时」的人层成分,两条独立的路给出同一个数吗')
for kind in ['S','z','a']:
    sub=T[T.kind==kind]
    g.asserted(f'{kind} 的人特异种植把信度推上去',bool((sub.sb_plant>sub.sb+0.05).all()),
               " ".join(f"k={int(r.k)}:{r.sb:.3f}->{r.sb_plant:.3f}" for _,r in sub.iterrows()))
g.asserted('三个量的人内置换零都在零附近',bool((T.sb_null.abs()<0.10).all()),
           " ".join(f"{k}{int(r.k)}:{r.sb_null:+.3f}" for k,(_,r) in zip(T.kind,T.iterrows())))
sd=float((P['a']-P['z']).std())
g.equivalent_within('a 与 z 是不是同一个数(边界 0.10)',float((P['a']-P['z']).mean()),sd,0.10)
g.require_resolvable_first('a 与 S 的差',float((P['S']-P['a']).mean()),sd)
g.offset_control('a vs z',float(P['a'].mean()),float(P['z'].mean()),sd,
                 null_kind='同一台机器、同一个 k 下 z 的分半信度(不是零假设,是被比较的对象)')
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
E01 A17 R01 -- "去衰减 <= 0.362 所以是三条近乎独立的轴" —— 这个 0.4 的阈值是刀尖吗?

#140b 的 NEXT 写着「0.362 与 0.349 几乎相等」。**先纠正我自己的 NEXT:那是两个不同的量。**
  0.349 = 族内 PCA 成分 vs **随机载荷**全类别坐标的留出 |相关|
  0.362 = A02 手工命名量表之间的 r_obs / sqrt(alpha_a * alpha_b)
数字接近是巧合,不是等价。**但底下的担忧合法,只是要用 A02 自己的量来问。**

跑一遍 A02/R10 拿到它的实际数字:
  POWER alpha = 0.686(3 题)· GAZE alpha = **0.163**(2 题)· COORD4 分半 = 0.337
  原始相关 P-G +0.079 · P-C +0.112 · G-C +0.085   -> 全部 <= 0.112
  去衰减    P-G +0.236 · P-C +0.233 · G-C **+0.362**
  判定写死在源码里:`mx < 0.4 -> "still three axes"`

**G-C 的 0.362 = 0.085 x 4.27,而那个 4.27 来自 alpha_G = 0.163。**
一个 2 题量表的 alpha 0.163 不是信度,是噪声。所以问题是:

  KNIFE  去衰减值的自助分布**跨过 0.4** -> A02 的判定是刀尖,阈值又是选的,
         那么 README 上那一行必须改成引用**原始相关**而不是去衰减最大值
  ROBUST 自助分布整段在 0.4 以下 -> 判定稳,只是表述可以更强

ESTIMAND        max 去衰减 |r| 的按人自助分布,以及它跨过 0.4 的比例。
                加上一条**信度匹配**的随机地板(alpha 落在与 GAZE/COORD4 同一区间的随机量表)。
IDENTIFICATION  与 A02/R10 **逐字相同**的公式与数据;只把"人"换成自助重抽。
SCOPE           A02/R10 的原口径。
WORLDS          KNIFE / ROBUST
KILL            条件式:自助必须先在一个已知量(原始相关)上给出合理窄的区间,才读去衰减的区间。
POSITIVE CTRL   同一构念的两半(把 POWER 的 3 题劈成 2+1)去衰减必须接近 1。
NEGATIVE CTRL   合成两组互相独立、内部一致的数据,去衰减必须接近 0。
                ⚠ 两个对照都必须**先把题目定向**(按与临时合成分的相关取符号),
                因为**没定向的 alpha 不是信度** —— 第一版随机翻符号,两个对照都塌了。
NOISE FLOOR     500 次按人自助。
MULTIPLICITY    3 对量表 x {原始, 去衰减},整格发表。
IMPOSSIBLE      alpha 本身是 tau-等价假设下的下界;2 题量表的 alpha 等于其相关的
                Spearman-Brown 值,没有更多信息。本轮不修这个,只量它的后果。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore'); rng=np.random.default_rng(89)
from lib.gates import Gate

# 与 A02/R10 逐字相同的加载与量表构造
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts()
pool=np.array(sorted(nblk[nblk>=8].index)); pm={p:i for i,p in enumerate(pool)}
G=np.load('data/derived/gcca_G.npy')
A=pd.read_csv('data/derived/agent_patient.csv')
def z(s): return (s-s.mean())/(s.std()+1e-9)
def alpha_np(Z):
    Z=Z[np.isfinite(Z).all(1)]; k=Z.shape[1]
    if len(Z)<200 or k<2: return np.nan
    return k/(k-1)*(1-Z.var(0,ddof=1).sum()/max(Z.sum(1).var(ddof=1),1e-12))
pc=[c for c in A.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
Pm=pd.concat([z(pd.to_numeric(A[c],errors='coerce'))*sg[c] for c in pc],axis=1).values
ex=[c for c in A.columns if 'exhibition' in c][0]; vo=[c for c in A.columns if 'voyeur' in c][0]
Gm=pd.concat([z(pd.to_numeric(A[ex],errors='coerce')),-z(pd.to_numeric(A[vo],errors='coerce'))],axis=1).values
Ps=np.nanmean(Pm,axis=1); Gs=np.nanmean(Gm,axis=1)
Cs=np.full(len(A),np.nan); Cs[pool]=G[:,3]
REL_C=0.337                                   # A02/R10 实测的 coord4 分半(块分半,不随人自助)
print(f"POWER {Pm.shape[1]} 题  GAZE {Gm.shape[1]} 题  coord4 分半 {REL_C}",flush=True)

def one(idx):
    aP=alpha_np(Pm[idx]); aG=alpha_np(Gm[idx])
    out={}
    for nm,(x,y,ra,rb) in {'P-G':(Ps,Gs,aP,aG),'P-C':(Ps,Cs,aP,REL_C),'G-C':(Gs,Cs,aG,REL_C)}.items():
        m=np.isfinite(x[idx])&np.isfinite(y[idx])
        if m.sum()<500: out[nm]=(np.nan,np.nan); continue
        r=float(np.corrcoef(x[idx][m],y[idx][m])[0,1])
        out[nm]=(abs(r),abs(r)/np.sqrt(max(ra,1e-6)*max(rb,1e-6)))
    return out,aP,aG

allidx=np.arange(len(A))
pt,aP0,aG0=one(allidx)
print(f"\n点估计:alpha_P {aP0:.3f}  alpha_G {aG0:.3f}")
print(f"  {'对':<6}{'原始|r|':>10}{'去衰减|r|':>12}")
for k,(r,d) in pt.items(): print(f"  {k:<6}{r:>10.3f}{d:>12.3f}")

rb=np.random.default_rng(zlib.crc32(b'A17R01'))
BS={k:[] for k in pt}; BSr={k:[] for k in pt}; BSa=[]
for _ in range(500):
    s=rb.integers(0,len(A),len(A))
    o,_,ag=one(s); BSa.append(ag)
    for k,(r,d) in o.items(): BSr[k].append(r); BS[k].append(d)
print(f"\n按人自助 500 次:")
print(f"  {'对':<6}{'原始|r| [2.5,97.5]':>26}{'去衰减|r| [2.5,97.5]':>28}")
rows=[]
for k in pt:
    r=np.array(BSr[k]); d=np.array(BS[k]); d=d[np.isfinite(d)]
    rows.append(dict(pair=k,r=pt[k][0],r_lo=np.percentile(r,2.5),r_hi=np.percentile(r,97.5),
                     d=pt[k][1],d_lo=np.percentile(d,2.5),d_hi=np.percentile(d,97.5),
                     over40=float((d>0.4).mean())))
    print(f"  {k:<6}{f'{pt[k][0]:.3f} [{np.percentile(r,2.5):.3f},{np.percentile(r,97.5):.3f}]':>26}"
          f"{f'{pt[k][1]:.3f} [{np.percentile(d,2.5):.3f},{np.percentile(d,97.5):.3f}]':>28}")
mx=np.array([max(BS[k][i] for k in pt if np.isfinite(BS[k][i])) for i in range(len(BS['P-G']))])
mx=mx[np.isfinite(mx)]
over=float((mx>0.4).mean())
aga=np.array(BSa)
print(f"\n  alpha_G 的自助 [2.5,97.5] = [{np.percentile(aga,2.5):.3f},{np.percentile(aga,97.5):.3f}]"
      f"  中位 {np.median(aga):.3f}")
print(f"  **max 去衰减|r| 的自助 [2.5,97.5] = [{np.percentile(mx,2.5):.3f},{np.percentile(mx,97.5):.3f}]**")
print(f"  **跨过 A02 写死的 0.4 的比例 = {over:.1%}**")

# 正 / 负对照(题目先定向 —— 没定向的 alpha 不是信度)
def orient(Z):
    s=np.nanmean(Z,axis=1)
    sgn=np.array([1. if np.corrcoef(Z[np.isfinite(Z[:,j])&np.isfinite(s),j],
                                    s[np.isfinite(Z[:,j])&np.isfinite(s)])[0,1]>=0 else -1.
                  for j in range(Z.shape[1])])
    return Z*sgn[None,:]
h1,h2=orient(Pm[:,:2]),orient(Pm[:,2:])
s1,s2=np.nanmean(h1,axis=1),np.nanmean(h2,axis=1)
m=np.isfinite(s1)&np.isfinite(s2)
a1=alpha_np(h1); r12=abs(np.corrcoef(s1[m],s2[m])[0,1])
pos=r12/np.sqrt(max(a1,1e-6)*max(a1,1e-6)) if np.isfinite(a1) and a1>0.05 else np.nan
print(f"\n正对照(POWER 的 3 题劈成 2+1,题目已定向):原始 {r12:.3f}  alpha_half {a1:.3f}  去衰减 {pos:.3f}")
rgs=np.random.default_rng(11); f1=rgs.normal(0,1,len(A)); f2=rgs.normal(0,1,len(A))
S1=orient(np.c_[tuple(f1+rgs.normal(0,1,len(A)) for _ in range(3))])
S2=orient(np.c_[tuple(f2+rgs.normal(0,1,len(A)) for _ in range(2))])
n1,n2=np.nanmean(S1,axis=1),np.nanmean(S2,axis=1)
neg=abs(np.corrcoef(n1,n2)[0,1])/np.sqrt(alpha_np(S1)*alpha_np(S2))
print(f"负对照(合成独立数据,内部一致):去衰减 {neg:.3f}")

D=pd.DataFrame(rows); D.to_csv(pathlib.Path(__file__).parent/'results'/'bootstrap.csv',index=False)
g=Gate('"去衰减<=0.362 所以三条轴"这个判定是刀尖吗')
g.asserted('正对照:同一构念的两半去衰减到高值',pos>0.7,f"{pos:.3f}")
g.asserted('负对照:合成独立数据去衰减到 ~0',neg<0.15,f"{neg:.3f}")
g.asserted('原始相关的自助区间窄且全部远低于任何阈值',
           all(r['r_hi']<0.20 for r in rows),
           " ".join(f"{r['pair']}[{r['r_lo']:.3f},{r['r_hi']:.3f}]" for r in rows))
g.asserted('⚠ A02 的 0.4 是**选定**阈值,从未对过任何地板',True,
           'A02/R10 源码 `mx<0.4 -> "still three axes"`;frontier §2 第九条')
g.asserted('max 去衰减|r| 的自助区间是否跨过那个选定阈值',over>0.05,
           f"[{np.percentile(mx,2.5):.3f},{np.percentile(mx,97.5):.3f}],跨过 0.4 的比例 {over:.1%}")
g.require_resolvable_first('去衰减值相对它自己的自助展布',float(pt['G-C'][1]),float(np.std(BS['G-C'])))
g.asserted('而不确定性的来源是一个 2 题量表的 alpha',True,
           f"alpha_G = {aG0:.3f},自助 [{np.percentile(aga,2.5):.3f},{np.percentile(aga,97.5):.3f}];"
           f"去衰减因子 1/sqrt(alpha_G*{REL_C}) = {1/np.sqrt(aG0*REL_C):.2f}x")
print(g)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A120 R377 -- `起始越早,羞耻越多` 当作一条独立声明

`#331b` 的 −0.1024 是**顺带**报的,而 `#324b` 的教训是:**顺带的数最容易被高估**。

ESTIMAND        ① **逐类别** `起始年龄 ↔ 羞耻`(31 个类别)—— 普遍,还是几个类别扛的;
                ② 控制该人的**类别数**与 **`S`**;
                ③ **去衰减 + 自助 95% 区间**(信度由**类别分半**给)。
KILL            **若逐类别普遍为负且控制后保留 -> 一条可发布的声明;
                若集中在少数类别 -> 那是那几个类别的性质,不是「起始年龄」的。**
POSITIVE CTRL   合成一个**只由起始年龄驱动**的羞耻结局 -> 全流程必须复原。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠⚠ 删失        起始年龄分档到 **26yo+**,当前年龄档到 **29-32** ——
                **晚起始被压缩**。必须报**删失比例**并做一次**只用未删失者**的敏感性。
⚠ guard 21     若某一步判零,交出三件套。
IMPOSSIBLE      起始年龄是**回溯自报**,而此刻的羞耻会污染这个回忆 —— 关联,不是因果。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
NC=ONS.shape[1]; nc=np.isfinite(ONS).sum(1)
MO=np.where(nc>=5,np.nanmean(ONS,1),np.nan)
CENS=np.where(np.isfinite(ONS),ONS==28,np.nan)
cens_rate=float(np.nanmean(CENS))
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]
m0=np.isfinite(MO)&np.isfinite(sh)&np.isfinite(S)&ok
zz=lambda v,m:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
def cor(u,v,m=None):
    k=np.isfinite(u)&np.isfinite(v)&(m0 if m is None else m)
    return float(np.corrcoef(u[k],v[k])[0,1]) if k.sum()>150 else np.nan
def partial(u,y,ctrl,mm=None):
    m=(m0 if mm is None else mm).copy()
    for c in ctrl: m&=np.isfinite(c)
    X=np.column_stack([np.ones(m.sum())]+[zz(c,m) for c in ctrl])
    ru=zz(u,m)-X@np.linalg.lstsq(X,zz(u,m),rcond=None)[0]
    ry=zz(y,m)-X@np.linalg.lstsq(X,zz(y,m),rcond=None)[0]
    return float(np.corrcoef(ru,ry)[0,1]),int(m.sum())
n=int(m0.sum())
print(f"n={n:,} · 类别 {NC} · **「26yo+」删失比例 {100*cens_rate:.1f}%**")
r0=cor(MO,sh)
r1,_=partial(MO,sh,[nc.astype(float)]); r2,_=partial(MO,sh,[nc.astype(float),S])
print(f"\n② `起始均值 ↔ 羞耻`:原样 **{r0:+.4f}** · 控类别数 **{r1:+.4f}** · +控 `S` **{r2:+.4f}**")
per=[]
for j in range(NC):
    v=ONS[:,j]; k=np.isfinite(v)&np.isfinite(sh)&ok
    if k.sum()>=300: per.append(float(np.corrcoef(v[k],sh[k])[0,1]))
per=np.array(per)
print(f"① 逐类别({len(per)} 个可算):中位 **{np.median(per):+.4f}** · "
      f"为负的 **{int((per<0).sum())}/{len(per)}** · 范围 [{per.min():+.4f}, {per.max():+.4f}]")
rg=np.random.default_rng(4)
o=rg.permutation(NC); h=NC//2
half=lambda idx:np.where(np.isfinite(ONS[:,idx]).sum(1)>=3,np.nanmean(ONS[:,idx],1),np.nan)
a,b=half(o[:h]),half(o[h:2*h])
rr=cor(a,b,np.isfinite(a)&np.isfinite(b)&ok); rel=2*rr/(1+rr)
print(f"③ 类别分半信度:`corr(半A, 半B)` = **{rr:+.4f}** -> SB **{rel:.4f}** · "
      f"去衰减 **{r0/np.sqrt(rel):+.4f}**")
idx=np.flatnonzero(m0); NBOOT=400; rgb=np.random.default_rng(707)
bs=[]
for _ in range(NBOOT):
    ii=idx[rgb.integers(0,len(idx),len(idx))]
    x,y=MO[ii],sh[ii]; aa,bb=a[ii],b[ii]
    g=np.isfinite(aa)&np.isfinite(bb)
    r_=float(np.corrcoef(aa[g],bb[g])[0,1]); rl=min(2*r_/(1+r_),1.0) if r_>0 else np.nan
    ro=float(np.corrcoef(x,y)[0,1])
    bs.append(ro/np.sqrt(rl) if np.isfinite(rl) and rl>0 else np.nan)
q=np.nanpercentile(bs,[2.5,50,97.5])
print(f"   **自助 95% 区间 [{q[0]:+.4f}, {q[2]:+.4f}]** 中位 {q[1]:+.4f}")
uncens=np.where(np.isfinite(ONS),ONS<28,False).all(1)|(np.nansum(CENS,1)==0)
mu=m0&uncens
r_u,n_u=partial(MO,sh,[nc.astype(float),S],mm=mu)
print(f"\n⚠⚠ 删失敏感性:只用**一个 26yo+ 也没有**的人(n={n_u:,},占 {100*n_u/n:.0f}%):"
      f"控类别数+`S` 后 **{r_u:+.4f}**(全样本 {r2:+.4f})")
def perm_finite(v,seed):
    z2=v.copy(); j2=np.flatnonzero(np.isfinite(z2))
    z2[j2]=z2[np.random.default_rng(seed).permutation(j2)]; return z2
nul=[partial(perm_finite(MO,500+i),sh,[nc.astype(float),S])[0] for i in range(25)]
qq=float(np.mean([abs(x)>=abs(r2) for x in nul]))
print(f"负对照(打乱人):**{np.mean(nul):+.4f} ± {np.std(nul):.4f}** · |零| ≥ |观测| **{qq:.3f}**")
rgp=np.random.default_rng(88)
y=np.full(NN,np.nan); y[m0]=-0.15*zz(MO,m0)+rgp.standard_normal(n)
p_,_=partial(MO,y,[nc.astype(float),S])
print(f"正对照(只由起始年龄驱动,真 −0.15):控制后读出 **{p_:+.4f}**")
mde=2.8/np.sqrt(max(n,1))
# ⚠ #300a:上页面前发明一个能弄坏它的旋钮 —— **分档中点是我编的**(19-25yo -> 22,26yo+ -> 28)。
#    换成**序号**(1..10,不假设任何间距)与**另一套中点**,看结论动不动。
RANK={k:i+1.0 for i,k in enumerate(BIN)}
ALT={'0-4yo':2.5,'5-6yo':6,'7-8yo':8,'9-10yo':10,'11-12yo':12,'13-14yo':14,
     '15-16yo':16,'17-18yo':18,'19-25yo':21,'26yo+':32}
KN=[]
for tag,MAP in (('我的中点',BIN),('序号 1..10',RANK),('另一套中点(26yo+ -> 32)',ALT)):
    O2=np.column_stack([d[c].map(MAP).values.astype(float) for c in ons])
    M2=np.where(np.isfinite(O2).sum(1)>=5,np.nanmean(O2,1),np.nan)
    v,_=partial(M2,sh,[nc.astype(float),S],mm=np.isfinite(M2)&np.isfinite(sh)&np.isfinite(S)&ok)
    KN.append((tag,v)); print(f"   旋钮 · {tag:<24} 控类别数+`S` 后 **{v:+.4f}**")
kv=[v for _,v in KN]
print(f"   -> 跨三套编码:**{min(kv):+.4f} … {max(kv):+.4f}**(极差 {max(kv)-min(kv):.4f};2×MDE {2*mde:.4f})")

T=pd.DataFrame([dict(v_arm='原样',v_r=r0),dict(v_arm='控类别数',v_r=r1),dict(v_arm='控类别数+S',v_r=r2),
                dict(v_arm='去衰减',v_r=float(q[1])),dict(v_arm='未删失',v_r=r_u)])
check_columns(T,'R377'); T.to_csv(pathlib.Path(__file__).parent/'results'/'onset.csv',index=False)
mde=2.8/np.sqrt(max(n,1))
print(f"\n发明的旋钮(起始年龄的数值编码):")
gg=Gate('起始越早羞耻越多,当独立声明')
gg.asserted('★ 发明的旋钮:换掉我编的分档中点,结论动不动',(max(kv)-min(kv))<2*mde,
            ' · '.join(f"{t} {v:+.4f}" for t,v in KN)+f" -> 极差 {max(kv)-min(kv):.4f}")
gg.asserted('★ 正对照:真 −0.15 必须被复原(|读出| > 0.08)',abs(p_)>0.08,f"读出 {p_:+.4f}")
gg.negative_control('★ 负对照:打乱人',float(np.mean(nul)),r2,null_spread=float(np.std(nul)),
    null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill ①:逐类别是否**普遍**为负(> 80% 为负)',
            (per<0).mean()>0.8,f"{int((per<0).sum())}/{len(per)} 为负,中位 {np.median(per):+.4f}")
gg.asserted('★ 注册的 kill ②:控类别数与 `S` 后是否保留(> 50% 且 > 2×MDE)',
            abs(r2)/max(abs(r0),1e-9)>0.5 and abs(r2)>2*mde,
            f"{r2:+.4f} / {r0:+.4f} = {100*abs(r2)/max(abs(r0),1e-9):.0f}% · 2×MDE {2*mde:.4f}")
gg.asserted('⚠⚠ 删失敏感性:只用未删失者时符号与量级是否保持',
            (r_u*r2>0) and abs(r_u)>abs(r2)*0.5,
            f"未删失 {r_u:+.4f}(n={n_u:,},{100*n_u/n:.0f}%)vs 全样本 {r2:+.4f};删失率 {100*cens_rate:.1f}%")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

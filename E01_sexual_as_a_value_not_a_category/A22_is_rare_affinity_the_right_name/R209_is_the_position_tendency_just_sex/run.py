import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A22 R209 -- 「偏爱冷门」是不是性别的影子?

#163a:「偏爱冷门」是一个可迁移的**位置**倾向 —— 把块分到内容最不像的两端仍保留 85%,
而纯位置种植保留 100%、纯内容种植塌到 −2%。

但它有一个没被问过的对手:`#101`/`#154` 说它**唯一挂得住的外部锚是性别**(+0.093),
**而性别本身也是一个"跨话题都成立"的位置变量** —— 男性若在**每个**块上都挑更冷门的,
就会产生一个完美可迁移的位置倾向,而它不是"这个人的口味",是"这个人的性别"。

    INDEPENDENT  把性别从**每个块层分数**里回归掉之后,信度与 85% 基本不变
                 -> 位置倾向是自己的东西,`#100` 的 0.46 不变
    SHADOW       塌下去 -> 「偏爱冷门」很大程度上是"男性在每个块上都挑更冷门的"的另一种说法,
                 **`#100` 的 0.46 要重新定价**

ESTIMAND        S 的分半信度(随机分半 与 内容最不相似分半),在**每块回归掉性别之后**,
                与 `#163` 的同一台机器、同一个 k、同一个 SB 对照。
IDENTIFICATION  性别按**块**回归(不是按人),所以它移除的正是"性别在这一块上的位置偏移"。
SCOPE           有性别、且 >=2k 个可算块的人。
WORLDS          INDEPENDENT / SHADOW
KILL            条件式:**种一个纯由性别驱动的位置倾向,残差化后必须塌掉** ——
                否则"残差化有效"没被证明,后面的读数不可读。
POSITIVE CTRL   见上。另加一个纯**个人**位置种植,残差化后必须**存活**(证明残差化没有误伤)。
NEGATIVE CTRL   每块独立跨人置换(`#163c` 修好的那个零)。
NOISE FLOOR     3 个分半种子。
MULTIPLICITY    k ∈ {4,5,6} x {随机, 不相似} x {去性别前, 去性别后} x 4 个臂,整格发表。
IMPOSSIBLE      性别在本 release 是二值自报;它与其它未测的位置变量(如社会经验)混在一起,
                所以"不是性别的影子"不等于"不是任何人口学变量的影子"。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A22_is_rare_affinity_the_right_name'
          /'R208_position_or_content'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('KS=[4,5,6]')[0])   # 跨轮依赖显式声明(P16)

male=pd.to_numeric(df.get('biomale'),errors='coerce').values
hasg=np.isfinite(male)
print(f"有性别的人 {hasg.sum():,}",flush=True)
# 每块回归掉性别
Br=Bc.copy(); expl=[]
for j in range(Bc.shape[1]):
    m=ok[:,j]&hasg
    if m.sum()<500: continue
    X=np.c_[np.ones(m.sum()),male[m]]
    b=np.linalg.lstsq(X,Bc[m,j],rcond=None)[0]
    r=Bc[m,j]-X@b
    expl.append(1-np.var(r)/max(np.var(Bc[m,j]),1e-12))
    Br[m,j]=r
print(f"性别在块层分数里解释掉的方差:中位 {100*np.median(expl):.2f}%  "
      f"最大 {100*np.max(expl):.2f}%",flush=True)

def rel2(base,kind,k,seed,mode='random',plant=0.,sexplant=0.):
    rg=np.random.default_rng(seed); u=rg.standard_normal(len(B))
    fam=set(order[:len(order)//2].tolist())
    P=make_perm(seed*7+11)
    A=[];Bb=[]
    for i in np.flatnonzero(nb>=2*k):
        if not hasg[i]: continue
        av=np.flatnonzero(ok[i])
        if len(av)<2*k: continue
        if mode=='random':
            p=rg.permutation(av); h1,h2=p[:k],p[k:2*k]
        else:
            lo=[j for j in av if j in fam]; hi=[j for j in av if j not in fam]
            if len(lo)<k or len(hi)<k: continue
            lo=sorted(lo,key=lambda j:pc1[j]); hi=sorted(hi,key=lambda j:-pc1[j])
            h1,h2=np.array(lo[:k]),np.array(hi[:k])
        x=base[i].copy() if kind=='real' else P[i].copy()
        if plant:    x=P[i]+plant*u[i]                       # 纯**个人**位置种植
        if sexplant: x=P[i]+sexplant*male[i]                 # 纯**性别驱动**的位置种植
        a,b2=np.nanmean(x[h1]),np.nanmean(x[h2])
        if np.isfinite(a) and np.isfinite(b2): A.append(a); Bb.append(b2)
    A=np.array(A); Bb=np.array(Bb)
    if len(A)<300: return np.nan,len(A)
    r=float(np.corrcoef(A,Bb)[0,1])
    return (2*r/(1+r) if r>-0.99 else np.nan), len(A)

KS=[4,5,6]; rows=[]
print(f"\n{'k':<4}{'分半':<11}{'去性别前':>10}{'去性别后':>10}{'保留':>7}{'性别种植前':>11}{'性别种植后':>11}{'个人种植后':>11}")
for k in KS:
    for mode in ['random','dissimilar']:
        sd_=zlib.crc32(f'{mode}{k}s'.encode())%9973
        pre,_=rel2(Bc,'real',k,sd_,mode); post,n=rel2(Br,'real',k,sd_,mode)
        sp_pre,_=rel2(Bc,'real',k,sd_+1,mode,sexplant=0.8)
        # 性别种植后:先种,再对同一批数据按块回归掉性别
        rgp=np.random.default_rng(sd_+1)
        Bs=make_perm(sd_*7+11)+0.8*np.where(hasg,male,0.)[:,None]
        Bsr=Bs.copy()
        for j in range(Bs.shape[1]):
            m=ok[:,j]&hasg
            if m.sum()<500: continue
            X=np.c_[np.ones(m.sum()),male[m]]
            Bsr[m,j]=Bs[m,j]-X@np.linalg.lstsq(X,Bs[m,j],rcond=None)[0]
        sp_post,_=rel2(Bsr,'real',k,sd_+1,mode)
        pp_post,_=rel2(Br,'real',k,sd_+2,mode,plant=0.8)
        rows.append(dict(k=k,mode=mode,n=n,pre=pre,post=post,
                         sexplant_pre=sp_pre,sexplant_post=sp_post,plant_post=pp_post))
        print(f"{k:<4}{mode:<11}{pre:>+10.4f}{post:>+10.4f}{100*post/pre:>6.0f}%"
              f"{sp_pre:>+11.4f}{sp_post:>+11.4f}{pp_post:>+11.4f}",flush=True)

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'sex_control.csv',index=False)
Dz=T[T['mode']=='dissimilar']; Rr=T[T['mode']=='random']
print(f"\n  不相似分半:去性别前 {Dz.pre.mean():+.4f} -> 后 {Dz.post.mean():+.4f}"
      f"(保留 {100*Dz.post.mean()/Dz.pre.mean():.0f}%)")
print(f"  随机分半  :去性别前 {Rr.pre.mean():+.4f} -> 后 {Rr.post.mean():+.4f}"
      f"(保留 {100*Rr.post.mean()/Rr.pre.mean():.0f}%)")
print(f"  对照:纯性别驱动的位置种植 {T.sexplant_pre.mean():+.4f} -> {T.sexplant_post.mean():+.4f}"
      f"(保留 {100*T.sexplant_post.mean()/T.sexplant_pre.mean():.0f}%)")
print(f"        纯个人位置种植(去性别后) {T.plant_post.mean():+.4f}")

sd=float((Dz.pre-Dz.post).std())
g=Gate('「偏爱冷门」是不是性别的影子')
g.asserted('正对照一:纯性别驱动的位置种植,残差化后塌掉',
           T.sexplant_post.mean()/max(T.sexplant_pre.mean(),1e-9)<0.3,
           f"{T.sexplant_pre.mean():+.4f} -> {T.sexplant_post.mean():+.4f}"
           f"(保留 {100*T.sexplant_post.mean()/T.sexplant_pre.mean():.0f}%)")
g.asserted('正对照二:纯个人位置种植,残差化后存活(残差化没误伤)',
           T.plant_post.mean()>0.8,f"{T.plant_post.mean():+.4f}")
g.require_resolvable_first('去性别前后的差',float((Dz.pre-Dz.post).mean()),sd)
g.offset_control('去性别后的信度 vs 去性别前',float(Dz.post.mean()),float(Dz.pre.mean()),sd,
                 null_kind='同一批人、同一个 k、同一种分半下未去性别的 SB 信度(不是零假设,是被解释的基准)')
g.no_sign_crossing('每个 k 上去性别后仍为正',list(Dz.post.values))
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A127 R391 -- 放宽纳入口径:n 换来多少,`c3⁻` 的信度掉多少

`#345c`:`cov>=8` 把 n 从 15,503 砍到 6,717,而这正是 `#390` 功效不够的来源。
**`cov>=6` 是 `CALIBER.md` 里已登记的旋钮 ①**(它对联合 R² 只值 0.17pp)。

ESTIMAND        `cov >= 8 / 6 / 4` 三档:各自的 **n** · `c3` 的**分半 |cos|**(`#349` 的做法,8 次劈半)·
                以及 `c3⁻ ↔ 羞耻` 与 `S ↔ 羞耻`。
KILL            **若 n 明显增加而 `c3` 的信度只掉一点 -> 宽口径可用,`#390` 可以重跑;
                若信度塌了 -> 宽口径不可用,而「这个问题在这份数据上功效不够」就是最终答案。**
POSITIVE CTRL   两半都用**全样本** C -> |cos| = 1(验证比较代码,`#349` 同款)。
NEGATIVE CTRL   随机 32 维向量 -> |cos| ≈ 0.14。
⚠ 换口径就是换对象(`#320b`)-> **同时报三档下的 `c3⁻ ↔ 羞耻` 与 `S ↔ 羞耻`**,
                确认那两个数没被口径改掉。
IMPOSSIBLE      `cov` 低的人本来就只答了少数块,他们的 `c3` 是从更少信息里估的 ——
                信度下降是**预期**的,本轮量的是**降多少**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
def loadings(rows,ref=None):
    return load_of(rows,ref=ref)
def scores(v,rows):
    m=np.zeros(NN,bool); m[rows]=True
    def prof_(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo_=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo_,np.nan); R[b]=R[b]-np.nanmean(np.where(m,R[b],np.nan))
        return R
    Ra,Rb=prof_(A),prof_(B)
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(v[:,None]*Zm).sum(0); den=(Fm*np.abs(v)[:,None]).sum(0)
    s=np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
    s=np.where(m,s,np.nan); return s
def Spos(rows):
    m=np.zeros(NN,bool); m[rows]=True
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nn=M.sum(1)
        v=np.where(nn>0,(M@rr)/np.maximum(nn,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(m&(cv>=1),ps/np.maximum(cv,1),np.nan)
def cor(u,v,m):
    k=np.isfinite(u)&np.isfinite(v)&m
    return float(np.corrcoef(u[k],v[k])[0,1]) if k.sum()>200 else np.nan
rows=[]
rg=np.random.default_rng(2929)
for TH in (8,6,4):
    mk=cov>=TH; R_=np.flatnonzero(mk&np.isfinite(sh))
    v=loadings(R_); c3=scores(v,R_); Sv=Spos(R_)
    if cor(c3,sh,mk)<0: v=-v; c3=-c3
    cs=[]
    for t in range(8):
        p=rg.permutation(R_); h=len(p)//2
        va=loadings(p[:h]); vb=loadings(p[h:2*h])
        cs.append(abs(float(va[:,]@vb)) if va.ndim==1 else np.nan)
    cs=np.array([x for x in cs if np.isfinite(x)])
    rows.append(dict(v_th=TH,n=int((mk&np.isfinite(sh)).sum()),cos=float(cs.mean()),cos_sd=float(cs.std()),
                     c3sh=cor(c3,sh,mk),Ssh=cor(Sv,sh,mk)))
    print(f"cov>={TH}: n=**{rows[-1]['n']:,}** · 分半 |cos| **{cs.mean():.4f} ± {cs.std():.4f}** · "
          f"`c3⁻↔羞耻` **{rows[-1]['c3sh']:+.4f}** · `S↔羞耻` **{rows[-1]['Ssh']:+.4f}**")
T=pd.DataFrame(rows); check_columns(T,'R391')
T.to_csv(pathlib.Path(__file__).parent/'results'/'cov.csv',index=False)
R8=np.flatnonzero((cov>=8)&np.isfinite(sh))
vF=loadings(R8)
print(f"\n正对照(两半都用全样本 C):|cos| **{abs(float(vF@vF)):.4f}**")
rgn=np.random.default_rng(7); nc=[]
for _ in range(200):
    a=rgn.standard_normal(NB); b=rgn.standard_normal(NB)
    nc.append(abs(float(a@b))/np.linalg.norm(a)/np.linalg.norm(b))
print(f"负对照(随机 32 维):|cos| **{np.mean(nc):.4f} ± {np.std(nc):.4f}**")
g8,g4=T.iloc[0],T.iloc[2]
gain=np.sqrt(g4.n/g8.n)
print(f"\n★ 从 cov>=8 到 cov>=4:n **{g8.n:,} -> {g4.n:,}**(√n 增益 **{gain:.2f}×**,"
      f"MDE 降到 **{100/gain:.0f}%**)· 分半 |cos| **{g8.cos:.4f} -> {g4.cos:.4f}** "
      f"(掉 **{100*(g8.cos-g4.cos)/g8.cos:.1f}%**)")
gg=Gate('放宽纳入口径的代价')
gg.asserted('★ 正对照:两半都用全样本 C 时 |cos| = 1',abs(float(vF@vF))>0.999,f"{abs(float(vF@vF)):.4f}")
gg.asserted('★ 负对照:随机 32 维 |cos| ≈ 0.14',abs(np.mean(nc)-0.14)<0.05,f"{np.mean(nc):.4f}")
gg.asserted('★ 注册的 kill:n 增加而信度只掉一点(掉幅 < 15%)',
            (g8.cos-g4.cos)/g8.cos<0.15,
            f"n {g8.n:,} -> {g4.n:,}(√n {gain:.2f}×)· |cos| {g8.cos:.4f} -> {g4.cos:.4f} "
            f"(掉 {100*(g8.cos-g4.cos)/g8.cos:.1f}%)")
gg.asserted('⚠ 换口径就是换对象:三档下 `c3⁻↔羞耻` 与 `S↔羞耻` 稳不稳',
            (T.c3sh.max()-T.c3sh.min())<0.03 and (T.Ssh.max()-T.Ssh.min())<0.03,
            ' · '.join(f"cov>={int(r.v_th)}: c3 {r.c3sh:+.4f} / S {r.Ssh:+.4f}" for _,r in T.iterrows()))
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

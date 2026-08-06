import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A111 R366 -- `r_true(POWER, SUBSTANCE) = +0.605` 的区间,跨不跨那两个门槛

`#320` 的 NEXT。`#143`/`R043` 的判决是:**预注册的 kill 是 `>0.70 -> 一个构念`、
`<0.45 且漂移 <30% -> 确认不同`、中间 -> UNVERIFIED**;实测 **+0.605** 落在中间带,判 UNVERIFIED。
**而那个判决是拿一个点估计去比两个门槛的。**

⚠ 两句话强度完全不同,而页面写的是前一句:
- 「**落在中间带**」= 我知道它在 0.45 与 0.70 之间;
- 「**这个设计分不开三种可能**」= 区间同时跨过两个门槛,连「在中间带」都不知道。

ESTIMAND        `r_true` 的**人层自助** 95% 区间(≥300 次),统计量与原判决同构
                (梯子上各级 `r_obs/√(rel_x·rel_y)` 的中位)。
KILL            **若区间落在 (0.45, 0.70) 内 -> 「中间带」成立,UNVERIFIED 的措辞正确;
                若区间跨过任一门槛 -> 措辞要改成「这个设计分不开」,那是一句更弱的话。**
POSITIVE CTRL   **`SUBSTANCE` 对自己**(不相交两半)的 `r_true` 区间必须**覆盖 1.0** ——
                一个连自己对自己都复原不了 1.0 的流程,它报的 0.605 不可信(`R043` 的原话)。
NEGATIVE CTRL   噪声对噪声:区间必须覆盖 0。
⚠ 报           梯子**每一级**的区间,不只报终点。
IMPOSSIBLE      自助传播抽样不确定性,不传播「指标是否平行」——
                后者由 `R043` 的漂移检验与 `#317a` 的不变性单独支持。
"""
import numpy as np, pandas as pd, warnings, hashlib, itertools
from scipy import stats
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
Ax=pd.read_csv('data/derived/agent_patient.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')

def z(s): return (s-s.mean())/(s.std()+1e-9)
pc=[c for c in Ax.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
POWER=[z(pd.to_numeric(Ax[c],errors='coerce')).values for c in pc]
SUB=[]
for qi in [7,8,9,11,83,6,10]:
    s=lg[lg.qi==qi]
    if not len(s): continue
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    lo=pd.Series(opt).str.lower()
    a=np.flatnonzero(lo.str.contains(r'(myself|\bmy\b)',regex=True).values)
    b=np.flatnonzero(lo.str.contains(r'(others|other )',regex=True).values)
    if len(a) and len(b):
        SUB.append(pd.Series(M[:,a].mean(1)-M[:,b].mean(1),index=ppl).reindex(df.index).values)
POWER=[p*sg[c] for p,c in zip(POWER,pc)]
rg0=np.random.default_rng(7); N=len(df)
NOISE=[rg0.normal(size=N) for _ in range(7)]; NOISE2=[rg0.normal(size=N) for _ in range(7)]
print(f"指标数 —— POWER {len(POWER)} · SUBSTANCE {len(SUB)}")
def comp(items,idx):
    M=np.column_stack([it[idx] for it in items]); return np.nanmean(M,1)
def rel_of(items,idx,seed):
    if len(items)<2: return np.nan
    r=np.random.default_rng(seed); o=r.permutation(len(items)); h=len(o)//2
    a=comp([items[i] for i in o[:h]],idx); b=comp([items[i] for i in o[h:2*h]],idx)
    m=np.isfinite(a)&np.isfinite(b)
    if m.sum()<300: return np.nan
    rr=float(np.corrcoef(a[m],b[m])[0,1]); return 2*rr/(1+rr) if rr>-1 else np.nan
def r_obs(X,Y,idx):
    a,b=comp(X,idx),comp(Y,idx); m=np.isfinite(a)&np.isfinite(b)
    return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>300 else np.nan
def ladder(X,Y,idx,seeds=(5,15,25)):
    v=[]
    for kx,ky,sd in itertools.product(range(2,len(X)+1),range(2,len(Y)+1),seeds):
        r=np.random.default_rng(sd*31+kx*7+ky)
        xs=[X[i] for i in r.choice(len(X),kx,replace=False)]
        ys=[Y[i] for i in r.choice(len(Y),ky,replace=False)]
        rx,ry=rel_of(xs,idx,sd),rel_of(ys,idx,sd)
        if not(np.isfinite(rx) and np.isfinite(ry)) or rx<=0.02 or ry<=0.02: continue
        ro=r_obs(xs,ys,idx)
        if np.isfinite(ro): v.append(ro/np.sqrt(rx*ry))
    return float(np.median(v)) if v else np.nan
ALL=np.arange(N)
pt=ladder(POWER,SUB,ALL)
print(f"点估计 `r_true(POWER, SUBSTANCE)` = **{pt:+.4f}**(`#143` 报 +0.605)")
B=300; rgb=np.random.default_rng(2718)
bs=np.array([ladder(POWER,SUB,rgb.integers(0,N,N)) for _ in range(B)])
q=np.nanpercentile(bs,[2.5,50,97.5])
print(f"\n★ 自助 {B} 次:**2.5% {q[0]:+.4f} · 中位 {q[1]:+.4f} · 97.5% {q[2]:+.4f}** "
      f"(宽度 **{q[2]-q[0]:.4f}**)")
cross45=q[0]<0.45<q[2]; cross70=q[0]<0.70<q[2]
print(f"   跨 **0.45**:{'是' if cross45 else '否'} · 跨 **0.70**:{'是' if cross70 else '否'}")
pcv=np.array([ladder(SUB[:3],SUB[3:],rgb.integers(0,N,N)) for _ in range(150)])
qp=np.nanpercentile(pcv,[2.5,50,97.5])
ngv=np.array([ladder(NOISE,NOISE2,rgb.integers(0,N,N)) for _ in range(120)])
qn=np.nanpercentile(ngv,[2.5,50,97.5])
print(f"\n正对照 `SUBSTANCE` 对自己:**[{qp[0]:+.4f}, {qp[2]:+.4f}]** 中位 {qp[1]:+.4f} · "
      f"覆盖 1.0 **{'是' if qp[0]<=1.0<=qp[2] else '否'}**")
print(f"负对照 噪声对噪声:**[{qn[0]:+.4f}, {qn[2]:+.4f}]** · 覆盖 0 **{'是' if qn[0]<=0<=qn[2] else '否'}**")
T=pd.DataFrame([dict(v_pair='POWER-SUBSTANCE',lo=q[0],mid=q[1],hi=q[2]),
                dict(v_pair='SUBSTANCE-自己(正)',lo=qp[0],mid=qp[1],hi=qp[2]),
                dict(v_pair='噪声-噪声(负)',lo=qn[0],mid=qn[1],hi=qn[2])])
check_columns(T,'R366'); T.to_csv(pathlib.Path(__file__).parent/'results'/'r_true_ci.csv',index=False)
gg=Gate('`r_true` 的区间跨不跨那两个门槛')
gg.asserted('★ 正对照:`SUBSTANCE` 对自己的区间必须覆盖 1.0',qp[0]<=1.0<=qp[2],
            f"[{qp[0]:+.4f}, {qp[2]:+.4f}] 中位 {qp[1]:+.4f}")
gg.asserted('★ 负对照:噪声对噪声的区间必须覆盖 0',qn[0]<=0<=qn[2],f"[{qn[0]:+.4f}, {qn[2]:+.4f}]")
gg.asserted('★ 注册的 kill:区间是否落在 (0.45, 0.70) 内',
            (not cross45) and (not cross70) and 0.45<q[0] and q[2]<0.70,
            f"[{q[0]:+.4f}, {q[2]:+.4f}] —— 跨 0.45 {'是' if cross45 else '否'} · "
            f"跨 0.70 {'是' if cross70 else '否'}")
gg.asserted('⚠ 边界:自助传播抽样不确定性,不传播「指标是否平行」',True,
            '后者由 `R043` 的漂移检验与 `#317a` 的不变性单独支持')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A50 R261 -- 陡峭部分是不是"晚半每人条目少"造的

`#215` 的 NEXT:Δ 是一个**人内相关**,条目越少 |rho| 的期望绝对值越大。
**而这个对手有一个已知的代数形式**:双向去均值让人内残差和 ≈ 0,`x` 又是中心化的,
所以 `corr(y,x)` 带一个约 **−1/(k−1)** 的负偏 ——
**`#214c` 里 `c=18` 那个 −0.0792 的置换零,很可能就是它。**

ESTIMAND        ① 按 `>c` 半的条目数 k 分层(6–8 · 9–12 · 13+),各层报 Δ **与该层自己的置换零**;
                ② **条目数匹配**:从 `≤16` 半抽出 k 分布与 `>18` 半相同的子样本,报它的 Δ。
KILL            **若匹配后 `≤16` 的净 Δ 逼近 `>18` 的净 Δ(差 < 2× 展布)->
                陡峭部分是条目数造的,`#214` 要整个重估。**
预测(跑之前)   若纯粹是 k 的代数偏,则**每层的净 Δ(扣掉该层自己的零)应当彼此接近**,
                且 `>18` 与 `≤16` 在同一 k 层上应当没有差别。
NEGATIVE CTRL   每层各自的题内跨人置换零(这正是 k 偏的直接度量)。
IMPOSSIBLE      k 与"晚"在数据里天然相关(晚半条目少)——
                所以匹配后的 `≤16` 子样本不是随机子集,只能用来**证伪**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N=len(V0); M=V0.shape[1]

def demean_np(A,iters=200,tol=1e-10):
    D=np.where(np.isfinite(A),A,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
def rho_vec(D,need=6):
    W=np.isfinite(D).astype(float); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    X=W*(rar0[None,:]-rb[:,None]); Y0=np.where(np.isfinite(D),D,0.0)
    yb=np.where(k>0,Y0.sum(1)/np.maximum(k,1),0.0); Yc=W*(Y0-yb[:,None])
    num=(Yc*X).sum(1); den=np.sqrt((X*X).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(len(D),np.nan); good=(k>=need)&(den>1e-12); out[good]=num[good]/den[good]
    return out, k

def half(Vm,c,late=True,shuffle=False,rng=None):
    A=np.where(np.isfinite(Vm)&((Vm>c) if late else (Vm<=c)),Vm,np.nan)
    if shuffle:
        for i in range(len(A)):
            idx=np.flatnonzero(np.isfinite(A[i]))
            if len(idx)>1: A[i,idx]=A[i,rng.permutation(idx)]
    return rho_vec(demean_np(A))

rng=np.random.default_rng(20260803)
BINS=[(6,8),(9,12),(13,99)]
rows=[]
for lab,c,late in (('>18',18,True),('≤16',16,False)):
    r,k=half(V0,c,late); rn,_=half(V0,c,late,shuffle=True,rng=np.random.default_rng(88))
    for lo,hi in BINS:
        m=np.isfinite(r)&(k>=lo)&(k<=hi); mn=np.isfinite(rn)&(k>=lo)&(k<=hi)
        if m.sum()<150: continue
        rows.append(dict(half=lab,kbin=f'{lo}–{hi if hi<99 else "+"}',n=int(m.sum()),
                         delta=float(np.mean(r[m])),null=float(np.mean(rn[mn])),
                         net=float(np.mean(r[m])-np.mean(rn[mn]))))
T=pd.DataFrame(rows); check_columns(T,'R261'); T.to_csv(pathlib.Path(__file__).parent/'results'/'strata.csv',index=False)
print(f"{'半':<6}{'k 层':<8}{'n':>7}{'Δ':>10}{'该层置换零':>11}{'净 Δ':>10}")
for _,q in T.iterrows():
    print(f"{q.half:<6}{q.kbin:<8}{q.n:>7,}{q.delta:>+10.4f}{q.null:>+11.4f}{q.net:>+10.4f}")

# 条目数匹配:从 ≤16 半抽出 k 分布与 >18 半相同的子样本
r18,k18=half(V0,18,True); r16,k16=half(V0,16,False)
n18,_=half(V0,18,True,shuffle=True,rng=np.random.default_rng(5))
n16,_=half(V0,16,False,shuffle=True,rng=np.random.default_rng(5))
m18=np.isfinite(r18); m16=np.isfinite(r16)
tgt=pd.Series(k18[m18].astype(int)).value_counts()
sel=[]
for kk,cnt in tgt.items():
    pool=np.flatnonzero(m16&(k16==kk))
    if len(pool)==0: continue
    sel.append(rng.choice(pool,min(cnt,len(pool)),replace=False))
sel=np.concatenate(sel) if sel else np.array([],int)
d16m=float(np.mean(r16[sel])); z16m=float(np.mean(n16[np.isfinite(n16)&(np.isin(np.arange(N),sel))]))
d18=float(np.mean(r18[m18])); z18=float(np.mean(n18[np.isfinite(n18)]))
print(f"\n条目数匹配:从 ≤16 半抽 {len(sel):,} 人(k 分布对齐 >18 半的 {int(m18.sum()):,} 人)")
print(f"  >18   Δ {d18:+.4f} · 零 {z18:+.4f} · **净 {d18-z18:+.4f}**")
print(f"  ≤16m  Δ {d16m:+.4f} · 零 {z16m:+.4f} · **净 {d16m-z16m:+.4f}**")
bs=[float(np.mean(r16[rng.choice(sel,len(sel),replace=True)])) for _ in range(200)]
sd=float(np.hypot(np.std(bs),0.0147))
print(f"  差(>18 净 − ≤16m 净)= {(d18-z18)-(d16m-z16m):+.4f} ± {sd:.4f}")

g=Gate('陡峭部分是不是条目数造的')
g.asserted('可判前提:分层里确实有 k 的梯度',len(T)>=4,f"{len(T)} 层")
g.asserted('k 偏确实存在且随 k 减小而变负(代数预测)',
           bool(T[T.half=='>18'].sort_values('kbin').null.iloc[0] <
                T[T.half=='>18'].sort_values('kbin').null.iloc[-1]),
           ' · '.join(f"{q.half}{q.kbin}:{q.null:+.4f}" for _,q in T.iterrows()))
g.offset_control('>18 净 Δ vs 条目数匹配的 ≤16 净 Δ',float(d18-z18),float(d16m-z16m),sd,
                 null_kind='从 ≤16 半抽出的、**条目数分布相同**的子样本 —— 不是零假设,'
                           '是"若只有 k 在起作用,>18 该落在哪"')
g.asserted('注册的 kill:匹配后两者净 Δ 不可分辨 -> 陡峭部分是条目数造的',
           abs((d18-z18)-(d16m-z16m))<=2*sd,
           f"差 {(d18-z18)-(d16m-z16m):+.4f} ± {sd:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

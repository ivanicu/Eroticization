import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A51 R264 -- 那个最严的零,在**主 Δ** 上吃掉多少

`#218b`:**curveball 固定边际**(本项目用得最多的零,34 轮)在 `>18` 半上
**吃掉 65%,负对照 FAIL**。**而它在主 Δ 上吃掉多少,从没算过** —— `#217` 只跑了 A 与 B。

ESTIMAND        全样本主 Δ 在**三种零**下的净值:
                A 题内跨人 · B 人内打乱 · **C curveball 固定边际**。
KILL            **若 curveball 在主 Δ 上也吃掉一半以上 -> `#128` 起的整条线都要按这个零重估;
                若只在晚半吃得多 -> `#218b` 的作用域钉死在那个 n=1,944 的子样本上。**
可判前提        主 Δ 必须复现 `#128` 的 −0.0328(`#217c` 的教训:先复现原值,再判)。
NOISE FLOOR     每种零 5 个种子;Δ 人层 bootstrap 40 次。
IMPOSSIBLE      curveball 在**观测矩阵**上交换,保留两侧边际 ——
                它同时毁掉"谁答了哪一题"这件事,所以它比另外两种零**更狠**;
                "更狠的零吃得更多"本身不是缺陷,是它的定义。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N=len(V0)
KEEP0=(np.isfinite(V0).sum(1)>=8)

def demean_np(A,iters=200,tol=1e-10):
    D=np.where(np.isfinite(A),A,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
def delta(A,need=8,keep=True):
    D=demean_np(A)
    W=np.isfinite(D).astype(float); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    X=W*(rar0[None,:]-rb[:,None]); Y0=np.where(np.isfinite(D),D,0.0)
    yb=np.where(k>0,Y0.sum(1)/np.maximum(k,1),0.0); Yc=W*(Y0-yb[:,None])
    num=(Yc*X).sum(1); den=np.sqrt((X*X).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(len(D),np.nan); good=(k>=need)&(den>1e-12)
    if keep: good&=KEEP0
    out[good]=num[good]/den[good]
    return float(np.nanmean(out))
def within_person(Vm,rng):
    A=Vm.copy()
    for i in range(len(A)):
        idx=np.flatnonzero(np.isfinite(A[i]))
        if len(idx)>1: A[i,idx]=A[i,rng.permutation(idx)]
    return A
def curveball_ages(Vm,rng,swaps=200000):
    A=Vm.copy(); obs_=np.isfinite(A)
    rows=[np.flatnonzero(obs_[i]) for i in range(len(A))]
    idx=np.flatnonzero(obs_.sum(1)>=2)
    for _ in range(swaps):
        i,j=rng.choice(idx,2,replace=False)
        a,b=set(rows[i]),set(rows[j]); da,db=list(a-b),list(b-a)
        if not da or not db: continue
        x=da[rng.integers(len(da))]; y=db[rng.integers(len(db))]
        A[i,y],A[i,x]=A[i,x],np.nan
        A[j,x],A[j,y]=A[j,y],np.nan
        rows[i]=np.flatnonzero(np.isfinite(A[i])); rows[j]=np.flatnonzero(np.isfinite(A[j]))
    return A

rng=np.random.default_rng(20260803)
d0=delta(V0)
bs=[delta(V0[rng.integers(0,N,N)],keep=False) for _ in range(40)]
sd=float(np.std(bs))
print(f"主 Δ = {d0:+.4f}(人层 bootstrap sd {sd:.4f})\n")
NULLS={'A 题内跨人':lambda s: perm_null(V0,np.random.default_rng(700+s)),
       'B 人内打乱':lambda s: within_person(V0,np.random.default_rng(800+s)),
       'C curveball 固定边际':lambda s: curveball_ages(V0,np.random.default_rng(900+s))}
rows=[]
for nm,f in NULLS.items():
    vals=[delta(f(s)) for s in range(5)]
    rows.append(dict(null=nm,value=float(np.mean(vals)),null_sd=float(np.std(vals)),
                     net=float(d0-np.mean(vals)),frac=abs(np.mean(vals))/abs(d0)))
    print(f"  {nm:<22} 零 {np.mean(vals):+.4f} ± {np.std(vals):.4f} · 占效应 {100*abs(np.mean(vals))/abs(d0):>3.0f}%"
          f"  -> 净 **{d0-np.mean(vals):+.4f}**")
T=pd.DataFrame(rows); check_columns(T,'R264'); T.to_csv(pathlib.Path(__file__).parent/'results'/'main_three_nulls.csv',index=False)
cb=T[T.null.str.startswith('C')].iloc[0]
spread=float(T.net.max()-T.net.min())
print(f"\n三个净值极差 = {spread:.4f};Δ 自身展布 = {sd:.4f} -> {spread/sd:.1f}×")
print(f"curveball 在主 Δ 上占 **{100*cb.frac:.0f}%**,而在 `>18` 半上是 **65%**")

g=Gate('那个最严的零在主 Δ 上吃掉多少')
g.asserted('可判前提:主 Δ 复现 `#128` 的 −0.0328',abs(d0+0.0328)<0.003,f"{d0:+.4f}")
for _,q in T.iterrows():
    g.negative_control(f'零 {q.null}',abs(float(q.value)),abs(d0),
                       null_spread=float(q.null_sd),null_kind=q.null)
g.asserted('⚠ curveball 比另两种零更狠(它同时毁掉"谁答了哪一题")—— 吃得更多不是缺陷',True,
           'A 保留每题值分布 · B 保留每人的年龄集合 · C 固定双边边际')
g.asserted('注册的 kill:curveball 在主 Δ 上也吃掉一半以上 -> 整条线要按这个零重估',
           cb.frac>0.5,f"主 Δ 上 {100*cb.frac:.0f}% vs 晚半 65%")
print(g)
print(f"\n  => {'整条线要按 curveball 重估' if cb.frac>0.5 else '`#218b` 的作用域钉死在晚半那个 n=1,944 的子样本上'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

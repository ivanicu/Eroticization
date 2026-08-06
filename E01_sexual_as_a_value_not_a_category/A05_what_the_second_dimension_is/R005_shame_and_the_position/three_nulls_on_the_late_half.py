import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A51 R263 -- 晚半上,三种零

`#217b`:零方案的分歧**只在受限的 `>18` 半上** —— 而那正是 `#214` 的 −0.224
与 `#216` 的 −0.4344 所在的地方。**两个数各自对着不同的零,而它们描述的是同一个现象。**

ESTIMAND        在 `>18` 半上用**三种零**跑同一个 Δ:
                A 题内跨人(`perm_null`)· B 人内打乱 · C **curveball 固定边际**
                (本项目用得最多的一种,34 次)。
KILL            **若三个净值彼此差异超过任一展布的 2 倍 ->「晚半有多强」在本项目里
                没有一个唯一的数,只能报区间,`#214`/`#216` 的措辞都要改成区间。**
NOISE FLOOR     每种零 5 个种子;Δ 本身人层 bootstrap 40 次。
IMPOSSIBLE      三种零**问的是三个不同的问题**,所以"它们该一致"本身不是先验 ——
                本轮判的是**它们给出的净值能不能被当作同一个数引用**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N=len(V0); Mc=V0.shape[1]
KEEP0=(np.isfinite(V0).sum(1)>=8)

def demean_np(A,iters=200,tol=1e-10):
    D=np.where(np.isfinite(A),A,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
def rho_late(A,need=6):
    D=demean_np(np.where(np.isfinite(A)&(A>18),A,np.nan))
    W=np.isfinite(D).astype(float); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    X=W*(rar0[None,:]-rb[:,None]); Y0=np.where(np.isfinite(D),D,0.0)
    yb=np.where(k>0,Y0.sum(1)/np.maximum(k,1),0.0); Yc=W*(Y0-yb[:,None])
    num=(Yc*X).sum(1); den=np.sqrt((X*X).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(len(D),np.nan); good=(k>=need)&(den>1e-12); out[good]=num[good]/den[good]
    return float(np.nanmean(out))

def within_person(Vm,rng):
    A=Vm.copy()
    for i in range(len(A)):
        idx=np.flatnonzero(np.isfinite(A[i]))
        if len(idx)>1: A[i,idx]=A[i,rng.permutation(idx)]
    return A
def curveball_ages(Vm,rng,swaps=200000):
    """固定边际:在**观测矩阵**上做 curveball 交换,再把年龄按新位置重排 ——
    保留每人条目数与每题被报人数(边际),毁掉具体配对。"""
    A=Vm.copy(); obs_=np.isfinite(A)
    rows=[np.flatnonzero(obs_[i]) for i in range(len(A))]
    idx=np.flatnonzero(obs_.sum(1)>=2)
    for _ in range(swaps):
        i,j=rng.choice(idx,2,replace=False)
        a,b=set(rows[i]),set(rows[j])
        da,db=list(a-b),list(b-a)
        if not da or not db: continue
        x=da[rng.integers(len(da))]; y=db[rng.integers(len(db))]
        A[i,y],A[i,x]=A[i,x],np.nan
        A[j,x],A[j,y]=A[j,y],np.nan
        rows[i]=np.flatnonzero(np.isfinite(A[i])); rows[j]=np.flatnonzero(np.isfinite(A[j]))
    return A

rng=np.random.default_rng(20260803)
d0=rho_late(V0)
bs=[rho_late(V0[rng.integers(0,N,N)]) for _ in range(40)]
sd=float(np.std(bs))
print(f">18 半 Δ = {d0:+.4f}(人层 bootstrap sd {sd:.4f})\n")
NULLS={'A 题内跨人':lambda s: perm_null(V0,np.random.default_rng(700+s)),
       'B 人内打乱':lambda s: within_person(V0,np.random.default_rng(800+s)),
       'C curveball 固定边际':lambda s: curveball_ages(V0,np.random.default_rng(900+s))}
rows=[]
for nm,f in NULLS.items():
    vals=[rho_late(f(s)) for s in range(5)]
    rows.append(dict(null=nm,value=float(np.mean(vals)),null_sd=float(np.std(vals)),
                     net=float(d0-np.mean(vals))))
    print(f"  {nm:<22} 零 {np.mean(vals):+.4f} ± {np.std(vals):.4f}  -> 净 **{d0-np.mean(vals):+.4f}**")
T=pd.DataFrame(rows); check_columns(T,'R263'); T.to_csv(pathlib.Path(__file__).parent/'results'/'three_nulls.csv',index=False)
spread=float(T.net.max()-T.net.min())
print(f"\n三个净值:{' · '.join(f'{v:+.4f}' for v in T.net)}")
print(f"  极差 = {spread:.4f};Δ 自身展布 = {sd:.4f} -> {spread/sd:.1f}×")

g=Gate('晚半有多强,有没有唯一的数')
g.asserted('可判前提:>18 半的 Δ 复现 `#214` 的 −0.3034',abs(d0+0.3034)<0.02,f"{d0:+.4f}")
for _,q in T.iterrows():
    g.negative_control(f'零 {q.null}',abs(float(q.value)),abs(d0),
                       null_spread=float(q.null_sd),null_kind=q.null)
g.asserted('⚠ 三种零问的是三个不同的问题 —— "它们该一致"不是先验',True,
           'A 毁人层配对保留每题值分布 · B 毁"这个人哪个类别配到哪个年龄" · C 固定双边边际')
g.asserted('注册的 kill:三个净值极差 > 2× Δ 自身展布 -> 只能报区间',spread>2*sd,
           f"极差 {spread:.4f} vs 2×sd {2*sd:.4f} = {spread/sd:.1f}×")
print(g)
print(f"\n  => {'没有唯一的数 —— 只能报区间 [' + f'{T.net.min():+.3f}, {T.net.max():+.3f}' + ']' if spread>2*sd else '三种零给出同一个数'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

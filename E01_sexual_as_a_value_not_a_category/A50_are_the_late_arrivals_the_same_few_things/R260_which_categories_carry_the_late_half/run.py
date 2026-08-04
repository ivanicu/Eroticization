import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A50 R260 -- 「晚来的」是不是同一批东西

`#214` 把曲线钉住了,但整条线一直在问**同一个量**在不同时间段上的大小。
**从没问过:那些"晚来的"兴趣,是不是同一批东西。**
`c=18` 那一半只有 1,944 人参与 —— 若他们贡献的条目集中在少数几个类别上,
那么 `−0.224` 描述的是**那几个类别**,不是"晚来的东西"。

ESTIMAND        ① 每个切点上 `>c` 半的**类别分布**(前 3 类占多少人-类别观测);
                ② **更直接的判据**:把 `>18` 半贡献最多的前 3 个类别**删掉**,Δ 还在不在。
KILL            **若 `c=18` 半的观测在前 3 类上超过一半 -> `#214` 的 −0.224
                要改写成关于那几个类别的陈述;
                或:删掉前 3 类后 Δ 塌到不可分辨 -> 同上。**
基线            31 个类别均匀分布时,前 3 类占 3/31 = **9.7%**。
NEGATIVE CTRL   打乱"哪些观测算晚"(人内打乱起始年龄)-> 集中度应回到基线附近。
POSITIVE CTRL   合成一个"晚来的只发生在 3 个类别里"的世界 -> 集中度必须接近 100%。
IMPOSSIBLE      删掉 3 个类别同时也减少了每人的条目数 -> Δ 变化里混着"更少 item"的效应;
                所以必须并排跑一个**随机删 3 类**的对照。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N=len(V0); M=V0.shape[1]
CUTS=[10,12,14,16,18]

def demean_np(A,iters=200,tol=1e-10):
    D=np.where(np.isfinite(A),A,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
def rho_vec(D,keep_cols=None,need=6):
    rr=rar0 if keep_cols is None else rar0[keep_cols]
    Dm=D if keep_cols is None else D[:,keep_cols]
    W=np.isfinite(Dm).astype(float); k=W.sum(1)
    rb=np.where(k>0,(W*rr[None,:]).sum(1)/np.maximum(k,1),0.0)
    X=W*(rr[None,:]-rb[:,None]); Y0=np.where(np.isfinite(Dm),Dm,0.0)
    yb=np.where(k>0,Y0.sum(1)/np.maximum(k,1),0.0); Yc=W*(Y0-yb[:,None])
    num=(Yc*X).sum(1); den=np.sqrt((X*X).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(len(Dm),np.nan); good=(k>=need)&(den>1e-12); out[good]=num[good]/den[good]
    return out

def counts(Vm,c,late=True):
    msk=(Vm>c) if late else (Vm<=c)
    return (np.isfinite(Vm)&msk).sum(0).astype(float)

def top3(cnt):
    s=cnt.sum(); return float(np.sort(cnt)[::-1][:3].sum()/s) if s>0 else np.nan

print(f"基线:31 个类别均匀时前 3 类占 {3/M:.1%}\n")
print(f"{'c':>4}{'>c 前3类占比':>13}{'≤c 前3类占比':>13}{'>c 观测数':>11}")
rows=[]
for c in CUTS:
    a=counts(V0,c,True); b=counts(V0,c,False)
    rows.append(dict(cut=c,late_top3=top3(a),early_top3=top3(b),late_obs=int(a.sum())))
    print(f"{c:>4}{top3(a):>13.1%}{top3(b):>13.1%}{int(a.sum()):>11,}")
T=pd.DataFrame(rows); check_columns(T,'R260'); T.to_csv(pathlib.Path(__file__).parent/'results'/'conc.csv',index=False)

rng=np.random.default_rng(20260803)
Vn=V0.copy()
for i in range(N):
    idx=np.flatnonzero(np.isfinite(Vn[i]))
    if len(idx)>1: Vn[i,idx]=Vn[i,rng.permutation(idx)]
print(f"\n负对照(人内打乱后 c=18):前 3 类占 {top3(counts(Vn,18,True)):.1%}")
Vp=V0.copy(); keep3=np.argsort(counts(V0,18,True))[::-1][:3]
mask=np.ones(M,bool); mask[keep3]=False
Vp[:,mask]=np.where(Vp[:,mask]>18,np.nan,Vp[:,mask])   # 晚条目只留在那 3 类里
print(f"正对照(晚来的只发生在 3 个类别里):前 3 类占 {top3(counts(Vp,18,True)):.1%}")

# ---- 更直接的判据:删掉前 3 类,Δ 还在不在 -----------------------------------
def delta_at(c, drop=None):
    cols=np.array([j for j in range(M) if drop is None or j not in set(drop)])
    D=demean_np(np.where(V0>c,V0,np.nan))
    return float(np.nanmean(rho_vec(D,keep_cols=cols)))
d_full=delta_at(18); d_drop=delta_at(18,drop=keep3)
rand=[delta_at(18,drop=rng.choice(M,3,replace=False)) for _ in range(20)]
print(f"\nc=18 的 Δ:全部类别 {d_full:+.4f} · 删掉前 3 贡献类 {d_drop:+.4f} · "
      f"随机删 3 类 {np.mean(rand):+.4f} ± {np.std(rand):.4f}")
bs=[]
for _ in range(60):
    i=rng.integers(0,N,N)
    D=demean_np(np.where(V0[i]>18,V0[i],np.nan))
    bs.append(float(np.nanmean(rho_vec(D,keep_cols=np.array([j for j in range(M) if j not in set(keep3)])))))
sd_drop=float(np.std(bs))
print(f"删掉前 3 类后的人层 bootstrap sd = {sd_drop:.4f}")

c18=T[T.cut==18].iloc[0]
g=Gate('晚来的是不是同一批东西')
g.asserted('正对照:合成的"晚来只发生在 3 类"必须读到接近 100%',
           top3(counts(Vp,18,True))>0.95,f"{top3(counts(Vp,18,True)):.1%}")
g.asserted('负对照:人内打乱后集中度回到基线附近',
           abs(top3(counts(Vn,18,True))-c18.late_top3)>0.05 or top3(counts(Vn,18,True))<0.35,
           f"打乱后 {top3(counts(Vn,18,True)):.1%} vs 真实 {c18.late_top3:.1%} vs 基线 {3/M:.1%}")
g.asserted('注册的 kill 之一:c=18 半的观测在前 3 类上超过一半',c18.late_top3>0.50,
           f"{c18.late_top3:.1%}")
g.resolvable('删掉前 3 贡献类之后的 Δ',d_drop,sd_drop)
g.offset_control('删前 3 贡献类 vs 随机删 3 类',d_drop,float(np.mean(rand)),sd_drop,
                 null_kind='随机删 3 个类别 —— 不是零假设,是"只因少了 3 个 item 该掉多少"')
g.asserted('注册的 kill 之二:删掉前 3 类后 Δ 塌到不可分辨',abs(d_drop)<=2*sd_drop,
           f"{d_drop:+.4f} ± {sd_drop:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

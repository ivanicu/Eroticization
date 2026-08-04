import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A51 R262 -- 账本里的「净值」,减掉的是哪一个零

`#216c`:这个项目至少有**两种置换零**(题内跨人 · 人内),
而它们在**同一个量**上给出**符号相反**的值(−0.0792 vs +0.1309)。
**而账本里没有一处写明用的是哪一个。**

ESTIMAND        ① 扫全部轮次的 `run.py`,统计各种零方案各出现多少次;
                ② 对 `#128` 的主 Δ **同时报两种零下的净值**。
KILL            **若两种零给出的净值差异超过任一自身展布的 2 倍 ->
                账本里所有"净 Δ"都必须补上零的种类。**
库改动          `negative_control` 已加 `null_kind`(标注式,不报错)——
                二百多个既有轮次不会一次全红(`#197a` 收紧名单时的同一个权衡)。
IMPOSSIBLE      静态扫描只看得见**代码里出现的名字**,看不见一个轮次实际把哪个零用在哪一句上;
                所以计数只在**存在**方向可读。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

# ---- ① 扫描 ------------------------------------------------------------------
PAT={'题内跨人(perm_null)':re.compile(r'\bperm_null\b'),
     '人内打乱(按人重排自己的值)':re.compile(r'人内打乱|rng\.permutation\(idx\)|permutation\(len\(a\)\)|shuffle=True'),
     '跨人打乱结局(y[bi]=permutation)':re.compile(r'yp\[bi\]=.*permutation|rb\.permutation\(y\[|rng\.permutation\(sh\['),
     'curveball / 固定边际':re.compile(r'curveball'),
     '标签打乱(组/稀有度)':re.compile(r'permutation\(rar|permutation\(b\)|打乱组标签|打乱稀有度')}
paths=sorted(pathlib.Path('.').glob('E01_*/A*/R*/run.py'))
rows=[]
for p in paths:
    t=p.read_text()
    hit=[k for k,rx in PAT.items() if rx.search(t)]
    if hit: # ⚠ `round` 是 pandas 方法名 —— `#197a` 收紧后的 check_columns 第三次在我自己的新轮次上抓到
        rows.append(dict(rnd=p.parents[0].name[:44],kinds='; '.join(hit),n_kinds=len(hit)))
S=pd.DataFrame(rows); check_columns(S,'R262 扫描')
S.to_csv(pathlib.Path(__file__).parent/'results'/'null_scan.csv',index=False)
print(f"扫 {len(paths)} 个 run.py,{len(S)} 个用到可识别的零方案\n")
from collections import Counter
cnt=Counter(k for r in rows for k in r['kinds'].split('; '))
for k,v in cnt.most_common(): print(f"  {v:>4}  {k}")
print(f"\n同一轮里用到 ≥2 种零方案的轮次:{int((S.n_kinds>=2).sum())}/{len(S)}")

# ---- ② 同一个量,两种零 -------------------------------------------------------
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N=len(V0)
def demean_np(A,iters=200,tol=1e-10):
    D=np.where(np.isfinite(A),A,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
# ⚠ #217b:第一版用 `need=6` 且不套 `KEEP`,算的是**另一批人**(Δ = −0.0401 而不是 −0.0328)。
#   `R173` 的 `betas` 只在 `KEEP`(≥8 个起始年龄)上取人。**same_scale 家族:同一个名字,不同的人群。**
KEEP0=(np.isfinite(V0).sum(1)>=8)
def rho_vec(D,need=8):
    W=np.isfinite(D).astype(float); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    X=W*(rar0[None,:]-rb[:,None]); Y0=np.where(np.isfinite(D),D,0.0)
    yb=np.where(k>0,Y0.sum(1)/np.maximum(k,1),0.0); Yc=W*(Y0-yb[:,None])
    num=(Yc*X).sum(1); den=np.sqrt((X*X).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(len(D),np.nan); good=(k>=need)&(den>1e-12)&KEEP0; out[good]=num[good]/den[good]
    return out
def within_person(Vm,rng):
    A=Vm.copy()
    for i in range(len(A)):
        idx=np.flatnonzero(np.isfinite(A[i]))
        if len(idx)>1: A[i,idx]=A[i,rng.permutation(idx)]
    return A
rng=np.random.default_rng(20260803)
d0=float(np.nanmean(rho_vec(demean_np(V0))))
n_item=[float(np.nanmean(rho_vec(demean_np(perm_null(V0,np.random.default_rng(700+s)))))) for s in range(5)]
n_pers=[float(np.nanmean(rho_vec(demean_np(within_person(V0,np.random.default_rng(800+s)))))) for s in range(5)]
bs=[float(np.nanmean(rho_vec(demean_np(V0[rng.integers(0,N,N)])))) for _ in range(40)]
sd=float(np.std(bs))
print(f"\n`#128` 主 Δ = {d0:+.4f}(bootstrap sd {sd:.4f})")
print(f"  零 A 题内跨人:{np.mean(n_item):+.4f} ± {np.std(n_item):.4f}  -> 净 **{d0-np.mean(n_item):+.4f}**")
print(f"  零 B 人内打乱:{np.mean(n_pers):+.4f} ± {np.std(n_pers):.4f}  -> 净 **{d0-np.mean(n_pers):+.4f}**")
diff=abs((d0-np.mean(n_item))-(d0-np.mean(n_pers)))
print(f"  两个净值之差 = {diff:.4f};任一自身展布 = {sd:.4f} -> {diff/sd:.1f}×")

T=pd.DataFrame([dict(null='A 题内跨人',value=float(np.mean(n_item)),net=float(d0-np.mean(n_item))),
                dict(null='B 人内打乱',value=float(np.mean(n_pers)),net=float(d0-np.mean(n_pers)))])
check_columns(T,'R262'); T.to_csv(pathlib.Path(__file__).parent/'results'/'two_nulls.csv',index=False)
g=Gate('账本里的净值减掉的是哪一个零')
g.asserted('可判前提:主 Δ 复现 `#128` 的 −0.0328',abs(d0+0.0328)<0.003,f"{d0:+.4f}")
g.asserted('扫描确认项目里存在 ≥2 种零方案',len(cnt)>=2,
           ' · '.join(f"{k.split('(')[0]}×{v}" for k,v in cnt.most_common()))
g.negative_control('零 A(题内跨人)',abs(float(np.mean(n_item))),abs(d0),
                   null_spread=float(np.std(n_item)),null_kind='题内跨人置换')
g.negative_control('零 B(人内打乱)',abs(float(np.mean(n_pers))),abs(d0),
                   null_spread=float(np.std(n_pers)),null_kind='人内打乱')
g.asserted('注册的 kill:两种零的净值差 > 2× 自身展布 -> 所有"净 Δ"都要补零的种类',
           diff>2*sd,f"差 {diff:.4f} vs 2×sd {2*sd:.4f} = {diff/sd:.1f}×")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

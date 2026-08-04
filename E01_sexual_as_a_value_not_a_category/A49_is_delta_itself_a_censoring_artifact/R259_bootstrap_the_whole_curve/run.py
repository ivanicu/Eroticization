import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A49 R259 -- 报曲线族,不是一条线

`#213c`:陡峭度几乎全靠 `c=18` 那个 n=1,944 的点。
`#213b`:真实切点只有约 5 个(起始年龄两年一档)。
**两条合起来说明这条曲线的形状本身没被稳住。**

ESTIMAND        5 个真实切点(10 · 12 · 14 · 16 · 18)上 `>c 半的 Δ`,
                **每次重抽人、重算整条曲线**,报**曲线族**。
KILL            **若 `c=18` 在重抽下的展布覆盖 `c=16` 的点值 ->
                「越晚越强」的陡峭部分不可分辨,`#213a` 的措辞要再退一步到
                「较晚的那一段更强,而强多少本设计定不了」。**
前置(硬)       逐人循环太慢,必须向量化;**向量化版必须与循环版逐位相同**,否则整轮作废。
NEGATIVE CTRL   人内打乱后的曲线族(同样重抽)。
IMPOSSIBLE      重抽的是**人**,而 `c` 越大参与的人越少 —— 右端的重抽展布天然更大,
                这正是本轮要量的东西,不是缺陷。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N=len(V0)
CUTS=[10,12,14,16,18]

def demean_np(A, iters=200, tol=1e-10):
    D=np.where(np.isfinite(A),A,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D

def rho_vec(D, need=6):
    """向量化:每人 corr(残差, 中心化稀有度)。"""
    W=np.isfinite(D).astype(float); k=W.sum(1)
    ok=k>=need
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    X=W*(rar0[None,:]-rb[:,None])
    Y0=np.where(np.isfinite(D),D,0.0)
    yb=np.where(k>0,Y0.sum(1)/np.maximum(k,1),0.0)
    Yc=W*(Y0-yb[:,None])
    num=(Yc*X).sum(1); den=np.sqrt((X*X).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(len(D),np.nan)
    good=ok&(den>1e-12)
    out[good]=num[good]/den[good]
    return out

def rho_loop(D, need=6):
    out=np.full(len(D),np.nan)
    for i in range(len(D)):
        idx=np.flatnonzero(np.isfinite(D[i]))
        if len(idx)<need: continue
        x=rar0[idx]-rar0[idx].mean(); y=D[i,idx]
        if x.std()<1e-9 or np.nanstd(y)<1e-9: continue
        out[i]=float(np.corrcoef(y,x)[0,1])
    return out

# ---- 硬前置:向量化必须与循环逐位相同 ----------------------------------------
Dt=demean_np(np.where(V0>14,V0,np.nan))
a,b=rho_vec(Dt),rho_loop(Dt)
m=np.isfinite(a)&np.isfinite(b)
maxdiff=float(np.nanmax(np.abs(a[m]-b[m]))); nmis=int((np.isfinite(a)!=np.isfinite(b)).sum())
print(f"向量化 vs 循环:最大差 {maxdiff:.2e} · 有效性不一致 {nmis} 人 · 共 {int(m.sum()):,}")
assert maxdiff<1e-9 and nmis==0, "向量化没有复现循环版 —— 整轮作废"
print("  ✅ 逐位相同")

def curve(idx, Vm=None):
    Vs=(V0 if Vm is None else Vm)[idx]
    return [float(np.nanmean(rho_vec(demean_np(np.where(Vs>c,Vs,np.nan))))) for c in CUTS]

base=curve(np.arange(N))
print(f"\n原曲线(>c 的 Δ):" + ' '.join(f"c={c}:{v:+.4f}" for c,v in zip(CUTS,base)))

rng=np.random.default_rng(20260803); B=60
boot=np.array([curve(rng.integers(0,N,N)) for _ in range(B)])
lo,hi=np.percentile(boot,[2.5,97.5],axis=0); md=np.median(boot,axis=0)
print(f"\n{'c':>4}{'点值':>10}{'重抽中位':>10}{'2.5%':>10}{'97.5%':>10}{'宽度':>9}")
for j,c in enumerate(CUTS):
    print(f"{c:>4}{base[j]:>+10.4f}{md[j]:>+10.4f}{lo[j]:>+10.4f}{hi[j]:>+10.4f}{hi[j]-lo[j]:>9.4f}")

Vn=perm_null(V0,np.random.default_rng(41))
bn=np.array([curve(rng.integers(0,N,N),Vn) for _ in range(20)])
print(f"\n负对照(人内打乱)曲线中位:" + ' '.join(f"{v:+.4f}" for v in np.median(bn,axis=0)))

i16,i18=CUTS.index(16),CUTS.index(18)
covers = (lo[i18] <= base[i16] <= hi[i18])
T=pd.DataFrame(dict(cut=CUTS,point=base,boot_median=md,lo=lo,hi=hi,width=hi-lo))
check_columns(T,'R259'); T.to_csv(pathlib.Path(__file__).parent/'results'/'curve_family.csv',index=False)
g=Gate('曲线的形状稳不稳')
g.asserted('硬前置:向量化与循环逐位相同',maxdiff<1e-9 and nmis==0,f"最大差 {maxdiff:.1e}")
g.asserted('可判前提:原曲线复现 `#213a`(c=14 处 −0.0704)',abs(base[CUTS.index(14)]+0.0704)<0.004,
           f"{base[CUTS.index(14)]:+.4f}")
g.negative_control('人内打乱曲线(c=18 处)',abs(float(np.median(bn,axis=0)[i18])),abs(base[i18]))
g.no_sign_crossing('重抽中位曲线同号',[float(v) for v in md])
g.asserted('⚠ c 越大参与的人越少 -> 右端展布天然更大,这是要量的东西不是缺陷',True,
           f"宽度 c=10 {hi[0]-lo[0]:.4f} -> c=18 {hi[i18]-lo[i18]:.4f}")
g.asserted('注册的 kill:c=18 的重抽区间覆盖 c=16 的点值 -> 陡峭部分不可分辨',covers,
           f"c=18 区间 [{lo[i18]:+.4f}, {hi[i18]:+.4f}] {'覆盖' if covers else '不覆盖'} c=16 的 {base[i16]:+.4f}")
print(g)
print(f"\n  => {'陡峭部分不可分辨 —— 措辞要再退一步' if covers else '陡峭部分站得住'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- c=18 的置换零是整条曲线里唯一非零的,必须扣掉再读 ------------------------
nullc=np.median(bn,axis=0)
print("\n---- 扣掉各自的置换零之后 ----")
print(f"{'c':>4}{'点值':>10}{'置换零':>10}{'零占效应':>10}{'净值':>10}")
net=[]
for j,c in enumerate(CUTS):
    frac=abs(nullc[j])/abs(base[j]); net.append(base[j]-nullc[j])
    print(f"{c:>4}{base[j]:>+10.4f}{nullc[j]:>+10.4f}{100*frac:>9.0f}%{base[j]-nullc[j]:>+10.4f}")
i16,i18=CUTS.index(16),CUTS.index(18)
print(f"\n净值下 c=18 / c=16 = {net[i18]/net[i16]:.1f}×(原始 {base[i18]/base[i16]:.1f}×)")
g2=Gate('扣掉零之后陡峭部分还在不在')
g2.asserted('c=18 的置换零是整条曲线里唯一明显非零的',
            abs(nullc[i18])>3*max(abs(v) for k,v in enumerate(nullc) if k!=i18),
            ' · '.join(f'c={c}:{v:+.4f}' for c,v in zip(CUTS,nullc)))
g2.offset_control('c=18 的点值 vs 它自己的置换零',float(base[i18]),float(nullc[i18]),
                  float((hi[i18]-lo[i18])/4),
                  null_kind='同一切点上的题内跨人置换 —— 保留缺失模式与每题值分布,只毁配对')
g2.no_sign_crossing('净值曲线同号',[float(v) for v in net])
g2.asserted('净值下陡峭部分仍在(c=18 净值 > 2× c=16 净值)',abs(net[i18])>2*abs(net[i16]),
            f"{net[i18]:+.4f} vs {net[i16]:+.4f} = {net[i18]/net[i16]:.1f}×")
print(g2)

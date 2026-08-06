import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A67 R290 -- Δ 是不是「第 8 个兴趣出现在几岁」的剂量-反应

`#244a`(D7):早期积累组 Δ = −0.14,同样多兴趣但来得晚的 Δ = +0.04。**机制完全空着。**

WORLDS          ① **窗口/剂量-反应**:Δ 随「第 8 个起始年龄」连续变化 ->
                   `#243a` 的二分只是这条连续曲线的两端
                ② **二分**:只有「早/晚」这个类别有意义,曲线在组内是平的
                ③ **都不是**:分组变量不是年龄,而是别的东西(如获得**速率**)
ESTIMAND        对所有 ≥8 个起始年龄的人算「第 8 个起始年龄」`A8`,按它分 5 层,各层算 Δ;
                并把 `A8` 作为协变量放进「早期组指示变量」的回归,看指示变量还剩多少。
KILL            **若 Δ 随 A8 单调上升(从负到正)且早期组指示变量在控制 A8 后塌掉
                -> 世界①,二分是人为的,应当改用连续量;
                若曲线在组内平、只有跨组有差 -> 世界②;
                若既不单调、控制后指示变量也仍在 -> 世界③,分组变量另有其物。**
POSITIVE CTRL   两端:① 种一个**随 A8 变化**的人效应 -> 曲线必须出现;
                ② 种一个**与 A8 无关**的均匀人效应 -> 曲线必须保持平。
NEGATIVE CTRL   人内跨人置换。
IMPOSSIBLE      `A8` 由**回忆的**起始年龄算出,`#289` 只扣掉了「爱得深 → 记得早」这一条通路。
                另外 `A8` 与总类别数相关(报得多的人更早到第 8 个),作为协变量同报。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
def rho_of(Vm):
    D=np.where(np.isfinite(Vm),Vm,np.nan)
    for _ in range(300):
        a=np.nanmean(D,0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
    W=np.isfinite(D); Z=np.where(W,D,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); ok=(k>=8)&(den>1e-12); out[ok]=num[ok]/den[ok]; return out
rho=rho_of(V0); base=np.isfinite(rho); NCAT=np.isfinite(V0).sum(1)
Vs=np.sort(np.where(np.isfinite(V0),V0,np.inf),axis=1)
A8=np.where(base,Vs[:,7],np.nan)                      # 第 8 个起始年龄
EARLY=base&(A8<=14)
print(f"n = {int(base.sum()):,};A8 范围 {np.nanmin(A8):.1f}–{np.nanmax(A8):.1f};"
      f"中位 {np.nanmedian(A8):.1f};corr(A8, 总类别数) = "
      f"{np.corrcoef(A8[base],NCAT[base])[0,1]:+.4f}")
print(f"⚠ 早期组(A8 ≤ 14)= {int(EARLY.sum()):,} 人 —— 与 `#243a` 的 2,806 "
      f"{'一致' if abs(int(EARLY.sum())-2806)<60 else '不一致(定义略有差别)'}")
rng=np.random.default_rng(20260804)
def D_of(mask, r=None):
    r=rho if r is None else r
    m=mask&np.isfinite(r); v=r[m]
    if len(v)<100: return np.nan,np.nan,int(m.sum())
    return (float(np.mean(v)),
            float(np.std([np.mean(v[i]) for i in (rng.choice(len(v),len(v),True) for _ in range(300))])),
            int(m.sum()))
qs=np.nanpercentile(A8[base],[20,40,60,80])
def strata(r=None,tag=''):
    out=[]
    edges=[-np.inf]+list(qs)+[np.inf]
    for lo,hi in zip(edges[:-1],edges[1:]):
        m=base&(A8>lo)&(A8<=hi); d=D_of(m,r)
        out.append((float(np.nanmedian(A8[m])),d[0],d[1],d[2]))
    if tag:
        print(f"  {tag:<26}"+' · '.join(f"A8≈{a:.0f}: {b:+.4f}" for a,b,_,_ in out))
    return out
print(f"\nΔ 按「第 8 个起始年龄」分 5 层:")
S=strata(tag='观测')
for a,b,s,n in S: print(f"    A8 中位 {a:>5.1f}  Δ = {b:+.4f} ± {s:.4f}  (n={n:,})")
mono=all(S[i][1]<=S[i+1][1]+0.005 for i in range(len(S)-1))
print(f"  单调不降 = **{mono}**;首末差 **{S[-1][1]-S[0][1]:+.4f}** vs 2×展布 "
      f"{2*np.hypot(S[0][2],S[-1][2]):.4f}")

# A8 作为协变量:早期组指示变量还剩多少
m=base&np.isfinite(A8)
zA=(A8[m]-A8[m].mean())/A8[m].std(); zE=EARLY[m].astype(float); zE=(zE-zE.mean())/zE.std()
zr=(rho[m]-rho[m].mean())/rho[m].std()
def fit(cols):
    X=np.column_stack([np.ones(len(zr))]+cols); b=np.linalg.lstsq(X,zr,rcond=None)[0]
    sd=[float(np.std([np.linalg.lstsq(X[i],zr[i],rcond=None)[0][j]
        for i in (rng.choice(len(zr),len(zr),True) for _ in range(200))])) for j in range(1,len(cols)+1)]
    return b,sd
bE,sE=fit([zE]); bA,sA=fit([zA]); bB,sB=fit([zE,zA])
print(f"\n  早期组指示单独  beta = {bE[1]:+.4f} ± {sE[0]:.4f}")
print(f"  A8 单独        beta = {bA[1]:+.4f} ± {sA[0]:.4f}")
print(f"  **同时:早期组 {bB[1]:+.4f} ± {sB[0]:.4f} · A8 {bB[2]:+.4f} ± {sB[1]:.4f}**(n = {int(m.sum()):,})")

# 正对照两端
zA_full=np.full(N,np.nan); zA_full[m]=zA
u1=np.where(base,np.nan_to_num(zA_full)*0.15+rng.standard_normal(N)*0.05,np.nan)
u2=np.where(base,rng.standard_normal(N)*0.15,np.nan)
p1=strata(r=np.where(base,rho+np.nan_to_num(u1),np.nan),tag='正对照①(随 A8 变化)')
p2=strata(r=np.where(base,rho+np.nan_to_num(u2),np.nan),tag='正对照②(与 A8 无关)')
g1=p1[-1][1]-p1[0][1]-(S[-1][1]-S[0][1]); g2=p2[-1][1]-p2[0][1]-(S[-1][1]-S[0][1])
print(f"  正对照增益:① {g1:+.4f}(必须 >0)· ② {g2:+.4f}(必须 ≈0)")
nul=[]
for _ in range(15):
    Vp=V0.copy()
    for j in range(Mc):
        idx=np.flatnonzero(np.isfinite(Vp[:,j])); Vp[idx,j]=Vp[rng.permutation(idx),j]
    r2=rho_of(Vp); nul.append(strata(r=r2)[-1][1]-strata(r=r2)[0][1])
print(f"  置换零的首末差 {np.mean(nul):+.4f} ± {np.std(nul):.4f}")

T=pd.DataFrame([dict(a8_median=a,delta=b,sd=s,n=n) for a,b,s,n in S])
check_columns(T,'R290'); T.to_csv(pathlib.Path(__file__).parent/'results'/'delta_by_a8.csv',index=False)

g=Gate('Δ 是不是「第 8 个起始年龄」的剂量-反应')
g.asserted('正对照两端:随 A8 变化的种入必须造出曲线,与 A8 无关的必须不造',
           g1>2*np.std(nul) and abs(g2)<abs(g1)/2, f"① {g1:+.4f} · ② {g2:+.4f} · 零展布 {np.std(nul):.4f}")
g.negative_control('置换零的首末差',abs(float(np.mean(nul))),abs(S[-1][1]-S[0][1]),
                   null_spread=float(np.std(nul)),null_kind='题内跨人置换起始年龄 —— 只打掉配对')
g.offset_control('★ 早期组指示在控制 A8 后 vs 单独',float(bB[1]),float(bE[1]),float(np.hypot(sB[0],sE[0])),
                 null_kind='它单独时的 beta —— 不是零假设,是「若二分本身有意义,控制 A8 后它该落在哪」')
g.asserted('★ 注册的 kill:Δ 随 A8 单调 且 早期组指示控制后塌掉 -> 世界①(连续剂量-反应)',
           mono and abs(bB[1])<2*sB[0],
           f"单调={mono};首末差 {S[-1][1]-S[0][1]:+.4f};早期组指示 {bE[1]:+.4f} -> {bB[1]:+.4f} ± {sB[0]:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

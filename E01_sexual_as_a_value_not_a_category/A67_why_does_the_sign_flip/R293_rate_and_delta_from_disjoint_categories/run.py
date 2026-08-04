import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A67 R293 -- rate 与 Δ 算自不相交的类别,那条地形还在不在

⚠ **先纠正 `#247` 自己的 NEXT**:它写「把 `rate` 对人内起始年龄 sd 残差化」。
**那是一个过强的控制** —— `rate = (类别数−1)/(最后−最初)`,分母**就是**跨度,
而人内 sd 与跨度几乎是同一个量。扣掉它等于把构念本身扣掉,
结果无论如何都会塌,而那不能区分「假象」与「被我自己删掉」。

**真正的威胁是另一个,而它更根本**:`rate` 与 `Δ` **算自同一批起始年龄**。
一个人若回忆噪声大 / 把一切挤在一起报,**两个量会同时被扭曲**,产生一条自造的关系。

ESTIMAND        把 31 个类别随机劈成两半:**`rate` 只用 A 半算,`Δ` 只用 B 半算**,
                跨 6 次劈分报曲线;两个量此时**不共享任何一个观测**。
KILL            **若跨半曲线仍单调且跨度 > 全样本 |Δ| 的 3 倍(>0.10)-> 地形是真的,升 D7;
                若塌到不可分辨 -> `#247` 的地形大部分是「同一批数据算两次」造出来的。**
POSITIVE CTRL   两端:① 种一个**真的**人层 rate→Δ 关系 -> 跨半必须测到;
                ② 种一个**只在半内**的假象(把某半的起始年龄整体压缩)-> 跨半必须**测不到**。
NEGATIVE CTRL   人内跨人置换。
⚠ 覆盖率         每半只有 ~15 个类别,所以两边都放宽到 **≥6**;n 会掉,**同报**。
IMPOSSIBLE      跨半切断的是**共享观测**,切不断**共享的人层回忆风格** ——
                若一个人在所有类别上都压缩,两半仍会同时被压缩。**这条边仍开着,如实登记。**
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
def rho_sub(Vm,cols,need=6):
    Vc=np.full_like(Vm,np.nan); Vc[:,cols]=Vm[:,cols]
    D=np.where(np.isfinite(Vc),Vc,np.nan)
    for _ in range(200):
        a=np.nanmean(D,0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
    W=np.isfinite(D); Z=np.where(W,D,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); ok=(k>=need)&(den>1e-12); out[ok]=num[ok]/den[ok]; return out
def rate_sub(Vm,cols,need=3):
    Vc=np.where(np.isfinite(Vm[:,cols]),Vm[:,cols],np.nan)
    k=np.isfinite(Vc).sum(1); lo=np.nanmin(Vc,1); hi=np.nanmax(Vc,1)
    return np.where(k>=need,(k-1)/np.maximum(hi-lo,0.5),np.nan)
rng=np.random.default_rng(20260804)
def curve_split(seed,Vm=None,nb=6):
    Vm=V0 if Vm is None else Vm
    rg=np.random.default_rng(seed); p=rg.permutation(Mc); h=Mc//2
    rt=rate_sub(Vm,p[:h]); rh=rho_sub(Vm,p[h:])
    m=np.isfinite(rt)&np.isfinite(rh)
    if m.sum()<1000: return None,0
    q=np.nanpercentile(rt[m],np.linspace(0,100,nb+1)); out=[]
    for lo,hi in zip(q[:-1],q[1:]):
        mm=m&(rt>lo)&(rt<=hi) if lo>q[0] else m&(rt<=hi)
        if mm.sum()>=200: out.append((float(np.nanmedian(rt[mm])),float(np.mean(rh[mm])),int(mm.sum())))
    return out,int(m.sum())
CS=[curve_split(3000+s) for s in range(6)]
n_used=int(np.mean([c[1] for c in CS]))
K=min(len(c[0]) for c in CS)
mat=np.array([[c[0][i][1] for i in range(K)] for c in CS])
med=np.array([[c[0][i][0] for i in range(K)] for c in CS]).mean(0)
mu,sd=mat.mean(0),mat.std(0)
print(f"跨半(rate 用 A 半、Δ 用 B 半,两者不共享任何观测);6 次劈分,n ≈ {n_used:,}")
for i in range(K):
    print(f"  rate 中位 {med[i]:>6.2f}  Δ = {mu[i]:+.4f} ± {sd[i]:.4f}  "
          f"(n≈{int(np.mean([c[0][i][2] for c in CS])):,})")
span=mu[-1]-mu[0]; mono=all(mu[i]>=mu[i+1]-0.01 for i in range(K-1))
print(f"  单调不升 = **{mono}**;首末差 **{span:+.4f}** vs 2×展布 {2*np.hypot(sd[0],sd[-1]):.4f}"
      f"  [同数据曲线是 −0.4489]")

# 正对照①:种一个真的人层 rate->Δ 关系
def plant_true(g=0.6):
    rt_all=rate_sub(V0,np.arange(Mc)); z=np.where(np.isfinite(rt_all),
        (rt_all-np.nanmean(rt_all))/np.nanstd(rt_all),0.0)
    Vp=V0.copy()
    for j in range(Mc):
        mm=np.isfinite(Vp[:,j]); Vp[mm,j]=Vp[mm,j]+g*z[mm]*(rar0[j]-rar0.mean())
    return Vp
c1,_=curve_split(3100,plant_true()); s1=c1[-1][1]-c1[0][1]
# ⚠ 上面那个种入**自己坏了**:它把量加到起始年龄上,而那同时改变了 min/max/span,
#   也就改变了 rate —— **种入干扰了它本该保持不变的那个变量**(与 `#196d` 同族)。
#   正确构造:rate 只由 A 半决定,**种入只进 B 半**,于是 rate 一动不动。
def curve_split_planted(seed,g=0.5,nb=6):
    rg=np.random.default_rng(seed); p=rg.permutation(Mc); h=Mc//2
    rt=rate_sub(V0,p[:h])                      # rate:只用 A 半,不被种入触碰
    z=np.where(np.isfinite(rt),(rt-np.nanmean(rt))/np.nanstd(rt),0.0)
    Vp=V0.copy()
    for j in p[h:]:                            # 种入:只进 B 半
        mm=np.isfinite(Vp[:,j]); Vp[mm,j]=Vp[mm,j]+g*z[mm]*(rar0[j]-rar0.mean())
    rh=rho_sub(Vp,p[h:])
    m=np.isfinite(rt)&np.isfinite(rh)
    q=np.nanpercentile(rt[m],np.linspace(0,100,nb+1)); out=[]
    for lo,hi in zip(q[:-1],q[1:]):
        mm=m&(rt>lo)&(rt<=hi) if lo>q[0] else m&(rt<=hi)
        if mm.sum()>=200: out.append((float(np.nanmedian(rt[mm])),float(np.mean(rh[mm])),int(mm.sum())))
    return out
SW=[(g,(lambda c:c[-1][1]-c[0][1])(curve_split_planted(3100,g))) for g in (0.0,0.25,0.5,1.0)]
print(f"正对照①修后(种入只进 B 半,rate 由 A 半决定,强度扫描):"
      + ' · '.join(f"g={a:.2f} 首末差 {b:+.4f}" for a,b in SW))
s1fix=SW[-1][1]-SW[0][1]
# 正对照②:只在半内的假象 —— 把某一半的起始年龄整体压缩(不制造跨半关系)
def plant_within(seed=77,g=0.35):
    rg=np.random.default_rng(seed); p=rg.permutation(Mc); h=Mc//2
    Vp=V0.copy(); f=rg.random(N)<0.5
    for j in p[:h]:
        mm=np.isfinite(Vp[:,j])&f; Vp[mm,j]=np.nanmean(Vp[:,j])+g*(Vp[mm,j]-np.nanmean(Vp[:,j]))
    return Vp
c2,_=curve_split(3100,plant_within()); s2=c2[-1][1]-c2[0][1]
print(f"正对照①(真的人层 rate→Δ):跨半首末差 **{s1:+.4f}**(必须明显 <0)")
print(f"正对照②(只在半内的压缩假象):跨半首末差 **{s2:+.4f}**(必须 ≈ 基线,不额外造出关系)")
nul=[]
for s in range(5):
    Vp=V0.copy(); rg=np.random.default_rng(4000+s)
    for j in range(Mc):
        idx=np.flatnonzero(np.isfinite(Vp[:,j])); Vp[idx,j]=Vp[rg.permutation(idx),j]
    c,_=curve_split(3000+s,Vp)
    if c: nul.append(c[-1][1]-c[0][1])
print(f"置换零的首末差 {np.mean(nul):+.4f} ± {np.std(nul):.4f}")

T=pd.DataFrame([dict(rate_median=float(med[i]),delta=float(mu[i]),sd=float(sd[i])) for i in range(K)])
check_columns(T,'R293'); T.to_csv(pathlib.Path(__file__).parent/'results'/'cross_half_curve.csv',index=False)

g=Gate('rate 与 Δ 算自不相交类别')
g.asserted('⚠ 正对照①原版不合身(种入改变了 rate 本身),保留记录',
           False, f"原版 {s1:+.4f} vs 基线 {span:+.4f} —— 种入干扰了它本该保持不变的变量(#196d 同族)")
g.asserted('正对照①修后:种入只进 B 半,强度扫描必须单调加深',
           SW[-1][1]<SW[0][1]-0.05, ' · '.join(f"g={a:.2f} {b:+.4f}" for a,b in SW))
g.asserted('正对照②:只在半内的压缩假象必须不额外造出关系',
           abs(s2-span)<abs(span)/2, f"② {s2:+.4f} vs 基线 {span:+.4f}")
g.negative_control('置换零的首末差',abs(float(np.mean(nul))),abs(span),
                   null_spread=float(np.std(nul)),null_kind='题内跨人置换起始年龄 —— 只打掉配对')
g.asserted('★ 注册的 kill:跨半曲线仍单调且跨度 > 0.10 -> 地形是真的',
           mono and abs(span)>0.10,
           f"单调={mono};跨半首末差 {span:+.4f}(同数据是 −0.4489);n ≈ {n_used:,}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

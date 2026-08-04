import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A109 R362 -- 去衰减是一个校正,还是一次观察

`#316a` 判 Ⓐ 赢,理由是**去衰减后**信号出现(+0.1293)。
**但去衰减是一个校正:它假设「若信度更高,相关就会是那么大」,而这个假设本身没被测过。**
这个项目里**每一个**去衰减数字都靠它。

⚠ **先排除一个陷阱**:「合成分与它自己两个题目的相关」是一个**代数恒等式**
(`r_合成 = (r1+r2)/√(2(1+r_12))`,`#267b` / guard 14 的坑)——**那不是检验。**

**真正能检验它的**:**人为把信度压下去**(给两个指标各加独立噪声),
信度与原始相关都会掉,而**去衰减值应当不变**。
**曲线平 -> 去衰减在这份数据上是一次观察;曲线不平 -> 这个校正在这里不成立。**

ESTIMAND        噪声强度 λ ∈ {0, 0.5, 1.0, 1.5, 2.0}:
                实测 Spearman–Brown 信度 · 原始 `corr(form_λ, 羞耻)` · **去衰减值**。
KILL            **若去衰减值跨 λ 的展布 < 它自身的种子展布 -> 校正成立,`#316a` 从校正变成观察;
                若系统性漂移 -> 校正在这里不成立,而这会影响项目里每一个去衰减数字。**
POSITIVE CTRL   **完全合成**的一对指标(真相关已知)-> 同一条流程必须把去衰减值稳定在真值上。
NEGATIVE CTRL   `perm_finite` 打乱人 -> 去衰减值必须在零附近(且**不随 λ 系统漂移**)。
⚠ 多种子       每个 λ 用 **8 个噪声种子**,报均值 ± 展布(`#314b`:单种子会造出假的非单调)。
IMPOSSIBLE      加噪只压**信度**,不改真相关 —— 所以本轮检的是校正对**已知方向**的响应,
                不是「真信度更高时会怎样」。后者需要更多指标,这份问卷没有。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ani=pd.to_numeric(d['animated'],errors='coerce').values.astype(float)
wri=pd.to_numeric(d['written'],errors='coerce').values.astype(float)

m0=np.isfinite(ani)&np.isfinite(wri)&np.isfinite(sh)
zz=lambda v:(v[m0]-v[m0].mean())/v[m0].std()
A0,W0,Y=zz(ani),zz(wri),zz(sh); n=len(Y)
print(f"n={n:,} · `corr(animated, written)` = **{np.corrcoef(A0,W0)[0,1]:+.4f}**")
def one(lam,seed,a=None,w=None,y=None):
    rg=np.random.default_rng(seed)
    a=(A0 if a is None else a)+lam*rg.standard_normal(n)
    w=(W0 if w is None else w)+lam*rg.standard_normal(n)
    y=Y if y is None else y
    r=float(np.corrcoef(a,w)[0,1]); rel=2*r/(1+r) if r>0 else np.nan
    f=(a-a.mean())/a.std()+(w-w.mean())/w.std()
    raw=float(np.corrcoef(f,y)[0,1])
    return rel,raw,(raw/np.sqrt(rel) if np.isfinite(rel) and rel>0 else np.nan)
LAMS=(0.0,0.5,1.0,1.5,2.0); NS=8
print(f"\n{'λ':>5}{'信度':>10}{'原始 r':>12}{'去衰减':>12}{'(± 种子展布)':>14}")
ROWS=[]
for lam in LAMS:
    v=np.array([one(lam,3000+100*i) for i in range(NS)])
    ROWS.append(dict(lam=lam,rel=v[:,0].mean(),raw=v[:,1].mean(),dis=v[:,2].mean(),
                     dsd=v[:,2].std(),rawsd=v[:,1].std()))
    print(f"{lam:>5.1f}{v[:,0].mean():>10.4f}{v[:,1].mean():>+12.4f}{v[:,2].mean():>+12.4f}"
          f"{v[:,2].std():>14.4f}")
T=pd.DataFrame(ROWS); check_columns(T,'R362')
T.to_csv(pathlib.Path(__file__).parent/'results'/'reliability_sweep.csv',index=False)
spread=float(T.dis.max()-T.dis.min()); seed_sd=float(T.dsd.mean())
raw_drop=float(T.raw.iloc[0]-T.raw.iloc[-1])
print(f"\n★ 去衰减值跨 λ 的展布 **{spread:.4f}** vs 各 λ 自身的种子展布均值 **{seed_sd:.4f}** "
      f"-> **{spread/max(seed_sd,1e-9):.2f}×**")
print(f"   同期原始 r 掉了 **{raw_drop:.4f}**(从 {T.raw.iloc[0]:+.4f} 到 {T.raw.iloc[-1]:+.4f}),"
      f"信度从 {T.rel.iloc[0]:.4f} 掉到 {T.rel.iloc[-1]:.4f}")
rgS=np.random.default_rng(11); LAT=rgS.standard_normal(n); TRUE=0.20
Ys=TRUE*LAT+np.sqrt(1-TRUE**2)*rgS.standard_normal(n)
As=LAT+0.8*rgS.standard_normal(n); Ws=LAT+0.8*rgS.standard_normal(n)
As=(As-As.mean())/As.std(); Ws=(Ws-Ws.mean())/Ws.std(); Ys=(Ys-Ys.mean())/Ys.std()
print(f"\n正对照:完全合成的一对指标,真相关 **{TRUE:.2f}**")
PC=[]
for lam in LAMS:
    v=np.array([one(lam,5000+100*i,a=As,w=Ws,y=Ys) for i in range(NS)])
    PC.append(v[:,2].mean()); print(f"   λ={lam:.1f}: 信度 {v[:,0].mean():.4f} · "
          f"原始 {v[:,1].mean():+.4f} · **去衰减 {v[:,2].mean():+.4f}**")
pc_sp=float(max(PC)-min(PC)); pc_bias=float(np.mean(PC)-TRUE)
print(f"   -> 跨 λ 展布 **{pc_sp:.4f}** · 相对真值的偏差 **{pc_bias:+.4f}**")
def perm_finite(v,seed):
    z2=v.copy(); z2=z2[np.random.default_rng(seed).permutation(len(z2))]; return z2
NG=[]
for lam in LAMS:
    v=np.array([one(lam,7000+100*i,y=perm_finite(Y,900+i)) for i in range(NS)])
    NG.append(v[:,2].mean())
print(f"负对照(打乱人):去衰减值 " + ' · '.join(f"{x:+.4f}" for x in NG) +
      f" -> 展布 **{max(NG)-min(NG):.4f}**,均值 **{np.mean(NG):+.4f}**")
gg=Gate('去衰减是一个校正,还是一次观察')
gg.asserted('★ 正对照:合成一对指标(真相关 0.20)-> 去衰减值必须稳定在真值上',
            pc_sp<0.05 and abs(pc_bias)<0.05,
            f"跨 λ 展布 {pc_sp:.4f} · 偏差 {pc_bias:+.4f}")
gg.asserted('★ 负对照:打乱人后去衰减值必须在零附近且不随 λ 漂',
            abs(np.mean(NG))<0.03 and (max(NG)-min(NG))<0.05,
            f"均值 {np.mean(NG):+.4f} · 展布 {max(NG)-min(NG):.4f}")
gg.asserted('⚠ 原始 r 确实随信度掉了(否则这一轮什么也没压到)',raw_drop>0.02,
            f"原始 r 从 {T.raw.iloc[0]:+.4f} 掉到 {T.raw.iloc[-1]:+.4f}(信度 "
            f"{T.rel.iloc[0]:.4f} -> {T.rel.iloc[-1]:.4f})")
gg.asserted('★ 注册的 kill:去衰减值跨 λ 的展布 < 自身种子展布',spread<seed_sd,
            f"跨 λ 展布 **{spread:.4f}** vs 种子展布 **{seed_sd:.4f}** "
            f"({spread/max(seed_sd,1e-9):.2f}×)")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

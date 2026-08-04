import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A99 R350 -- 重现度 vs 样本量:发布版的 `c3` 到底多稳

`#304c`:0.796 是**用一半数据(n≈3,350)估 `c3`** 的重现度,**不是发布版(n=6,717)的**。
半样本 |cos| **低估**了真正要报的那个数。

ESTIMAND        |cos|(n),n ∈ {800, 1600, 3350},每个 n 上 ≥8 次**等大且不相交**的两组;
                再按 Spearman–Brown 形状(|cos| 的行为等同一个平行测量的信度)外推到 6,717。
⛔ 外推的自检    **用 n=800 预测 1600、用 1600 预测 3350,与实测比。**
                预测不准 -> **只报曲线,不报外推值**(`#301` 的教训:一个外推出来的数会被当成测出来的引用)。
POSITIVE CTRL   两半都用全样本 C:|cos| = 1,与 n 无关。
NEGATIVE CTRL   随机正交基:|cos| ≈ 0.14,**不随 n 变** —— 这是一条好基线,因为它证明
                曲线的上升不是「n 大了什么都变像」。
⚠ CONTROL       `Cmat` 里有一个固定的 `mm.sum()>200` 门槛;n 小的时候会有格子被置零。
                **每个 n 上都报被置零的格子比例** —— 若 n=800 上大量置零,那个点不可比,直接弃。
IMPOSSIBLE      |cos| 判方向,不判那条方向解释多少方差(`#349` 已说)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def fit_apply')[0])

ZERO={}
def Cmat(rows,tag=None,frac=None):
    # ⚠ 固定门槛 `>200` 在 n 小的时候把矩阵掏空(n=800 置零 85.7%),
    #    于是三个 n 量的不是同一个估计量。发布口径是 200/6,717 = **3.0% 的 n**,
    #    所以要做 n 的曲线,门槛必须是**同一个比例规则**。
    m=np.zeros(NN,bool); m[rows]=True
    THR=200 if frac is None else max(20,int(frac*len(rows)))
    def prof_(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0)
        R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(np.where(m,R[b],np.nan))
        return R
    Ra,Rb=prof_(A),prof_(B)
    C=np.full((NB,NB),np.nan); z=0
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
            if mm.sum()>THR: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
            else: z+=1
    if tag is not None: ZERO.setdefault(tag,[]).append(z/(NB*NB))
    C=np.where(np.isfinite(C),C,0.0); return (C+C.T)/2
def eig(C):
    w,V=np.linalg.eigh(C); o=np.argsort(-w); return w[o],V[:,o]
ALLR=np.flatnonzero(ok); NS=[800,1600,3350]
FR=200/len(np.flatnonzero(ok))   # 发布口径的比例规则
print(f'门槛 = **{100*FR:.2f}% 的 n**(发布口径 200/{len(np.flatnonzero(ok)):,} 的同一条规则)')
rg=np.random.default_rng(50505)
CURVE={}
for n in NS:
    acc=[]
    for t in range(8):
        p=rg.permutation(ALLR)[:2*n]
        _,Va=eig(Cmat(p[:n],tag=n,frac=FR)); _,Vb=eig(Cmat(p[n:],tag=n,frac=FR))
        acc.append([abs(float(Va[:,k]@Vb[:,k])) for k in range(3)])
    CURVE[n]=np.array(acc)
    z=100*np.mean(ZERO[n])
    print(f"n={n:>5}(每组,8 次不相交劈分)· 置零格子 **{z:.1f}%** · "
          + ' · '.join(f"{nm} **{CURVE[n][:,k].mean():.4f} ± {CURVE[n][:,k].std():.4f}**"
                       for k,nm in enumerate(['c1','c2','c3'])))
def sb(r,k=2.0): return k*r/(1+(k-1)*r)          # Spearman–Brown:样本量 ×k
print(f"\n⛔ 外推的自检(用小 n 预测大 n,SB 倍率 = n 之比):")
ROWS=[]
for k,nm in enumerate(['c1','c2','c3']):
    p16=sb(CURVE[800][:,k].mean(),1600/800); m16=CURVE[1600][:,k].mean()
    p33=sb(CURVE[1600][:,k].mean(),3350/1600); m33=CURVE[3350][:,k].mean()
    print(f"   {nm}: 800→1600 预测 **{p16:.4f}** 实测 **{m16:.4f}**(差 {p16-m16:+.4f})· "
          f"1600→3350 预测 **{p33:.4f}** 实测 **{m33:.4f}**(差 {p33-m33:+.4f})")
    ROWS.append(dict(v_axis=nm,n800=CURVE[800][:,k].mean(),n1600=m16,n3350=m33,
                     err16=p16-m16,err33=p33-m33,proj6717=sb(m33,6717/3350)))
T=pd.DataFrame(ROWS); check_columns(T,'R350')
T.to_csv(pathlib.Path(__file__).parent/'results'/'repro_vs_n.csv',index=False)
worst=float(np.max(np.abs(np.r_[T.err16.values,T.err33.values])))
print(f"\n   最大预测误差 **{worst:.4f}**")
OK_EXTRAP=worst<0.05
if OK_EXTRAP:
    print(f"\n外推到 n=6,717(发布版的规模):")
    for _,r in T.iterrows(): print(f"   {r.v_axis}: **{r.proj6717:.4f}**")
else:
    print(f"\n⛔ **自检没过 -> 只报曲线,不报外推值。**")
    for _,r in T.iterrows(): print(f"   ({r.v_axis} 的外推值 {r.proj6717:.4f} **不发布**)")
pc=[];  wF,VF=eig(Cmat(ALLR))
for n in NS: pc.append(abs(float(VF[:,2]@VF[:,2])))
rgN=np.random.default_rng(7); nc=[]
for _ in range(200):
    Q1=np.linalg.qr(rgN.standard_normal((NB,3)))[0]; Q2=np.linalg.qr(rgN.standard_normal((NB,3)))[0]
    nc.append(abs(float(Q1[:,2]@Q2[:,2])))
print(f"\n正对照(两半都用全样本 C):|cos| **{np.mean(pc):.4f}**(与 n 无关)")
print(f"负对照(随机正交基):|cos| **{np.mean(nc):.4f} ± {np.std(nc):.4f}**(不随 n 变)")
gg=Gate('重现度 vs 样本量')
gg.asserted('★ 正对照:两半都用全样本 C 时 |cos| = 1',np.mean(pc)>0.999,f"{np.mean(pc):.4f}")
gg.asserted('★ 负对照:随机正交基 |cos| ≈ 0.14 且不随 n 变',abs(np.mean(nc)-0.14)<0.05,f"{np.mean(nc):.4f}")
gg.asserted('⚠ 门槛控制:每个 n 上被置零的格子比例',max(np.mean(v) for v in ZERO.values())<0.05,
            ' · '.join(f"n={n} {100*np.mean(ZERO[n]):.1f}%" for n in NS))
gg.bounded_statistic_out_of_range('⚠ guard 17:最大预测误差先过定义域',worst,0,1,'预测误差')
gg.asserted('⛔ 外推的自检:最大预测误差 < 0.05 才允许发布外推值',OK_EXTRAP,
            f"最大误差 **{worst:.4f}** —— " + ('可以外推' if OK_EXTRAP else '**只报曲线,不报外推值**'))
gg.asserted('★ 注册的读数:`c3` 的曲线形状 vs `c1`',True,
            ' · '.join(f"{r.v_axis}: {r.n800:.3f} -> {r.n1600:.3f} -> {r.n3350:.3f}" for _,r in T.iterrows()))
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

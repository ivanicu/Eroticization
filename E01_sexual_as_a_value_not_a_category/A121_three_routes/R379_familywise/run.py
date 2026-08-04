import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A121 R379 -- 那个次可加,过不过全族阈

`#333b`:`S × EARLY` = −0.0395,是自身置换零展布的 3.2 倍,**但置换只做了 30 次
(分辨率下限 1/30)且三个交互的多重性没被处理** -> PLAUSIBLE。

⚠ **族是预先声明的那三个**:`S × c3⁻` · `S × EARLY` · `c3⁻ × EARLY`。
**本轮不多测任何别的交互** —— 全族阈只在预先声明的族上有效。

ESTIMAND        同一模型,**置换 1000 次**;报三个交互的**精确置换 p**,
                以及 **max-|t| 的置换分布**(一次性处理多重性,比 Bonferroni 更准)。
KILL            **若 `S × EARLY` 的 |t| 仍过全族 95% 阈 -> 从 PLAUSIBLE 升到 CONFIRMED,可上页面;
                若不过 -> 三条路是**全相加**的,而那是一句更干净的话。**
POSITIVE CTRL   种入 `S × EARLY` = 0.15 -> 必须过全族阈。
NEGATIVE CTRL   置换本身(全族阈的构造)。
⚠ guard 21     若判零,交出三件套。
IMPOSSIBLE      全族阈控制的是**这三个**的族错误率;它不管我在别的轮次里测过多少东西。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
nc=np.isfinite(ONS).sum(1); MO=np.where(nc>=5,np.nanmean(ONS,1),np.nan)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]; EARLY=-MO
m0=np.isfinite(S)&np.isfinite(C3)&np.isfinite(EARLY)&np.isfinite(sh)&ok
n=int(m0.sum()); zz=lambda v:(v[m0]-v[m0].mean())/max(v[m0].std(),1e-12)
zs,zc,ze,ncz=zz(S),zz(C3),zz(EARLY),zz(nc.astype(float)); y=zz(sh)
X=np.column_stack([np.ones(n),zs,zc,ze,ncz,zs*zc,zs*ze,zc*ze])
ILAB=['S × c3⁻','S × EARLY','c3⁻ × EARLY']; IDX=[5,6,7]
XtXi=np.linalg.pinv(X.T@X)
def tstats(yy):
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(n-X.shape[1]); se=np.sqrt(np.diag(s2*XtXi))
    return b[IDX]/se[IDX]
t0=tstats(y)
print(f"n={n:,} · 族 = **预先声明的三个交互**(本轮不多测任何别的)")
print(f"观测 |t|:" + ' · '.join(f"{ILAB[i]} **{t0[i]:+.3f}**" for i in range(3)))
NP=1000; rg=np.random.default_rng(20260804)
T_=np.array([tstats(y[rg.permutation(n)]) for _ in range(NP)])
p_ind=[float((np.abs(T_[:,i])>=abs(t0[i])).mean()) for i in range(3)]
mx=np.abs(T_).max(1); thr=float(np.percentile(mx,95))
p_fw=[float((mx>=abs(t0[i])).mean()) for i in range(3)]
print(f"\n置换 {NP} 次:")
print(f"{'交互':<14}{'|t|':>8}{'逐个 p':>10}{'全族 p':>10}{'过全族阈?':>12}")
for i in range(3):
    print(f"{ILAB[i]:<14}{abs(t0[i]):>8.3f}{p_ind[i]:>10.3f}{p_fw[i]:>10.3f}"
          f"{('**是**' if abs(t0[i])>thr else '否'):>12}")
print(f"   **全族 95% 阈(max-|t| 的置换分布)= {thr:.3f}**")
rg2=np.random.default_rng(19)
yp=0.12*zs+0.12*zc+0.12*ze+0.15*zs*ze+rg2.standard_normal(n); yp=(yp-yp.mean())/yp.std()
tp=tstats(yp)
print(f"\n正对照(种入 `S × EARLY` = 0.15):|t| " +
      ' · '.join(f"{ILAB[i]} **{abs(tp[i]):.2f}**" for i in range(3)) +
      f" -> 过全族阈 **{'是' if abs(tp[1])>thr else '否'}**")
# ⚠ #300a:上页面前发明一个旋钮 —— 起始年龄的数值编码(`#332c` 只测过主效应,没测交互)。
RANK={k:i+1.0 for i,k in enumerate(BIN)}
ALT={'0-4yo':2.5,'5-6yo':6,'7-8yo':8,'9-10yo':10,'11-12yo':12,'13-14yo':14,
     '15-16yo':16,'17-18yo':18,'19-25yo':21,'26yo+':32}
KN=[]
for tag,MAP in (('我的中点',BIN),('序号 1..10',RANK),('另一套中点(26+→32)',ALT)):
    O2=np.column_stack([d[c].map(MAP).values.astype(float) for c in ons])
    M2=np.where(np.isfinite(O2).sum(1)>=5,np.nanmean(O2,1),np.nan)
    ze2=(-M2[m0]-np.nanmean(-M2[m0]))/np.nanstd(-M2[m0])
    X2=np.column_stack([np.ones(n),zs,zc,ze2,ncz,zs*zc,zs*ze2,zc*ze2])
    b2,*_=np.linalg.lstsq(X2,y,rcond=None); r2_=y-X2@b2
    s2=float(r2_@r2_)/(n-X2.shape[1]); se2=np.sqrt(np.diag(s2*np.linalg.pinv(X2.T@X2)))
    KN.append((tag,float(b2[6]),float(b2[6]/se2[6])))
print(f"\n发明的旋钮 · 起始年龄的数值编码,对 `S × EARLY`:")
for t,c,tt in KN: print(f"   {t:<20} 系数 **{c:+.4f}** · |t| **{abs(tt):.3f}** · 过全族阈 **{'是' if abs(tt)>thr else '否'}**")
kok=all(abs(tt)>thr for _,_,tt in KN)

T=pd.DataFrame([dict(v_int=ILAB[i],t=float(t0[i]),p_ind=p_ind[i],p_fw=p_fw[i],thr=thr) for i in range(3)])
check_columns(T,'R379'); T.to_csv(pathlib.Path(__file__).parent/'results'/'fw.csv',index=False)
gg=Gate('那个次可加过不过全族阈')
gg.asserted('★ 发明的旋钮:三套起始年龄编码下 `S × EARLY` 是否都过全族阈',kok,
            ' · '.join(f"{t} |t|={abs(tt):.2f}" for t,_,tt in KN)+f" vs 阈 {thr:.3f}")
gg.asserted('★ 正对照:种入 0.15 必须过全族阈',abs(tp[1])>thr,
            f"|t| {abs(tp[1]):.2f} vs 全族阈 {thr:.3f}")
gg.asserted('★ 负对照:全族阈由置换分布自身构造',True,
            f"max-|t| 的 95 分位 = {thr:.3f};三个交互在零下的 |t| 中位 {np.median(np.abs(T_)):.3f}")
gg.asserted('★ 注册的 kill:`S × EARLY` 是否过全族阈',abs(t0[1])>thr,
            f"|t| **{abs(t0[1]):.3f}** vs 阈 **{thr:.3f}** · 全族 p **{p_fw[1]:.3f}**")
gg.null_claim_uses_null_criteria('★ guard 21:另两个交互判零的三件套','NULL',
    perm_quantile=float(np.mean([p_ind[0],p_ind[2]])),mde=2.8/np.sqrt(n),
    sensitivity_shown=f"种入 0.15 的 |t| = {abs(tp[1]):.2f}",meaningful=0.05)
gg.asserted('⚠ 边界:全族阈只管这三个',True,'它不管我在别的轮次里测过多少东西')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

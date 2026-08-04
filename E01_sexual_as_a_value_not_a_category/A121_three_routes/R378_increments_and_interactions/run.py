import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A121 R378 -- 羞耻的第三条路,与前两条是什么关系

`#308a`:`S` 与 `c3⁻` **相加**(交互 −0.0150,置换零 +0.0051 ± 0.0135,植入 0.15 抓到 +0.1502)。
`#332a`:**起始年龄**是第三个相关物(−0.0810,控类别数与 `S` 之后)。
**三条一起是什么关系?**

⚠ **关键控制**:`#130` 说人群共享一张时间表(内容类早、关系类晚),
所以「起始早」与「兴趣冷门」在**人群层**就相关 —— 先报 `corr(起始均值, S)`,再做偏分析。

ESTIMAND        三者一起进模型:各自的**增量 ΔR²**(去掉它之后 R² 掉多少)+ **三个两两交互**;
                交互对 `perm_finite` 置换零。
KILL            **若起始年龄的增量与 `S`/`c3` 相当且交互为零 -> 羞耻有**三条相加**的路;
                若它被 `S` 吃掉 -> 它是「早来的兴趣更冷门」的影子。**
POSITIVE CTRL   合成一个**带已知交互**的结局 -> 必须被抓到(`#308b` 同款)。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ guard 21     若交互判零,交出三件套。
IMPOSSIBLE      三者都是**同一份自报问卷**里的量;共享方法方差不可分离。
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
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
EARLY=-MO                                       # ⚠ 约定:三条都指向「更多羞耻」
m0=np.isfinite(S)&np.isfinite(C3)&np.isfinite(EARLY)&np.isfinite(sh)&ok
n=int(m0.sum()); zz=lambda v:(v[m0]-v[m0].mean())/max(v[m0].std(),1e-12)
zs,zc,ze=zz(S),zz(C3),zz(EARLY); y=zz(sh); ncz=zz(nc.astype(float))
print(f"n={n:,} · 符号约定:三条都指向更多羞耻(`EARLY = −起始均值`)")
print(f"⚠ 关键控制,先报:`corr(起始均值, S)` = **{np.corrcoef(MO[m0],S[m0])[0,1]:+.4f}** · "
      f"`corr(EARLY, c3⁻)` = **{np.corrcoef(ze,zc)[0,1]:+.4f}** · "
      f"`corr(S, c3⁻)` = **{np.corrcoef(zs,zc)[0,1]:+.4f}**")
def r2(cols):
    X=np.column_stack([np.ones(n)]+cols); b,*_=np.linalg.lstsq(X,y,rcond=None)
    r=y-X@b; return 1-float(r@r)/float(((y-y.mean())**2).sum())
MAIN=[zs,zc,ze,ncz]
full=r2(MAIN)
LAB=['S 位置','c3⁻ 广度型','EARLY 起始早','(控制)类别数']
print(f"\n主效应模型 R² = **{100*full:.3f}%**(含类别数控制)")
for i,l in enumerate(LAB):
    inc=full-r2([c for j,c in enumerate(MAIN) if j!=i])
    print(f"   {l:<14} 增量 ΔR² = **{100*inc:.3f}pp**")
INT=[zs*zc,zs*ze,zc*ze]; ILAB=['S × c3⁻','S × EARLY','c3⁻ × EARLY']
fullI=r2(MAIN+INT)
print(f"\n加三个两两交互后 R² = **{100*fullI:.3f}%**(交互合计增量 **{100*(fullI-full):.3f}pp**)")
X=np.column_stack([np.ones(n)]+MAIN+INT); b,*_=np.linalg.lstsq(X,y,rcond=None)
for i,l in enumerate(ILAB): print(f"   {l:<12} 系数 **{b[5+i]:+.4f}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
def coefs(yy):
    bb,*_=np.linalg.lstsq(X,yy,rcond=None); return bb[5:8]
nul=np.array([coefs((lambda v:(v[m0]-v[m0].mean())/v[m0].std())(perm_finite(sh,300+i))) for i in range(30)])
print(f"负对照(打乱人)三个交互:" +
      ' · '.join(f"{ILAB[i]} {nul[:,i].mean():+.4f} ± {nul[:,i].std():.4f}" for i in range(3)))
q=[float(np.mean(np.abs(nul[:,i])>=abs(b[5+i]))) for i in range(3)]
print(f"   |零| ≥ |观测| 的比例:" + ' · '.join(f"{ILAB[i]} **{q[i]:.3f}**" for i in range(3)))
rg=np.random.default_rng(19)
yp=0.12*zs+0.12*zc+0.12*ze+0.15*zs*ze+rg.standard_normal(n)
yp=(yp-yp.mean())/yp.std(); bp=coefs(yp)
print(f"\n正对照(种入 `S × EARLY` = 0.15):读出 " +
      ' · '.join(f"{ILAB[i]} **{bp[i]:+.4f}**" for i in range(3)))
T=pd.DataFrame([dict(v_term=l,v_inc=100*(full-r2([c for j,c in enumerate(MAIN) if j!=i]))) for i,l in enumerate(LAB)]
               +[dict(v_term=ILAB[i],v_inc=float(b[5+i])) for i in range(3)])
check_columns(T,'R378'); T.to_csv(pathlib.Path(__file__).parent/'results'/'three.csv',index=False)
incE=full-r2([c for j,c in enumerate(MAIN) if j!=2]); incS=full-r2([c for j,c in enumerate(MAIN) if j!=0])
gg=Gate('羞耻的三条路是什么关系')
gg.asserted('★ 正对照:种入 `S × EARLY` = 0.15 必须被抓到且另两个落零',
            abs(bp[1])>0.08 and abs(bp[0])<0.06 and abs(bp[2])<0.06,
            ' · '.join(f"{ILAB[i]} {bp[i]:+.4f}" for i in range(3)))
gg.negative_control('★ 负对照:打乱人后的 `S × c3⁻`',float(nul[:,0].mean()),float(b[5]),
    null_spread=float(nul[:,0].std()),null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill ①:起始年龄的增量是否与 `S` 相当(> 50%)',
            incE/max(incS,1e-9)>0.5,f"EARLY {100*incE:.3f}pp vs S {100*incS:.3f}pp = {100*incE/max(incS,1e-9):.0f}%")
gg.asserted('★ 注册的 kill ②:三个两两交互是否都落在置换零里',
            all(abs(b[5+i])<2*nul[:,i].std() for i in range(3)),
            ' · '.join(f"{ILAB[i]} {b[5+i]:+.4f} vs 2×{nul[:,i].std():.4f}" for i in range(3)))
gg.null_claim_uses_null_criteria('★ guard 21:若交互判零,三件套在不在','NULL',
    perm_quantile=float(np.mean(q)),mde=2.8*float(np.mean(nul.std(0))),
    sensitivity_shown=f"种入 0.15 读出 {bp[1]:+.4f}",meaningful=0.05)
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

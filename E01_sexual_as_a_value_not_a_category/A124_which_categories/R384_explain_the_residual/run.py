import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A124 R384 -- 类别层羞耻残差里的那四成,能不能被**已有的量**解释

`#338b`:残差 sd 是抽样噪声的 **4.94 倍** —— 里面有真东西。
**而命名是这个项目死过五次的地方,所以换一个不需要命名的问法。**

⚠ **注册时列的 ② 不可算,如实报**:`#327` 的 `c3⁻` 载荷是**块**仪器上的,
起始类别与块**题目不相交**(`#303` 的 `assert`),**没有一一对应** ——
**不硬凑映射。** 改用**类别层可算的已有量**。

ESTIMAND        把 `#338b` 的残差回归到四个**已有**的类别层量:
                ① 该类别的**平均起始年龄**(`#332`:早 -> 羞耻多)
                ② 该类别起始年龄的 **sd**
                ③ 报了这一类的人的**平均 `S`**
                ④ 报了这一类的人的**平均 `c3⁻`**
                报 R² 与偏相关,**并报置换零下的 R² 分布**。
KILL            **若某一个吃掉大部分 -> 那四成有名字了,而且名字是已有的量;
                若都吃不掉 -> 类别层有一个这个项目还没有的量,那是一个诚实的开放问题。**
POSITIVE CTRL   合成一个**只由 ① 驱动**的残差 -> 回归必须抓到。
NEGATIVE CTRL   随机残差 -> R² 必须落在置换零里(31 点、4 预测量的期望 ≈ 4/31)。
IMPOSSIBLE      **31 个点**,功效极低;所有结论只能读成「有没有一个大到看得见的解释」。
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
HAS=np.isfinite(ONS); NC=HAS.shape[1]
PREV=HAS.mean(0); RAR=-np.log(np.clip(PREV,1e-4,1.))
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
okS=np.isfinite(sh)
mu=np.array([sh[okS&HAS[:,j]].mean() for j in range(NC)])
z=lambda v:(v-np.nanmean(v))/np.nanstd(v)
X0=np.column_stack([np.ones(NC),z(RAR)])
res=mu-X0@np.linalg.lstsq(X0,mu,rcond=None)[0]
F={'①平均起始年龄':np.array([np.nanmean(ONS[HAS[:,j],j]) for j in range(NC)]),
   '②起始年龄 sd':np.array([np.nanstd(ONS[HAS[:,j],j]) for j in range(NC)]),
   '③拥有者平均 S':np.array([np.nanmean(S[okS&HAS[:,j]&ok]) for j in range(NC)]),
   '④拥有者平均 c3⁻':np.array([np.nanmean(C3[okS&HAS[:,j]&ok]) for j in range(NC)])}
print(f"⚠ 注册时列的「`c3⁻` 载荷对应」**不可算** —— 起始类别与块题目不相交,没有一一对应。")
print(f"   改用类别层可算的四个已有量。31 个点,4 个预测量。\n")
X=np.column_stack([np.ones(NC)]+[z(v) for v in F.values()])
def r2(y,XX=X):
    b,*_=np.linalg.lstsq(XX,y,rcond=None); r=y-XX@b
    return 1-float(r@r)/float(((y-y.mean())**2).sum()),b
R2,b=r2(res)
print(f"R² = **{R2:.4f}**")
for i,l in enumerate(F,1):
    cols=[0]+[k for k in range(1,5) if k!=i]
    Xo=X[:,cols]
    rv=res-Xo@np.linalg.lstsq(Xo,res,rcond=None)[0]
    rx=X[:,i]-Xo@np.linalg.lstsq(Xo,X[:,i],rcond=None)[0]
    print(f"   {l:<16} 系数 **{b[i]:+.4f}** · 偏相关 **{float(np.corrcoef(rv,rx)[0,1]):+.4f}**")
rg=np.random.default_rng(818); NP=2000
nul=np.array([r2(rg.permutation(res))[0] for _ in range(NP)])
q=float((nul>=R2).mean())
print(f"\n置换零(打乱残差 {NP} 次):R² **{nul.mean():.4f} ± {nul.std():.4f}** "
      f"(解析期望 4/31 = {4/NC:.3f})· **零里 ≥ 观测 {q:.3f}**")
rgp=np.random.default_rng(9)
syn=z(F['①平均起始年龄'])+0.5*rgp.standard_normal(NC)
print(f"正对照(只由 ① 驱动的合成残差):R² **{r2(syn)[0]:.4f}**")
ngv=[r2(rgp.standard_normal(NC))[0] for _ in range(300)]
print(f"负对照(随机残差):R² **{np.mean(ngv):.4f} ± {np.std(ngv):.4f}**")
T=pd.DataFrame([dict(v_term='R²',v_val=R2),dict(v_term='零均值',v_val=float(nul.mean())),
                dict(v_term='零分位',v_val=q)]
               +[dict(v_term=l,v_val=float(b[i])) for i,l in enumerate(F,1)])
check_columns(T,'R384'); T.to_csv(pathlib.Path(__file__).parent/'results'/'explain.csv',index=False)
gg=Gate('那四成能不能被已有的量解释')
gg.asserted('★ 正对照:只由 ① 驱动的合成残差必须被抓到(R² > 0.5)',r2(syn)[0]>0.5,f"R² {r2(syn)[0]:.4f}")
gg.asserted('★ 负对照:随机残差的 R² ≈ 4/31',abs(np.mean(ngv)-4/NC)<0.06,
            f"{np.mean(ngv):.4f} ± {np.std(ngv):.4f} vs {4/NC:.3f}")
gg.asserted('★ 注册的 kill:四个已有量能不能吃掉那四成(置换分位 < 0.05)',q<0.05,
            f"R² **{R2:.4f}** · 零 {nul.mean():.4f} ± {nul.std():.4f} · 分位 **{q:.3f}**")
gg.null_claim_uses_null_criteria('★ guard 21:这个零可不可发布','NULL',
    perm_quantile=q,mde=float(np.percentile(nul,95)),
    sensitivity_shown=f"只由 ① 驱动的合成残差 R² = {r2(syn)[0]:.3f}",meaningful=0.30)
gg.asserted('⚠ 如实报:注册时列的「`c3⁻` 载荷对应」不可算',True,
            '起始类别与块题目不相交(`#303` 的 assert),没有一一对应 —— **不硬凑映射**')
gg.asserted('⚠ 边界:31 个点,功效极低',True,'只能读成「有没有一个大到看得见的解释」')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

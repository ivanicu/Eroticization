import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A125 R387 -- 那 42% 是不是「能不能被实施」

`#339b` 的开放问题:类别层羞耻的残差(42%)现有的量一个也解释不了。
人眼里有一个显而易见的读法(更羞耻的一端有真实对象与法律后果,更不羞耻的一端不可能发生)——
**但那是一个名字,而这里已经死过五次命名。所以不起名,问一个可测的代理。**

ESTIMAND        用问卷里已有的 **「我已经实践/试验过所有让我兴奋的东西」**(`41kpfir`)
                做类别层的量:逐类别算「报了这一类的人在这道题上的平均分」;
                问它能不能吃掉 `#338b` 的残差(单独,以及在四个已有量之上的增量)。
KILL            **若能 -> 那 42% 有一个已有量的名字(可实施性),而它本来就在数据里;
                若不能 -> 那 42% 连这个最显而易见的代理都不是。**
POSITIVE CTRL   合成一个只由该代理驱动的残差 -> 必须被抓到。
NEGATIVE CTRL   打乱残差 2000 次的 R² 分布。
⚠⚠ 共享方法方差 该题与羞耻**同源**(同一份自报问卷)-> **先报人层 `corr(该题, 羞耻)`**,
                并把它写进读法:**一个正相关可以只是「答题的人一致地在两道题上都说是」。**
IMPOSSIBLE      31 个点,功效极低;而且「已实践」混着**机会**与**意愿**,不能读成「可实施性」本身。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
ACT=next(c for c in d.columns if '41kpfir' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
act=pd.to_numeric(d[ACT],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
HAS=np.isfinite(ONS); NCA=HAS.shape[1]
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
okS=np.isfinite(sh); okA=np.isfinite(act)
kk=okS&okA
print(f"⚠⚠ 共享方法方差,先报:人层 `corr(已实践, 羞耻)` = "
      f"**{np.corrcoef(act[kk],sh[kk])[0,1]:+.4f}**(n={int(kk.sum()):,})")
print(f"   **一个正相关可以只是「答题的人一致地在两道题上都说是」—— 读法里必须带着这句。**\n")
PREV=HAS.mean(0); RAR=-np.log(np.clip(PREV,1e-4,1.))
mu=np.array([sh[okS&HAS[:,j]].mean() for j in range(NCA)])
z=lambda v:(v-np.nanmean(v))/np.nanstd(v)
X0=np.column_stack([np.ones(NCA),z(RAR)])
res=mu-X0@np.linalg.lstsq(X0,mu,rcond=None)[0]
ACTj=np.array([np.nanmean(act[okA&HAS[:,j]]) for j in range(NCA)])
F4=[np.array([np.nanmean(ONS[HAS[:,j],j]) for j in range(NCA)]),
    np.array([np.nanstd(ONS[HAS[:,j],j]) for j in range(NCA)]),
    np.array([np.nanmean(S[okS&HAS[:,j]&ok]) for j in range(NCA)]),
    np.array([np.nanmean(C3[okS&HAS[:,j]&ok]) for j in range(NCA)])]
def r2(y,XX):
    b,*_=np.linalg.lstsq(XX,y,rcond=None); r=y-XX@b
    return 1-float(r@r)/float(((y-y.mean())**2).sum())
Xa=np.column_stack([np.ones(NCA),z(ACTj)])
X4=np.column_stack([np.ones(NCA)]+[z(v) for v in F4])
X5=np.column_stack([np.ones(NCA)]+[z(v) for v in F4]+[z(ACTj)])
ra,r4,r5=r2(res,Xa),r2(res,X4),r2(res,X5)
print(f"① 单独:`corr(残差, 类别层已实践)` = **{np.corrcoef(res,ACTj)[0,1]:+.4f}** · R² **{ra:.4f}**")
print(f"② 四个已有量 R² **{r4:.4f}** -> 加上它 **{r5:.4f}**(**增量 {r5-r4:+.4f}**)")
rg=np.random.default_rng(313); NP=2000
nula=np.array([r2(rg.permutation(res),Xa) for _ in range(NP)])
nul5=np.array([r2(rg.permutation(res),X5) for _ in range(NP)])
qa=float((nula>=ra).mean()); q5=float((nul5>=r5).mean())
print(f"   置换零:单独 **{nula.mean():.4f} ± {nula.std():.4f}**(分位 **{qa:.3f}**)· "
      f"五量 **{nul5.mean():.4f} ± {nul5.std():.4f}**(分位 **{q5:.3f}**)")
rgp=np.random.default_rng(11)
syn=z(ACTj)+0.5*rgp.standard_normal(NCA)
print(f"\n正对照(只由该代理驱动的合成残差):单独 R² **{r2(syn,Xa):.4f}**")
ngv=[r2(rgp.standard_normal(NCA),Xa) for _ in range(300)]
print(f"负对照(随机残差):单独 R² **{np.mean(ngv):.4f} ± {np.std(ngv):.4f}**(解析 1/31 = {1/NCA:.3f})")
T=pd.DataFrame([dict(v_arm='单独',v_r2=ra,v_q=qa),dict(v_arm='四量',v_r2=r4,v_q=np.nan),
                dict(v_arm='五量',v_r2=r5,v_q=q5)])
check_columns(T,'R387'); T.to_csv(pathlib.Path(__file__).parent/'results'/'act.csv',index=False)
gg=Gate('那 42% 是不是「能不能被实施」')
gg.asserted('★ 正对照:只由该代理驱动的合成残差必须被抓到(R² > 0.5)',r2(syn,Xa)>0.5,
            f"R² {r2(syn,Xa):.4f}")
gg.asserted('★ 负对照:随机残差的单独 R² ≈ 1/31',abs(np.mean(ngv)-1/NCA)<0.05,
            f"{np.mean(ngv):.4f} ± {np.std(ngv):.4f} vs {1/NCA:.3f}")
gg.asserted('★ 注册的 kill:这个代理能不能吃掉残差(单独或增量,置换分位 < 0.05)',
            qa<0.05 or q5<0.05,
            f"单独 R² {ra:.4f}(分位 {qa:.3f})· 五量 {r5:.4f}(分位 {q5:.3f})· 增量 {r5-r4:+.4f}")
gg.null_claim_uses_null_criteria('★ guard 21:若判零,三件套在不在',
    'NULL' if (qa>=0.05 and q5>=0.05) else 'EFFECT',
    perm_quantile=qa,mde=float(np.percentile(nula,95)),
    sensitivity_shown=f"只由代理驱动的合成残差 R² = {r2(syn,Xa):.3f}",meaningful=0.30)
gg.asserted('⚠⚠ 共享方法方差已先报,并写进读法',True,
            f"人层 `corr(已实践, 羞耻)` = {np.corrcoef(act[kk],sh[kk])[0,1]:+.4f} —— "
            f"一个正相关可以只是答题一致性")
gg.asserted('⚠ 边界:「已实践」混着机会与意愿',True,'不能读成「可实施性」本身;31 个点功效极低')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

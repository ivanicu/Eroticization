import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A124 R383 -- 哪些类别的拥有者平均更羞耻

`#337a` 建起了一张 147,387 行的**人×类别**长表,而这个项目此前**全部**在人层做。

⚠⚠ **跑之前写下,免得结果一出来就读错**:
这仍然是**人层羞耻的重排**,**不是**「他为**这一类**羞耻」——
后者需要逐类别的羞耻题,而问卷**没有**。**本轮读到的一切都是「报了这一类的人整体更羞耻」。**

ESTIMAND        逐类别「报了这一类的人的平均羞耻」(31 个数);
                报它与 `rar_j` 的相关(**已声明会高**),再看**残差**两端是哪些类别;
                **加权(按 n_j)与未加权两版都报。**
KILL            **若残差两端有可读的图样 -> 类别层有 `rar` 之外的东西;
                若残差是噪声 -> 类别层的羞耻剖面就是稀有度,没有别的。**
POSITIVE CTRL   合成一个**只由某一个类别驱动**的羞耻 -> 残差必须把那一类挑出来。
NEGATIVE CTRL   打乱**人**(保持类别结构)-> 残差图样必须消失。
⚠ 每类别的均值有自己的 se(n_j 差异很大)-> **报每类别的 se,并做加权版**。
IMPOSSIBLE      31 个点的回归,功效极低;残差的排序**不是**一个检验,是一张描述性的表。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')

BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
HAS=np.column_stack([np.isfinite(d[c].map(BIN).values.astype(float)) for c in ons])
NC=HAS.shape[1]; PREV=HAS.mean(0); RAR=-np.log(np.clip(PREV,1e-4,1.))
NAM=[str(c) for c in ons]
okS=np.isfinite(sh)
mu=np.array([sh[okS&HAS[:,j]].mean() for j in range(NC)])
se=np.array([sh[okS&HAS[:,j]].std()/np.sqrt(max((okS&HAS[:,j]).sum(),1)) for j in range(NC)])
nj=np.array([(okS&HAS[:,j]).sum() for j in range(NC)])
print(f"类别 {NC} · 每类别 n 从 **{nj.min():,}** 到 **{nj.max():,}** · "
      f"平均羞耻从 **{mu.min():+.3f}** 到 **{mu.max():+.3f}**")
z=lambda v:(v-v.mean())/v.std()
def fit(w=None):
    X=np.column_stack([np.ones(NC),z(RAR)])
    W=np.ones(NC) if w is None else w/w.mean()
    b=np.linalg.lstsq(X*W[:,None],mu*W,rcond=None)[0]
    r=mu-X@b
    ss=1-float((W*r**2).sum())/float((W*(mu-np.average(mu,weights=W))**2).sum())
    return b,r,ss
bU,rU,r2U=fit(); bW,rW,r2W=fit(w=1/np.maximum(se,1e-9)**2)
print(f"\n① 与 `rar_j` 的相关(**跑前已声明会高**):"
      f"未加权 **{np.corrcoef(mu,RAR)[0,1]:+.4f}**(R² {r2U:.3f})· "
      f"加权 R² **{r2W:.3f}**")
o=np.argsort(-rU)
print(f"\n② 残差两端(未加权;⚠ 这是描述性的表,不是检验):")
for i in o[:5]: print(f"   {rU[i]:+.3f}  (n={nj[i]:>5,}, se {se[i]:.3f})  {NAM[i][:60]}")
print("   ---")
for i in o[::-1][:5]: print(f"   {rU[i]:+.3f}  (n={nj[i]:>5,}, se {se[i]:.3f})  {NAM[i][:60]}")
print(f"\n   残差 sd **{rU.std():.4f}** vs 各类别 se 中位 **{np.median(se):.4f}** -> "
      f"**{rU.std()/np.median(se):.2f}×**(> 1 才说明残差里有比抽样噪声更多的东西)")
print(f"   加权残差与未加权残差的相关 **{np.corrcoef(rU,rW)[0,1]:+.4f}**")
rg=np.random.default_rng(505); NP=300
def resid_of(shv):
    m2=np.array([shv[okS&HAS[:,j]].mean() for j in range(NC)])
    X=np.column_stack([np.ones(NC),z(RAR)])
    return m2-X@np.linalg.lstsq(X,m2,rcond=None)[0]
NUL=np.array([resid_of(sh[okS][rg.permutation(int(okS.sum()))].astype(float)
              if False else (lambda v: v)(np.where(okS,sh[rg.permutation(NN)],np.nan)))
              for _ in range(NP)])
print(f"负对照(打乱人,保持类别结构):残差 sd **{np.nanmean(NUL.std(1)):.4f} ± {np.nanstd(NUL.std(1)):.4f}**")
J=int(np.argmax(nj*0+ (nj>2000)*RAR))
syn=sh.copy().astype(float); syn[HAS[:,J]]+=1.0
rs=resid_of(syn); k=int(np.argmax(rs))
print(f"正对照(只给类别 #{J} 的人 +1.0):残差最大的是 **#{k}**(应 #{J}),"
      f"残差 **{rs[k]:+.3f}** vs 其余中位 {np.median(np.delete(rs,k)):+.3f}")
T=pd.DataFrame([dict(v_cat=NAM[j][:44],v_mu=float(mu[j]),v_se=float(se[j]),v_n=int(nj[j]),
                     v_rar=float(RAR[j]),v_res=float(rU[j])) for j in range(NC)])
check_columns(T,'R383'); T.to_csv(pathlib.Path(__file__).parent/'results'/'cats.csv',index=False)
gg=Gate('哪些类别的拥有者平均更羞耻')
gg.asserted('★ 正对照:只给一个类别加 +1.0,残差必须把它挑出来',k==J,f"挑出 #{k},应 #{J}")
gg.asserted('★ 负对照:打乱人后残差 sd 必须掉到抽样水平',
            np.nanmean(NUL.std(1))<rU.std(),
            f"零 {np.nanmean(NUL.std(1)):.4f} vs 观测 {rU.std():.4f}")
gg.asserted('★ 注册的 kill:残差里有没有比抽样噪声更多的东西',rU.std()>np.median(se),
            f"残差 sd {rU.std():.4f} vs 各类别 se 中位 {np.median(se):.4f} "
            f"({rU.std()/np.median(se):.2f}×)")
gg.asserted('⚠⚠ 跑前声明:这是人层羞耻的重排,不是「他为这一类羞耻」',True,
            '后者需要逐类别的羞耻题,而问卷没有')
gg.asserted('⚠ 边界:31 个点,残差排序不是检验',True,'它是一张描述性的表')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

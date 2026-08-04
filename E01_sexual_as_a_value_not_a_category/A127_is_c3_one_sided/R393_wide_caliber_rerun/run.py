import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A127 R393 -- 用宽口径重跑 `#390`,窄口径在同一次运行里作参照臂

`#390`:`c3⁻` 上半斜率约是下半的两倍,差 **+0.0741** vs 2×se **0.0792** —— 刚好没过;
guard 21 说这个零也不可发布(MDE 0.111 ≥ 有意义量 0.10)。**两个方向都读不了。**
`#346a`:宽口径(`cov>=4`,n = 12,720)对 `c3` 几乎免费(信度只掉 10%)。

⚠ **口径已登记**(`CALIBER.md` ⑩)**且页面已注明**(`#348b` 验证过)—— **现在可以用它做新推断。**
⚠ **门槛由同一次运行的参照臂给**(`#329b`:从别轮拿数当门槛已经错过四次)。

ESTIMAND        `cov>=4` 与 `cov>=8` **同一次运行**:以中位为节点的分段斜率差 · 各自的 se 与 MDE。
KILL            **若差保持而 MDE 降到它之下 -> 「`c3⁻` 上半更陡」成立,是一个新的形状;
                若差随 n 缩小 -> 它本来就是噪声,这条线关掉。**
POSITIVE CTRL   两个口径各自种入**已知单向**的关系 -> 都必须被抓到。
NEGATIVE CTRL   打乱人。
⚠ 衰减         `c3` 在宽口径下信度掉 10% -> 斜率会被**额外衰减**,所以**参照臂是必须的**。
IMPOSSIBLE      两个口径的**人群不同**(宽口径多了 6,003 个低覆盖的人,而 `#346b` 说他们系统上不一样)
                —— 所以这不是「同一个问题、更多数据」,而是「一个更大也更杂的人群」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
def build(TH):
    R_=np.flatnonzero((cov>=TH)&np.isfinite(sh))
    v=load_of(R_)
    s=score_of(v)
    s=np.where(cov>=TH,s,np.nan)
    k=np.isfinite(s)&np.isfinite(sh)
    if float(np.corrcoef(s[k],sh[k])[0,1])<0: s=-s
    m=np.isfinite(s)&np.isfinite(sh)&(cov>=TH)
    return s,m
def seg(x,yv,m,q=0.5):
    n=int(m.sum())
    z=lambda v:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
    xx,yy=z(x),(yv[m]-yv[m].mean())/max(yv[m].std(),1e-12)
    kn=np.quantile(xx,q); hi=(xx>kn).astype(float)
    X=np.column_stack([np.ones(n),xx,(xx-kn)*hi,hi])
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(n-X.shape[1]); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[1]),float(b[1]+b[2]),float(b[2]),float(se[2]),n
res={}
for TH in (8,4):
    s,m=build(TH)
    lo,hi,dif,sed,n=seg(s,sh,m)
    res[TH]=dict(lo=lo,hi=hi,dif=dif,se=sed,n=n,s=s,m=m)
    print(f"cov>={TH}(n={n:,}):下半 **{lo:+.4f}** · 上半 **{hi:+.4f}** · "
          f"**差 {dif:+.4f}**(se {sed:.4f},|t| {abs(dif/sed):.2f},**MDE {2.8*sed:.4f}**)")
rg=np.random.default_rng(31)
print(f"\n正对照(两个口径各自种入**已知单向**):")
for TH in (8,4):
    s,m=res[TH]['s'],res[TH]['m']; n=res[TH]['n']
    z=(s[m]-s[m].mean())/s[m].std()
    y=np.full(NN,np.nan); y[m]=np.where(z>0,0.25*z,0)+rg.standard_normal(n)
    _,_,dd,ss,_=seg(s,y,m)
    print(f"   cov>={TH}: 差 **{dd:+.4f}**(|t| {abs(dd/ss):.2f})")
    res[TH]['pc']=(dd,ss)
def perm(v,seed):
    z2=v.copy(); z2=z2[np.random.default_rng(seed).permutation(len(z2))]; return z2
for TH in (8,4):
    s,m=res[TH]['s'],res[TH]['m']
    nul=np.array([seg(s,np.where(m,perm(sh,900+i),np.nan),m)[2] for i in range(120)])
    res[TH]['nul']=(float(nul.mean()),float(nul.std()),
                    float((np.abs(nul)>=abs(res[TH]['dif'])).mean()))
    print(f"负对照 cov>={TH}:差 **{nul.mean():+.4f} ± {nul.std():.4f}** · "
          f"|零| ≥ |观测| **{res[TH]['nul'][2]:.3f}**")
r8,r4=res[8],res[4]
print(f"\n★ 参照臂对比:窄 差 **{r8['dif']:+.4f}**(MDE {2.8*r8['se']:.4f})-> "
      f"宽 差 **{r4['dif']:+.4f}**(MDE {2.8*r4['se']:.4f})")
T=pd.DataFrame([dict(v_th=8,dif=r8['dif'],se=r8['se'],n=r8['n']),
                dict(v_th=4,dif=r4['dif'],se=r4['se'],n=r4['n'])])
check_columns(T,'R393'); T.to_csv(pathlib.Path(__file__).parent/'results'/'wide.csv',index=False)
gg=Gate('宽口径重跑:`c3⁻` 上半更不更陡')
gg.asserted('★ 正对照:两个口径都必须抓到已知单向',
            abs(r8['pc'][0]/r8['pc'][1])>3 and abs(r4['pc'][0]/r4['pc'][1])>3,
            f"窄 |t| {abs(r8['pc'][0]/r8['pc'][1]):.2f} · 宽 |t| {abs(r4['pc'][0]/r4['pc'][1]):.2f}")
gg.negative_control('★ 负对照(宽口径):打乱人',r4['nul'][0],r4['dif'],
    null_spread=r4['nul'][1],null_kind='打乱人 —— 保持 `c3⁻` 分布,只打散羞耻')
gg.asserted('★ 注册的 kill:宽口径下差是否越过它自己的 2×se',abs(r4['dif'])>2*r4['se'],
            f"宽 差 **{r4['dif']:+.4f}** vs 2×se **{2*r4['se']:.4f}**(|零| ≥ |观测| {r4['nul'][2]:.3f});"
            f"参照臂 窄 {r8['dif']:+.4f} vs 2×se {2*r8['se']:.4f}")
gg.null_claim_uses_null_criteria('★ guard 21:若判对称(零),三件套在不在',
    'NULL' if abs(r4['dif'])<=2*r4['se'] else 'EFFECT',
    perm_quantile=r4['nul'][2],mde=2.8*r4['se'],
    sensitivity_shown=f"宽口径单向合成 |t| {abs(r4['pc'][0]/r4['pc'][1]):.1f}",meaningful=0.10)
gg.asserted('⚠ 边界:两个口径的**人群不同**',True,
            '宽口径多了 6,003 个低覆盖的人,而 `#346b` 说他们系统上不一样 —— '
            '这不是「同一个问题、更多数据」,是「一个更大也更杂的人群」')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

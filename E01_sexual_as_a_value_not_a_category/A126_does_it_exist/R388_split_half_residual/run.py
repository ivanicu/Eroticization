import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A126 R388 -- 那 42% 到底是不是**类别**的性质

`#339`(四个已有量)· `#341`(概念重叠)· `#342`(可实施性)——**三轮各排除一个候选**,
而这是一条会无限延伸的线。**换一个能一次性定界的问法。**

ESTIMAND        **分半复现**:把**人**随机劈两半,各自算 31 个类别的残差
                (类别平均羞耻 − 稀有度拟合),判两条 **31 维残差向量**的相关;重复 **20** 次;
                再按 Spearman–Brown 外推到全样本。
KILL            **若两半的残差高度相关 -> 它是**类别**的稳定性质,值得继续找,而且给出了信度上界;
                若不相关 -> 那个「4.94×」是抽样噪声在类别层的堆积,`#339b` 要重新解释,
                这条线可以关掉。**
POSITIVE CTRL   合成一个**类别层真信号** -> 分半相关必须高。
NEGATIVE CTRL   合成**纯人层噪声** -> 分半相关必须 ≈ 0。
⚠ 每半 n 减半 -> 各类别 se 变大;**同时报每半的 se 中位**。
IMPOSSIBLE      分半相关是**信度**,不是**有效性**:一个稳定的类别层量也可以是别的东西的影子。
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
NCA=HAS.shape[1]; okS=np.isfinite(sh); POOL=np.flatnonzero(okS)
z=lambda v:(v-np.nanmean(v))/np.nanstd(v)
def resid(rows,y=None):
    yy=sh if y is None else y
    m=np.zeros(NN,bool); m[rows]=True
    mu=np.array([yy[m&HAS[:,j]].mean() if (m&HAS[:,j]).sum()>=30 else np.nan for j in range(NCA)])
    se=np.array([yy[m&HAS[:,j]].std()/np.sqrt(max((m&HAS[:,j]).sum(),1)) for j in range(NCA)])
    P=HAS[rows].mean(0); RAR=-np.log(np.clip(P,1e-4,1.))
    g=np.isfinite(mu)
    X=np.column_stack([np.ones(g.sum()),z(RAR[g])])
    r=np.full(NCA,np.nan); r[g]=mu[g]-X@np.linalg.lstsq(X,mu[g],rcond=None)[0]
    return r,se
def sh_split(y=None,seed=0,T=20):
    rg=np.random.default_rng(seed); out=[];ses=[]
    for _ in range(T):
        p=rg.permutation(POOL); h=len(p)//2
        a,sa=resid(p[:h],y); b,sb=resid(p[h:2*h],y)
        k=np.isfinite(a)&np.isfinite(b)
        if k.sum()>=10: out.append(float(np.corrcoef(a[k],b[k])[0,1])); ses.append(float(np.median(np.r_[sa,sb])))
    return np.array(out),np.array(ses)
r,ses=sh_split(seed=4242)
sb=lambda x:2*x/(1+x)
print(f"n={len(POOL):,} · 20 次随机人劈半,每半 n≈{len(POOL)//2:,}")
print(f"★ **分半相关(31 维残差向量)= {r.mean():+.4f} ± {r.std():.4f}** · "
      f"范围 [{r.min():+.4f}, {r.max():+.4f}]")
print(f"   Spearman–Brown 外推到全样本:**{sb(r.mean()):.4f}**")
print(f"   每半各类别 se 中位 **{ses.mean():.4f}**(全样本是 0.0294)")
rg=np.random.default_rng(7)
TRUE=rg.standard_normal(NCA)
per=np.zeros(NN)
cnt=HAS.sum(1)
per=np.where(cnt>0,(HAS*TRUE[None,:]).sum(1)/np.maximum(cnt,1),0.0)
ypos=np.where(okS,0.6*z(per)+rg.standard_normal(NN),np.nan)
rp,_=sh_split(y=ypos,seed=99)
yneg=np.where(okS,rg.standard_normal(NN),np.nan)
rn,_=sh_split(y=yneg,seed=98)
print(f"\n正对照(合成**类别层真信号**):分半相关 **{rp.mean():+.4f} ± {rp.std():.4f}**")
print(f"负对照(纯人层噪声):分半相关 **{rn.mean():+.4f} ± {rn.std():.4f}**")
T=pd.DataFrame([dict(v_arm='观测',v_r=float(r.mean()),v_sd=float(r.std())),
                dict(v_arm='正对照',v_r=float(rp.mean()),v_sd=float(rp.std())),
                dict(v_arm='负对照',v_r=float(rn.mean()),v_sd=float(rn.std()))])
check_columns(T,'R388'); T.to_csv(pathlib.Path(__file__).parent/'results'/'sh.csv',index=False)
gg=Gate('那 42% 是不是类别的性质')
gg.asserted('★ 正对照:合成类别层真信号 -> 分半相关必须高(> 0.6)',rp.mean()>0.6,
            f"{rp.mean():+.4f} ± {rp.std():.4f}")
gg.asserted('★ 负对照:纯人层噪声 -> 分半相关必须 ≈ 0',abs(rn.mean())<0.15,
            f"{rn.mean():+.4f} ± {rn.std():.4f}")
gg.asserted('★★ 注册的 kill(这条线的存亡判据):两半残差是否高度相关(> 0.5)',
            r.mean()>0.5,
            f"**{r.mean():+.4f} ± {r.std():.4f}**(SB 外推 {sb(r.mean()):.4f})· "
            f"正对照 {rp.mean():+.4f} · 负对照 {rn.mean():+.4f}")
gg.asserted('⚠ 边界:分半相关是**信度**,不是**有效性**',True,
            '一个稳定的类别层量也可以是别的东西的影子')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

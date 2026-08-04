import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A104 R355 -- 一个性开放的成长环境,最多能挡掉多少羞耻

`#309c`:压抑 +0.1091 / 开放 +0.0789,差 0.030,各格 se 0.017 —— **读不出来**。
**「挡不掉」和「挡掉的量小于我能测的」是两句不同的话,而页面只该出现后一句。**
所以本轮报的是**界**和 **MDE**,不是零。

ESTIMAND        成长期性开放度是**三档**(压抑/中性/开放)。剂量-反应:
                羞耻 ~ zRS + zRC + **zRS·U + zRC·U**(U 为 −1/0/+1),
                报两个交互系数的 **95% CI**;并报**这个 n 下的 MDE**。
KILL            **若 CI 排除零 -> 成长环境确实缓冲,报缓冲比例;
                若 CI 含零 -> 报**上界**:「最多挡掉 X%」,并同时报 MDE ——
                只有 MDE 明显小于一个有意义的缓冲量,这个上界才有内容。**
POSITIVE CTRL   合成一个**已知按成长环境缓冲 30%** 的羞耻结局 -> 剂量-反应必须抓到。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ 混淆         成长环境与**当前年龄**、**关系风格**相关 -> **同一轮里放进模型**,不事后加。
IMPOSSIBLE      成长环境是**回溯自报**,而报告它的人此刻的羞耻会污染这个回忆;
                本轮测的是**关联**,不是「成长环境**造成**了什么」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
UP=d['How "sexually liberated" was your upbringing? (fs700v2)'].map(
    {'Repressed':-1.,'Neutral':0.,'Liberated':1.}).values.astype(float)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
AG=d['age'].map(AGE).values.astype(float)
RL=d['Personally, your preferred relationship style is: (4jib23m)'].map(
    {'Monogamous':0.,'Not monogamous':1.}).values.astype(float)
m0=np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(UP)&np.isfinite(AG)&np.isfinite(RL)&ok
def resid(a,b,m):
    out=np.full(NN,np.nan); x=b[m]; x=(x-x.mean())/x.std()
    out[m]=a[m]-np.polyval(np.polyfit(x,a[m],1),x); return out
RS=resid(S,C3,m0); RC=resid(C3,S,m0)
n=int(m0.sum())
print(f"n={n:,} · 成长环境三档 " + ' · '.join(f"{k}: {int((UP[m0]==v).sum()):,}"
      for k,v in (('压抑',-1),('中性',0),('开放',1))))
print(f"⚠ 混淆(同一轮放进模型):corr(U, 年龄) = **{np.corrcoef(UP[m0],AG[m0])[0,1]:+.4f}** · "
      f"corr(U, 关系风格) = **{np.corrcoef(UP[m0],RL[m0])[0,1]:+.4f}**")
def design(y,m):
    z=lambda v:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
    u=z(UP)
    X=np.column_stack([np.ones(m.sum()),z(RS),z(RC),u,z(AG),z(RL),
                       z(RS)*u,z(RC)*u,                       # ★ 两个剂量-反应交互
                       z(RS)*z(AG),z(RC)*z(AG),z(RS)*z(RL),z(RC)*z(RL)])  # ⚠ 混淆的交互也放进去
    return X,(y[m]-y[m].mean())/y[m].std()
def fit(y,m=None):
    m=m0 if m is None else m
    X,yy=design(y,m); b,*_=np.linalg.lstsq(X,yy,rcond=None)
    r=yy-X@b; s2=float(r@r)/(len(yy)-X.shape[1])
    cov=s2*np.linalg.pinv(X.T@X); se=np.sqrt(np.diag(cov))
    return b,se
b,se=fit(sh)
LAB={6:'zRS × U',7:'zRC × U'}
print(f"\n剂量-反应交互(U = −1 压抑 / 0 / +1 开放,已同时控制年龄与关系风格的交互):")
for i,l in LAB.items():
    lo,hi=b[i]-1.96*se[i],b[i]+1.96*se[i]
    print(f"   {l:<10} **{b[i]:+.4f}** · 95% CI **[{lo:+.4f}, {hi:+.4f}]** · se {se[i]:.4f}")
MAIN={6:b[1],7:b[2]}
print(f"\n主效应(U=0 时):zRS **{b[1]:+.4f}** · zRC **{b[2]:+.4f}**")
print(f"⇒ 从压抑(U=−1)到开放(U=+1),系数变化 = 2×交互:")
for i,l in LAB.items():
    ch=2*b[i]; lo,hi=2*(b[i]-1.96*se[i]),2*(b[i]+1.96*se[i])
    base=MAIN[i]
    print(f"   {l.split(' ')[0]:<6} 变化 **{ch:+.4f}**(相对主效应 **{100*ch/base:+.1f}%**)· "
          f"95% CI 的缓冲比例 **[{100*hi/base:+.1f}%, {100*lo/base:+.1f}%]**")
MDE={i:2.8*se[i] for i in LAB}
print(f"\n**MDE**(80% 功效、α=.05 双侧 ≈ 2.8×se):")
for i,l in LAB.items():
    print(f"   {l:<10} 交互的 MDE **{MDE[i]:.4f}** -> 换算成缓冲比例 **{200*MDE[i]/MAIN[i]:.1f}%**")
rg=np.random.default_rng(31)
def perm_finite(v,seed):
    z=v.copy(); j=np.flatnonzero(np.isfinite(z))
    z[j]=z[np.random.default_rng(seed).permutation(j)]; return z
nul=np.array([[fit(perm_finite(sh,400+i))[0][k] for k in LAB] for i in range(15)])
print(f"负对照(打乱人):zRS×U **{nul[:,0].mean():+.4f} ± {nul[:,0].std():.4f}** · "
      f"zRC×U **{nul[:,1].mean():+.4f} ± {nul[:,1].std():.4f}**")
z=lambda v,m:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
print(f"\n正对照:合成一个**已知按成长环境缓冲 30%** 的羞耻结局")
CT={}
for buf in (0.0,0.30):
    y=np.full(NN,np.nan)
    zr,zc,u=z(RS,m0),z(RC,m0),z(UP,m0)
    y[m0]=0.10*zr*(1-buf*u)+0.10*zc*(1-buf*u)+rg.standard_normal(n)
    bb,ss=fit(y); CT[buf]=(bb[6],bb[1])
    print(f"   缓冲 {100*buf:.0f}%: zRS×U **{bb[6]:+.4f}** · 主效应 {bb[1]:+.4f} -> "
          f"读出的缓冲 **{-200*bb[6]/max(bb[1],1e-9):.1f}%**")
T=pd.DataFrame([dict(v_term=l,coef=b[i],se=se[i],lo=b[i]-1.96*se[i],hi=b[i]+1.96*se[i],
                     main=MAIN[i],mde=MDE[i]) for i,l in LAB.items()])
check_columns(T,'R355'); T.to_csv(pathlib.Path(__file__).parent/'results'/'dose_mde.csv',index=False)
gg=Gate('性开放的成长环境最多能挡掉多少羞耻')
gg.asserted('★ 正对照:已知缓冲 30% 必须被抓到(读出 > 15%)',
            -200*CT[0.30][0]/max(CT[0.30][1],1e-9)>15,
            f"缓冲 30% -> 读出 {-200*CT[0.30][0]/max(CT[0.30][1],1e-9):.1f}%;"
            f"缓冲 0% -> 读出 {-200*CT[0.0][0]/max(CT[0.0][1],1e-9):.1f}%")
gg.negative_control('★ 负对照:打乱人后的 zRS×U',float(nul[:,0].mean()),b[6],
    null_spread=float(nul[:,0].std()),null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill:CI 排不排除零',
            (b[6]-1.96*se[6])*(b[6]+1.96*se[6])>0 or (b[7]-1.96*se[7])*(b[7]+1.96*se[7])>0,
            f"zRS×U CI [{b[6]-1.96*se[6]:+.4f}, {b[6]+1.96*se[6]:+.4f}] · "
            f"zRC×U CI [{b[7]-1.96*se[7]:+.4f}, {b[7]+1.96*se[7]:+.4f}]")
gg.asserted('★ 上界要有内容:MDE 必须明显小于一个有意义的缓冲量(取 30%)',
            max(200*MDE[i]/MAIN[i] for i in LAB)<30,
            ' · '.join(f"{l} MDE={200*MDE[i]/MAIN[i]:.1f}%" for i,l in LAB.items())+" vs 30%")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

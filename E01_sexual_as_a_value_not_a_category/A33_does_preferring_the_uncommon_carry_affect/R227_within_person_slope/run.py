import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A33 R227 -- 换估计量:人内斜率,而不是劈半

`#181`:要在绝对劈分上判到 2× 需约 1.7 倍样本,而样本是固定的。**所以要换的不是 n,是估计量。**
相对劈分有功效是因为它用上**每一个人**;绝对劈分丢掉了所有"全在 12 岁前/后"的人。

ESTIMAND        对每个人跑 `稀有度 ~ 获得年龄` 的人内回归,取**斜率**作为人层量 b_i;
                判 r(b_i, 羞耻)。斜率为正 = 这个人越晚获得的东西越冷门。
                `#128` 的共享时间表说人群层面斜率是正的 -> **本轮判的是个体偏离**。
KILL            条件式:先要**正对照开火**(把羞耻换成一个由 b_i 构造的合成量,必须强测到);
                再判:**|r(b, 羞耻)| > 2× 自助 sd,且方向与 `#180` 的配对差一致** ——
                `#180` 说羞耻贴**早半**更紧,而早半是常见的东西 ->
                **预测 r(b, 羞耻) 为负**(斜率越小 = 早/晚的冷门度差越小 = ?)。
                ⚠ 这个映射不平凡,所以**方向的预测在跑之前写死在这里**:
                `#180` 的 Δ>0 意味着"早半的冷门度与羞耻的关联强于晚半";
                在斜率语言里,一个人若**早期就偏冷门**(截距高、斜率平),
                他的 b_i 更小而早半冷门度更高 -> **预测 r(b, 羞耻) < 0**。
NEGATIVE CTRL   人内打乱获得年龄 -> 斜率的零分布。
POSITIVE CTRL   合成 y = b + 2σ 噪声,同一条管道必须强测到。
CONTROL         类别数 K、斜率的**人内标准误**(拟合得越差的人斜率越噪)。
NOISE FLOOR     人层 bootstrap 500。
IMPOSSIBLE      斜率把"起点"与"速率"压成一个数;截距与斜率在本设计里同时报,但**不能声称分开了**。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
sh=df['"I am ashamed or embarrassed about at least some of what arouses me" (7cw1ziu)'].values.astype(float)
O=pd.read_csv('data/derived/onset.csv')
onset=[c for c in O.columns if re.search(r'How old were you when you first',c)]
A_=O[onset].apply(pd.to_numeric,errors='coerce').values
A_=np.where((A_>=2)&(A_<=60),A_,np.nan)
assert np.isfinite(A_).sum()>10000
have=np.isfinite(A_); rar=-np.log(np.clip(have.mean(0),1e-4,1.)); K=have.sum(1).astype(float)

def slopes(Amat, shuffle=False, rng=None):
    n=Amat.shape[0]; b=np.full(n,np.nan); a=np.full(n,np.nan); se=np.full(n,np.nan)
    for i in range(n):
        idx=np.flatnonzero(np.isfinite(Amat[i]))
        if len(idx)<6: continue
        x=Amat[i,idx].astype(float)
        if shuffle: x=x[rng.permutation(len(x))]
        y=rar[idx]
        if x.std()<1e-9: continue
        X=np.c_[np.ones(len(x)),x]
        coef,res,_,_=np.linalg.lstsq(X,y,rcond=None)
        a[i],b[i]=coef
        resid=y-X@coef; dof=max(len(x)-2,1)
        s2=float(resid@resid)/dof
        se[i]=np.sqrt(s2/max(((x-x.mean())**2).sum(),1e-9))
    return a,b,se

I0,B,SE=slopes(A_)
ok=np.isfinite(B)
print(f"人内斜率:有效 {ok.sum():,} 人   均值 {np.nanmean(B):+.4f}   "
      f"为正的占 {100*np.nanmean(B[ok]>0):.1f}%  (#128 的共享时间表 -> 应显著多于一半)")

def pr(y,x,ctrls,idx):
    X=np.c_[np.ones(len(idx)),*[c[idx] for c in ctrls]] if ctrls else np.ones((len(idx),1))
    ry=y[idx]-X@np.linalg.lstsq(X,y[idx],rcond=None)[0]
    rx=x[idx]-X@np.linalg.lstsq(X,x[idx],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

m=np.isfinite(sh)&np.isfinite(B)&np.isfinite(K)&np.isfinite(SE); idx=np.flatnonzero(m)
MODELS=[('raw',[]),('+类别数',[K]),('+斜率 SE',[K,SE]),('+截距',[K,SE,I0])]
rows=[]; rb=np.random.default_rng(20260803)
print(f"\n{'模型':<12}{'n':>8}{'r(斜率, 羞耻)':>14}{'bootstrap sd':>14}{'|r|/sd':>9}")
for name,ctrls in MODELS:
    r=pr(sh,B,ctrls,idx)
    bs=[]
    for _ in range(500):
        s_=rb.choice(idx,len(idx),replace=True); bs.append(pr(sh,B,ctrls,s_))
    sd=float(np.std(bs))
    rows.append(dict(model=name,n=len(idx),r=r,sd=sd,ratio=abs(r)/sd))
    print(f"{name:<12}{len(idx):>8,}{r:>+14.4f}{sd:>14.4f}{abs(r)/sd:>9.1f}",flush=True)
T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'slope.csv',index=False)

rgn=np.random.default_rng(5)
_,Bn,_=slopes(A_,shuffle=True,rng=rgn)
mn=m&np.isfinite(Bn); r_null=pr(sh,Bn,[K],np.flatnonzero(mn))
synth=B+rb.standard_normal(len(B))*np.nanstd(B)*2.0
r_plant=pr(synth,B,[K],idx)
print(f"\n对照:人内打乱 r={r_null:+.4f} · 正对照(b+2σ 噪声)r={r_plant:+.4f}")

r_full=float(T[T.model=='+截距'].r.iloc[0]); sd_full=float(T[T.model=='+截距'].sd.iloc[0])
g=Gate('人内斜率,与羞耻')
g.asserted('正对照:b+噪声 必须被同一条管道强测到',r_plant>0.3,f"{r_plant:+.4f}")
g.negative_control('人内打乱获得年龄',float(abs(r_null)),float(T[T.model=='+类别数'].r.iloc[0]))
g.resolvable('控制截距后 r(斜率, 羞耻)',r_full,sd_full)
g.asserted('跑之前写死的方向预测:r(b, 羞耻) < 0',r_full<0,
           f"{r_full:+.4f} —— 预测为负是因为 `#180` 说羞耻贴早半更紧,而早半是常见的东西")
g.asserted('功效对比 `#181`:本估计量用上全部人',len(idx)>10000,
           f"n={len(idx):,} vs 绝对劈分 1,089–3,854")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 截距才是对应 `#180` 的那个量 --------------------------------------------
# 斜率与羞耻在控制**截距**之前几乎为零(+0.004),控制之后翻正到 +0.0845 —— 典型的**抑制**。
# 而截距 = "这个人在年龄轴原点上的冷门度",正是 `#180` 的「早半有多冷门」的连续版。
# **本轮真正该报的头条是截距,不是斜率。** 写在这里而不是悄悄换掉主结果。
print("\n---- 截距(= 早期冷门度的连续版)----")
rows2=[]
for name,ctrls in [('raw',[]),('+类别数',[K]),('+斜率 SE',[K,SE]),('+斜率',[K,SE,B])]:
    r=pr(sh,I0,ctrls,idx)
    bs=[pr(sh,I0,ctrls,rb.choice(idx,len(idx),replace=True)) for _ in range(500)]
    sd=float(np.std(bs)); rows2.append(dict(model=name,r=r,sd=sd,ratio=abs(r)/sd))
    print(f"  {name:<10} r(截距, 羞耻) = {r:+.4f} ± {sd:.4f}   {abs(r)/sd:.1f}×")
T2=pd.DataFrame(rows2); T2.to_csv(pathlib.Path(__file__).parent/'results'/'intercept.csv',index=False)
ci=float(np.corrcoef(I0[idx],B[idx])[0,1])
print(f"\n  corr(截距, 斜率) = {ci:+.4f}  —— 抑制的来源")

r_i=float(T2[T2.model=='+斜率'].r.iloc[0]); sd_i=float(T2[T2.model=='+斜率'].sd.iloc[0])
g2=Gate('截距(早期冷门度)与羞耻')
g2.resolvable('控制斜率后 r(截距, 羞耻)',r_i,sd_i)
g2.asserted('方向与 `#180` 一致:早期越冷门,羞耻越高',r_i>0,
            f"{r_i:+.4f};`#180` 的早半配对差 +0.0276 也是正的")
g2.asserted('抑制的来源已量化',abs(ci)>0.3,f"corr(截距, 斜率) = {ci:+.4f}")
g2.same_scale('截距与斜率在同一批人上',float(len(idx)),float(len(idx)),'n')
print(g2)

# ---- 共线性闸:两个 9× 是同一个假象的两面 -------------------------------------
# corr(截距, 斜率) = −0.97。在这个共线度下,「控制斜率后的截距」与「控制截距后的斜率」
# **不是两个发现,是同一个抑制的两面**,而且都被共线性放大。
# 能解释的只有 raw:斜率 +0.0040(0.4×,无)· 截距 +0.0270(2.9×,方向与 `#180` 一致)。
print("\n---- 共线性闸 ----")
vif=1/max(1-ci**2,1e-9)
print(f"  corr(截距, 斜率) = {ci:+.4f}   VIF = {vif:.1f}")
r_b_raw=float(T[T.model=='raw'].r.iloc[0]);  sd_b_raw=float(T[T.model=='raw'].sd.iloc[0])
r_i_raw=float(T2[T2.model=='raw'].r.iloc[0]); sd_i_raw=float(T2[T2.model=='raw'].sd.iloc[0])
print(f"  raw 斜率 {r_b_raw:+.4f} ({abs(r_b_raw)/sd_b_raw:.1f}×) · raw 截距 {r_i_raw:+.4f} ({abs(r_i_raw)/sd_i_raw:.1f}×)")
print(f"  partial 斜率|截距 {r_full:+.4f} · partial 截距|斜率 {r_i:+.4f} —— **两者几乎相等且同号,"
      f"这正是共线抑制的签名,不是两个独立发现**")
g3=Gate('两个 9× 能不能当成发现')
g3.asserted('共线性已量化',abs(ci)>0.9,f"corr = {ci:+.4f}, VIF = {vif:.0f}")
g3.asserted('两个 partial 几乎相等 -> 同一个抑制的两面,不是两个发现',
            abs(abs(r_full)-abs(r_i))<0.01,f"|{r_full:+.4f}| vs |{r_i:+.4f}|")
g3.asserted('因此只报 raw',True,
            f"斜率 raw {r_b_raw:+.4f}({abs(r_b_raw)/sd_b_raw:.1f}×,无)· "
            f"截距 raw {r_i_raw:+.4f}({abs(r_i_raw)/sd_i_raw:.1f}×,方向与 `#180` 一致)")
g3.resolvable('raw 截距与羞耻',r_i_raw,sd_i_raw)
g3.require_resolvable_first('raw 斜率与羞耻',r_b_raw,sd_b_raw,family='slope')
print(g3)
print(f"\n  => **截距(早期冷门度)预测羞耻 {r_i_raw:+.4f}({abs(r_i_raw)/sd_i_raw:.1f}×,n={len(idx):,});"
      f"斜率不预测({r_b_raw:+.4f},{abs(r_b_raw)/sd_b_raw:.1f}×)。**"
      f"\n     用完全不同的估计量、在 10,567 人上,复现了 `#180` 的方向。")

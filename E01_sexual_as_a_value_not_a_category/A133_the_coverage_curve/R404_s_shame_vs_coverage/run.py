import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A133 R404 -- `S ↔ 羞耻` 沿覆盖阈的整条曲线

⚠ **本轮的起点是一次 L3 闸口的拦截**(第三次):我以为页面上有「第二道羞耻题」,
读了才发现那是 **therapeutic 判别项**,不是第二个测量 —— 页面早写着。
**而读页面时看到一条从没被连起来的线。**

三块已有的证据,分别在三个地方,从没并排过:
- 页面自己:**最严覆盖阈(≥16 块,n=950)羞耻衰减到 +0.0156**,6/7 个阈同号;
- `#346b`:**放宽到 ≥4 时是 +0.1715**(≥8 是 +0.1185);
- `#357a`:**`S` 与广度共线**(`corr(S, 勾选数) = +0.608`,页面已带零 +0.719)。

**⇒ 合起来:`S ↔ 羞耻` 可能沿覆盖阈**单调下降** ——
即这个关系在 `S` 被**测得最差**的地方**最强**。那与测量误差的方向相反。**

ESTIMAND        覆盖阈 ∈ {4,6,8,10,12,16}:各自的 n · `corr(S, 羞耻)` · **自助 95% 区间** ·
                **控制勾选数后**的偏相关;并判**单调性**。
KILL            **若单调下降且区间在两端不重叠 -> 这是一条关于「谁被纳入」的强信号,
                页面上的 +0.1185 必须带上整条曲线;
                若曲线平或非单调 -> 那两个端点只是抽样,页面照旧。**
POSITIVE CTRL   合成一个**真值固定**的 `S↔y` 关系 -> 曲线必须**平**(证明曲线不是纳入本身造出来的)。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ 嵌套样本     六个阈的样本是**嵌套**的 -> 相邻点高度相关,**不能当独立点做趋势检验**;
                本轮报**区间是否重叠**与**两端差**,不报 p。
IMPOSSIBLE      覆盖阈同时改变**样本**与 `S` 的**估计质量**;本轮分不开这两者。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

cov=np.zeros(NN); PK=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1; PK[ppl]+=M.sum(1)
def Spos(mask):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nn=M.sum(1)
        v=np.where(nn>0,(M@rr)/np.maximum(nn,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(mask&(cv>=1),ps/np.maximum(cv,1),np.nan)
def arm(TH,y=None):
    mk=cov>=TH; S=Spos(mk); yy=sh if y is None else y
    m=mk&np.isfinite(S)&np.isfinite(yy)&np.isfinite(PK)
    n=int(m.sum())
    if n<200: return None
    a,b=S[m],yy[m]; pk=PK[m]
    r=float(np.corrcoef(a,b)[0,1])
    z=lambda v:(v-v.mean())/max(v.std(),1e-12)
    X=np.column_stack([np.ones(n),z(pk)])
    ra=z(a)-X@np.linalg.lstsq(X,z(a),rcond=None)[0]
    rb=z(b)-X@np.linalg.lstsq(X,z(b),rcond=None)[0]
    rp=float(np.corrcoef(ra,rb)[0,1])
    rg=np.random.default_rng(500+TH)
    bs=np.array([float(np.corrcoef(a[i],b[i])[0,1]) for i in (rg.integers(0,n,n) for _ in range(400))])
    q=np.percentile(bs,[2.5,97.5])
    return dict(th=TH,n=n,r=r,lo=q[0],hi=q[1],rp=rp)
THS=[4,6,8,10,12,16]
print(f"{'阈':>4}{'n':>9}{'corr(S,羞耻)':>14}{'95% 区间':>22}{'控勾选数后':>12}")
R=[]
for TH in THS:
    a=arm(TH)
    if a is None: continue
    R.append(a)
    print(f"{TH:>4}{a['n']:>9,}{a['r']:>+14.4f}   [{a['lo']:+.4f}, {a['hi']:+.4f}]{a['rp']:>+12.4f}")
lo_end,hi_end=R[0],R[-1]
mono=all(R[i]['r']>=R[i+1]['r']-1e-9 for i in range(len(R)-1))
ovl=(lo_end['lo']<hi_end['hi'] and hi_end['lo']<lo_end['hi'])
print(f"\n★ 两端:阈 {lo_end['th']} **{lo_end['r']:+.4f}** [{lo_end['lo']:+.4f}, {lo_end['hi']:+.4f}] "
      f"-> 阈 {hi_end['th']} **{hi_end['r']:+.4f}** [{hi_end['lo']:+.4f}, {hi_end['hi']:+.4f}]")
print(f"   单调下降 **{'是' if mono else '否'}** · 两端区间重叠 **{'是' if ovl else '否'}** · "
      f"差 **{lo_end['r']-hi_end['r']:+.4f}**")
rg=np.random.default_rng(31)
S8=Spos(cov>=4); m8=np.isfinite(S8)
ysyn=np.full(NN,np.nan)
zz=(S8[m8]-np.nanmean(S8[m8]))/np.nanstd(S8[m8])
ysyn[m8]=0.15*zz+rg.standard_normal(int(m8.sum()))
print(f"\n正对照(真值固定 0.15 的合成结局):")
PC=[]
for TH in THS:
    a=arm(TH,y=ysyn)
    if a: PC.append(a['r']); print(f"   阈 {TH:>2}: **{a['r']:+.4f}**")
print(f"   曲线全距 **{max(PC)-min(PC):.4f}**(平 = 曲线不是纳入本身造出来的)")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
NG=[arm(8,y=perm_finite(sh,700+i))['r'] for i in range(20)]
print(f"负对照(阈 8,打乱人):**{np.mean(NG):+.4f} ± {np.std(NG):.4f}**")
T=pd.DataFrame(R); check_columns(T,'R404')
T.to_csv(pathlib.Path(__file__).parent/'results'/'curve.csv',index=False)
gg=Gate('`S ↔ 羞耻` 沿覆盖阈的曲线')
gg.asserted('★ 正对照:真值固定的合成结局 -> 曲线必须平(全距 < 0.05)',max(PC)-min(PC)<0.05,
            f"全距 {max(PC)-min(PC):.4f}")
gg.negative_control('★ 负对照:阈 8 打乱人',float(np.mean(NG)),float([x for x in R if x['th']==8][0]['r']),
    null_spread=float(np.std(NG)),null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill:曲线是否单调下降且两端区间不重叠',mono and not ovl,
            f"单调 {'是' if mono else '否'} · 两端重叠 {'是' if ovl else '否'} · "
            f"{lo_end['r']:+.4f} -> {hi_end['r']:+.4f}")
gg.asserted('⚠ 嵌套样本:相邻点高度相关,不做趋势 p 值',True,
            '本轮报区间是否重叠与两端差,不报 p')
gg.asserted('⚠ 边界:覆盖阈同时改变样本与 `S` 的估计质量',True,'本轮分不开这两者')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

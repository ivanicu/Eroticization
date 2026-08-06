import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A130 R399 -- 「类别数」不是控制项,是一个预测量

`#354c`:「一个人报了多少个起始类别」解释了 `animated` **增量**份额的 **17.6%**,而对羞耻是 **0.0%**。
**一个控制项贡献了六分之一,它就不再是控制项。**

⚠⚠ **跑之前写下**:**增量份额**(**唯一**贡献)与**「扣掉后掉多少」**(**总**贡献)
**不是同一个量** —— 共线时前者小后者大。**两个都报,分开读。**

ESTIMAND        ① `corr(animated, 类别数)`,**去衰减**(类别分半信度 `#332` = 0.9242)+ **自助 95% 区间**;
                ② 六坐标对 `animated` 的联合 R²,**有/无**类别数;
                ③ 同样两项对**羞耻**作参照。
KILL            **若扣掉它之后 `animated` 的 R² 掉幅大 -> 页面上「最好解释的一格 7.79%」
                要带一句「其中约 X 是报告广度」;若掉幅小 -> 增量份额高只是共线的记账方式。**
POSITIVE CTRL   合成一个**只由类别数驱动**的结局 -> 全流程必须复原。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ 窄口径       绝对量比较(`CALIBER.md` ⑩)。
IMPOSSIBLE      「报了多少类别」既是**兴趣广度**也是**答题行为**;本轮分不开这两者。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ani=pd.to_numeric(d['animated'],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
onsc=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in onsc])
HASc=np.isfinite(ONS); ncat=HASc.sum(1).astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); CO=[Q[0],Q[1],Q[2],Q[3],-Q[4],Q[5]]
m0=ok.copy()
for q_ in CO: m0&=np.isfinite(q_)
m0&=np.isfinite(sh)&np.isfinite(ani)
n=int(m0.sum()); z=lambda v:(v[m0]-v[m0].mean())/max(v[m0].std(),1e-12)
Z=[z(q_) for q_ in CO]; zn=z(ncat)
print(f"n={n:,} · 窄口径")
print(f"⚠⚠ 增量份额(唯一贡献)与「扣掉后掉多少」(总贡献)**不是同一个量** —— 两个都报,分开读。\n")
rg=np.random.default_rng(9)
o=rg.permutation(HASc.shape[1]); h=HASc.shape[1]//2
ha,hb=HASc[:,o[:h]].sum(1).astype(float),HASc[:,o[h:2*h]].sum(1).astype(float)
rr=float(np.corrcoef(ha[m0],hb[m0])[0,1]); rel=2*rr/(1+rr)
for nm,y in (('animated',ani),('羞耻',sh)):
    r0=float(np.corrcoef(zn,z(y))[0,1]); dis=r0/np.sqrt(max(rel,1e-9))
    idx=np.flatnonzero(m0); rgb=np.random.default_rng(404)
    bs=[]
    for _ in range(400):
        ii=idx[rgb.integers(0,len(idx),len(idx))]
        a1,b1=ha[ii],hb[ii]; rr2=float(np.corrcoef(a1,b1)[0,1]); rl=min(2*rr2/(1+rr2),1.0)
        ro=float(np.corrcoef(ncat[ii],y[ii])[0,1]); bs.append(ro/np.sqrt(max(rl,1e-9)))
    q=np.nanpercentile(bs,[2.5,97.5])
    print(f"① {nm:<10} `corr(结局, 类别数)` **{r0:+.4f}** · 信度 {rel:.4f} · "
          f"**去衰减 {dis:+.4f}** · 95% 区间 **[{q[0]:+.4f}, {q[1]:+.4f}]**")
def r2(cols,y):
    yy=z(y); X=np.column_stack([np.ones(n)]+cols); b,*_=np.linalg.lstsq(X,yy,rcond=None)
    r=yy-X@b; return 1-float(r@r)/float(((yy-yy.mean())**2).sum())
print()
for nm,y in (('animated',ani),('羞耻',sh)):
    w=r2(Z+[zn],y); wo=r2(Z,y); only=r2([zn],y)
    print(f"② {nm:<10} 六坐标+类别数 **{100*w:.3f}%** · **只六坐标 {100*wo:.3f}%** · "
          f"只类别数 **{100*only:.3f}%** · **扣掉类别数掉 {100*(w-wo):.3f}pp({100*(w-wo)/w:.1f}%)**")
rgp=np.random.default_rng(77)
ysyn=np.full(NN,np.nan); ysyn[m0]=0.4*zn+rgp.standard_normal(n)
w=r2(Z+[zn],ysyn); wo=r2(Z,ysyn)
print(f"\n正对照(只由类别数驱动):六+类别数 **{100*w:.3f}%** · 只六坐标 **{100*wo:.3f}%** · "
      f"掉 **{100*(w-wo):.3f}pp**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[float(np.corrcoef(z(perm_finite(ncat,300+i)),z(ani))[0,1]) for i in range(20)]
print(f"负对照(打乱人):`corr(animated, 类别数)` **{np.mean(nul):+.4f} ± {np.std(nul):.4f}**")
wA=r2(Z+[zn],ani); woA=r2(Z,ani)
T=pd.DataFrame([dict(v_arm='animated 六+n',v_r2=100*wA),dict(v_arm='animated 只六',v_r2=100*woA),
                dict(v_arm='羞耻 六+n',v_r2=100*r2(Z+[zn],sh)),dict(v_arm='羞耻 只六',v_r2=100*r2(Z,sh))])
check_columns(T,'R399'); T.to_csv(pathlib.Path(__file__).parent/'results'/'ncat.csv',index=False)
gg=Gate('「类别数」不是控制项,是一个预测量')
gg.asserted('★ 正对照:只由类别数驱动的结局必须复原(掉幅 > 5pp)',100*(w-wo)>5,
            f"掉 {100*(w-wo):.3f}pp")
gg.negative_control('★ 负对照:打乱人后 `corr(animated, 类别数)`',float(np.mean(nul)),
    float(np.corrcoef(zn,z(ani))[0,1]),null_spread=float(np.std(nul)),
    null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill:扣掉类别数后 `animated` 的 R² 掉幅是否大(> 10%)',
            (wA-woA)/wA>0.10,
            f"六+类别数 {100*wA:.3f}% -> 只六坐标 {100*woA:.3f}%,掉 {100*(wA-woA):.3f}pp"
            f"({100*(wA-woA)/wA:.1f}%)")
gg.asserted('⚠⚠ 增量份额 ≠ 扣掉后掉多少',True,
            f"`#354a` 的增量份额是 17.6%(**唯一**贡献),本轮的掉幅是 {100*(wA-woA)/wA:.1f}%(**总**贡献)")
gg.asserted('⚠ 边界:「报了多少类别」既是兴趣广度也是答题行为',True,'本轮分不开这两者')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A114 R369 -- `form` 是不是「反事实性」

`#316`:`form`(动画与文字共有的东西)去衰减后与羞耻 **+0.1293**;
它**不是媒介**(与 `pornhabit`/付费/开始年龄/暴力比例的相关全 ≤0.052),
**也不是答题风格**(扣掉后仍保留 61%)。**那它是什么?**

⚠ 这份问卷里**没有**独立的「虚构 vs 真人」题(唯二的媒介题就是 `animated`/`written`,
它们**构成** `form`)。**但有一组更好的、而且我从没在 `form` 里用过的**:
四道「想象**成为**另一性别」的题(**存在为** / **独自自慰为** × 生理女 / 生理男)。

**判据规则,跑之前钉死(确定性规则,没有研究者自由度)**:
- **反事实半** = 与本人 `biomale` **相反**的那两道(男 -> 女的两题;女 -> 男的两题);
- **同构半** = 与本人性别**一致**的那两道。
**「情欲化一个不是你的身体」正是反事实性;而同构半是同一题族、同一格式的对照。**

ESTIMAND        `corr(form, 反事实半)` 与 `corr(form, 同构半)`,并报两者之差;
                同时扣掉 `biomale` 与答题风格。
KILL            **若反事实半明显强于同构半 -> `form` 里有「反事实性」这个成分;
                若两半相当 -> 不是反事实性,是这一族题共有的别的东西(那也要说出来);
                若同构半更强 -> 反过来,而那会杀掉这个假设。**
POSITIVE CTRL   合成一个**只由反事实半驱动**的量 -> 判据必须抓到,且同构半必须落零。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ guard 15     若用剖面,必须同时给分数层相关(本轮直接给分数层,不用剖面)。
IMPOSSIBLE      「反事实」在这里只有**性别**这一个维度的操作化;
                它测不到「虚构角色 / 非人 / 不可能场景」那几种反事实。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)

g=lambda k:next(c for c in d.columns if k in str(c))
EXF=g('existing (in *nonsexual* situations) as a biological *female*')
EXM=g('existing (in *nonsexual* situations) as a biological *male*')
MBF=g('masturbating alone as a biological female')
MBM=g('masturbating alone as a biological male')
SHAME=next(c for c in d.columns if 'ashamed' in str(c))
num=lambda c: pd.to_numeric(d[c],errors='coerce').values.astype(float)
exf,exm,mbf,mbm=num(EXF),num(EXM),num(MBF),num(MBM)
sex=num('biomale'); sh=num(SHAME)
ani,wri=num('animated'),num('written')
zz=lambda v,m:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
m0=np.isfinite(ani)&np.isfinite(wri)&np.isfinite(sex)&np.isfinite(exf)&np.isfinite(exm)&np.isfinite(mbf)&np.isfinite(mbm)
FORM=np.full(NN,np.nan); FORM[m0]=(zz(ani,m0)+zz(wri,m0))/2
CF=np.full(NN,np.nan); CG=np.full(NN,np.nan)      # ⚠ 规则在跑之前钉死:反事实 = 与本人性别相反
male=(sex==1)
CF[m0]=np.where(male[m0],(zz(exf,m0)+zz(mbf,m0))/2,(zz(exm,m0)+zz(mbm,m0))/2)
CG[m0]=np.where(male[m0],(zz(exm,m0)+zz(mbm,m0))/2,(zz(exf,m0)+zz(mbf,m0))/2)
LKc=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
LK=np.column_stack([num(c) for c in LKc])
STY=[np.nanmean(LK,1),np.nanstd(LK,1),np.isfinite(LK).sum(1).astype(float)]
def cor(u,v,m=None):
    k=np.isfinite(u)&np.isfinite(v)&(m0 if m is None else m)
    return (float(np.corrcoef(u[k],v[k])[0,1]),int(k.sum())) if k.sum()>200 else (np.nan,0)
def partial(u,y,ctrl):
    m=np.isfinite(u)&np.isfinite(y)&m0
    for c in ctrl: m&=np.isfinite(c)
    X=np.column_stack([np.ones(m.sum())]+[zz(c,m) for c in ctrl])
    ru=zz(u,m)-X@np.linalg.lstsq(X,zz(u,m),rcond=None)[0]
    ry=zz(y,m)-X@np.linalg.lstsq(X,zz(y,m),rcond=None)[0]
    return float(np.corrcoef(ru,ry)[0,1]),int(m.sum())
print(f"n={int(m0.sum()):,} · 男 {int(male[m0].sum()):,} / 女 {int((~male[m0]).sum()):,}")
r_cf,n_=cor(FORM,CF); r_cg,_=cor(FORM,CG)
p_cf,_=partial(FORM,CF,[sex]+STY); p_cg,_=partial(FORM,CG,[sex]+STY)
print(f"\n★ `form` ↔ **反事实半** **{r_cf:+.4f}** · ↔ **同构半** **{r_cg:+.4f}** · "
      f"差 **{r_cf-r_cg:+.4f}**(n={n_:,})")
print(f"   扣掉 `biomale` + 答题风格后:反事实 **{p_cf:+.4f}** · 同构 **{p_cg:+.4f}** · "
      f"差 **{p_cf-p_cg:+.4f}**")
print(f"⚠ 两半彼此相关 `corr(CF, CG)` = **{cor(CF,CG)[0]:+.4f}**")
print(f"⚠ 各自与 `biomale`:CF **{cor(CF,sex)[0]:+.4f}** · CG **{cor(CG,sex)[0]:+.4f}** · "
      f"form **{cor(FORM,sex)[0]:+.4f}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[cor(perm_finite(FORM,600+i),CF)[0] for i in range(20)]
print(f"负对照(打乱人)↔反事实半:**{np.mean(nul):+.4f} ± {np.std(nul):.4f}**")
rg=np.random.default_rng(13)
SYN=np.full(NN,np.nan); SYN[m0]=zz(CF,m0)+rg.standard_normal(int(m0.sum()))
s_cf,_=cor(SYN,CF); s_cg,_=cor(SYN,CG); sp_cf,_=partial(SYN,CF,[sex]+STY); sp_cg,_=partial(SYN,CG,[sex]+STY)
print(f"\n正对照(只由反事实半驱动的合成量):↔反事实 **{s_cf:+.4f}** · ↔同构 **{s_cg:+.4f}** · "
      f"扣控制后 {sp_cf:+.4f} / {sp_cg:+.4f}")
T=pd.DataFrame([dict(v_arm='原样',cf=r_cf,cg=r_cg,v_diff=r_cf-r_cg),
                dict(v_arm='扣sex+风格',cf=p_cf,cg=p_cg,v_diff=p_cf-p_cg)])
check_columns(T,'R369'); T.to_csv(pathlib.Path(__file__).parent/'results'/'counterfactual.csv',index=False)
mde=2.8/np.sqrt(max(n_,1))
gg=Gate('`form` 是不是反事实性')
gg.asserted('★ 正对照:只由反事实半驱动的合成量 -> 反事实强、同构弱',
            s_cf>0.5 and abs(s_cg)<0.3 and (sp_cf-sp_cg)>0.3,
            f"↔反事实 {s_cf:+.4f} · ↔同构 {s_cg:+.4f} · 扣控制后差 {sp_cf-sp_cg:+.4f}")
gg.negative_control('★ 负对照:打乱人后 ↔ 反事实半',float(np.mean(nul)),r_cf,
    null_spread=float(np.std(nul)),null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill:反事实半是否明显强于同构半(差 > 2×MDE)',
            (r_cf-r_cg)>2*mde,
            f"原样差 **{r_cf-r_cg:+.4f}** · 扣控制后 **{p_cf-p_cg:+.4f}** vs 2×MDE **{2*mde:.4f}**")
gg.asserted('⚠ 边界:「反事实」只操作化到性别这一个维度',True,
            '测不到「虚构角色 / 非人 / 不可能场景」那几种反事实')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

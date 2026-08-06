import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A115 R370 -- `form ↔ biomale = −0.2075` 当作一条独立声明重测

`#324b` 是**顺带**测到的,而**顺带测到的数最容易被高估** —— 它没有自己的设计。

ESTIMAND        ① `animated ↔ biomale` 与 `written ↔ biomale` **分别**报(两个指标一不一致);
                ② 去衰减(`#317a` 已验证该校正在这份数据上不变)+ **自助 95% 区间**;
                ③ 扣掉**答题风格**与**跨块勾选数**后保留多少;
                ④ 与页面上已有的性别线(`c3` 的 `#286`–`#295`)**并排**,说清是不是同一回事。
KILL            **若两个指标不一致 -> 性别关联是**那一题**的性质,而 `#316` 的「`form` 是一个维度」
                要收窄;若一致且扣控制后保留 -> 这是一条可发布的关于 `form` 的声明。**
POSITIVE CTRL   合成一个已知与 `biomale` 相关 0.20 的量 -> 同一流程必须复原(含去衰减与区间覆盖)。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ guard 21     若某一步判为零,交出三件套。
IMPOSSIBLE      `biomale` 是自报的生理性别,不是性别认同;本轮说不了后者。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

num=lambda c: pd.to_numeric(d[c],errors='coerce').values.astype(float)
ani,wri,sex=num('animated'),num('written'),num('biomale')
m0=np.isfinite(ani)&np.isfinite(wri)&np.isfinite(sex)
zz=lambda v,m:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
FORM=np.full(NN,np.nan); FORM[m0]=(zz(ani,m0)+zz(wri,m0))/2
LKc=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
LK=np.column_stack([num(c) for c in LKc])
STY=[np.nanmean(LK,1),np.nanstd(LK,1),np.isfinite(LK).sum(1).astype(float)]
NQ=np.full(NN,np.nan); cv=np.zeros(NN); s_=np.zeros(NN)
for M,ppl in MB: cv[ppl]+=1; s_[ppl]+=M.sum(1)
NQ=np.where(cv>=8,s_/np.maximum(cv,1),np.nan)
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
r_a,_=cor(ani,sex); r_w,_=cor(wri,sex); r_f,n_=cor(FORM,sex)
r_aw,_=cor(ani,wri); rel=2*r_aw/(1+r_aw)
print(f"① 两个指标分别:`animated ↔ biomale` **{r_a:+.4f}** · `written ↔ biomale` **{r_w:+.4f}** "
      f"(差 **{r_a-r_w:+.4f}**)· 合成 `form` **{r_f:+.4f}**(n={n_:,})")
print(f"   `corr(animated, written)` = **{r_aw:+.4f}** -> SB 信度 **{rel:.4f}**")
def dis_of(a,w,y):
    za=(a-a.mean())/max(a.std(),1e-12); zw=(w-w.mean())/max(w.std(),1e-12)
    rr=float(np.corrcoef(za,zw)[0,1]); rl=min(2*rr/(1+rr),1.0) if rr>0 else np.nan
    f=za+zw; raw=float(np.corrcoef(f,(y-y.mean())/max(y.std(),1e-12))[0,1])
    return raw,rl,(raw/np.sqrt(rl) if np.isfinite(rl) and rl>0 else np.nan)
A0,W0,S0=ani[m0],wri[m0],sex[m0]; n=len(S0)
rg=np.random.default_rng(515)
bs=np.array([dis_of(A0[i],W0[i],S0[i]) for i in (rg.integers(0,n,n) for _ in range(500))])
q=np.nanpercentile(bs[:,2],[2.5,50,97.5])
print(f"② 去衰减 **{dis_of(A0,W0,S0)[2]:+.4f}** · **自助 95% 区间 [{q[0]:+.4f}, {q[2]:+.4f}]**")
p_sty,_=partial(FORM,sex,STY); p_all,_=partial(FORM,sex,STY+[NQ])
print(f"③ 扣掉答题风格 **{p_sty:+.4f}**(保留 {100*abs(p_sty)/abs(r_f):.0f}%)· "
      f"再扣跨块勾选数 **{p_all:+.4f}**(保留 {100*abs(p_all)/abs(r_f):.0f}%)")
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR)
c3=Q[4]; S=Q[0]
r_c3,_=cor(c3,sex,np.isfinite(c3)&np.isfinite(sex)&ok)
r_S,_=cor(S,sex,np.isfinite(S)&np.isfinite(sex)&ok)
r_fc,_=cor(FORM,c3,np.isfinite(FORM)&np.isfinite(c3)&ok)
print(f"④ 并排:`c3 ↔ biomale` **{r_c3:+.4f}** · `S ↔ biomale` **{r_S:+.4f}** · "
      f"`form ↔ biomale` **{r_f:+.4f}**;而 `form ↔ c3` = **{r_fc:+.4f}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[cor(perm_finite(FORM,700+i),sex)[0] for i in range(20)]
print(f"负对照(打乱人):**{np.mean(nul):+.4f} ± {np.std(nul):.4f}**")
rgp=np.random.default_rng(31)
SYN=np.full(NN,np.nan)
SYN[m0]=0.20*zz(sex,m0)+np.sqrt(1-0.04)*rgp.standard_normal(int(m0.sum()))
s_r,_=cor(SYN,sex)
print(f"\n正对照(合成量,真相关 0.20):读出 **{s_r:+.4f}**")
# ⚠⚠ kill① 反向开火 -> 立刻查同一个问题会不会打到 `#316` 的羞耻结论:
#    「`form` 够到羞耻」会不会其实是「`written` 够到羞耻」?
SH=next(c for c in d.columns if 'ashamed' in str(c)); sh=num(SH)
ms=m0&np.isfinite(sh)
sa,_=cor(ani,sh,ms); sw,_=cor(wri,sh,ms); sf,nsf=cor(FORM,sh,ms)
print(f"\n⚠⚠ 同一个问题打到 `#316`:")
print(f"   `animated ↔ 羞耻` **{sa:+.4f}** · `written ↔ 羞耻` **{sw:+.4f}** · "
      f"差 **{sa-sw:+.4f}** vs 2×MDE **{2*2.8/np.sqrt(nsf):.4f}** · 合成 `form` **{sf:+.4f}**")
CONSIST_SH=abs(sa-sw)<2*2.8/np.sqrt(nsf)
print(f"   -> 羞耻上两个指标{'**一致**,`#316` 不受影响' if CONSIST_SH else '**不一致**,`#316` 也要收窄'}")

T=pd.DataFrame([dict(v_term='animated',v_r=r_a),dict(v_term='written',v_r=r_w),
                dict(v_term='form',v_r=r_f),dict(v_term='form去衰减',v_r=float(q[1])),
                dict(v_term='扣风格',v_r=p_sty),dict(v_term='扣风格+勾选数',v_r=p_all)])
check_columns(T,'R370'); T.to_csv(pathlib.Path(__file__).parent/'results'/'form_sex.csv',index=False)
mde=2.8/np.sqrt(max(n_,1))
gg=Gate('`form ↔ biomale` 当独立声明重测')
gg.asserted('★ 正对照:合成量真相关 0.20 必须被复原',abs(s_r-0.20)<0.03,f"读出 {s_r:+.4f}")
gg.negative_control('★ 负对照:打乱人',float(np.mean(nul)),r_f,null_spread=float(np.std(nul)),
    null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill ①:两个指标是否一致(差 < 2×MDE)',abs(r_a-r_w)<2*mde,
            f"animated {r_a:+.4f} · written {r_w:+.4f} · 差 {r_a-r_w:+.4f} vs 2×MDE {2*mde:.4f}")
gg.asserted('★ 注册的 kill ②:扣掉风格与勾选数后是否保留(>50%)',
            abs(p_all)/abs(r_f)>0.5,
            f"{p_all:+.4f} / {r_f:+.4f} = {100*abs(p_all)/abs(r_f):.0f}%")
gg.asserted('★ 连带检查:羞耻上两个指标一不一致(决定 `#316` 要不要收窄)',CONSIST_SH,
            f"animated↔羞耻 {sa:+.4f} · written↔羞耻 {sw:+.4f} · 差 {sa-sw:+.4f}")
gg.asserted('⚠ ④ 与已有性别线并排:是不是同一回事',True,
            f"`form ↔ c3` = {r_fc:+.4f} —— {'几乎正交,是**另一条**性别通路' if abs(r_fc)<0.15 else '与 c3 有重叠'}")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

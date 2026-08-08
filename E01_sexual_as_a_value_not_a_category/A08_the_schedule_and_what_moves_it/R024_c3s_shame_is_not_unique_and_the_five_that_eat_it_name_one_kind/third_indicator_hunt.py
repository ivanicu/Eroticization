import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A139 R419 -- 问卷里还有没有第三个指向 `form` 的题

页面上最脆的一步:`#316`/`#325` —— `form`(`animated`/`written`)的 Spearman–Brown 可靠性只有 **0.380**,
而 disattenuation **拿这个 0.380 当分母**。分母越小,放大越猛,而 0.380 已经很小。
**若能找到第三个同构念的题,α 会升,分母会稳,页面上那一步就不再是最脆的。**

ESTIMAND        ① `form` 的残差(去掉六坐标)对**所有**数值型人层量求相关,报**最大值的分布**;
                ② 对最高的几个,逐个问「**它是不是同一个构念**」——
                   判据 = (a) 与 `animated`、`written` **都**正相关且两者之比在 [0.5, 2](载荷大致均衡);
                          (b) α({a, w, X}) > α({a, w})。
KILL(条件式)  仅当判据在正对照上**开火**、在负对照上**不开火** -> 判:是否存在任何一个候选同时满足 (a)(b)。
                有 -> 页面上 0.380 那一步可以被修;没有 -> **这份问卷里没有第三个指标**,
                那条「测得太粗」的话就是最终结论,而不是一个待办。
POSITIVE CTRL   合成一个**真的第三个平行指标**(共同潜变量 + 噪声)-> 判据必须开火。
                ⚠ `#374b` 的教训:**一个从未开火过的判据,它的每一个「不是」都是沉默,不是无罪。**
NEGATIVE CTRL   合成纯噪声列 -> 判据必须不开火。
⚠ 零的种类     `offset_control`:**最大相关的零绝不是零** —— 在数百个候选里取最大值,最大值天然为正。
                零 = 对**同样多的合成噪声列**取最大相关的分布。
⚠ 多重性       报分布不报单格。
⚠ `#357b`      「相关高」可能只是**它已经在页面上了**(同一个量的另一种写法)-> 命中的要逐个人工看名字。
IMPOSSIBLE      「同一构念」在只有两个指标时**无法**从数据本身判定(两点定不出一条曲线);
                本轮的判据是**必要条件**,不是充分条件。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_R410=(ROOT/'E01_sexual_as_a_value_not_a_category/A136_is_c3_shame_its_own/R410_commonality_vs_person_variables/run.py').read_text()
exec(_R410.split('"""',2)[2].split('def num(c)')[0])

AN=pd.to_numeric(d['animated'],errors='coerce').values.astype(float)
WR=pd.to_numeric(d['written'],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR)
base=ok.copy()
for q_ in Q: base&=np.isfinite(q_)
m0=base&np.isfinite(AN)&np.isfinite(WR)
n0=int(m0.sum())
def sb(a,b,g):
    r=float(np.corrcoef(a[g],b[g])[0,1]); return 2*r/(1+r),r
SB2,R_AW=sb(AN,WR,m0)
print(f"n={n0:,} · corr(animated, written) = **{R_AW:+.4f}** · Spearman–Brown α₂ = **{SB2:.4f}**"
      f"(页面上写的是 0.380)\n")

# form 的残差:去掉六坐标
def resid6(v,g):
    X=np.column_stack([np.ones(int(g.sum()))]+[ (q_[g]-q_[g].mean())/q_[g].std() for q_ in Q])
    y=v[g]; b,*_=np.linalg.lstsq(X,y,rcond=None); out=np.full(NN,np.nan); out[g]=y-X@b; return out
FORM=np.where(m0,(AN+WR)/2,np.nan); FR=resid6(FORM,m0)

# 候选:所有能转成数值、且在 m0 上有 ≥2000 个有限值的列(排除两个指标本身)
CAND=[]
for c in d.columns:
    if c in ('animated','written'): continue
    v=pd.to_numeric(d[c],errors='coerce').values.astype(float)
    if np.isfinite(v[m0]).sum()>=2000 and np.nanstd(v[m0])>1e-9: CAND.append((c,v))
print(f"① 候选 **{len(CAND)}** 个数值型人层量;对 `form` 残差求相关:")
cor=[]
for c,v in CAND:
    g=m0&np.isfinite(v)&np.isfinite(FR)
    if g.sum()<2000: continue
    cor.append((abs(float(np.corrcoef(FR[g],v[g])[0,1])),c,v))
cor.sort(reverse=True,key=lambda t:t[0])
print(f"   最大 |r| = **{cor[0][0]:.4f}**({cor[0][1][:56]})")
for k,(r_,c,_) in enumerate(cor[:6]): print(f"   {k+1}. {r_:.4f}  {c[:70]}")

rgO=np.random.default_rng(1234)
NOFF=200; K=len(cor)
offmax=[]
for s in range(NOFF):
    rg=np.random.default_rng(2000+s)
    mx=0.
    for _ in range(K):
        v=rg.standard_normal(NN); g=m0&np.isfinite(FR)
        mx=max(mx,abs(float(np.corrcoef(FR[g],v[g])[0,1])))
    offmax.append(mx)
offmax=np.array(offmax)
print(f"⚠ offset 零(**同样多({K})的合成噪声列取最大相关**,{NOFF} 次):"
      f"**{offmax.mean():.4f} ± {offmax.std():.4f}** · 95 分位 **{np.percentile(offmax,95):.4f}**")
print(f"   -> 实测最大 {cor[0][0]:.4f} "
      f"{'**越阈**' if cor[0][0]>np.percentile(offmax,95) else '**未越阈**'}\n")

# ② 同一构念的判据
def alpha3(a,b,c,g):
    Z=np.column_stack([(x[g]-x[g].mean())/x[g].std() for x in (a,b,c)])
    C=np.corrcoef(Z.T); k=3; rbar=(C.sum()-k)/(k*(k-1))
    return k*rbar/(1+(k-1)*rbar)
A2=2*R_AW/(1+R_AW)
def same_construct(v):
    g=m0&np.isfinite(v)
    if g.sum()<2000: return False,{}
    ra=float(np.corrcoef(v[g],AN[g])[0,1]); rw=float(np.corrcoef(v[g],WR[g])[0,1])
    bal=(ra>0 and rw>0 and 0.5<=ra/max(rw,1e-9)<=2.0)
    a3=alpha3(AN,WR,v,g)
    return bool(bal and a3>A2),dict(ra=ra,rw=rw,a3=a3)
rgC=np.random.default_rng(7)
lat=np.where(m0,(AN+WR)/2,np.nan)
vpos=np.where(m0,lat+rgC.standard_normal(NN)*np.nanstd(lat[m0]),np.nan)
vneg=np.where(m0,rgC.standard_normal(NN),np.nan)
fp,dp=same_construct(vpos); fn,dn=same_construct(vneg)
print(f"② 判据的对照(`#374b`:一个从未开火过的判据,它的每一个「不是」都是沉默):")
print(f"   正对照(真的第三个平行指标):开火 **{fp}** · r_a {dp.get('ra',0):+.3f} · r_w {dp.get('rw',0):+.3f} · "
      f"α₃ {dp.get('a3',0):.4f} vs α₂ {A2:.4f}")
print(f"   负对照(纯噪声):开火 **{fn}** · r_a {dn.get('ra',0):+.3f} · r_w {dn.get('rw',0):+.3f} · "
      f"α₃ {dn.get('a3',0):.4f}")

rows=[]
for r_,c,v in cor[:25]:
    f_,dd=same_construct(v)
    rows.append(dict(v_col=c[:70],v_absr=r_,v_ra=dd.get('ra',np.nan),v_rw=dd.get('rw',np.nan),
                     v_a3=dd.get('a3',np.nan),v_same=bool(f_)))
T=pd.DataFrame(rows); check_columns(T,'R419')
T.to_csv(pathlib.Path(__file__).parent/'results'/'candidates.csv',index=False)
HIT=T[T.v_same]
print(f"\n   前 25 个候选里,满足「同一构念」判据的:**{len(HIT)}**")
for r in HIT.itertuples():
    print(f"      ✅ {r.v_col[:60]:<62} r_a {r.v_ra:+.3f} · r_w {r.v_rw:+.3f} · α₃ **{r.v_a3:.4f}** (α₂ {A2:.4f})")

# ---- ③ ⚠ 页面上有一处**算术**自相矛盾,发现于本轮,不需要新数据 ----
# `#316` 写:「Spearman–Brown 可靠性只有 **0.380**(`animated`↔`written` = **+0.468**)」。
# 但 SB(r) = 2r/(1+r) 是一个恒等式:
inv=lambda a: a/(2-a)
print(f"\n③ ⚠ 页面算术自检(不需要新数据):")
print(f"   若 r = +0.468 -> SB = **{2*0.468/(1+0.468):.4f}**(页面写 0.380)")
print(f"   若 SB = 0.380 -> r = **{inv(0.380):+.4f}**(页面写 +0.468)")
print(f"   **两个数不能都对。** 而 2 × {inv(0.380):.4f} = **{2*inv(0.380):.4f}** ≈ 0.468 ——")
print(f"   **这正是把 SB 的分子 `2r` 当成了 `r`。**")
for lbl,g_ in (('本轮掩码(去六坐标后仍有限)',m0),('仅两题都有作答',np.isfinite(AN)&np.isfinite(WR)),
               ('两题 + ok',ok&np.isfinite(AN)&np.isfinite(WR))):
    rr=float(np.corrcoef(AN[g_],WR[g_])[0,1])
    print(f"   实测 [{lbl}] n={int(g_.sum()):,} · r = **{rr:+.4f}** · SB = **{2*rr/(1+rr):.4f}**")
print(f"   **⇒ `animated`↔`written` ≈ +0.23,不是 +0.468;而 0.380 那个可靠性是对的。**")

g=Gate('问卷里还有没有第三个指向 form 的题')
g.asserted('★ 正对照:判据必须在**真的第三个平行指标**上开火',fp,
           f"α₃ {dp.get('a3',0):.4f} > α₂ {A2:.4f}",kind='control')
g.asserted('★ 负对照:判据必须在纯噪声上不开火',not fn,f"开火 {fn}",kind='control')
g.asserted('★ offset 零非退化(数百个候选取最大值,天然为正)',offmax.std()>0,
           f"{offmax.mean():.4f} ± {offmax.std():.4f}",kind='control')
if fp and (not fn):
    g.asserted('★ 注册的 kill:存在至少一个同构念的第三指标',len(HIT)>0,
               f"命中 {len(HIT)} 个:{HIT.v_col.tolist()[:3]}")
else:
    g.asserted('★ 注册的 kill(判据未通过对照 -> 不判)',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

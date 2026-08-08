import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A108 R361 -- `form_i` 为什么是唯一一个够不到羞耻的维度

这一整段(`#308`–`#315`)都在问「什么东西碰羞耻」。**反过来问一次更省力。**
`form`(`animated` 与 `written` 共有的那个东西)是唯一一个**不通向羞耻**的常设维度。

⚠ **必须用未正交化的 `form`** —— `form_i` 与六坐标按构造正交,任何与它们的相关**代数上为零**
(`#277b` 的恒等式陷阱)。羞耻不是六坐标之一,但 `form_i` 的正交化会把与羞耻**经由六坐标的那部分**
也抹掉,所以两个都报。

**三个互斥的世界,预测不同:**
| | 去衰减后 ↔ 羞耻 | ↔ 消费/媒介类结局 | ↔ 答题风格 |
|---|---|---|---|
| **Ⓐ 它太弱**(信度 0.34) | **出现** | 中等 | ~0 |
| **Ⓑ 它是媒介不是内容** | **仍 ~0** | **强** | ~0 |
| **Ⓒ 它是测量伪影** | ~0 | ~0(扣掉风格后) | **强** |

ESTIMAND        ① `corr(form, 羞耻)` 的**去衰减**值(用 `animated`↔`written` 的分半信度);
                ② `form` 与消费/媒介类结局(`pornhabit` · 付费 · 开始年龄 · 暴力比例)的相关;
                ③ `form` 与**答题风格**(跨 Likert 的均值 · sd · 作答数)的相关。
KILL            **三条预测互斥,一轮可判。若三个都不成立 -> 第四个世界,而那要另写。**
POSITIVE CTRL   合成一个**已知只碰媒介类结局**的量 -> ② 的判据必须抓到、① 与 ③ 必须落零。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ guard 21     若结论是「它确实不碰羞耻」,必须交出置换分位数 · MDE · 正对照灵敏度。
IMPOSSIBLE      羞耻只有**一道题**,信度不可估 -> 去衰减只用 `form` 那一侧,
                所以报的是**真相关的上界**,不是点估计。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR)
ani=pd.to_numeric(d['animated'],errors='coerce').values.astype(float)
wri=pd.to_numeric(d['written'],errors='coerce').values.astype(float)
m0=np.isfinite(ani)&np.isfinite(wri)&ok
zz=lambda v,m:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
FORM=np.full(NN,np.nan); FORM[m0]=(zz(ani,m0)+zz(wri,m0))/2
base=ok.copy()
for q_ in Q: base&=np.isfinite(q_)
X6=np.column_stack([np.ones(int(base.sum()))]+[zz(q_,base) for q_ in Q])
def resid_on6(v):
    y=v[base]; f=np.isfinite(y); out=np.full(NN,np.nan); idx=np.flatnonzero(base)
    b=np.linalg.lstsq(X6[f],(y[f]-y[f].mean())/y[f].std(),rcond=None)[0]
    out[idx[f]]=(y[f]-y[f].mean())/y[f].std()-X6[f]@b; return out
FORM_I=np.full(NN,np.nan)
ra,rw=resid_on6(ani),resid_on6(wri)
g=np.isfinite(ra)&np.isfinite(rw); FORM_I[g]=(ra[g]+rw[g])/2
def cor(u,v,m=None):
    k=np.isfinite(u)&np.isfinite(v)&(ok if m is None else m)
    return (float(np.corrcoef(u[k],v[k])[0,1]),int(k.sum())) if k.sum()>200 else (np.nan,0)
r_aw,_=cor(ani,wri); rel=2*r_aw/(1+r_aw)
r_sh,n_sh=cor(FORM,sh); r_shi,_=cor(FORM_I,sh)
dis=r_sh/np.sqrt(max(rel,1e-9))
mde=2.8/np.sqrt(max(n_sh,1))
print(f"`corr(animated, written)` = **{r_aw:+.4f}** -> Spearman–Brown 信度 **{rel:.4f}**")
print(f"① ↔羞耻:`form`(未正交化)**{r_sh:+.4f}**(n={n_sh:,})· "
      f"`form_i`(正交化)**{r_shi:+.4f}** · **去衰减上界 {dis:+.4f}** · MDE **{mde:.4f}**")
MEDIA={'pornhabit':'pornhabit','付费':'Have you paid for porn/erotic content? (pnltks8)',
       '开始年龄':'At what age did you begin watching porn or reading erotic content at least semiregularly? (ugf1hyy)',
       '暴力比例':'How much of the porn or erotica you watch is violent? (joz3g52)'}
print(f"\n② ↔消费/媒介类结局:")
med=[]
for nm,c in MEDIA.items():
    v=pd.to_numeric(d[c],errors='coerce')
    if v.notna().sum()<1000:
        v=d[c].astype('category').cat.codes.astype(float).replace(-1,np.nan)
    r_,n_=cor(FORM,v.values.astype(float)); med.append(abs(r_) if np.isfinite(r_) else np.nan)
    print(f"   {nm:<10} **{r_:+.4f}**(n={n_:,})")
LK=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in
                    [c for c in d.columns if d[c].dtype!=object and
                     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]])
STY={'风格·均值':np.nanmean(LK,1),'风格·sd':np.nanstd(LK,1),'风格·作答数':np.isfinite(LK).sum(1).astype(float)}
print(f"\n③ ↔答题风格({LK.shape[1]} 道 Likert 上建的三个人层量):")
sty=[]
for nm,v in STY.items():
    r_,_=cor(FORM,v); sty.append(abs(r_)); print(f"   {nm:<10} **{r_:+.4f}**")
print(f"\n判决表:① 去衰减 {abs(dis):.4f}(MDE {mde:.4f})· ② max|r| **{np.nanmax(med):.4f}** · "
      f"③ max|r| **{max(sty):.4f}**")
# ⚠ ③ 的 0.1824 与 ① 的去衰减 0.1293 同量级 -> **答题风格是个活的混淆,必须扣掉再看**。
def partial(u,y,ctrl):
    m=np.isfinite(u)&np.isfinite(y)&ok
    for c in ctrl: m&=np.isfinite(c)
    X=np.column_stack([np.ones(m.sum())]+[zz(c,m) for c in ctrl])
    ru=zz(u,m)-X@np.linalg.lstsq(X,zz(u,m),rcond=None)[0]
    ry=zz(y,m)-X@np.linalg.lstsq(X,zz(y,m),rcond=None)[0]
    return float(np.corrcoef(ru,ry)[0,1]),int(m.sum())
CT=list(STY.values())
r_p,n_p=partial(FORM,sh,CT); dis_p=r_p/np.sqrt(max(rel,1e-9))
print(f"\n⚠ 扣掉答题风格后:`form` ↔ 羞耻 **{r_p:+.4f}**(原 {r_sh:+.4f})· "
      f"**去衰减 {dis_p:+.4f}**(原 {dis:+.4f})· 保留 **{100*abs(dis_p)/max(abs(dis),1e-9):.1f}%**")
r_pi,_=partial(FORM_I,sh,CT)
print(f"   `form_i`(已扣六坐标)再扣风格:**{r_pi:+.4f}**")

def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[cor(perm_finite(FORM,700+i),sh)[0] for i in range(20)]
q=float(np.mean([abs(x)>=abs(r_sh) for x in nul]))
print(f"负对照(打乱人)↔羞耻:{np.mean(nul):+.4f} ± {np.std(nul):.4f} · |零| ≥ |观测| 的比例 **{q:.3f}**")
rg=np.random.default_rng(9)
pv=pd.to_numeric(d[MEDIA['pornhabit']],errors='coerce')
if pv.notna().sum()<1000: pv=d[MEDIA['pornhabit']].astype('category').cat.codes.astype(float).replace(-1,np.nan)
pv=pv.values.astype(float); mp=np.isfinite(pv)&ok
SYN=np.full(NN,np.nan); SYN[mp]=zz(pv,mp)+rg.standard_normal(int(mp.sum()))
s1,_=cor(SYN,sh); s2=max(abs(cor(SYN,pv)[0]),0); s3=max(abs(cor(SYN,v)[0]) for v in STY.values())
print(f"\n正对照(只碰媒介类的合成量):① ↔羞耻 **{s1:+.4f}** · ② ↔pornhabit **{s2:.4f}** · ③ ↔风格 **{s3:.4f}**")
T=pd.DataFrame([dict(v_test='①去衰减↔羞耻',v_val=dis),dict(v_test='②媒介 max|r|',v_val=float(np.nanmax(med))),
                dict(v_test='③风格 max|r|',v_val=float(max(sty))),dict(v_test='MDE',v_val=mde)])
check_columns(T,'R361'); T.to_csv(pathlib.Path(__file__).parent/'results'/'three_worlds.csv',index=False)
gg=Gate('`form` 为什么够不到羞耻')
gg.asserted('★ 正对照:只碰媒介类的合成量 -> ② 抓到、① 与 ③ 落零',
            s2>0.30 and abs(s1)<0.06 and s3<0.10,
            f"① {s1:+.4f} · ② {s2:.4f} · ③ {s3:.4f}")
gg.negative_control('★ 负对照:打乱人后 `form` ↔ 羞耻',float(np.mean(nul)),r_sh,
    null_spread=float(np.std(nul)),null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ Ⓐ 它太弱:去衰减后羞耻相关是否出现(> 2×MDE)',abs(dis)>2*mde,
            f"去衰减 {dis:+.4f} vs 2×MDE {2*mde:.4f}(信度 {rel:.4f})")
gg.asserted('★ Ⓑ 它是媒介不是内容:② 是否明显强于 ①',np.nanmax(med)>2*abs(dis),
            f"② {np.nanmax(med):.4f} vs 2×|①| {2*abs(dis):.4f}")
gg.asserted('★ Ⓒ 它是测量伪影:③ 是否强',max(sty)>0.20,f"③ {max(sty):.4f}")
gg.asserted('★ 关键控制:扣掉答题风格后 Ⓐ 还成不成立(去衰减 > 2×MDE)',
            abs(dis_p)>2*mde,
            f"扣风格后去衰减 **{dis_p:+.4f}**(原 {dis:+.4f},保留 {100*abs(dis_p)/max(abs(dis),1e-9):.1f}%)"
            f" vs 2×MDE {2*mde:.4f}")
gg.null_claim_uses_null_criteria('★ guard 21:若判「它确实不碰羞耻」','NULL',
    perm_quantile=q,mde=mde,sensitivity_shown=f"合成媒介量 ② {s2:.2f} 抓到",meaningful=0.05)
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

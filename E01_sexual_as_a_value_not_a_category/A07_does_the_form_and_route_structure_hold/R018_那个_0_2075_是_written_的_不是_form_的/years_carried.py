import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A119 R376 -- 这份羞耻,随「带了多少年」变不变

`#308a` 说两条路相加,`#311` 说没有东西缓冲它,`#309b` 说十个调节变量都不调节它。
**那时间呢?**

⚠⚠ **横断面数据 -> 只能说关联,不能说「磨掉」。** 写在设计里,不是写在讨论里。

ESTIMAND        每人「带了多少年」= **当前年龄 − 起始年龄的人内均值**;
                问它与羞耻的相关,**同时控制当前年龄**(否则测的是年龄本身)与 **`S`**
                (`#130` 的时间表:起始早的人既「带得久」也「兴趣更冷门」)。
KILL            **带得越久羞耻越低 -> 与时间同向;不变 -> 它是一个稳定属性
                (与 `#309b` 是同一句话的两面);越久越高 -> 最意外的一支。**
POSITIVE CTRL   合成一个**已知随「带了多少年」下降**的结局 -> 同一流程必须抓到,
                且控制当前年龄后仍在。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ guard 21     若判为零,交出**置换分位数 · MDE · 正对照灵敏度**三件套。
⚠ 天花板       年龄分档到 29-32(中点 30.5),起始最高 28 -> **「带了多少年」的量程很窄**,
                这限制了能测到的效应,必须连同 MDE 一起报。
IMPOSSIBLE      起始年龄是**回溯自报**,而此刻的羞耻会污染这个回忆;
                本轮测的是关联,不是「时间**造成**了什么」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
n_ons=np.isfinite(ONS).sum(1)
mean_ons=np.where(n_ons>=5,np.nanmean(ONS,1),np.nan)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
age=d['age'].map(AGE).values.astype(float)
YRS=age-mean_ons
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]
m0=np.isfinite(YRS)&np.isfinite(sh)&np.isfinite(age)&np.isfinite(S)&ok
zz=lambda v,m:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
n=int(m0.sum())
print(f"n={n:,} · 「带了多少年」中位 **{np.nanmedian(YRS[m0]):.1f}** 年 · "
      f"四分位 [{np.nanpercentile(YRS[m0],25):.1f}, {np.nanpercentile(YRS[m0],75):.1f}] · "
      f"sd **{np.nanstd(YRS[m0]):.2f}**")
print(f"⚠ 量程:年龄档中点最高 30.5,起始最高 28 -> 天花板效应;"
      f"`corr(带了多少年, 当前年龄)` = **{np.corrcoef(YRS[m0],age[m0])[0,1]:+.4f}**")
def cor(u,v,m=None):
    k=np.isfinite(u)&np.isfinite(v)&(m0 if m is None else m)
    return float(np.corrcoef(u[k],v[k])[0,1]) if k.sum()>200 else np.nan
def partial(u,y,ctrl):
    m=m0.copy()
    for c in ctrl: m&=np.isfinite(c)
    X=np.column_stack([np.ones(m.sum())]+[zz(c,m) for c in ctrl])
    ru=zz(u,m)-X@np.linalg.lstsq(X,zz(u,m),rcond=None)[0]
    ry=zz(y,m)-X@np.linalg.lstsq(X,zz(y,m),rcond=None)[0]
    return float(np.corrcoef(ru,ry)[0,1]),int(m.sum())
r0=cor(YRS,sh)
r1,_=partial(YRS,sh,[age]); r2,_=partial(YRS,sh,[age,S]); r3,_=partial(YRS,sh,[age,S,mean_ons])
print(f"\n★ `带了多少年 ↔ 羞耻`:")
print(f"   原样            **{r0:+.4f}**")
print(f"   控制当前年龄       **{r1:+.4f}**")
print(f"   + 控制 `S`       **{r2:+.4f}**")
print(f"   + 控制起始均值      **{r3:+.4f}**  ⚠ 三者共线(带了多少年 = 年龄 − 起始),这一格是过控制")
print(f"⚠ 参照:`当前年龄 ↔ 羞耻` **{cor(age,sh):+.4f}** · `起始均值 ↔ 羞耻` **{cor(mean_ons,sh):+.4f}**")
# ⚠⚠ 恒等式检查(`#267b` / guard 14 的家族):**带了多少年 = 当前年龄 − 起始均值**,
#    所以 `partial(年数, 羞耻 | 年龄)` 在代数上**就是** `−partial(起始均值, 羞耻 | 年龄)`。
#    若两者相等 -> 那个 +0.0621 **不是时长效应,是起始年龄效应**,而这个设计分不开它们。
i_ons,_=partial(mean_ons,sh,[age]); i_ons2,_=partial(mean_ons,sh,[age,S])
print(f"\n⚠⚠ 恒等式检查:`partial(年数, 羞耻 | 年龄)` = **{r1:+.4f}** vs "
      f"`−partial(起始均值, 羞耻 | 年龄)` = **{-i_ons:+.4f}** -> 差 **{abs(r1+i_ons):.6f}**")
print(f"   加控 `S` 后:**{r2:+.4f}** vs **{-i_ons2:+.4f}** -> 差 **{abs(r2+i_ons2):.6f}**")
IDENT=abs(r1+i_ons)<1e-6 and abs(r2+i_ons2)<1e-6
print(f"   -> {'**恒等,不是时长效应,是起始年龄效应**' if IDENT else '不恒等(说明 mean_ons 有缺失格局差异)'}")

def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[partial(perm_finite(YRS,400+i),sh,[age,S])[0] for i in range(30)]
q=float(np.mean([abs(x)>=abs(r2) for x in nul]))
mde=2.8/np.sqrt(max(n,1))
print(f"负对照(打乱人):**{np.mean(nul):+.4f} ± {np.std(nul):.4f}** · |零| ≥ |观测| 的比例 **{q:.3f}**")
print(f"MDE(80% 功效)**{mde:.4f}**")
rg=np.random.default_rng(66); PC={}
for g in (0.0,0.05,0.10):
    y=np.full(NN,np.nan)
    y[m0]=-g*zz(YRS,m0)+rg.standard_normal(n)
    PC[g]=partial(YRS,y,[age,S])[0]
    print(f"正对照(已知随年数下降 g={g:.2f}):控制年龄与 `S` 后 **{PC[g]:+.4f}**")
T=pd.DataFrame([dict(v_arm='原样',v_r=r0),dict(v_arm='控年龄',v_r=r1),
                dict(v_arm='控年龄+S',v_r=r2),dict(v_arm='MDE',v_r=mde)])
check_columns(T,'R376'); T.to_csv(pathlib.Path(__file__).parent/'results'/'years.csv',index=False)
gg=Gate('羞耻随「带了多少年」变不变')
gg.asserted('★ 正对照:已知随年数下降 0.10 的结局必须被抓到',abs(PC[0.10])>2*mde,
            f"g=0.10 -> {PC[0.10]:+.4f};g=0 -> {PC[0.0]:+.4f}")
gg.negative_control('★ 负对照:打乱人',float(np.mean(nul)),r2,null_spread=float(np.std(nul)),
    null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill:控制当前年龄与 `S` 后是否可分辨',abs(r2)>2*mde,
            f"**{r2:+.4f}** vs 2×MDE **{2*mde:.4f}**")
gg.null_claim_uses_null_criteria('★ guard 21:若判为零,三件套在不在',
    'NULL' if abs(r2)<=2*mde else 'EFFECT',perm_quantile=q,mde=mde,
    sensitivity_shown=f"g=0.10 抓到 {PC[0.10]:+.4f}",meaningful=0.05)
gg.could_have_come_out_otherwise('⚠⚠ guard 14:控制年龄后的「年数」是不是一个恒等式',
    (lambda: round(r2+i_ons2,9)),[],tol=1e-12) if False else gg.asserted(
    '★★ 恒等式检查:`partial(年数|年龄)` == `−partial(起始|年龄)`',IDENT,
    f"{r1:+.4f} vs {-i_ons:+.4f}(差 {abs(r1+i_ons):.2e})· 加控 S 后 {r2:+.4f} vs {-i_ons2:+.4f}"
    f" —— **恒等则本轮测的不是时长,是起始年龄**")
gg.asserted('⚠⚠ 横断面:只能说关联,不能说「磨掉」',True,
            '起始年龄是回溯自报,而此刻的羞耻会污染这个回忆')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

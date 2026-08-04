import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A131 R400 -- 六个坐标带着多少「报告广度」,遍及 29 个结局

`#355b`:类别数单独能解释 `animated` 的 2.648%,而它的**独有**增量只有 1.137pp ——
**剩下约 1.51pp 已经藏在六个坐标里。而对羞耻是 0.0pp。**
**那对其余 27 个结局呢?**

ESTIMAND        对全部 29 个结局各算:`R²(只类别数)` · `独有增量` = R²(六+n) − R²(六) ·
                **「广度渗漏」= 只类别数 R² − 独有增量**(= 已被六坐标吸收的那部分);
                报**分布**,以及渗漏与该结局**联合 R²** 的相关。
KILL            **若渗漏普遍且与联合 R² 正相关 -> 「六个坐标解释了 1.4%」里有一部分是报告广度,
                要写进页面口径注;若只有 `animated` 一个 -> 那是那一题的性质。**
POSITIVE CTRL   合成一个只由类别数驱动的结局 -> 渗漏必须接近它的全部 R²。
NEGATIVE CTRL   合成一个与类别数**正交**的结局 -> 渗漏必须 ≈ 0。
⚠ 多重性       29 个结局 -> **报分布,不挑最大的那个**(`#309c` 的陷阱)。
⚠ 窄口径       绝对量比较(`CALIBER.md` ⑩)。
IMPOSSIBLE      「渗漏」是一个**分解**的记账项,不是一个因果量;它只说「这部分共线」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A89_where_is_the_non_invariance/R333_gender_referential_split/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def curve(rows')[0])

inv=pd.read_csv('data/derived/inventory.csv')
BINo={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
      '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
onsc=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BINo).notna().sum()>300]
ncat=np.column_stack([np.isfinite(d[c].map(BINo).values.astype(float)) for c in onsc]).sum(1).astype(float)
ALLR=np.flatnonzero(ok); CO=coords(ALLR); CO=[CO[0],CO[1],CO[2],CO[3],-CO[4],CO[5]]
base=ok.copy()
for q_ in CO: base&=np.isfinite(q_)
rows=[]
for nm,y in OUT:
    m=base&np.isfinite(y)&np.isfinite(ncat)
    if m.sum()<300: continue
    n=int(m.sum()); z=lambda v:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
    Z=[z(q_) for q_ in CO]; zn=z(ncat); yy=z(y)
    def r2(cols):
        X=np.column_stack([np.ones(n)]+cols); b,*_=np.linalg.lstsq(X,yy,rcond=None)
        r=yy-X@b; return 1-float(r@r)/float(((yy-yy.mean())**2).sum())
    only=r2([zn]); six=r2(Z); both=r2(Z+[zn])
    rows.append(dict(v_out=str(nm)[:40],only=100*only,uniq=100*(both-six),
                     leak=100*(only-(both-six)),six=100*six,n=n))
T=pd.DataFrame(rows); check_columns(T,'R400')
T.to_csv(pathlib.Path(__file__).parent/'results'/'leak.csv',index=False)
print(f"{len(T)} 个结局 · 窄口径 · ⚠ 报分布,不挑最大的那个")
print(f"\n**广度渗漏**(= 只类别数 R² − 独有增量,单位 pp):")
print(f"   中位 **{T.leak.median():.3f}** · 均值 **{T.leak.mean():.3f}** · "
      f"范围 [{T.leak.min():.3f}, {T.leak.max():.3f}]")
print(f"   为正的 **{int((T.leak>0).sum())}/{len(T)}** · > 0.5pp 的 **{int((T.leak>0.5).sum())}**")
print(f"   相对该结局的六坐标 R²:中位 **{(100*T.leak/T.six.clip(lower=1e-9)).median():.1f}%**")
r=float(np.corrcoef(T.leak,T.six)[0,1])
print(f"\n★ `corr(渗漏, 六坐标联合 R²)` = **{r:+.4f}** —— "
      f"{'解释力越大的结局,渗漏越多' if r>0.3 else '与联合 R² 无关'}")
print(f"\n渗漏最大的四个(⚠ **列出来是为了说明分布的形状,不是挑它们**):")
for _,x in T.nlargest(4,'leak').iterrows():
    print(f"   {x.leak:>6.3f}pp(占其六坐标 R² 的 {100*x.leak/max(x.six,1e-9):>5.1f}%)  {x.v_out}")
rg=np.random.default_rng(88)
m=base&np.isfinite(ncat); n=int(m.sum()); z=lambda v:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
Z=[z(q_) for q_ in CO]; zn=z(ncat)
def leak_of(yv):
    yy=(yv-yv.mean())/yv.std()
    def r2(cols):
        X=np.column_stack([np.ones(n)]+cols); b,*_=np.linalg.lstsq(X,yy,rcond=None)
        r=yy-X@b; return 1-float(r@r)/float(((yy-yy.mean())**2).sum())
    only=r2([zn]); six=r2(Z); both=r2(Z+[zn]); return 100*only,100*(only-(both-six))
# ⚠ 这个统计量是标准的**共同性**(commonality = 只A + 只B − 合并),它**可以为负**(抑制)。
#    第一版的两个对照都是我设计错了:
#    ① 正对照的门槛(渗漏 > 只类别数 R² 的 50%)问的是一个**不由它决定**的比例 ——
#       种入越强,独有增量越大,共同性占比越小;**正确的问法是共同性 > 0 且明显。**
#    ② 负对照把结局建在**与六坐标正交的 ncat 残差**上 -> **必然**产生抑制(共同性为负),
#       那检的是我的构造,不是统计量。**正确的负对照:结局与 `ncat` 本身无关。**
o1,l1=leak_of(0.4*zn+rg.standard_normal(n))
print(f"\n正对照(只由类别数驱动):只类别数 R² **{o1:.3f}%** · **共同性 {l1:.3f}pp**")
NG=[]
for t in range(20):
    NG.append(leak_of(rg.standard_normal(n))[1])
NG=np.array(NG)
o2,l2=float('nan'),float(NG.mean())
print(f"负对照(**与 `ncat` 无关**的纯噪声结局,20 次):共同性 **{NG.mean():+.4f} ± {NG.std():.4f}pp**")
gg=Gate('六个坐标带着多少报告广度')
gg.asserted('★ 正对照(改正后):只由类别数驱动 -> 共同性必须明显为正(> 1pp)',l1>1.0,
            f"共同性 {l1:.3f}pp(只类别数 R² {o1:.3f}%)—— "
            f"⚠ 第一版的门槛「占比 > 50%」问的是一个不由它决定的比例,是我设计错了")
gg.asserted('★ 负对照(改正后):与 `ncat` **无关**的纯噪声结局 -> 共同性 ≈ 0',
            abs(NG.mean())<3*max(NG.std(),1e-9) and abs(NG.mean())<0.2,
            f"{NG.mean():+.4f} ± {NG.std():.4f}pp —— "
            f"⚠ 第一版把结局建在与六坐标正交的 ncat 残差上,那必然产生抑制,检的是我的构造")
gg.asserted('★ 注册的 kill:渗漏是否普遍(> 半数为正)且与联合 R² 正相关',
            (T.leak>0).mean()>0.5 and r>0.3,
            f"为正 {int((T.leak>0).sum())}/{len(T)} · 中位 {T.leak.median():.3f}pp · "
            f"`corr(渗漏, 六坐标 R²)` {r:+.4f}")
gg.asserted('⚠ 多重性:报分布,不挑最大的那个',True,'29 个结局 —— 上面列出的四个是形状,不是结论')
gg.asserted('⚠ 边界:「渗漏」是分解的记账项,不是因果量',True,'它只说「这部分共线」')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

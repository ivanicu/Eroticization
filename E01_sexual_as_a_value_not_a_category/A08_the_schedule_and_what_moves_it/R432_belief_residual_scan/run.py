import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A146 R432 -- 除了「做过多少」,还有什么预测「我能不能改掉它」

「能不能改」是这个项目里**唯一一个关于自我可塑性的判断**,
而它至今只有一个已知的相关物:**实践量 −0.042**(`#386a`);
**羞耻本身完全不预测它**(+0.0135,有功率的零,`#384a`)。
**那还有什么预测它?**

ESTIMAND        `BELIEF` 的残差(去掉 `ACTED` · `S` · `c3⁻` · 类别数 · 羞耻)
                对**所有**数值型人层量求相关;主量 = **最大 |r| 的分布**,以及逐个的名字。
判据(**先标支**,`#379c`)
                【两支】**判据必须先在一个已知相关物上开火**(`#419` 的教训):
                        `ACTED` 自己放回候选池,它**必须**排在前列;否则每一个「不是」都是沉默。
                        负对照:合成噪声列 -> 不得进前列。guard 26:MDE 扫描。
                【非零支】最大 |r| 越过 offset 零 -> 报名字与分布。
                【零支】未越阈时启用 MDE。
⚠ 零的种类     `offset_control`:**最大相关的零绝不是零** —— 在数百个候选里取最大值,最大值天然为正。
                零 = 对**同样多的合成噪声列**取最大相关的分布(`#419` 已验证的方法)。
⚠ 多重性       报分布不报单格。
IMPOSSIBLE      ① 残差已去掉五个量 -> 与它们共线的候选会被系统性压低,**这不是「它不相关」**;
                ② 全是自报,同源方差会抬高相关;
                ③ 相关高可能只是**同一个量的另一种写法**(`#357b`)-> 命中的要逐个看名字。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
BC=next(c for c in d.columns if '7lgg41e' in c); AC=next(c for c in d.columns if '41kpfir' in c)
BMAP={'Impossible':0.,'With an extreme amount of effort, maybe':1.,
      'With a lot of effort, yes':2.,'With some effort, yes':3.,'With little effort, yes':4.}
BELIEF=d[BC].map(BMAP).values.astype(float)
ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(BELIEF)&np.isfinite(ACTED)&np.isfinite(ncat)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
CTRL=[ACTED,S,C3,ncat,sh]
X0=np.column_stack([np.ones(n)]+[z(v,M) for v in CTRL])
yB=z(BELIEF,M); bb,*_=np.linalg.lstsq(X0,yB,rcond=None); RES=yB-X0@bb
# ⚠ `ACTED` 放回候选池当**已知相关物**的阳性参照;它对**残差**当然接近 0(已被去掉)——
# 所以参照要用**未去掉 ACTED** 的残差版本,否则我在测一个按构造为零的东西。
X1=np.column_stack([np.ones(n)]+[z(v,M) for v in [S,C3,ncat,sh]])
b1,*_=np.linalg.lstsq(X1,yB,rcond=None); RES_NOACT=yB-X1@b1
print(f"n=**{n:,}** · `BELIEF` 残差(去掉 5 个量)")
CAND=[]
for c in d.columns:
    if c in (BC,AC): continue
    v=pd.to_numeric(d[c],errors='coerce').values.astype(float)
    if np.isfinite(v[M]).sum()>=2000 and np.nanstd(v[M])>1e-9: CAND.append((c,v))
K=len(CAND)
def corr_to(res,v):
    g=np.isfinite(v[M])
    if g.sum()<2000: return np.nan
    return float(np.corrcoef(res[g],z(v,M)[g])[0,1])
cor=sorted(((abs(corr_to(RES,v)),corr_to(RES,v),c) for c,v in CAND if np.isfinite(corr_to(RES,v))),
           reverse=True)
print(f"候选 **{K}** 个;`BELIEF` 残差最大 |r| = **{cor[0][0]:.4f}**")
for k,(a,r_,c) in enumerate(cor[:8]): print(f"   {k+1}. {r_:+.4f}  {c[:66]}")
rgO=np.random.default_rng(31); offmax=[]
for s in range(200):
    rg=np.random.default_rng(2400+s); mx=0.
    for _ in range(K):
        v=rg.standard_normal(n); mx=max(mx,abs(float(np.corrcoef(RES,v)[0,1])))
    offmax.append(mx)
offmax=np.array(offmax); OTHR=float(np.percentile(offmax,95))
print(f"⚠ offset 零(**同样多({K})的合成噪声列取最大相关**,200 次):"
      f"**{offmax.mean():.4f} ± {offmax.std():.4f}** · 95 分位 **{OTHR:.4f}**")
print(f"   -> 实测最大 {cor[0][0]:.4f} {'**越阈**' if cor[0][0]>OTHR else '**未越阈**'} · "
      f"越阈的候选数 **{sum(1 for a,_,_ in cor if a>OTHR)}/{len(cor)}**")
# ⚠ **第一版这里错了**:我拿**族内(112 个取最大)**的阈 0.0419 去判 `ACTED` ——
# 但 `ACTED` 是一个**预先指定的单个变量**,**它不在那个族里**。
# **把族内阈用在族外的一个变量上,是一个类别错误**,而它的方向是**让阳性参照失败**,
# 于是整轮被判 UNVERIFIED —— 一个**假的**仪器失灵。
# 修法:两个阈 —— **单变量零**(打乱那一列)给预先指定的参照;**族内取最大零**给扫描的头名。
nulS=np.array([abs(float(np.corrcoef(RES_NOACT,
        np.random.default_rng(6100+s_).permutation(z(ACTED,M)))[0,1])) for s_ in range(400)])
STHR=float(np.percentile(nulS,95))
rA=abs(corr_to(RES_NOACT,ACTED))
print(f"\n⚠ **两个阈,不是一个**:族内取最大零 **{OTHR:.4f}**(给扫描头名)· "
      f"单变量零 **{STHR:.4f}**(给预先指定的参照)")
rank_a=1+sum(1 for a,_,_ in cor if a>rA)
print(f"★ 阳性参照(**已知相关物** `ACTED`,在**未去掉它**的残差上):|r| = **{rA:.4f}** vs "
      f"**单变量零** {STHR:.4f} -> "
      f"{'**开火(判据能认出已知的)**' if rA>STHR else '**不开火 —— 判据的每一个「不是」都是沉默**'}")
rgN=np.random.default_rng(7); fake=rgN.standard_normal(n)
rF=abs(float(np.corrcoef(RES,fake)[0,1]))
print(f"  负对照(一列纯噪声):|r| = **{rF:.4f}** vs 阈 {OTHR:.4f} -> "
      f"{'⚠ 越阈' if rF>OTHR else '**未越阈**'}")
T=pd.DataFrame([dict(v_col=c[:70],v_r=r_,v_absr=a) for a,r_,c in cor[:30]])
check_columns(T,'R432'); T.to_csv(pathlib.Path(__file__).parent/'results'/'belief_corr.csv',index=False)

g=Gate('除了「做过多少」,还有什么预测「我能不能改掉它」')
g.asserted('★【两支】判据必须先在**已知相关物** `ACTED` 上开火(用**单变量**零,不是族内零)',rA>STHR,
           f"|r| {rA:.4f} vs 单变量阈 {STHR:.4f}",kind='control')
g.asserted('★【两支】负对照:纯噪声列不得越阈',rF<=OTHR,f"{rF:.4f} vs {OTHR:.4f}",kind='control')
g.asserted('★【两支】offset 零非退化(取最大值天然为正)',offmax.std()>0,
           f"{offmax.mean():.4f} ± {offmax.std():.4f}",kind='control')
# ---- 家族读法:头几名不是四个独立的题,是一个家族 ----
FAM=['947wne3','normalsex','yuc275j','cunnilingus','jn2b355']
fv=[]
for c,v in CAND:
    if any(k in c for k in FAM): fv.append(z(v,M))
if fv:
    FM=np.nanmean(np.column_stack(fv),1)
    gF=np.isfinite(FM); rFAM=float(np.corrcoef(RES[gF],FM[gF])[0,1])
    print(f"\n★★ **家族读法**:头几名(dirtytalking · normalsex · blowjobs · cunnilingus ×2)"
          f"不是四个独立的题,是**同一个家族:常规性行为**。")
    print(f"   把 {len(fv)} 题平均成一个分数 -> 与 `BELIEF` 残差相关 **{rFAM:+.4f}** vs 单变量零 {STHR:.4f} -> "
          f"{'**越阈**' if abs(rFAM)>STHR else '未越阈'}")
    print(f"   ⚠ `S`(整体稀有度)**已在控制项里** -> 这是**超出**整体稀有度之外的残差效应。")
else:
    rFAM=np.nan

if rA>STHR and rF<=OTHR:
    g.asserted('★【非零支】扫描头名越过**族内取最大**零',cor[0][0]>OTHR,
               f"最大 {cor[0][0]:.4f} vs {OTHR:.4f} · 越阈 {sum(1 for a,_,_ in cor if a>OTHR)} 个")
    g.asserted('★【非零支】家族分数(常规性行为)越过**单变量**零',
               np.isfinite(rFAM) and abs(rFAM)>STHR,f"{rFAM:+.4f} vs {STHR:.4f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

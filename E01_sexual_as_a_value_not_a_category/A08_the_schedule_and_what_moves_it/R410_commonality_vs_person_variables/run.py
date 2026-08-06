import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A136 R410 -- `c3⁻` 的那份羞耻,能不能被任何**非坐标**的人层量解释

页面上羞耻有三条路:`c3⁻` 0.949pp · `S` 0.491pp · `EARLY` 0.484pp。
`#308a` 说 `c3⁻` 与 `S` **相加**;`#334` 说 `c3⁻` 与 `EARLY` 也**相加**。
**`c3⁻` 是三条里唯一与另外两条都不重叠的。** 那它自己那份,是不是被别的什么东西带着?

⚠ **先写下一件必须分清的事**:`#309b` 已经证明这些变量**不调节**羞耻的两条路。
**「不调节」与「不共线」是两件事** —— 一个变量可以完全不改变斜率,却把整条关系吃掉。
本轮问的是后者。

ESTIMAND        对每一个非坐标人层量 V,在**羞耻**上做共同性分解:
                `唯一(c3⁻|V) = R²(c3⁻,V) − R²(V)` · `共享 = R²(c3⁻) + R²(V) − R²(c3⁻,V)`
                主量 = **保留率 = 唯一(c3⁻|V) / R²(c3⁻)**;并报**全 panel 联合**的保留率。
KILL(条件式)  仅当正/负对照都过 -> 判:**是否有任何 V 的共享成分越过它自己的置换零的族内阈**。
                有 -> 那一个必须写进 `c3⁻` 的读法;没有 -> `c3⁻` 的那份羞耻是**独有**的。
POSITIVE CTRL   合成 V = c3⁻ + 噪声(真共线)-> 共享必须大、保留率必须塌。
NEGATIVE CTRL   合成 V = 纯噪声 -> 共享 ≈ 0、保留率 ≈ 1。
⚠ 零的种类     `offset_control`:**共享成分的零不该是零** —— 任何两个变量在有限样本上都共享一点。
                所以阈值取**该 V 自己的置换零的分位**(打乱 V 的人,保住缺失格局),不取 0。
⚠ 多重性       panel 里每个 V 一次 -> 报分布 + 族内(Bonferroni)分位。
⚠ 联合的零     打乱**整行 panel**(保住 V 之间的相关与维数)-> 零吸收同样的过拟合。
IMPOSSIBLE      共同性是记账项,不定因果方向;而「没有一个 V 吃掉它」不等于「不存在这样的量」——
                只等于**这个问卷里没有**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]

def num(c):  return pd.to_numeric(d[c],errors='coerce').values.astype(float)
def mp(c,m):
    v=d[c].map(m).values.astype(float)
    if np.isfinite(v).sum()<300:  # 有些列是数字与档位混排
        v=np.where(np.isfinite(v),v,pd.to_numeric(d[c],errors='coerce').values.astype(float))
    return v
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
# ⚠ 编码假设,写在跑之前:`TotalMentalIllness` 与 `childhood_adversity` 在公开发布里
# **只标出阳性**('Any'),其余是 NaN。把 NaN 读成 0 会把「没有」与「没作答」并成一档。
# 保留它们(它们在心理上是这一批里最相关的两个),但**标注这条假设**,并在联合里报有/无它们两种。
ASSUMED_ZERO={'心理疾病(Any)','童年逆境(Any)'}
PANEL=[
 ('年龄',              mp('age',AGE)),
 ('生理男',            num('biomale')),
 ('开放性',            num('opennessvariable')),
 ('尽责性',            num('consciensiousnessvariable')),
 ('外向性',            num('extroversionvariable')),
 ('神经质',            num('neuroticismvariable')),
 ('宜人性',            num('agreeablenessvariable')),
 ('无力感',            num('powerlessnessvariable')),
 ('心理疾病(Any)',    np.where(d['TotalMentalIllness'].eq('Any').values,1.,0.)),
 ('童年逆境(Any)',    np.where(d['childhood_adversity'].eq('Any').values,1.,0.)),
 ('成长环境开放度',    mp('How "sexually liberated" was your upbringing? (fs700v2)',
                        {'Repressed':-1.,'Neutral':0.,'Liberated':1.})),
 ('关系风格(非单偶)', mp('Personally, your preferred relationship style is: (4jib23m)',
                        {'Monogamous':0.,'Not monogamous':1.})),
 ('童年性别容忍度',    mp('childhood_gender_tolerance',{'Tolerant':1.,'Medium':0.,'Intolerant':-1.})),
 ('政治(保守←→自由)',mp('politics',{'Liberal':-1.,'Moderate':0.,'Conservative':1.})),
 ('异性恋',            mp('straightness',{'Straight':1.,'Not straight':0.})),
 ('超重',              mp('bmi',{'Not overweight':0.,'Overweight+':1.})),
 ('性伴数',            mp('sexcount',{'0':0.,'1-2':1.5,'3-7':5.,'8-20':14.,'21+':25.})),
 ('勾选的类别总数',    num('totalfetishcategory')),
 ('色情消费习惯',      num('pornhabit')),
 ('作答诚实度',        mp('How honest were you when answering this survey? (g1vao1y)',
                        {'Totally honest':1.,'Mostly honest':0.})),
]
base=ok&np.isfinite(C3)&np.isfinite(sh)
print(f"羞耻列 = {SHAME[:60]}… · 基础 n={int(base.sum()):,}\n")
print(f"⚠ 先分清:`#309b` 证明的是这些变量**不调节**羞耻;本轮问的是**共线**,是两件事。\n")

def decomp(v,y,c3,m):
    """返回 (R2_c3, R2_v, R2_both, 唯一c3, 共享, 保留率, n)。v 可以是 (n,k) 矩阵。"""
    V=np.atleast_2d(v.T).T if v.ndim==1 else v
    mm=m&np.isfinite(y)&np.isfinite(c3)&np.isfinite(V).all(1)
    n=int(mm.sum())
    if n<300: return None
    yy=y[mm]; yy=(yy-yy.mean())/max(yy.std(),1e-12); sst=float(((yy-yy.mean())**2).sum())
    def r2(cols):
        X=np.column_stack([np.ones(n)]+cols); b,*_=np.linalg.lstsq(X,yy,rcond=None)
        r=yy-X@b; return 100*(1-float(r@r)/sst)
    zc=(c3[mm]-c3[mm].mean())/max(c3[mm].std(),1e-12)
    Vz=np.column_stack([(V[mm,j]-V[mm,j].mean())/max(V[mm,j].std(),1e-12) for j in range(V.shape[1])])
    a=r2([zc]); b_=r2([Vz]); ab=r2([zc,Vz])
    uniq=ab-b_; shar=a+b_-ab
    return a,b_,ab,uniq,shar,uniq/max(a,1e-12),n

def perm_finite(v,seed):
    z=v.copy(); j=np.flatnonzero(np.isfinite(z))
    z[j]=z[np.random.default_rng(seed).permutation(j)]; return z

NP_=200; K=len(PANEL); FW=100*(1-0.05/K)          # 族内 Bonferroni 分位
print(f"逐个 V(每个 V 的阈是**它自己**的置换零的 {FW:.2f} 分位,族内 Bonferroni,K={K}):")
rows=[]
for nm,v in PANEL:
    r=decomp(v,sh,C3,base)
    if r is None: print(f"   {nm:<16} n 不足,跳过"); continue
    a,b_,ab,uq,sh_,ret,n=r
    nul=np.array([decomp(perm_finite(v,1000+s),sh,C3,base)[4] for s in range(NP_)])
    thr=float(np.percentile(nul,FW)); hit=sh_>thr
    rows.append(dict(v_name=nm,v_r2c3=a,v_r2v=b_,v_uniq=uq,v_shar=sh_,v_ret=ret,
                     v_thr=thr,v_hit=bool(hit),v_n=n))
    print(f"   {nm:<16} R²(V) {b_:5.2f}pp · 共享 **{sh_:+.4f}pp** (阈 {thr:+.4f}) · "
          f"保留 **{ret:6.1%}** {'⚠ 越阈' if hit else ''}")
T=pd.DataFrame(rows); check_columns(T,'R410')
T.to_csv(pathlib.Path(__file__).parent/'results'/'percoord.csv',index=False)

R2C3=float(T.v_r2c3.median())
print(f"\n`c3⁻` 单独对羞耻的 R² 中位 = **{R2C3:.3f}pp**")
print(f"保留率:最低 **{T.v_ret.min():.1%}**({T.loc[T.v_ret.idxmin(),'v_name']})· "
      f"中位 {T.v_ret.median():.1%} · 越阈的 **{int(T.v_hit.sum())}/{len(T)}**")

# ---- 全 panel 联合 ----
KEEP=[(n,v) for n,v in PANEL if n in set(T.v_name) and np.nanstd(v)>1e-9]
DROP=[n for n,v in PANEL if (n,None) not in [(k,None) for k,_ in KEEP]]
print(f"\n⚠ 联合里**没有**放进去的(个体 n<300 或方差为 0,不静默):{DROP if DROP else '无'}")
NAMES=[n for n,_ in KEEP]; MAT=np.column_stack([v for _,v in KEEP]); K2=len(KEEP)
rj=decomp(MAT,sh,C3,base)
aJ,bJ,abJ,uqJ,shJ,retJ,nJ=rj
print(f"\n全 panel 联合(k={K2},n={nJ:,}):R²(panel) **{bJ:.3f}pp** · R²(两者) {abJ:.3f}pp")
print(f"   `c3⁻` 的唯一成分 **{uqJ:.3f}pp** / 单独 {aJ:.3f}pp -> **保留 {retJ:.1%}**")
nulJ=[]
for s in range(NP_):
    rg=np.random.default_rng(5000+s); mm=np.isfinite(MAT).all(1)&base&np.isfinite(sh)&np.isfinite(C3)
    idx=np.flatnonzero(mm); P=MAT.copy(); P[idx]=MAT[rg.permutation(idx)]   # ⚠ 整行打乱,保住 V 间相关与维数
    rr=decomp(P,sh,C3,base)
    if rr: nulJ.append(rr[4])
nulJ=np.array(nulJ); thrJ=float(np.percentile(nulJ,95))
print(f"   联合共享 **{shJ:+.3f}pp** vs 整行置换零 {nulJ.mean():+.3f}±{nulJ.std():.3f} · 95 分位 {thrJ:+.3f} · "
      f"{'⚠ 越阈' if shJ>thrJ else '**未越阈**'}")

# ---- 对照 ----
rg=np.random.default_rng(77); nn=NN
vpos=np.where(np.isfinite(C3),C3+rg.standard_normal(nn)*np.nanstd(C3),np.nan)
vneg=np.where(np.isfinite(C3),rg.standard_normal(nn),np.nan)
rp=decomp(vpos,sh,C3,base); rn=decomp(vneg,sh,C3,base)
print(f"\n正对照(V = c3⁻ + 同幅噪声,真共线):共享 **{rp[4]:+.4f}pp** · 保留 **{rp[5]:.1%}**")
print(f"负对照(V = 纯噪声):共享 **{rn[4]:+.4f}pp** · 保留 **{rn[5]:.1%}**")

g=Gate('c3⁻ 的那份羞耻能否被非坐标的人层量解释')
POS=(rp[4]>0.05)and(rp[5]<0.60); NEG=(abs(rn[4])<0.05)and(rn[5]>0.95)
g.asserted('★ 正对照:真共线的 V -> 共享大且保留率塌到 60% 以下',POS,
           f"共享 {rp[4]:+.4f}pp · 保留 {rp[5]:.1%}",kind='control')
g.asserted('★ 负对照:纯噪声 V -> 共享 ≈0 且保留 ≥95%',NEG,
           f"共享 {rn[4]:+.4f}pp · 保留 {rn[5]:.1%}",kind='control')
if POS and NEG:
    g.asserted('★ 注册的 kill:没有任何单个 V 的共享越过它自己置换零的族内阈',
               int(T.v_hit.sum())==0,
               f"越阈 {int(T.v_hit.sum())}/{len(T)} —— " +
               (', '.join(T[T.v_hit].v_name.tolist()) if T.v_hit.any() else '无'))
    g.asserted('★ 全 panel 联合的共享未越过整行置换零的 95 分位',shJ<=thrJ,
               f"{shJ:+.3f}pp vs 阈 {thrJ:+.3f}")
else:
    g.asserted('★ 注册的 kill(对照未过 -> UNVERIFIED,不判)',False,
               'UNVERIFIED:正/负对照未双双通过')
print(g)
print(f"\n最大的三口:" + ' · '.join(
    f"**{r.v_name}** 吃掉 {100*(1-r.v_ret):.1f}%" for r in
    T.nsmallest(3,'v_ret').itertuples()))
print(f"⚠ 刻度(**不是检验**,k 不同):一个**真是 c3⁻ 加同幅噪声**的变量保留 {rp[5]:.1%};"
      f"20 个普通人层量合起来保留 **{retJ:.1%}**")
print(f"\n方向(共享成分是**哪一种人**——共享大不说方向,所以直接读相关):")
for r in T.nsmallest(5,'v_ret').itertuples():
    v=dict(PANEL)[r.v_name]
    m=base&np.isfinite(v)&np.isfinite(sh)&np.isfinite(C3)
    rc=np.corrcoef(v[m],C3[m])[0,1]; rs=np.corrcoef(v[m],sh[m])[0,1]
    print(f"   {r.v_name:<16} ↔c3⁻ **{rc:+.4f}** · ↔羞耻 **{rs:+.4f}** · "
          f"吃掉 {100*(1-r.v_ret):4.1f}% · {'同号(真共线)' if rc*rs>0 else '⚠ 异号 -> 抑制,不是共线'}")
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

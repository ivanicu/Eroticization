import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A153 R446 -- 那 1,790 人是**一个类型**,还是**四条羞耻路径的混合**

`#401e③`:「存在这群人」不等于「这群人是一个类型」。**本轮问那一句。**

比较的是**两格**:**「都高」**(羞耻高 + 疗愈高)vs **「只羞耻」**(羞耻高 + 疗愈低)。
**两格的羞耻都高** -> 差别只可能来自「疗愈」那一维,而不是「羞耻多少」。

两个活着的世界:
**A 一个类型** -> 已有的人层量里**至少一个**把两格分开,且方向可读;
**B 混合** -> **没有一个**分得开 -> 「都高」只是**羞耻的人里恰好也觉得疗愈的那些**,
不是一群有自己特征的人。

ESTIMAND        八个量(`S` · `D` · `c1` · `c2` · `c3⁻` · 清晰度 · 五题分数 · `EARLY`)
                在两格之间的标准化均值差;主量 = **族内 max-|t|**。
判据(**先标支**,`#379c`)
                【两支】guard 24 先给四个特征向量坐标定向 · 负对照用**越阈率** · guard 26 用 **MDE 扫描**。
                【非零支】族内 max-|t| 越阈 -> 世界 A,报是哪个量、什么方向。
                【零支】未越阈 -> 世界 B,启用 MDE。
⚠ 零的种类     `offset_control`:**两格之间任何量的差的零绝不是零**(任意两组人都有差)->
                零 = **随机等大小分组**(格大小照旧,人打乱)后的族内 max-|t| 分布。
⚠ 多重性       8 个量 -> **族内 max-|t| 阈**(`#334` 的做法),**不逐条判**。
IMPOSSIBLE      ① 六坐标需要块覆盖 ≥8 -> 两格的 n 会缩小,**同轮报缩小前后的格大小**;
                ② 「没有一个量分得开」只说**这八个量**分不开,不说不存在这样的量;
                ③ 两格都在「羞耻高」这一侧 -> 本轮**完全不测**羞耻本身的差异。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR)
THC=next(c for c in d.columns if 'vmq8jqw' in str(c))
th=pd.to_numeric(d[THC],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
INT=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i]) for i in range(5)]),1)
EARLY=np.where(np.isfinite(O).sum(1)>0,np.nanmean(np.where(np.isfinite(O),O,np.nan),1),np.nan)
FULL=np.isfinite(sh)&np.isfinite(th)
ms=float(np.median(sh[FULL])); mt=float(np.median(th[FULL]))
BOTH=FULL&(sh>ms)&(th>mt); SHONLY=FULL&(sh>ms)&(th<=mt)
print(f"全样本四格(`#401b`):都高 **{int(BOTH.sum()):,}** · 只羞耻 **{int(SHONLY.sum()):,}**")
NAMES=['S 位置','D 块间对比','c1','c2','c3⁻','清晰度','五题(常规不色)','EARLY 平均起始']
XS=[S,Q[1],Q[2],Q[3],-Q[4],Q[5],INT,EARLY]
gA=Gate('特征向量坐标进模型前先定向')
for nm,v in (('c1',Q[2]),('c2',Q[3]),('c3⁻',-Q[4])):
    gA.eigenvector_is_anchored(f'★ `{nm}` 对着羞耻定向',v,sh,'羞耻')
print(gA)
M=ok.copy()
for v in XS: M&=np.isfinite(v)
B=BOTH&M; Sx=SHONLY&M
print(f"\n⚠ IMPOSSIBLE ①:六坐标要块覆盖 ≥8 -> 缩小后:都高 **{int(B.sum()):,}** · 只羞耻 **{int(Sx.sum()):,}**")
n=int((B|Sx).sum())
def tstats(gb,gs):
    out=[]
    for v in XS:
        a=v[gb]; b=v[gs]
        se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
        sd=np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
        out.append(((a.mean()-b.mean())/max(sd,1e-12),(a.mean()-b.mean())/max(se,1e-12)))
    return np.array(out)
TS=tstats(B,Sx)
print(f"\n八个量(标准化均值差 d · |t|):")
for i,nm in enumerate(NAMES):
    print(f"   {nm:<16} d **{TS[i,0]:+.4f}** · |t| **{abs(TS[i,1]):.3f}**")
NP_=400; idx=np.flatnonzero(B|Sx); nb=int(B.sum()); mt_=[]
for s_ in range(NP_):
    rg=np.random.default_rng(6200+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:nb]]=True; gs=np.zeros(NN,bool); gs[p[nb:]]=True
    mt_.append(float(np.max(np.abs(tstats(ga,gs)[:,1]))))
mt_=np.array(mt_); THR=float(np.percentile(mt_,95))
MAXT=float(np.max(np.abs(TS[:,1]))); WHO=NAMES[int(np.argmax(np.abs(TS[:,1])))]
print(f"\n⚠ offset 零(**随机等大小分组** {NP_} 次;**任意两组人都有差,所以零不是零**):")
print(f"   族内 max-|t| = **{mt_.mean():.3f} ± {mt_.std():.3f}** · 95 分位 **{THR:.3f}**")
print(f"   实测 max-|t| = **{MAXT:.3f}**({WHO})-> "
      f"{'**越阈 -> 世界 A**' if MAXT>THR else '**未越阈 -> 世界 B(混合)**'}")
negs=np.array([float(np.max(np.abs(tstats(
    (lambda p: (lambda g: g)(np.isin(np.arange(NN),p[:nb])))(np.random.default_rng(95000+s).permutation(idx)),
    (lambda p: np.isin(np.arange(NN),p[nb:]))(np.random.default_rng(95000+s).permutation(idx)))[:,1])))
    for s in range(200)])
rate=float((negs>THR).mean())
print(f"\n负对照(**越阈率**,随机分组 200 次):**{100*rate:.1f}%**(合格 1–12%)")
# ⚠ **第一版的 MDE 扫描是空的,而它让 guard 26 通过了。**
# 我在**真实的两格**上种偏移,而真实的 max-|t| 已经是 6.803 —— **种什么都 100% 检出**,
# 因为最大值由**五题**给出,与我种的 `c3⁻` 无关。**它量的不是灵敏度,是「6.803 > 2.752」。**
# (与 `#387c`/`#372c②` 同族:**MDE 扫描种的必须能改变主量**;这里种的改不动最大值。)
# 修法:**在随机分组的基线上种** —— 那里 max-|t| 本来落在零里,种下去才看得出检出率。
print(f"\nguard 26 = **MDE 扫描**(⚠ 在**随机分组的基线**上种,不在真实两格上),每级 30 次:")
MDE=None
for gg in (0.05,0.10,0.15,0.25):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(400+int(gg*100)*79+s_); p2=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p2[:nb]]=True; gs=np.zeros(NN,bool); gs[p2[nb:]]=True
        bak=XS[4]; v2=XS[4].copy(); v2[ga]=v2[ga]+gg*np.nanstd(XS[4][B|Sx]); XS[4]=v2
        t2=float(np.max(np.abs(tstats(ga,gs)[:,1]))); XS[4]=bak
        if t2>THR: hit+=1
    print(f"   偏移 **{gg:.2f} sd** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.30
print(f"   **MDE = {MDE_:.2f} sd** · 有意义的组间差 = **0.20 sd**(小效应)")

pd.DataFrame([dict(v_name=NAMES[i],v_d=TS[i,0],v_t=TS[i,1]) for i in range(8)]+
             [dict(v_name='_thr',v_d=np.nan,v_t=THR),dict(v_name='_mde',v_d=MDE_,v_t=np.nan),
              dict(v_name='_n',v_d=int(B.sum()),v_t=int(Sx.sum()))]).to_csv(
    pathlib.Path(__file__).parent/'results'/'type_or_mix.csv',index=False)
g=Gate('那 1,790 人是一个类型,还是四条路径的混合')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
# ⚠ **第三次**把「零的争议幅度」传给 guard 26,而结果是**开火**的(`#372c①` · `#384d` · 本轮)。
# **一个零的争议幅度是「有意义的效应量」;一个非零结果的争议幅度是「实测到的效应」** ——
# 前者问「我本来想看见多小的东西」,后者问「我的仪器在**实际找到的那个大小**上工作吗」。
# 本轮 max-|t| = 6.803 >> 2.752,是**非零**支,所以争议幅度 = **实测的 d = 0.3081**。
# ⇒ 并且这一次不是靠记住:`guard 26` 现在**要求显式声明 branch**(`#402a`,`#383a` 的做法:改接口)。
OBS_D=float(np.max(np.abs(TS[:,0])))
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs **实测**效应(非零支)',MDE_,OBS_D,True,
    what='MDE 扫描 80% 检出',branch='non_null')
g.asserted('★【两支】offset 零非退化(任意两组人都有差)',mt_.std()>0,
           f"{mt_.mean():.3f} ± {mt_.std():.3f}",kind='control')
if 0.01<=rate<=0.12:
    if MAXT>THR:
        g.asserted('★【非零支】族内 max-|t| 越阈 -> 世界 A(是一个类型)',True,
                   f"{MAXT:.3f} vs {THR:.3f} · 最大的是 **{WHO}**(d {TS[int(np.argmax(np.abs(TS[:,1]))),0]:+.4f})")
    else:
        g.asserted('★【零支】未越阈且 MDE ≤ 0.20 sd -> 世界 B(混合)',MDE_<=0.20,
                   f"max-|t| {MAXT:.3f} vs {THR:.3f} · MDE {MDE_:.2f} sd")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

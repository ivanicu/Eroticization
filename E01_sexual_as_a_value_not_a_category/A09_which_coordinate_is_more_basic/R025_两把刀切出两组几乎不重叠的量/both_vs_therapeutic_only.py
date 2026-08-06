import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A161 R459 -- 在觉得疗愈的人里面,羞耻的那一维留下了什么

`#402a` 比的是**「都高」vs「只羞耻」** —— 两格都在**羞耻高**那一侧,所以差别只来自**疗愈**那一维。
**⇒ 另一刀从没做过:「都高」vs「只疗愈」** —— 两格都在**疗愈高**那一侧,差别只来自**羞耻**那一维。

两个活着的世界:
**A 羞耻在疗愈者里有自己的签名** -> 八个量里至少一个把两格分开;
**B 没有** -> 在觉得疗愈的人里,羞耻**不被这八个量中的任何一个标记** ——
   而那与 `#402a`(羞耻者里疗愈**有**签名)放在一起,是一个**不对称**的发现。

ESTIMAND        八个量(`S`·`D`·`c1`·`c2`·`c3⁻`·清晰度·五题·`EARLY`)在两格间的标准化均值差;
                主量 = **族内 max-|t|**。
判据(**先标支**,`#379c`)
                【两支】guard 24 先给特征向量坐标定向 · 负对照用**越阈率** ·
                        guard 26 **显式传 branch**,**MDE 扫描在随机分组基线上种**(`#402b` 的修法)。
                【非零支】族内 max-|t| 越阈 -> 世界 A,报是哪个量、什么方向;
                【零支】未越阈 -> 世界 B,启用 MDE。
⚠ 零的种类     `offset_control`:**两格之间任何量的差的零绝不是零**(任意两组人都有差)->
                零 = **随机等大小分组**(格大小照旧)的族内 max-|t| 分布。
⚠ 多重性       8 个量 -> **族内 max-|t| 阈**,不逐条判。
IMPOSSIBLE      ① 六坐标需块覆盖 ≥8 -> 格会缩小(**同轮报缩小前后**);
                ② 两格都在「疗愈高」一侧 -> **完全不测**疗愈本身的差异;
                ③ 「没有一个量分得开」只说**这八个量**分不开。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR)
THC=next(c for c in d.columns if 'vmq8jqw' in str(c)); th=pd.to_numeric(d[THC],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
INT=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i]) for i in range(5)]),1)
EARLY=np.where(np.isfinite(O).sum(1)>0,np.nanmean(np.where(np.isfinite(O),O,np.nan),1),np.nan)
FULL=np.isfinite(sh)&np.isfinite(th)
ms=float(np.median(sh[FULL])); mt=float(np.median(th[FULL]))
BOTH=FULL&(sh>ms)&(th>mt); THONLY=FULL&(sh<=ms)&(th>mt)
print(f"全样本:都高 **{int(BOTH.sum()):,}** · **只疗愈 {int(THONLY.sum()):,}**"
      f"(`#401b` 的 20.22%)")
NAMES=['S 位置','D 块间对比','c1','c2','c3⁻','清晰度','五题(常规不色)','EARLY 平均起始']
XS=[S,Q[1],Q[2],Q[3],-Q[4],Q[5],INT,EARLY]
gA=Gate('特征向量坐标进模型前先定向')
for nm,v in (('c1',Q[2]),('c2',Q[3]),('c3⁻',-Q[4])):
    gA.eigenvector_is_anchored(f'★ `{nm}` 对着羞耻定向',v,sh,'羞耻')
print(gA)
M=ok.copy()
for v in XS: M&=np.isfinite(v)
B=BOTH&M; T2=THONLY&M
print(f"\n⚠ IMPOSSIBLE ①:六坐标要块覆盖 ≥8 -> 缩小后:都高 **{int(B.sum()):,}** · 只疗愈 **{int(T2.sum()):,}**")
def tstats(gb,gs):
    out=[]
    for v in XS:
        a=v[gb]; b=v[gs]
        se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
        sd=np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
        out.append(((a.mean()-b.mean())/max(sd,1e-12),(a.mean()-b.mean())/max(se,1e-12)))
    return np.array(out)
TS=tstats(B,T2)
print(f"\n八个量(标准化均值差 d · |t|):")
for i,nm in enumerate(NAMES):
    print(f"   {nm:<16} d **{TS[i,0]:+.4f}** · |t| **{abs(TS[i,1]):.3f}**")
NP_=400; idx=np.flatnonzero(B|T2); nb=int(B.sum()); mt_=[]
for s_ in range(NP_):
    rg=np.random.default_rng(8700+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:nb]]=True; gs=np.zeros(NN,bool); gs[p[nb:]]=True
    mt_.append(float(np.max(np.abs(tstats(ga,gs)[:,1]))))
mt_=np.array(mt_); THR=float(np.percentile(mt_,95))
MAXT=float(np.max(np.abs(TS[:,1]))); WHO=NAMES[int(np.argmax(np.abs(TS[:,1])))]
print(f"\n⚠ offset 零(**随机等大小分组** {NP_} 次;**任意两组人都有差,所以零不是零**):")
print(f"   族内 max-|t| **{mt_.mean():.3f} ± {mt_.std():.3f}** · 95 分位 **{THR:.3f}**")
print(f"   实测 **{MAXT:.3f}**({WHO})-> "
      f"{'**越阈 -> 世界 A**' if MAXT>THR else '**未越阈 -> 世界 B**'}")
negs=[]
for s_ in range(200):
    rg=np.random.default_rng(99800+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:nb]]=True; gs=np.zeros(NN,bool); gs[p[nb:]]=True
    negs.append(float(np.max(np.abs(tstats(ga,gs)[:,1]))))
negs=np.array(negs); rate=float((negs>THR).mean())
print(f"\n负对照(**越阈率**,随机分组 200 次):**{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(⚠ **在随机分组基线上种**,`#402b` 的修法),每级 30 次:")
MDE=None
for gg in (0.05,0.10,0.15,0.20,0.25):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(3+int(gg*100)*139+s_); p=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p[:nb]]=True; gs=np.zeros(NN,bool); gs[p[nb:]]=True
        bak=XS[6]; v2=XS[6].copy(); v2[ga]=v2[ga]+gg*np.nanstd(XS[6][B|T2]); XS[6]=v2
        t2=float(np.max(np.abs(tstats(ga,gs)[:,1]))); XS[6]=bak
        if t2>THR: hit+=1
    print(f"   偏移 **{gg:.2f} sd** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.30
NONNULL=MAXT>THR
OBS=float(np.max(np.abs(TS[:,0])))
CONT=OBS if NONNULL else 0.20
print(f"   **MDE = {MDE_:.2f} sd** · 争议幅度 **{CONT:.4f}**({'实测' if NONNULL else '有意义(0.20 sd 小效应)'})")
pd.DataFrame([dict(v_name=NAMES[i],v_d=TS[i,0],v_t=TS[i,1]) for i in range(8)]+
             [dict(v_name='_thr',v_d=np.nan,v_t=THR),dict(v_name='_n',v_d=int(B.sum()),v_t=int(T2.sum()))]
            ).to_csv(pathlib.Path(__file__).parent/'results'/'other_cut.csv',index=False)
g=Gate('在觉得疗愈的人里面,羞耻的那一维留下了什么')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,what='随机分组基线上种',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
g.asserted('★【两支】offset 零非退化(任意两组人都有差)',mt_.std()>0,
           f"{mt_.mean():.3f} ± {mt_.std():.3f}",kind='control')
if 0.01<=rate<=0.12:
    if NONNULL:
        g.asserted('★【非零支】族内 max-|t| 越阈 -> 世界 A(羞耻在疗愈者里有签名)',True,
                   f"{MAXT:.3f} vs {THR:.3f} · 最大的是 **{WHO}**(d {TS[int(np.argmax(np.abs(TS[:,1]))),0]:+.4f})")
    else:
        g.asserted('★【零支】未越阈且 MDE ≤ 0.20 sd -> 世界 B',MDE_<=0.20,
                   f"max-|t| {MAXT:.3f} vs {THR:.3f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

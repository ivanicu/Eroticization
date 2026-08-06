import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A163 R461 -- 「两样都不强」的那 2,566 人是什么样的人

`#416c` 关掉了「那群人是独立类型」这条路(两维相加)。
而四格里**最大的那一格从没被描述过**:**都低**(不羞耻,也不觉得疗愈)。

两个活着的世界:
**A 什么都淡** -> 八个量的差**同号**,且集中在「更普通、更晚、更不冷门」——
   那意味着「两样都不强」是一个**缺失**;
**B 有自己的方向** -> 差**不同号** -> 「两样都不强」是一种**状态**,不是缺失。

ESTIMAND        ① 八个量在「都低」vs**其余三格合起来**之间的标准化均值差;
                   主量 = **族内 max-|t|**(**连续**主量);
                ② **符号一致性**(多数符号的比例)—— **离散**主量,判 A / B。
判据(**先标支**,`#379c`)
                【两支】guard 24 先定向 · 负对照用**越阈率** ·
                        guard 26 **对两个主量各调一次**,**显式传 `main_quantity`**(`#407c`),
                        **基线上种 + 网格一开始就加密**(`#402b` / `#403b` 三次教训)。
                【非零支】max-|t| 越阈 -> 这一格确实与众不同,报是哪些量;
                          再按**符号一致性**判 A(同号)/ B(不同号)。
                【零支】max-|t| 未越阈 -> 这一格与其余三格在这八个量上**没有可见差别**。
⚠ 零的种类     `offset_control`:**任意一格 vs 其余的差的零绝不是零**(任意两组人都有差)->
                零 = **随机等大小分组**(组大小照旧)的族内 max-|t| / 符号一致性分布。
IMPOSSIBLE      ① 「其余三格」不是一个同质的组 -> 这一刀只答「都低是否不同」,不答「不同于谁」;
                ② 中位数分割粗;③ 六坐标需覆盖 ≥8 -> 格缩小(同轮报)。
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
NAMES=['S 位置','D 块间对比','c1','c2','c3⁻','清晰度','五题(常规不色)','EARLY 平均起始']
XS=[S,Q[1],Q[2],Q[3],-Q[4],Q[5],INT,EARLY]
gA=Gate('特征向量坐标进模型前先定向')
for nm,v in (('c1',Q[2]),('c2',Q[3]),('c3⁻',-Q[4])):
    gA.eigenvector_is_anchored(f'★ `{nm}` 对着羞耻定向',v,sh,'羞耻')
print(gA)
M=ok.copy()
for v in XS: M&=np.isfinite(v)
LOW=FULL&M&(sh<=ms)&(th<=mt); REST=FULL&M&~LOW
print(f"\n**都低 {int(LOW.sum()):,}** vs 其余三格合起来 **{int(REST.sum()):,}**")
def tstats(gb,gs):
    out=[]
    for v in XS:
        a=v[gb]; b=v[gs]
        se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
        sd=np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
        out.append(((a.mean()-b.mean())/max(sd,1e-12),(a.mean()-b.mean())/max(se,1e-12)))
    return np.array(out)
TS=tstats(LOW,REST)
print(f"\n八个量(「都低」减「其余」· d · |t|):")
for i,nm in enumerate(NAMES):
    print(f"   {nm:<16} d **{TS[i,0]:+.4f}** · |t| **{abs(TS[i,1]):.3f}**")
npos=int((TS[:,0]>0).sum()); CONS=max(npos,8-npos)/8
print(f"\n   符号:为正 **{npos}/8** · **一致性 {CONS:.3f}**")
NP_=400; idx=np.flatnonzero(LOW|REST); nl=int(LOW.sum()); mt_=[]; cs_=[]
for s_ in range(NP_):
    rg=np.random.default_rng(9500+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:nl]]=True; gs=np.zeros(NN,bool); gs[p[nl:]]=True
    t2=tstats(ga,gs); mt_.append(float(np.max(np.abs(t2[:,1]))))
    q=int((t2[:,0]>0).sum()); cs_.append(max(q,8-q)/8)
mt_=np.array(mt_); cs_=np.array(cs_)
THR=float(np.percentile(mt_,95)); CTHR=float(np.percentile(cs_,95))
MAXT=float(np.max(np.abs(TS[:,1]))); WHO=NAMES[int(np.argmax(np.abs(TS[:,1])))]
print(f"\n⚠ offset 零(**随机等大小分组** {NP_} 次;**任意两组人都有差,所以零不是零**):")
print(f"   族内 max-|t| **{mt_.mean():.3f} ± {mt_.std():.3f}** · 95 分位 **{THR:.3f}**")
print(f"   符号一致性 **{cs_.mean():.3f} ± {cs_.std():.3f}** · 95 分位 **{CTHR:.3f}**")
print(f"   -> max-|t| **{MAXT:.3f}**({WHO}){'**越阈**' if MAXT>THR else '未越阈'} · "
      f"一致性 **{CONS:.3f}** {'**越阈**' if CONS>CTHR else '未越阈'}")
negs=[];cneg=[]
for s_ in range(200):
    rg=np.random.default_rng(99880+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:nl]]=True; gs=np.zeros(NN,bool); gs[p[nl:]]=True
    t2=tstats(ga,gs); negs.append(float(np.max(np.abs(t2[:,1]))))
    q=int((t2[:,0]>0).sum()); cneg.append(max(q,8-q)/8)
negs=np.array(negs); cneg=np.array(cneg)
rate=float((negs>THR).mean()); crate=float((cneg>CTHR).mean())
print(f"\n负对照(**越阈率**):max-|t| **{100*rate:.1f}%** · 一致性 **{100*crate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(⚠ 基线上种 + 网格一开始就加密),每级 30 次:")
MDE=None; det=[]
for gg in (0.05,0.10,0.15,0.18,0.20,0.25):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(1+int(gg*1000)*151+s_); p=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p[:nl]]=True; gs=np.zeros(NN,bool); gs[p[nl:]]=True
        bak=XS[0]; v2=XS[0].copy(); v2[ga]=v2[ga]+gg*np.nanstd(XS[0][M]); XS[0]=v2
        t2=float(np.max(np.abs(tstats(ga,gs)[:,1]))); XS[0]=bak
        if t2>THR: hit+=1
    det.append(hit/30); print(f"   偏移 **{gg:.2f} sd** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.30
NONNULL=MAXT>THR; OBS=float(np.max(np.abs(TS[:,0]))); CONT=OBS if NONNULL else 0.20
print(f"   **MDE = {MDE_:.2f} sd** · 争议幅度 **{CONT:.4f}**")
pd.DataFrame([dict(v_name=NAMES[i],v_d=TS[i,0],v_t=TS[i,1]) for i in range(8)]+
             [dict(v_name='_thr',v_d=np.nan,v_t=THR),dict(v_name='_cons',v_d=CONS,v_t=CTHR),
              dict(v_name='_n',v_d=int(LOW.sum()),v_t=int(REST.sum()))]).to_csv(
    pathlib.Path(__file__).parent/'results'/'quiet_cell.csv',index=False)
g=Gate('「两样都不强」的那 2,566 人是什么样的人')
g.asserted('★【两支】负对照(max-|t|):**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.asserted('★【两支】负对照(一致性):**越阈率** ≈5%',0.01<=crate<=0.15,f"{100*crate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26(**连续**主量 max-|t|)',MDE_,CONT,True,what='基线上种,网格已加密',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26(**离散**主量 符号一致性,`#407c` 的接口)',MDE_,CONT,True,
    what='同一扫描',branch='non_null' if NONNULL else 'null',
    main_quantity='discrete_count',sweep_detection=det)
g.asserted('★【两支】offset 零非退化',mt_.std()>0,f"{mt_.mean():.3f} ± {mt_.std():.3f}",kind='control')
if 0.01<=rate<=0.12:
    if NONNULL:
        g.asserted('★【非零支】max-|t| 越阈 -> 这一格确实与众不同',True,
                   f"{MAXT:.3f} vs {THR:.3f} · 最大 **{WHO}**(d {TS[int(np.argmax(np.abs(TS[:,1]))),0]:+.4f})")
        g.asserted('★【非零支/另判】符号一致性越阈 -> 世界 A(什么都淡);未越 -> 世界 B(有自己的方向)',
                   CONS>CTHR,f"一致性 {CONS:.3f} vs 阈 {CTHR:.3f} · 为正 {npos}/8")
    else:
        g.asserted('★【零支】max-|t| 未越阈且 MDE ≤ 0.20 -> 这一格没有可见差别',MDE_<=0.20,
                   f"{MAXT:.3f} vs {THR:.3f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

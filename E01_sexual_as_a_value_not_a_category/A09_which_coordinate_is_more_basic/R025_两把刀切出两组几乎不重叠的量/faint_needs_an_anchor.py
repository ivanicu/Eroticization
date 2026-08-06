import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A163 R462 -- 「什么都淡」要有一个锚才测得了,而那个锚不能是羞耻,也不能是我

`#417b`:符号一致性测不了「淡」,因为**八个量没有共同朝向**。
NEXT 提的是「按页面已有结论定向」。⚠ **但那条路本身有两个坑,写在跑之前:**

**① 按**羞耻**定向是**循环的**。** 页面上这些量的朝向都是**对着羞耻**记的
(`S` +0.119 · `c3⁻` +0.129 · `EARLY` −0.102 · 五题 +0.066)——
**而「都低」那一格**就是按羞耻定义的**。按羞耻定向再问它是否在低羞耻那侧,是被构造保证的**(`#396b` 同类)。
**② 按**我的读法**定向,正是 `#417b` 刚批评的那个错。**

**⇒ 需要一个锚:既不是羞耻、也不是我读出来的,而且方向由构造固定。**
`Totalsexacts` = **认可的性行为计数** —— **越大越卷入,方向由「它是一个计数」给出**;
**它不参与四格的定义**,所以不循环。

ESTIMAND        把八个量各自按 `sign(corr(量, Totalsexacts))` 定向到「**更卷入**」;
                主量 = 「都低」减「其余」在**定向后**八个量上的差 —— **有多少个为负**(离散计数)。
判据(**先标支**,`#379c`)
                【两支】**锚的方向必须先自证**:`corr(Totalsexacts, 勾选类别总数) > 0`(两个计数应同向);
                        负对照用**越阈率**;guard 26 **传 `main_quantity='discrete_count'` + 扫描**。
                【非零支】定向后**一致地为负**(越过 offset 零)-> 世界 A(**什么都淡**);
                【零支】未越阈 -> 世界 B(**有自己的方向**),启用 MDE。
⚠ 零的种类     `offset_control`:**定向后一致性的零绝不是零** ——
                八个量本身彼此相关,而且都被同一个锚定向过,**所以即使无关也会偏向一致**。
                零 = **随机等大小分组**(组大小照旧)后**用同一套朝向**重算一致性的分布。
IMPOSSIBLE      ① `Totalsexacts` 是「卷入」的**一个**操作化,不是它的全部;
                ② 定向用的是**全样本**相关 -> 与格无关,不循环,但也**不保证**它是心理上正确的朝向;
                ③ 中位数分割粗;格从 2,566 缩到 1,662。
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
ANC=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
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
M&=np.isfinite(ANC)
LOW=FULL&M&(sh<=ms)&(th<=mt); REST=FULL&M&~LOW
# ---- 锚先自证 ----
ncat_=pd.to_numeric(d['totalfetishcategory'],errors='coerce').values.astype(float)
ga=M&np.isfinite(ncat_)
rAnc=float(np.corrcoef(ANC[ga],ncat_[ga])[0,1])
print(f"\n⚠ **锚先自证**:`corr(Totalsexacts, 勾选类别总数)` = **{rAnc:+.4f}** "
      f"-> {'**同向,锚可用**' if rAnc>0 else '**⚠ 反向 —— 锚不可用**'}")
print(f"   ⚠ 锚**不参与**四格的定义(四格只用羞耻与治疗性)-> **不循环**")
SGN=[]
print(f"\n按 `sign(corr(量, Totalsexacts))` 定向到「更卷入」:")
for i,nm in enumerate(NAMES):
    g=M&np.isfinite(XS[i])
    r=float(np.corrcoef(XS[i][g],ANC[g])[0,1]); SGN.append(1.0 if r>0 else -1.0)
    print(f"   {nm:<16} corr(量, 锚) **{r:+.4f}** -> 朝向 **{'+' if r>0 else '−'}**")
XO=[SGN[i]*XS[i] for i in range(8)]
def diffs(gb,gs):
    out=[]
    for v in XO:
        a=v[gb]; b=v[gs]
        sd=np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
        out.append((a.mean()-b.mean())/max(sd,1e-12))
    return np.array(out)
DD=diffs(LOW,REST); nneg=int((DD<0).sum())
print(f"\n定向后「都低 − 其余」的八个差:")
for i,nm in enumerate(NAMES): print(f"   {nm:<16} **{DD[i]:+.4f}**")
print(f"   **为负(= 更不卷入)的:{nneg}/8**")
NP_=1000; idx=np.flatnonzero(LOW|REST); nl=int(LOW.sum()); nul=[]
for s_ in range(NP_):
    rg=np.random.default_rng(9900+s_); p=rg.permutation(idx)
    ga2=np.zeros(NN,bool); ga2[p[:nl]]=True; gs=np.zeros(NN,bool); gs[p[nl:]]=True
    nul.append(int((diffs(ga2,gs)<0).sum()))
nul=np.array(nul); HI=float(np.percentile(nul,95))
print(f"\n⚠ offset 零(**随机等大小分组** {NP_} 次,**用同一套朝向**;"
      f"**八个量彼此相关且被同一个锚定向 -> 即使无关也偏向一致**):")
print(f"   为负个数 **{nul.mean():.2f} ± {nul.std():.2f}** · 95 分位 **{HI:.1f}**")
print(f"   实测 **{nneg}** -> {'**越阈 -> 世界 A(什么都淡)**' if nneg>HI else '**未越阈 -> 世界 B**'}")
negs=[]
for s_ in range(300):
    rg=np.random.default_rng(99920+s_); p=rg.permutation(idx)
    ga2=np.zeros(NN,bool); ga2[p[:nl]]=True; gs=np.zeros(NN,bool); gs[p[nl:]]=True
    negs.append(int((diffs(ga2,gs)<0).sum()))
negs=np.array(negs); rate=float((negs>HI).mean())
print(f"\n负对照(**越阈率**,300 次):**{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(⚠ **离散主量**,基线上种,网格加密),每级 30 次:")
MDE=None; det=[]
for gg in (0.05,0.10,0.15,0.20,0.30):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(4+int(gg*1000)*157+s_); p=rg.permutation(idx)
        ga2=np.zeros(NN,bool); ga2[p[:nl]]=True; gs=np.zeros(NN,bool); gs[p[nl:]]=True
        XB=list(XO)
        for i in range(8):
            v2=XO[i].copy(); v2[ga2]=v2[ga2]-gg*np.nanstd(XO[i][M]); XB[i]=v2
        a=[];  # 用带偏移的版本重算
        for v in XB:
            aa=v[ga2]; bb=v[gs]
            sd=np.sqrt((aa.var(ddof=1)+bb.var(ddof=1))/2); a.append((aa.mean()-bb.mean())/max(sd,1e-12))
        if int((np.array(a)<0).sum())>HI: hit+=1
    det.append(hit/30); print(f"   八个量一起下移 **{gg:.2f} sd** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.40
NONNULL=nneg>HI
CONT=0.20
print(f"   **MDE = {MDE_:.2f} sd** · 有意义的「整体下移」= **{CONT:.2f} sd**")
pd.DataFrame([dict(v_name=NAMES[i],v_sign=SGN[i],v_d=DD[i]) for i in range(8)]+
             [dict(v_name='_nneg',v_sign=np.nan,v_d=nneg),
              dict(v_name='_thr',v_sign=np.nan,v_d=HI)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'oriented.csv',index=False)
g=Gate('「什么都淡」——用一个非循环、非我读出的锚')
g.asserted('★【两支】锚先自证:两个计数应同向',rAnc>0,f"{rAnc:+.4f}",kind='control')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.15,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26(**离散**主量 + 扫描,`#407c`)',MDE_,CONT,True,
    what='八个量一起下移',branch='non_null' if NONNULL else 'null',
    main_quantity='discrete_count',sweep_detection=det)
g.asserted('★【两支】offset 零非退化(同一套朝向下,随机分组也偏向一致)',nul.std()>0,
           f"{nul.mean():.2f} ± {nul.std():.2f}",kind='control')
if rAnc>0 and 0.01<=rate<=0.15:
    g.asserted('★【非零支】定向后一致地为负 -> 世界 A(什么都淡)',NONNULL,
               f"为负 {nneg}/8 vs 阈 {HI:.1f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

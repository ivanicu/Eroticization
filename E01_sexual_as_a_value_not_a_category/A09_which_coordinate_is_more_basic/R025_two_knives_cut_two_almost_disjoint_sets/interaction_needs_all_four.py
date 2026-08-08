import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A162 R460 -- 第三把刀是空的;真正的交互检验需要**第四格**

`#415` 的 NEXT 提的是**第三把刀**(只羞耻 vs 只疗愈),
并写着「**若 `S`/`EARLY` 在这一刀上接近零,说明它们标记的是「都高」这一格(交互)**」。

⚠⚠ **跑之前先算代数(`#387b` 的做法,第四次)—— 两处都错:**

**① 第三把刀的点估计是**恒等式**。** 三把刀都是同三个格均值的差:
`刀1 = m_都高 − m_只羞耻` · `刀2 = m_都高 − m_只疗愈` · `刀3 = m_只羞耻 − m_只疗愈`
**=> `刀3 ≡ 刀2 − 刀1`,精确成立,不需要任何可加性假设。**
**⇒ 点估计携带零比特新信息**(`#396b` 的同一类)。

**② 而那个读法本身也是错的。** `刀3 ≈ 0` 只等于 `刀1 ≈ 刀2` ——
**那意味着这个量沿两个维度**移动得一样多**,即**可加且对称**,恰恰**不是**交互。**
**真正的交互 = 四个格均值**不可加**,而它需要**第四格(都低)** —— 本项目**从没用过它**。

**⇒ 本轮做两件事:① 用数据验证那个恒等式(便宜、决定性);② 做**真正的**交互检验(四格)。**

ESTIMAND        ① `刀3 − (刀2 − 刀1)` 的最大绝对值(应为 0);
                ② 八个量各自的**交互项** `(m_都高 − m_只羞耻) − (m_只疗愈 − m_都低)`
                   —— 即四格的**不可加性**;主量 = **族内 max-|t|**。
判据(**先标支**,`#379c`)
                【两支】① 恒等式必须精确成立(否则是我的代数错)· guard 24 先定向 ·
                        负对照用**越阈率** · guard 26 **显式 branch**,**基线上种**。
                【非零支】② 族内 max-|t| 越阈 -> 存在真正的交互,报是哪个量;
                【零支】未越阈 -> **两个维度是可加的**,而那是页面上一句新话。
⚠ 零的种类     `offset_control`:**四格交互项的零绝不是零**(四组人各有采样噪声)->
                零 = **把四格的标签随机重排(各格大小照旧)**后的族内 max-|t| 分布。
IMPOSSIBLE      ① 中位数分割 -> 四格是粗的;② 六坐标需覆盖 ≥8 -> 四格都会缩小(同轮报);
                ③ 「可加」不等于「两个维度独立」。
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
BOTH=FULL&M&(sh>ms)&(th>mt); SHO=FULL&M&(sh>ms)&(th<=mt)
THO=FULL&M&(sh<=ms)&(th>mt); LOW=FULL&M&(sh<=ms)&(th<=mt)
print(f"\n四格(坐标掩码后):都高 **{int(BOTH.sum()):,}** · 只羞耻 **{int(SHO.sum()):,}** · "
      f"只疗愈 **{int(THO.sum()):,}** · **都低 {int(LOW.sum()):,}**(← **从没用过的那一格**)")
def mean_(v,g): return float(np.nanmean(v[g]))
print(f"\n① **恒等式验证**(`刀3 ≡ 刀2 − 刀1`,精确,不需要可加性):")
gaps=[]
for i,nm in enumerate(NAMES):
    c1=mean_(XS[i],BOTH)-mean_(XS[i],SHO); c2=mean_(XS[i],BOTH)-mean_(XS[i],THO)
    c3=mean_(XS[i],SHO)-mean_(XS[i],THO)
    gaps.append(abs(c3-(c2-c1)))
print(f"   八个量里 `|刀3 − (刀2 − 刀1)|` 的**最大值** = **{max(gaps):.2e}**")
print(f"   **⇒ 恒等式成立到浮点精度。第三把刀的点估计携带零比特新信息。**")
print(f"\n② **真正的交互**(四格的不可加性):`(都高 − 只羞耻) − (只疗愈 − 都低)`")
def inter(v,gB,gS,gT,gL):
    a=v[gB];b=v[gS];c=v[gT];e=v[gL]
    est=(a.mean()-b.mean())-(c.mean()-e.mean())
    se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)+c.var(ddof=1)/len(c)+e.var(ddof=1)/len(e))
    sd=np.sqrt((a.var(ddof=1)+b.var(ddof=1)+c.var(ddof=1)+e.var(ddof=1))/4)
    return est/max(sd,1e-12),est/max(se,1e-12)
TS=np.array([inter(v,BOTH,SHO,THO,LOW) for v in XS])
for i,nm in enumerate(NAMES):
    print(f"   {nm:<16} 交互 d **{TS[i,0]:+.4f}** · |t| **{abs(TS[i,1]):.3f}**")
sizes=[int(BOTH.sum()),int(SHO.sum()),int(THO.sum()),int(LOW.sum())]
idx=np.flatnonzero(BOTH|SHO|THO|LOW); NP_=400; mt_=[]
for s_ in range(NP_):
    rg=np.random.default_rng(9100+s_); p=rg.permutation(idx); c=0; gs=[]
    for k in sizes:
        g=np.zeros(NN,bool); g[p[c:c+k]]=True; gs.append(g); c+=k
    mt_.append(float(np.max(np.abs(np.array([inter(v,*gs) for v in XS])[:,1]))))
mt_=np.array(mt_); THR=float(np.percentile(mt_,95))
MAXT=float(np.max(np.abs(TS[:,1]))); WHO=NAMES[int(np.argmax(np.abs(TS[:,1])))]
print(f"\n⚠ offset 零(**四格标签随机重排(各格大小照旧)** {NP_} 次;"
      f"**四组人各有采样噪声 -> 零不是零**):")
print(f"   族内 max-|t| **{mt_.mean():.3f} ± {mt_.std():.3f}** · 95 分位 **{THR:.3f}**")
print(f"   实测 **{MAXT:.3f}**({WHO})-> "
      f"{'**越阈 -> 存在真正的交互**' if MAXT>THR else '**未越阈 -> 两个维度是可加的**'}")
negs=[]
for s_ in range(200):
    rg=np.random.default_rng(99850+s_); p=rg.permutation(idx); c=0; gs=[]
    for k in sizes:
        g=np.zeros(NN,bool); g[p[c:c+k]]=True; gs.append(g); c+=k
    negs.append(float(np.max(np.abs(np.array([inter(v,*gs) for v in XS])[:,1]))))
negs=np.array(negs); rate=float((negs>THR).mean())
print(f"\n负对照(**越阈率**,随机重排 200 次):**{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(⚠ **在随机重排基线上种**,`#402b` 的修法),每级 30 次:")
MDE=None
# ⚠ 第一版网格 (0.05,0.10,0.15,0.20,0.25):0.20 给 60%、0.25 给 90%,**80% 点在中间** ->
# MDE 报成 0.25 > 有意义 0.20 -> 门 FAIL。**`#403b`/`#411b` 的规矩第三次生效**:
# 网格分辨率把 MDE 系统性报**高**,而**加密网格改的是测量的分辨率,判据(0.20)一个字没动**。
for gg in (0.05,0.10,0.15,0.18,0.20,0.22,0.25):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(2+int(gg*1000)*149+s_); p=rg.permutation(idx); c=0; gs=[]
        for k in sizes:
            g=np.zeros(NN,bool); g[p[c:c+k]]=True; gs.append(g); c+=k
        bak=XS[7]; v2=XS[7].copy(); v2[gs[0]]=v2[gs[0]]+gg*np.nanstd(XS[7][M]); XS[7]=v2
        t2=float(np.max(np.abs(np.array([inter(v,*gs) for v in XS])[:,1]))); XS[7]=bak
        if t2>THR: hit+=1
    print(f"   都高格偏移 **{gg:.2f} sd** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.30
NONNULL=MAXT>THR
OBS=float(np.max(np.abs(TS[:,0])))
CONT=OBS if NONNULL else 0.20
print(f"   **MDE = {MDE_:.2f} sd** · 争议幅度 **{CONT:.4f}**({'实测' if NONNULL else '有意义(0.20 sd)'})")
pd.DataFrame([dict(v_name=NAMES[i],v_d=TS[i,0],v_t=TS[i,1],v_idgap=gaps[i]) for i in range(8)]+
             [dict(v_name='_thr',v_d=np.nan,v_t=THR,v_idgap=np.nan),
              dict(v_name='_sizes',v_d=np.nan,v_t=np.nan,v_idgap=np.nan)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'interaction4.csv',index=False)
g=Gate('第三把刀是空的;真正的交互需要第四格')
g.asserted('★【两支】① 恒等式 `刀3 ≡ 刀2 − 刀1` 精确成立(否则是我的代数错)',max(gaps)<1e-9,
           f"最大偏差 {max(gaps):.2e}",kind='control')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,what='随机重排基线上种',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
g.asserted('★【两支】offset 零非退化(四组人各有采样噪声)',mt_.std()>0,
           f"{mt_.mean():.3f} ± {mt_.std():.3f}",kind='control')
if max(gaps)<1e-9 and 0.01<=rate<=0.12:
    if NONNULL:
        g.asserted('★【非零支】② 族内 max-|t| 越阈 -> 存在真正的交互',True,
                   f"{MAXT:.3f} vs {THR:.3f} · 最大的是 **{WHO}**(d {TS[int(np.argmax(np.abs(TS[:,1]))),0]:+.4f})")
    else:
        g.asserted('★【零支】未越阈且 MDE ≤ 0.20 sd -> **两个维度是可加的**',MDE_<=0.20,
                   f"max-|t| {MAXT:.3f} vs {THR:.3f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

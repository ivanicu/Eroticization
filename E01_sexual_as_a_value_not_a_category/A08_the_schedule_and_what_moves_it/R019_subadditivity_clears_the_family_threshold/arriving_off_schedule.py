import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A160 R458 -- 一个类别「比它该来的时候更早到」,有意义吗

⚠ **先按 `#396b` 查一条腿**:`corr(类别中位起始, 稀有度) = +0.3406` **已经在 `#370b` 估过** ——
**冷门的类别在人口层确实来得更晚。那条腿不是新信息。**

**真正没被问过的是**:把中位起始**对稀有度回归**之后的**残差** ——
**「这个类别比它的稀有度所预测的**更早**还是**更晚**到」** —— 有没有意义。

两个活着的世界:
**A 提前到达有代价** -> 残差与该类别的**平均羞耻**相关(**更早到 -> 更羞耻**,残差为负 -> 羞耻高);
**B 只有稀有度算数** -> 残差什么也不带,而那意味着**「发育时间表」是一条一维的线**,
   偏离它本身不携带信息 —— 那也是页面上一句新话。

⚠⚠ **`n_eff` = 31 个**类别**,不是 15,503 个人。** 这一条写在最前面,不写在脚注里。

ESTIMAND        31 个起始类别:中位起始年龄 · 稀有度 · 该类别报告者的平均羞耻;
                `off = 中位起始 − 稀有度回归的预测`;主量 = `corr(off, 平均羞耻)`。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;guard 26 **显式传 branch**,网格一开始就加密。
                【非零支】越阈 -> 按符号判(**负 = 早到的更羞耻**);
                【零支】未越阈**且 MDE < 0.35**(31 点上的中等相关)-> 世界 B。
⚠ 零的种类     `offset_control`:**这个零不该是零** ——
                「中位起始」与「平均羞耻」**都是在重叠的人群上聚合的**,
                所以即使两者无关,共享的人也可能造出基线关联。
                零 = **在人层打乱羞耻**(`lib.nulls.perm_in`),**重算每个类别的平均羞耻**,再算相关
                —— 这**精确保住**了类别之间的人群重叠结构。
IMPOSSIBLE      ① 31 点 -> 相关的展布很大,**任何结论都只能是粗的**;
                ② 平均羞耻按报告者算 -> 报告人数少的类别更不稳(**同轮报每类别 n**);
                ③ 「早到有代价」与「早到的东西本来就更让人羞耻」在本设计里分不开。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
SHC=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHC],errors='coerce').values.astype(float)
BINo={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
      '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
inv=pd.read_csv('data/derived/inventory.csv')
ONS=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BINo).notna().sum()>300]
O=np.column_stack([d[c].map(BINo).values.astype(float) for c in ONS])
K=O.shape[1]
print(f"⚠⚠ **`n_eff` = **{K}** 个类别,不是 {NN:,} 个人。**(写在最前面,不写在脚注里)\n")
prev=np.isfinite(O).mean(0); RAR=-np.log(np.clip(prev,1e-4,1.))
med=np.array([np.nanmedian(O[np.isfinite(O[:,j]),j]) for j in range(K)])
def cat_shame(shv):
    return np.array([np.nanmean(shv[np.isfinite(O[:,j])&np.isfinite(shv)]) for j in range(K)])
CS=cat_shame(sh)
nrep=np.isfinite(O).sum(0)
print(f"⚠ **`#392e`:三个量各自先看清楚**")
print(f"   中位起始 [{med.min():.1f}, {med.max():.1f}] 岁 · 稀有度 [{RAR.min():.2f}, {RAR.max():.2f}] · "
      f"平均羞耻 [{CS.min():.3f}, {CS.max():.3f}]")
print(f"   每类别报告人数 [{nrep.min():,}, {nrep.max():,}] · 中位 {int(np.median(nrep)):,}")
print(f"⚠ 已估过的那条腿(`#370b`):corr(中位起始, 稀有度) = **{np.corrcoef(med,RAR)[0,1]:+.4f}**\n")
b=np.polyfit(RAR,med,1); OFF=med-np.polyval(b,RAR)
R=float(np.corrcoef(OFF,CS)[0,1])
print(f"残差 `off`(**正 = 比它的稀有度所预测的更晚到**)· "
      f"**corr(off, 平均羞耻) = {R:+.4f}**(n_eff = {K})")
order=np.argsort(OFF)
print(f"\n   **最「提前到」的 4 个类别**(off 最负):")
for j in order[:4]:
    print(f"      {str(ONS[j])[:52]:<54} off {OFF[j]:+5.2f} 岁 · 平均羞耻 {CS[j]:.3f} · n={nrep[j]:,}")
print(f"   **最「延后到」的 4 个**(off 最正):")
for j in order[-4:]:
    print(f"      {str(ONS[j])[:52]:<54} off {OFF[j]:+5.2f} 岁 · 平均羞耻 {CS[j]:.3f} · n={nrep[j]:,}")
NP_=1000; Mall=np.isfinite(sh)
nul=np.array([float(np.corrcoef(OFF,cat_shame(perm_in(sh,Mall,8300+s)))[0,1]) for s in range(NP_)])
THR=float(np.percentile(np.abs(nul),95))
print(f"\n⚠ offset 零(**在人层打乱羞耻,重算每个类别的平均羞耻** {NP_} 次;"
      f"**两个量都在重叠人群上聚合 -> 这个零不该是零**):")
print(f"   **{nul.mean():+.4f} ± {nul.std():.4f}** · |值| 95 分位 **{THR:.4f}**")
print(f"   实测 **{R:+.4f}** -> **{(R-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(R)>THR else '**未越阈**'}")
negs=np.array([float(np.corrcoef(OFF,cat_shame(perm_in(sh,Mall,88300+s)))[0,1]) for s in range(300)])
rate=float((np.abs(negs)>THR).mean())
print(f"\n负对照(**越阈率**,300 次):**{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(⚠ 网格一开始就加密),每级 30 次:")
MDE=None
for gg in (0.15,0.25,0.30,0.35,0.45,0.60):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(5+int(gg*100)*137+s_)
        y=gg*(OFF-OFF.mean())/OFF.std()+np.sqrt(max(1-gg*gg,1e-9))*rg.standard_normal(K)
        if abs(float(np.corrcoef(OFF,y)[0,1]))>THR: hit+=1
    print(f"   种植相关 **{gg:.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.70
NONNULL=abs(R)>THR
CONT=abs(R) if NONNULL else 0.35
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 **{CONT:.3f}**({'实测' if NONNULL else '有意义(31 点上的中等相关)'})")
pd.DataFrame([dict(v_cat=str(ONS[j])[:60],v_med=med[j],v_rar=RAR[j],v_off=OFF[j],
                   v_shame=CS[j],v_n=int(nrep[j])) for j in range(K)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'off_schedule.csv',index=False)
g=Gate('一个类别「比它该来的时候更早到」,有意义吗')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,what='网格一开始就加密',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
g.asserted('★【两支】offset 零非退化(两个量都在重叠人群上聚合)',nul.std()>0,
           f"{nul.mean():+.4f} ± {nul.std():.4f}",kind='control')
if 0.01<=rate<=0.12:
    if NONNULL:
        g.asserted('★【非零支】越阈 -> 按符号判(**负 = 早到的更羞耻**)',True,
                   f"{R:+.4f} vs 阈 {THR:.4f} -> **{'早到的更羞耻' if R<0 else '晚到的更羞耻'}**")
    else:
        g.asserted('★【零支】未越阈**且 MDE < 0.35** -> 世界 B(偏离时间表不携带信息)',MDE_<0.35,
                   f"{R:+.4f} vs {THR:.4f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

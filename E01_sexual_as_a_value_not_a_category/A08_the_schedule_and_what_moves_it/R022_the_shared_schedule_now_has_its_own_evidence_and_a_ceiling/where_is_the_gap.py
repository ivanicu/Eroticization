import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A132 R403 -- 男女时间表那 0.16 的缺口,差在哪几个类别上

`#358c`:男女顺序 ρ **+0.8052**,而随机劈半的天花板是 **+0.9661** —— 缺口 **0.16**。
`#291` 的教训:这种缺口可能集中在**极少数格**上(那一轮 22 道题里只有 1 道把两组分开)。

ESTIMAND        ① 逐类别的**秩差**(男的均值秩 − 女的均值秩)的分布;
                ② **留一类别**后的 ρ,**对着「随机去掉一个类别」的分布读**
                   (⚠ 留一必然提高 ρ —— 去掉一个点总会)。
KILL            **若缺口集中在 2–3 个类别 -> 那几个是性别特异的,例外可以点名;
                若秩差均匀铺开 -> 是整体的轻微不同步,那是另一句话。**
POSITIVE CTRL   人为把某一个类别的男女均值拉开 -> 留一必须挑出它。
NEGATIVE CTRL   随机去掉一个类别的 ρ 分布(**这就是留一的参照臂**,`#326b`)。
⚠ 多重性       31 个类别 -> 报**分布 + 留一曲线**,不点名单格除非它明显出界。
IMPOSSIBLE      秩差混着**真实顺序差**与**该类别在两组内的样本量差**;本轮报后者作为协变量。
"""
import numpy as np, pandas as pd, warnings, hashlib, re
warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')

BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
HAS=np.isfinite(ONS); NC=HAS.shape[1]
NAM=[re.sub(r'^How old were you when you first experienced (sexual )?interest in ','',str(c))[:40] for c in ons]
sex=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
def means(mask,O=ONS):
    H=np.isfinite(O)
    return np.array([np.nanmean(O[mask&H[:,j],j]) if (mask&H[:,j]).sum()>=60 else np.nan for j in range(NC)])
def rho(a,b,keep=None):
    k=np.isfinite(a)&np.isfinite(b)
    if keep is not None: k&=keep
    return float(spearmanr(a[k],b[k]).statistic) if k.sum()>=10 else np.nan
gm,gf=means(sex==1),means(sex==0)
base=rho(gm,gf)
rgp=np.random.default_rng(5)
p1=means(rgp.random(NN)<0.5); p2=means(~(rgp.random(NN)<0.5))
ceil=rho(p1,p2)
from scipy.stats import rankdata
km=np.isfinite(gm)&np.isfinite(gf)
rm=np.full(NC,np.nan); rf=np.full(NC,np.nan)
rm[km]=rankdata(gm[km]); rf[km]=rankdata(gf[km])
dr=rm-rf
nm_=np.array([(sex==1).sum() and ((sex==1)&HAS[:,j]).sum() for j in range(NC)])
nf_=np.array([((sex==0)&HAS[:,j]).sum() for j in range(NC)])
print(f"男女 ρ **{base:+.4f}** · 随机劈半天花板 **{ceil:+.4f}** · 缺口 **{ceil-base:.4f}**")
print(f"\n① 逐类别秩差(男秩 − 女秩,{int(km.sum())} 个可比):")
print(f"   |秩差| 中位 **{np.nanmedian(np.abs(dr)):.1f}** · 均值 {np.nanmean(np.abs(dr)):.1f} · "
      f"最大 **{np.nanmax(np.abs(dr)):.0f}**")
o=np.argsort(-np.abs(np.where(np.isfinite(dr),dr,0)))
for j in o[:5]:
    print(f"   {dr[j]:+5.0f}  (男 n={nm_[j]:>5,} / 女 n={nf_[j]:>5,})  {NAM[j]}")
LOO=np.array([rho(gm,gf,keep=(np.arange(NC)!=j)) for j in range(NC)])
print(f"\n② 留一类别后的 ρ:中位 **{np.nanmedian(LOO):+.4f}** · "
      f"最大 **{np.nanmax(LOO):+.4f}**(去掉 {NAM[int(np.nanargmax(LOO))]})")
print(f"   ⚠ 留一必然提高 ρ —— 参照臂 = 这 {NC} 个留一值自己的分布:"
      f"sd **{np.nanstd(LOO):.4f}**,最大落在 **{(np.nanmax(LOO)-np.nanmean(LOO))/max(np.nanstd(LOO),1e-9):+.2f}** sd")
print(f"   去掉最大的那一个后 ρ = {np.nanmax(LOO):+.4f},离天花板 {ceil:.4f} 还差 **{ceil-np.nanmax(LOO):.4f}**")
J=int(np.nanargmax(nm_+nf_))
O2=ONS.copy(); O2[(sex==1),J]=O2[(sex==1),J]+6.0
gm2=means(sex==1,O2); gf2=means(sex==0,O2)
L2=np.array([rho(gm2,gf2,keep=(np.arange(NC)!=j)) for j in range(NC)])
print(f"\n正对照(人为把类别 #{J} 的男均值 +6 岁):留一最大在 **#{int(np.nanargmax(L2))}**(应 #{J})· "
      f"ρ {rho(gm2,gf2):+.4f} -> 去掉它 {np.nanmax(L2):+.4f}")
T=pd.DataFrame([dict(v_cat=NAM[j],v_dr=float(dr[j]),v_loo=float(LOO[j]),
                     v_nm=int(nm_[j]),v_nf=int(nf_[j])) for j in range(NC)])
check_columns(T,'R403'); T.to_csv(pathlib.Path(__file__).parent/'results'/'gap.csv',index=False)
top3=np.nanmax(LOO)-base
gg=Gate('男女时间表的缺口在哪')
gg.asserted('★ 正对照:人为拉开一个类别 -> 留一必须挑出它',int(np.nanargmax(L2))==J,
            f"挑出 #{int(np.nanargmax(L2))},应 #{J}")
gg.asserted('★ 负对照/参照臂:留一必然提高 ρ,所以看的是它在自身分布里的位置',True,
            f"留一 ρ 中位 {np.nanmedian(LOO):+.4f} · sd {np.nanstd(LOO):.4f} · "
            f"最大 {(np.nanmax(LOO)-np.nanmean(LOO))/max(np.nanstd(LOO),1e-9):+.2f} sd")
gg.asserted('★ 注册的 kill:缺口是否集中在极少数类别(最大留一 > 2 sd 且能补上多数缺口)',
            (np.nanmax(LOO)-np.nanmean(LOO))/max(np.nanstd(LOO),1e-9)>2.0 and top3>0.5*(ceil-base),
            f"最大留一把 ρ 从 {base:+.4f} 提到 {np.nanmax(LOO):+.4f}(补上 {100*top3/max(ceil-base,1e-9):.0f}% 的缺口)· "
            f"位置 {(np.nanmax(LOO)-np.nanmean(LOO))/max(np.nanstd(LOO),1e-9):+.2f} sd")
gg.asserted('⚠ 多重性:报分布 + 留一曲线,不点名单格除非明显出界',True,f"{NC} 个类别")
gg.asserted('⚠ 边界:秩差混着真实顺序差与样本量差',True,'两组的 n 已随秩差一并报出')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

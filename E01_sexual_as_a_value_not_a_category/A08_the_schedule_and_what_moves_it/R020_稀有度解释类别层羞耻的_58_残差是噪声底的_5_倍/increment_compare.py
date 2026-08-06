import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A130 R398 -- `animated` 为什么是最好解释的那一格

`#286a` 的联合 R² 曲线里 **`animated` 是最高的一格(7.79%)**,是羞耻(2.8%)的 **2.8 倍**,
**而这个项目从没问过为什么。**

⚠ **`animated` 是 `form` 的两个指标之一**(`#325`)-> **先报它与预测量的相关**,
并把「它与预测量共享内容」写进读法。
⚠ **用窄口径**:这是**绝对量**的比较,`CALIBER.md` ⑩ 的「⛔ 不可用」栏禁止宽口径。

ESTIMAND        对 `animated` 与 `羞耻` 各做**六坐标 + 类别数**的增量分解 ΔR²,**并排**;
                并报每个坐标与两个结局的**单相关**。
KILL            **若同一批坐标以**同样的比例**在起作用 -> `animated` 只是「更容易被预测」;
                若比例明显不同(例如 `c3⁻` 占压倒性)-> 六个坐标里有一个是**媒介/形式**的,
                而那会回头照亮 `#316`/`#325` 的 `form`。**
POSITIVE CTRL   合成两个结局:一个由 `S` 主导、一个由 `c3⁻` 主导 -> 分解必须分别指出。
NEGATIVE CTRL   `perm_finite` 打乱人 -> 所有增量必须落零。
IMPOSSIBLE      `animated` 与六坐标同源(都来自同一份问卷的兴趣题)——
                共享方法方差不可分离;本轮报的是**比例**,不是「解释了什么」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ani=pd.to_numeric(d['animated'],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
onsc=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in onsc])
ncat=np.isfinite(ONS).sum(1)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR)
NAMES=['S 位置','D 块间对比','c1','c2','c3⁻','清晰度']
CO=[Q[0],Q[1],Q[2],Q[3],-Q[4],Q[5]]
m0=ok.copy()
for q_ in CO: m0&=np.isfinite(q_)
m0&=np.isfinite(sh)&np.isfinite(ani)
n=int(m0.sum()); z=lambda v:(v[m0]-v[m0].mean())/max(v[m0].std(),1e-12)
Z=[z(q_) for q_ in CO]+[z(ncat.astype(float))]
LAB=NAMES+['(控制)类别数']
print(f"n={n:,} · 窄口径(`cov>=8`)—— 绝对量比较,`CALIBER.md` ⑩ 禁宽口径")
print(f"⚠ `animated` 是 `form` 的指标之一,先报它与预测量的单相关:")
for l,zz_ in zip(LAB,Z):
    print(f"   {l:<12} ↔`animated` **{np.corrcoef(zz_,z(ani))[0,1]:+.4f}** · "
          f"↔羞耻 **{np.corrcoef(zz_,z(sh))[0,1]:+.4f}**")
def dec(y):
    yy=z(y)
    def r2(cols):
        X=np.column_stack([np.ones(n)]+cols); b,*_=np.linalg.lstsq(X,yy,rcond=None)
        r=yy-X@b; return 1-float(r@r)/float(((yy-yy.mean())**2).sum())
    full=r2(Z)
    return full,[full-r2([c for j,c in enumerate(Z) if j!=i]) for i in range(len(Z))]
fA,iA=dec(ani); fS,iS=dec(sh)
print(f"\n{'':<12}{'animated':>14}{'羞耻':>14}{'比':>10}")
print(f"{'联合 R²':<12}{100*fA:>13.3f}%{100*fS:>13.3f}%{fA/max(fS,1e-9):>10.2f}×")
for i,l in enumerate(LAB):
    print(f"{l:<12}{100*iA[i]:>13.3f}%{100*iS[i]:>13.3f}%"
          f"{(iA[i]/max(iS[i],1e-9)):>10.2f}×")
shA=np.array(iA)/max(np.sum(np.abs(iA)),1e-12); shS=np.array(iS)/max(np.sum(np.abs(iS)),1e-12)
print(f"\n★ **份额**(各坐标增量占该结局总增量的比例):")
for i,l in enumerate(LAB):
    print(f"   {l:<12} animated **{100*shA[i]:>5.1f}%** · 羞耻 **{100*shS[i]:>5.1f}%**")
print(f"   两条份额向量的相关 **{np.corrcoef(shA,shS)[0,1]:+.4f}**"
      f"(高 = 同一批坐标以同样比例在起作用)")
rg=np.random.default_rng(51)
print(f"\n正对照:")
for tag,yv in (('只由 `S` 主导',0.3*Z[0]+rg.standard_normal(n)),
               ('只由 `c3⁻` 主导',0.3*Z[4]+rg.standard_normal(n))):
    y2=np.full(NN,np.nan); y2[m0]=yv
    f2,i2=dec(y2); k=int(np.argmax(i2))
    print(f"   {tag:<16} 最大增量在 **{LAB[k]}**({100*i2[k]:.3f}pp)")
    if 'S' in tag and 'c3' not in tag: pcS=k
    else: pcC=k
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=np.array([dec(perm_finite(ani,600+i))[1] for i in range(30)])
print(f"负对照(打乱人):各增量 **{100*nul.mean():.4f} ± {100*nul.std():.4f}pp**")
T=pd.DataFrame([dict(v_term=LAB[i],v_ani=100*iA[i],v_sh=100*iS[i],
                     v_shareA=100*shA[i],v_shareS=100*shS[i]) for i in range(len(LAB))])
check_columns(T,'R398'); T.to_csv(pathlib.Path(__file__).parent/'results'/'cmp.csv',index=False)
gg=Gate('`animated` 为什么最好解释')
gg.asserted('★ 正对照:两个合成结局必须各自指出主导坐标',pcS==0 and pcC==4,
            f"S 主导 -> {LAB[pcS]} · c3⁻ 主导 -> {LAB[pcC]}")
gg.asserted('★ 负对照:打乱人后各增量落零',abs(100*nul.mean())<0.05,
            f"{100*nul.mean():.4f} ± {100*nul.std():.4f}pp")
gg.asserted('★ 注册的 kill:两条**份额**向量是否高度一致(> 0.8 = 只是更容易被预测)',
            float(np.corrcoef(shA,shS)[0,1])>0.8,
            f"份额相关 **{np.corrcoef(shA,shS)[0,1]:+.4f}** · "
            f"联合 R² {100*fA:.3f}% vs {100*fS:.3f}%({fA/max(fS,1e-9):.2f}×)")
gg.asserted('⚠⚠ `animated` 与六坐标同源',True,
            '都来自同一份问卷的兴趣题 —— 共享方法方差不可分离;本轮报的是**比例**,不是「解释了什么」')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

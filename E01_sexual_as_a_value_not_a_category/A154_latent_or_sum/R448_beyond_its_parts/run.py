import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A154 R448 -- `BOTH` 带着成分之外的东西吗:潜变量,还是只是加法

`#403d③`:「合起来更强」**不等于**「存在一个潜变量」。
检验:把 `BOTH` 与它的两个成分一起放进**「治疗性」**(一个**没有参与构造它**的结局),
看 `BOTH` 控制两个成分后**是否还在**。

⚠⚠ **先验证不是恒等式(`#387b` 的做法,写在 NEXT 里、跑之前):**
`BOTH_sum = z(S) − z(五题)` 是两个成分的**精确线性组合** ->
放进同一个设计矩阵是 **秩亏**,`BOTH_sum` 的系数**不可识别**。
**⇒ 那个版本的检验按构造不可能跑,而不是「跑出来是零」。本轮当场验证这一点并跳过它。**
`BOTH_min = min(z(S), −z(五题))` 是**非线性**的 -> **它可以携带成分之外的东西**。**只有它可跑。**

ESTIMAND        `治疗性 ~ BOTH_min + z(S) + (−z(五题)) + 羞耻 + c3⁻ + 类别数`;
                主量 = `BOTH_min` 的偏系数。
判据(**先标支**,`#379c`)
                【两支】**秩检查**必须先证明 `sum` 版秩亏、`min` 版满秩(否则我在读一个不可识别的系数)·
                        负对照用**越阈率** · guard 26 **显式传 branch**。
                【非零支】`BOTH_min` 控制两成分后**仍越阈** -> 它带着成分之外的东西;
                【零支】**塌到零** -> 它的全部内容就是两个成分,**「潜变量」这个词不能用**。
⚠ 零的种类     `offset_control`:`BOTH_min` 与两成分**构造性相关** -> **这个零绝不是零**;
                零 = **`lib.nulls.perm_in`** 在掩码内打乱结局后的分布(保住三者之间的构造关系)。
IMPOSSIBLE      ① `min` 的非线性可能只是**截断**,而截断本身能解释一点方差,**那不是「潜变量」**;
                ② 「治疗性」与两成分本来就有关(`#398b`)-> 这不是一个干净的外部结局,只是**没参与构造**的结局;
                ③ 「还在」不等于「是一个心理构念」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
THC=next(c for c in d.columns if 'vmq8jqw' in str(c))
th=pd.to_numeric(d[THC],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
INT=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i]) for i in range(5)]),1)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(th)&np.isfinite(INT)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
A=np.full(NN,np.nan); Bv=np.full(NN,np.nan)
A[M]=z(S,M); Bv[M]=-z(INT,M)                       # 两个成分:冷门 · 常规也管用
SUM=A+Bv; MIN=np.where(np.isfinite(A)&np.isfinite(Bv),np.minimum(A,Bv),np.nan)
print(f"n=**{n:,}**")
def rank_of(cols):
    X=np.column_stack([np.ones(n)]+[z(v,M) for v in cols])
    return int(np.linalg.matrix_rank(X)),X.shape[1]
r1,k1=rank_of([SUM,A,Bv]); r2,k2=rank_of([MIN,A,Bv])
print(f"\n⚠⚠ **秩检查(跑之前,`#387b` 的做法)**:")
print(f"   `[1, BOTH_sum, z(S), −z(五题)]` -> 秩 **{r1}** / 列数 **{k1}** -> "
      f"**{'秩亏 —— `BOTH_sum` 的系数不可识别,这个检验按构造不可能跑' if r1<k1 else '满秩'}**")
print(f"   `[1, BOTH_min, z(S), −z(五题)]` -> 秩 **{r2}** / 列数 **{k2}** -> "
      f"**{'满秩 —— 可跑' if r2==k2 else '秩亏'}**")
print(f"   ⇒ **只跑 `min` 版。`sum` 版不是「跑出来是零」,是「不可能跑」。**")
assert r1<k1 and r2==k2, "秩检查与预期不符 —— 停手"
CTRL=[A,Bv,sh,C3,ncat]
def fit(y,xs):
    X=np.column_stack([np.ones(n)]+[z(v,M) for v in xs]); yy=z(y,M)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(n-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
bAlone,sAlone=fit(th,[MIN,sh,C3,ncat])
bJ,sJ=fit(th,[MIN]+CTRL)
NP_=400
nul=np.array([fit(perm_in(th,M,8200+s),[MIN]+CTRL)[0][0] for s in range(NP_)])
THR=float(np.percentile(np.abs(nul),95))
print(f"\n① `BOTH_min` 不控制两成分:**{bAlone[0]:+.4f}** (se {sAlone[0]:.4f})")
print(f"② `BOTH_min` **控制两成分后**:**{bJ[0]:+.4f}** (se {sJ[0]:.4f})")
print(f"   两成分自己:z(S) **{bJ[1]:+.4f}** · −z(五题) **{bJ[2]:+.4f}**")
print(f"\n⚠ offset 零(**`lib.nulls.perm_in`** 掩码内打乱结局;`BOTH_min` 与两成分**构造性相关**,"
      f"所以这个零不该是零):**{nul.mean():+.5f} ± {nul.std():.5f}** · |值| 95 分位 **{THR:.5f}**")
print(f"   -> **{(bJ[0]-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**仍越阈 -> 带着成分之外的东西**' if abs(bJ[0])>THR else '**塌掉 -> 它的全部内容就是两个成分**'}")
negs=np.array([fit(perm_in(th,M,97000+s),[MIN]+CTRL)[0][0] for s in range(200)])
rate=float((np.abs(negs)>THR).mean())
print(f"\n负对照(**越阈率**):**{100*rate:.1f}%**")
MDE=None
print(f"guard 26 = **MDE 扫描**,每级 30 次:")
for gg in (0.02,0.03,0.04,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(200+int(gg*1000)*89+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(MIN,M)+rg.standard_normal(n)
        if abs(fit(y,[MIN]+CTRL)[0][0])>THR: hit+=1
    print(f"   种植 **{gg:+.3f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
NONNULL=abs(bJ[0])>THR
CONT=abs(float(bJ[0])) if NONNULL else 0.05
print(f"   **MDE = {MDE_:.3f}** · 争议幅度 = **{CONT:.4f}**({'实测(非零支)' if NONNULL else '有意义(零支)'})")
pd.DataFrame([dict(v_alone=bAlone[0],v_joint=bJ[0],v_se=sJ[0],v_thr=THR,v_mde=MDE_,
                   v_bS=bJ[1],v_bI=bJ[2],v_n=n,v_rank_sum=r1,v_rank_min=r2)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'latent.csv',index=False)
g=Gate('BOTH 带着成分之外的东西吗')
g.asserted('★【两支】秩检查:`sum` 版秩亏(不可识别)· `min` 版满秩(可跑)',r1<k1 and r2==k2,
           f"sum {r1}/{k1} · min {r2}/{k2}",kind='control')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,
    what='MDE 扫描 80% 检出',branch='non_null' if NONNULL else 'null')
g.asserted('★【两支】offset 零非退化(构造性相关)',nul.std()>0,
           f"{nul.mean():+.5f} ± {nul.std():.5f}",kind='control')
if 0.01<=rate<=0.12 and r1<k1 and r2==k2:
    if NONNULL:
        g.asserted('★【非零支】`BOTH_min` 控制两成分后仍越阈 -> 带着成分之外的东西',True,
                   f"{bJ[0]:+.4f} vs {THR:.4f}")
    else:
        g.asserted('★【零支】塌掉且 MDE < 0.05 -> 它的全部内容就是两个成分',MDE_<0.05,
                   f"{bJ[0]:+.4f} vs {THR:.4f} · MDE {MDE_:.3f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

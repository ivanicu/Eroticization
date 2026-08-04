import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A118 R373 -- 那条界线是不是**非内容**的

`#327c`:两件不相交的仪器给出同一条界线,而**四次命名都失败**。
**也许问题在于我一直在找**内容学**的名字,而这条界线可能是**结构**性的。**

ESTIMAND        把 `c3` 的 32 维**块载荷**回归到三个**与内容无关**的块属性上:
                ① **块内平均流行度**(多少人选)· ② **块内选项数** · ③ **块内流行度离散度**;
                报 R² 与各自偏相关,并给**人层自助**的 R² 区间。
KILL            **若某一个解释掉大部分 -> 这条界线是**结构**性质而非内容维度,
                那会一次性解释四次命名为什么都失败;
                若三个都解释不了 -> 它确实是内容的,而这份数据里没有能命名它的题。**
POSITIVE CTRL   合成一个**只由流行度驱动**的载荷向量 -> 回归必须抓到(R² > 0.8)。
NEGATIVE CTRL   随机载荷向量 -> R² 必须 ≈ 3/32 ≈ 0.09(三个自由度在 32 点上的期望)。
⚠⚠ 同义反复警告  **流行度与「越轨」在概念上高度重合**(越轨的东西本来就少人选)。
                **所以若 ① 解释掉大部分,那不是一个发现,是一个同义反复的度量** ——
                真正的问题会变成「**除了『少人选』之外还剩什么**」,
                那一步要用**控制流行度后的残差载荷**问,而本轮同时报它。
IMPOSSIBLE      n=32 个块,三个预测量 —— 这个回归的功效很低,R² 的区间会很宽。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

PREV=np.array([float(M.mean(0).mean()) for M,_ in MB])
NOPT=np.array([float(M.shape[1]) for M,_ in MB])
DISP=np.array([float(M.mean(0).std()) for M,_ in MB])
X=np.column_stack([np.ones(NB)]+[ (v-v.mean())/v.std() for v in (PREV,NOPT,DISP)])
def fit(v):
    b,*_=np.linalg.lstsq(X,v,rcond=None); r=v-X@b
    return 1-float(r@r)/float(((v-v.mean())**2).sum()), b
r2,b=fit(v0)
def part(j):
    cols=[0]+[i for i in (1,2,3) if i!=j]
    Xo=X[:,cols]
    rv=v0-Xo@np.linalg.lstsq(Xo,v0,rcond=None)[0]
    rx=X[:,j]-Xo@np.linalg.lstsq(Xo,X[:,j],rcond=None)[0]
    return float(np.corrcoef(rv,rx)[0,1])
LAB=['块内平均流行度','块内选项数','流行度离散度']
print(f"n = **{NB}** 块 · 三个与内容无关的属性")
print(f"**R² = {r2:.4f}**")
for j,l in enumerate(LAB,1):
    print(f"   {l:<14} 系数 **{b[j]:+.4f}** · 偏相关 **{part(j):+.4f}**")
NBOOT=200; rg=np.random.default_rng(3131)
bs=[]
for _ in range(NBOOT):
    vv=load_of(ALLR[rg.integers(0,len(ALLR),len(ALLR))],ref=v0)
    bs.append(fit(vv)[0])
q=np.percentile(bs,[2.5,50,97.5])
print(f"   人层自助 {NBOOT} 次的 R²:**[{q[0]:.4f}, {q[2]:.4f}]** 中位 {q[1]:.4f}")
rgs=np.random.default_rng(8)
syn=(PREV-PREV.mean())/PREV.std()+0.3*rgs.standard_normal(NB)
print(f"\n正对照(只由流行度驱动的合成载荷):R² **{fit(syn)[0]:.4f}**")
ng=[fit(rgs.standard_normal(NB))[0] for _ in range(400)]
print(f"负对照(随机载荷 400 次):R² **{np.mean(ng):.4f} ± {np.std(ng):.4f}**(解析期望 3/32 = {3/NB:.3f})")
# ⚠ 同义反复那一步:控制流行度之后,载荷还剩什么图样
Xp=np.column_stack([np.ones(NB),(PREV-PREV.mean())/PREV.std()])
res=v0-Xp@np.linalg.lstsq(Xp,v0,rcond=None)[0]
o=np.argsort(-res)
print(f"\n⚠ 控制流行度后的**残差载荷**,两端各 4 块(这是「除了少人选还剩什么」):")
for i in o[:4]: print(f"   {res[i]:+.3f}  {NAMES[i][:56]}")
print("   ---")
for i in o[::-1][:4]: print(f"   {res[i]:+.3f}  {NAMES[i][:56]}")
print(f"   残差与原载荷的相关 **{np.corrcoef(res,v0)[0,1]:+.4f}**"
      f"(=√(1−R²_流行度) 的量级,报出来是为了说明残差还剩多少)")
T=pd.DataFrame([dict(v_term='R²',v_val=r2)]+[dict(v_term=l,v_val=part(j)) for j,l in enumerate(LAB,1)]
               +[dict(v_term='R²自助lo',v_val=float(q[0])),dict(v_term='R²自助hi',v_val=float(q[2]))])
check_columns(T,'R373'); T.to_csv(pathlib.Path(__file__).parent/'results'/'structural.csv',index=False)
gg=Gate('那条界线是不是非内容的')
gg.asserted('★ 正对照:只由流行度驱动的合成载荷 -> R² > 0.8',fit(syn)[0]>0.8,f"R² {fit(syn)[0]:.4f}")
gg.asserted('★ 负对照:随机载荷的 R² ≈ 3/32',abs(np.mean(ng)-3/NB)<0.05,
            f"{np.mean(ng):.4f} ± {np.std(ng):.4f} vs {3/NB:.3f}")
gg.asserted('★ 注册的 kill:三个结构属性是否解释掉大部分(R² > 0.5)',r2>0.5,
            f"R² **{r2:.4f}** 自助 [{q[0]:.4f}, {q[2]:.4f}] · 随机基线 {np.mean(ng):.4f}")
gg.asserted('⚠⚠ 同义反复警告:流行度与「越轨」概念上高度重合',True,
            '若 ① 解释掉大部分,那不是发现是同义反复;真正的问题是「除了少人选还剩什么」,'
            '本轮同时报了控制流行度后的残差载荷')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

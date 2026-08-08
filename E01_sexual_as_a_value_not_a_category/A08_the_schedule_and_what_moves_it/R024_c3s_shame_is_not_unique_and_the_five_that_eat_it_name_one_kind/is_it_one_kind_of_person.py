import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A136 R411 -- 「同一种人」是我读出来的,还是数据里的一个对象

`#366c` 我写的是「最大的五口方向完全一致 —— 它们是同一种人」。
**「同一种人」是一个一维性主张,而我从来没测过它。**

ESTIMAND        ① 五个量(性伴数 · 外向性 · 年龄 · 神经质 · 无力感)按 `c3⁻` 对齐符号后的
                   **相关矩阵特征谱**,主量 = **PC1 份额**;
                ② 「五个同号」在 20 个变量里到底罕不罕见 —— 随机 5 元组全同号的频率;
                ③ 若一维成立:PC1 单独能吃掉 `c3⁻` 那份羞耻的多少(对比 `#366a` 的 40.4%)。
KILL(条件式)  仅当三个对照都过 -> 判:**PC1 份额是否越过 offset 零的 95 分位**。
                越过 -> 「同一种人」是数据里的对象,那 40% 配拿一个名字;
                没越过 -> 页面/账上的「同一种人」要退回成**五条并列的相关**。
POSITIVE CTRL   合成真一维五元组(单因子 + 噪声)-> PC1 份额必须远高于 offset 零。
NEGATIVE CTRL   五个独立噪声 -> 谱必须平(PC1 ≈ 逐个打乱的零)。
⚠ 零的种类     **`offset_control`,而且这个零绝不该是零,原因要点名**:
                (i) 五个变量的 PC1 天然 > 1/5(纯采样噪声就有谱展开);
                (ii) **更要命的是我按 `c3⁻` 对齐了符号** —— 让每个 `corr(V_i, c3⁻)>0`
                     会在 V 之间**诱导**一个 ≈ r_i·r_j 的相关下限。
                所以零 = **合成五个变量,每个与 `c3⁻` 的相关等于实测的 |r_i|,彼此其余部分独立**,
                走同一条对齐流程 -> 得到的 PC1 份额就是「只由共同指向 c3⁻ 造出来的一维性」。
⚠ 多重性       ② 里随机 5 元组是穷举/大样本抽取,报频率不报单例。
IMPOSSIBLE      PC1 是一个方向,不是一个人;一维性成立也只说这五个量共变,
                不说存在一个心理构念。命名仍然要独立的清白语境编码者(`#203c`)。
"""
import numpy as np, pandas as pd, warnings, hashlib, itertools
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A136_is_c3_shame_its_own/R410_commonality_vs_person_variables/run.py').read_text()
# ⚠ R410 的正文里自己也有一个 `_SRC`(指向 R347),exec 之后会把这个名字冲掉 ——
# 所以先把 R410 的源存到一个不会被覆盖的名字上,再 exec。
_R410=_SRC
exec(_R410.split('"""',2)[2].split('base=ok&np.isfinite(C3)')[0])
exec('def decomp'+_R410.split('def decomp',1)[1].split('def perm_finite')[0])  # 跨轮依赖显式化(P16)
base=ok&np.isfinite(C3)&np.isfinite(sh)

FIVE=['性伴数','外向性','年龄','神经质','无力感']
PD_=dict(PANEL)
m=base.copy()
for n_ in FIVE: m&=np.isfinite(PD_[n_])
n=int(m.sum()); print(f"n={n:,} · 五个量 = {' · '.join(FIVE)}\n")

def spectrum(cols):
    """cols: list of 1-d arrays already restricted to a mask. 返回 (PC1 份额, 特征值)."""
    Z=np.column_stack([(c-c.mean())/max(c.std(),1e-12) for c in cols])
    C=np.corrcoef(Z.T); w=np.linalg.eigvalsh(C)[::-1]
    return float(w[0]/w.sum()), w

R_C3=[float(np.corrcoef(PD_[n_][m],C3[m])[0,1]) for n_ in FIVE]
SGN=[1. if r>0 else -1. for r in R_C3]
cols=[SGN[i]*PD_[FIVE[i]][m] for i in range(5)]
pc1,w=spectrum(cols)
print("① 对齐后的相关矩阵(对齐 = 让每个量都指向更高的 `c3⁻`):")
Zc=np.column_stack([(c-c.mean())/c.std() for c in cols]); CM=np.corrcoef(Zc.T)
print("        "+"".join(f"{x[:5]:>8}" for x in FIVE))
for i,x in enumerate(FIVE):
    print(f"   {x[:5]:<6}"+"".join(f"{CM[i,j]:+8.3f}" for j in range(5)))
off=CM[np.triu_indices(5,1)]
print(f"\n   非对角:中位 **{np.median(off):+.4f}** · 范围 [{off.min():+.4f}, {off.max():+.4f}] · "
      f"为正 **{int((off>0).sum())}/10**")
print(f"   特征值 " + " · ".join(f"{x:.3f}" for x in w) + f"  -> **PC1 份额 {pc1:.1%}**")

# ---- offset 零:只由「共同指向 c3⁻」造出来的一维性 ----
def offset_draw(seed):
    rg=np.random.default_rng(seed); c=(C3[m]-C3[m].mean())/C3[m].std()
    cs=[]
    for r in R_C3:
        a=abs(r); v=a*c+np.sqrt(max(1-a*a,1e-9))*rg.standard_normal(n)
        cs.append(v)                      # 已按 c3⁻ 正向,等于走完同一条对齐流程
    return spectrum(cs)[0]
NP_=400
nul_off=np.array([offset_draw(3000+s) for s in range(NP_)])
thr=float(np.percentile(nul_off,95))
print(f"\n⚠ offset 零(五个量各自与 `c3⁻` 相关 |r| = {', '.join(f'{abs(r):.3f}' for r in R_C3)},"
      f"其余独立,走同一条对齐流程):")
print(f"   PC1 份额零 = **{nul_off.mean():.2%} ± {nul_off.std():.2%}** · 95 分位 **{thr:.2%}**")
print(f"   实测 **{pc1:.2%}** -> {'**越阈**' if pc1>thr else '**未越阈**'} · "
      f"距零 {(pc1-nul_off.mean())/max(nul_off.std(),1e-12):+.2f} sd")

# ---- 对照 ----
rgc=np.random.default_rng(11)
f0=rgc.standard_normal(n)
pos=[0.7*f0+np.sqrt(1-0.49)*rgc.standard_normal(n) for _ in range(5)]
neg=[rgc.standard_normal(n) for _ in range(5)]
p_pos=spectrum(pos)[0]; p_neg=spectrum(neg)[0]
print(f"\n正对照(真单因子,载荷 0.7):PC1 **{p_pos:.2%}**")
print(f"负对照(五个独立噪声):PC1 **{p_neg:.2%}**")

# ---- ② 「五个同号」在 20 个变量里罕不罕见 ----
NM=[x for x,_ in PANEL]
agree={}
for x in NM:
    v=PD_[x]; mm=base&np.isfinite(v)
    if mm.sum()<300: continue
    rc=np.corrcoef(v[mm],C3[mm])[0,1]; rs=np.corrcoef(v[mm],sh[mm])[0,1]
    agree[x]=rc*rs>0
K=len(agree); nagree=sum(agree.values())
print(f"\n② 「同号」= corr(V,c3⁻)·corr(V,羞耻) > 0。在 {K} 个可算的变量里,**{nagree} 个同号**"
      f"({nagree/K:.0%})。")
subs=list(itertools.combinations(sorted(agree),5))
allsame=sum(1 for s_ in subs if all(agree[x] for x in s_))
print(f"   随机 5 元组全同号的频率 = **{allsame}/{len(subs)} = {allsame/len(subs):.1%}**")
print(f"   ⇒ 「最大的五口全同号」" +
      (" **本身并不罕见** —— 这一半是我读出来的" if allsame/len(subs)>0.10 else " **是罕见的**"))

# ---- ③ PC1 单独吃掉多少 ----
Zc2=np.column_stack([(c-c.mean())/c.std() for c in cols])
Cm=np.corrcoef(Zc2.T); ev=np.linalg.eigh(Cm)[1][:,-1]
pc=np.full(NN,np.nan); pc[m]=Zc2@ev
r3=decomp(pc,sh,C3,base)
print(f"\n③ PC1 单独:R²(PC1) {r3[1]:.3f}pp · 共享 **{r3[4]:+.3f}pp** · 保留 **{r3[5]:.1%}** "
      f"(全 20 个 panel 是 59.6%,`#366a`)")

# ---- ④ ⚠ 最强混淆:这五个是**在同一份数据上按同一个判据挑出来的** ----
# 所以「PC1 拿走 92%」必须对着**随机 5 元组的 PC1** 比,而不是对着 0。
TAKE_FULL=100*(1-0.596)          # #366a:全 20 个 panel 拿走的比例
def pc1_take(names,seed=0):
    mm=base.copy()
    for x in names: mm&=np.isfinite(PD_[x])
    if mm.sum()<300: return None
    rr=[float(np.corrcoef(PD_[x][mm],C3[mm])[0,1]) for x in names]
    cs=[(1. if rr[i]>0 else -1.)*PD_[names[i]][mm] for i in range(len(names))]
    Z=np.column_stack([(c-c.mean())/max(c.std(),1e-12) for c in cs])
    e=np.linalg.eigh(np.corrcoef(Z.T))[1][:,-1]
    q=np.full(NN,np.nan); q[mm]=Z@e
    r=decomp(q,sh,C3,base)
    return None if r is None else 100*(1-r[5])
sel=pc1_take(FIVE)
rg4=np.random.default_rng(909); pool=sorted(agree)
rand=[pc1_take(list(rg4.choice(pool,5,replace=False))) for _ in range(300)]
rand=np.array([x for x in rand if x is not None])
print(f"\n④ ⚠ 最强混淆:这五个是**在同一份数据上按同一个判据挑出来的**。")
print(f"   选出的五元组 PC1 拿走 **{sel:.1f}%**(全 20 个 panel 拿走 {TAKE_FULL:.1f}%,占 **{100*sel/TAKE_FULL:.0f}%**)")
print(f"   随机 5 元组 PC1 拿走 **{rand.mean():.1f}% ± {rand.std():.1f}%** · 95 分位 **{np.percentile(rand,95):.1f}%** "
      f"· 最大 {rand.max():.1f}%")
print(f"   -> 选出的比随机高 **{(sel-rand.mean())/max(rand.std(),1e-12):+.1f} sd**,"
      f"超过随机的 **{100*(rand<sel).mean():.1f}%**")

g=Gate('「同一种人」是数据里的对象,还是我读出来的')
C_POS=p_pos>0.55; C_NEG=abs(p_neg-nul_off.mean())<0.10; C_OFF=nul_off.std()>0
g.asserted('★ 正对照:真单因子 -> PC1 份额 >55%',C_POS,f"{p_pos:.2%}",kind='control')
g.asserted('★ 负对照:五个独立噪声 -> PC1 贴着零',C_NEG,
           f"{p_neg:.2%} vs 零 {nul_off.mean():.2%}",kind='control')
g.asserted('★ offset 零非退化(它绝不该是零:对齐诱导了下限)',C_OFF,
           f"零 = {nul_off.mean():.2%} ± {nul_off.std():.2%},远高于 1/5",kind='control')
if C_POS and C_NEG and C_OFF:
    g.asserted('★ 注册的 kill:PC1 份额越过 offset 零的 95 分位',pc1>thr,
               f"实测 {pc1:.2%} vs 阈 {thr:.2%}")
    g.asserted('★ ④ 选出的五元组 PC1 越过**随机 5 元组**的 95 分位(选择偏差的零)',
               sel>float(np.percentile(rand,95)),
               f"{sel:.1f}% vs 阈 {np.percentile(rand,95):.1f}%(随机均值 {rand.mean():.1f}%)")
else:
    g.asserted('★ 注册的 kill(对照未过 -> 不判)',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A136 R413 -- 换估计量,不换问题:交互项能不能看见分段设计看不见的东西

`#368d`:分段设计的 MDE 是 **+0.20**,而有意义的集中是 +0.05 —— 差 4×,所以「均匀」报不了。
**但 0.20 是那个估计量的极限,不是数据的极限**:五段各 1,322 人,全距是**五个噪声量的极差**,
那是所有统计量里最不稳的一个。**同一个问题,换一个把全部 6,609 人用进一个系数的估计量。**

ESTIMAND        `羞耻 ~ c3⁻ + PC1 + c3⁻×PC1`(全标准化),主量 = **交互系数**。
                「92% 由某一端撑起」= 交互非零;「摊在整条线上」= 交互为零。
⚠ 顺序          **先算 MDE,再看系数**(`#312a`)。MDE 若仍 > 0.05,**就不要看那个系数**。
KILL(条件式)  仅当三个对照都过 **且 MDE < 0.05** -> 判:交互是否越过 offset 零的 95 分位。
POSITIVE CTRL   沿用 `#368` 的种植(只有最高五分位有关系)-> 交互必须亮。
NEGATIVE CTRL   五段真的一样 -> 交互必须是零。
⚠ 零的种类     `offset_control`,**这个零不该是零**:`c3⁻` 与 PC1 相关 +0.1632,
                乘积项与两个主效应共线,而 `c3⁻` 与羞耻本身有关系 -> 乘积项会**借到**主效应的一部分。
                零 = **只置换 PC1**(保住 `c3⁻ ↔ 羞耻`,只打断 PC1 的配对)后重算交互的分布。
⚠ guard 24     PC1 必须先定向,否则交互的符号无意义。
IMPOSSIBLE      交互是一个线性形状;若真实的集中是非单调的(只在中间),它照样看不见。
                本轮把这一条明说,而不是让「零」去承担它。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_R412=(ROOT/'E01_sexual_as_a_value_not_a_category/A136_is_c3_shame_its_own/R412_where_along_the_line/run.py').read_text()
exec(_R412.split('"""',2)[2].split('QS=np.quantile')[0])
QS=np.quantile(PC[mm],[0,.2,.4,.6,.8,1.]); QS[-1]+=1e-9
def seg(i):
    j=idx[(PC[idx]>=QS[i])&(PC[idx]<QS[i+1])]; return j

zc=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def inter(y,perm_pc=None):
    """返回 (交互系数, se, n)。perm_pc: 打乱后的 PC(offset 零用)。"""
    p=PC if perm_pc is None else perm_pc
    g=mm&np.isfinite(y)&np.isfinite(p)
    n=int(g.sum())
    a=zc(C3,g); b=zc(p,g); yy=zc(y,g)
    X=np.column_stack([np.ones(n),a,b,a*b])
    bb,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@bb
    s2=float(r@r)/(n-4); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(bb[3]),float(se[3]),n

def perm_pc(seed):
    rg=np.random.default_rng(seed); z=PC.copy()
    j=np.flatnonzero(np.isfinite(z)&mm); z[j]=z[rg.permutation(j)]; return z

def plant(g_,seed,shape='linear'):
    """在真数据的 c3⁻ 上种一个交互,噪声独立。shape: linear=沿 PC1 线性;q5=只有最高段。"""
    rg=np.random.default_rng(seed); y=np.full(NN,np.nan)
    j=np.flatnonzero(mm); a=zc(C3,mm); b=zc(PC,mm)
    w=b if shape=='linear' else (PC[j]>=QS[4]).astype(float)
    y[j]=0.11*a+g_*a*w+rg.standard_normal(len(j))
    return y

# ---------- ① 先算 MDE(在看系数之前) ----------
NUL0=np.array([inter(sh,perm_pc(6000+s))[0] for s in range(400)])
THR=float(np.percentile(np.abs(NUL0),95))
print(f"⚠ offset 零(**只置换 PC1**,保住 `c3⁻ ↔ 羞耻`;这个零不是零 —— "
      f"corr(c3⁻,PC1) = {np.corrcoef(C3[mm],PC[mm])[0,1]:+.4f},乘积项与主效应共线):")
print(f"   交互零 = **{NUL0.mean():+.5f} ± {NUL0.std():.5f}** · |值| 95 分位 **{THR:.5f}**\n")
print("① **先算 MDE,再看系数**(`#312a`)—— 沿 PC1 线性种植,每级 40 次:")
MDE=None
for g_ in (0.010,0.020,0.030,0.040,0.050,0.080):
    hit=sum(1 for s in range(40) if abs(inter(plant(g_,70000+s))[0])>THR)
    print(f"   种植 **{g_:+.3f}** -> 检出 **{hit}/40 = {hit*2.5:>5.1f}%**")
    if MDE is None and hit>=32: MDE=g_
MDE_=MDE if MDE else 0.10
print(f"   **MDE(80% 检出)= {MDE_:.3f}**(分段设计是 0.20,`#368d`)· "
      f"有意义的效应 0.05 -> {'**MDE < 0.05,可以看系数**' if MDE_<0.05 else '**MDE ≥ 0.05,不要看系数**'}")

# ---------- ② 对照 ----------
b_pos,se_pos,_=inter(plant(0.15,4242,'q5'))
b_neg,se_neg,_=inter(plant(0.0,4243,'linear'))
print(f"\n② 正对照(只有最高五分位有关系,g=0.15):交互 **{b_pos:+.5f}** (se {se_pos:.5f}) vs 阈 {THR:.5f}")
print(f"   负对照(五段真的一样,g=0):交互 **{b_neg:+.5f}** (se {se_neg:.5f}) vs 阈 {THR:.5f}")

# ---------- ③ 实测 ----------
b_obs,se_obs,n_obs=inter(sh)
print(f"\n③ 实测(n={n_obs:,}):交互 **{b_obs:+.5f}** · se {se_obs:.5f} · "
      f"95% CI [{b_obs-1.96*se_obs:+.5f}, {b_obs+1.96*se_obs:+.5f}]")
print(f"   vs offset 零 {NUL0.mean():+.5f} ± {NUL0.std():.5f} -> "
      f"**{(b_obs-NUL0.mean())/max(NUL0.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(b_obs)>THR else '**未越阈**'}")
pd.DataFrame([dict(v_b=b_obs,v_se=se_obs,v_n=n_obs,v_thr=THR,v_mde=MDE_,
                   v_bpos=b_pos,v_bneg=b_neg)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'interaction.csv',index=False)

# ---------- ④ ⚠ scope:交互是一个**线性**形状 ----------
# 正对照只是擦线过(+0.02295 vs 0.02243),而它种的是**阶跃**(只有最高段)。
# 所以这个零的 MDE 有**两个**,必须分开报,否则 scope 被悄悄扩大(§2 的 regime)。
print(f"\n④ ⚠ scope:交互是**线性**形状 —— 阶跃式集中它看不清。分开报两个 MDE:")
MDE_STEP=None
for g_ in (0.05,0.10,0.15,0.20,0.30):
    hit=sum(1 for s_ in range(40) if abs(inter(plant(g_,80000+s_,'q5'))[0])>THR)
    print(f"   **阶跃**(只有最高段)种植 {g_:+.3f} -> 检出 **{hit}/40 = {hit*2.5:>5.1f}%**")
    if MDE_STEP is None and hit>=32: MDE_STEP=g_
MDE_S=MDE_STEP if MDE_STEP else 0.40
print(f"   -> **线性 MDE {MDE_:.3f} · 阶跃 MDE {MDE_S:.3f}**")
print(f"   ⇒ 这个零覆盖的是**线性调节**;阶跃式集中只被 `#368e` 的界(≥+0.20)挡住,"
      f"{'两者之间没有缝' if MDE_S<=0.20 else '**两者之间有一条缝**'}。")

g=Gate('交互项能不能看见分段设计看不见的东西')
CP=abs(b_pos)>THR; CN=abs(b_neg)<=THR; CO=NUL0.std()>0
g.asserted('★ 正对照:只有最高段有关系 -> 交互必须越阈',CP,f"{b_pos:+.5f} vs {THR:.5f}",kind='control')
g.asserted('★ 负对照:五段真的一样 -> 交互必须是零',CN,f"{b_neg:+.5f} vs {THR:.5f}",kind='control')
g.asserted('★ offset 零非退化(乘积项与主效应共线)',CO,
           f"零 {NUL0.mean():+.5f} ± {NUL0.std():.5f}",kind='control')
g.eigenvector_is_anchored('★ guard 24:PC1 已定向',PC,C3,'c3⁻')
if CP and CN and CO and MDE_<0.05:
    g.asserted('★ 注册的 kill:交互越过 offset 零的 95 分位(= 由某一端撑起)',
               abs(b_obs)>THR,f"{b_obs:+.5f} vs {THR:.5f}")
    g.null_claim_uses_null_criteria('★ guard 21(线性):这个零可发布吗','NULL',
        perm_quantile=float((np.abs(NUL0)>abs(b_obs)).mean()),mde=MDE_,
        sensitivity_shown=True,meaningful=0.05)
    g.asserted('★ scope:阶跃 MDE 与 `#368e` 的界(0.20)之间没有缝',MDE_S<=0.20,
               f"阶跃 MDE {MDE_S:.3f} vs 界 0.20")
else:
    g.asserted('★ 注册的 kill(MDE ≥ 0.05 或对照未过 -> 不看系数)',False,
               f"MDE {MDE_:.3f} · 对照 {CP}/{CN}/{CO}")
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

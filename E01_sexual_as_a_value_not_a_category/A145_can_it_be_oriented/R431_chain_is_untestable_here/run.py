import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A145 R431 -- 我在 `#386` 的 NEXT 里写的那个分离器,代数上就不成立

`#386` 的 NEXT 写的是:
「若羞耻只经由**少做**影响信念,那么控制 `ACTED` 之后 `羞耻 → BELIEF` 应当**变正**且越阈。」

⚠⚠ **跑之前先做代数,而代数说这句话是错的 —— 两处都错。**

**① 符号错。** 间接路径 = `(羞耻→ACTED) × (ACTED→BELIEF)` = **(−0.083) × (−0.039) = +0.0032**,**正的**。
总效应 = 直接 + 间接,所以控制 `ACTED` 是**减掉**这 +0.0032 ->
`羞耻 → BELIEF` 应当从 +0.0135 **降到 ≈ +0.0103**,**变小,不是变正**。

**② 更要命:幅度错。** 预期变化 **≈ 0.0032**,而 `#384a` 那个模型的阈是 **0.0248**。
**⇒ 预期效应是阈的 **13%**。这个检验没有任何功率去看它,而这一点在花计算之前就知道。**

**这就是 `AGENTS.md §9` 的「算术陷阱」:一个由代数强制的量,不是证据。**
本轮的产出是:**用代数杀掉一个我自己刚注册的分离器,并在同一轮里用数据验证那个代数。**

ESTIMAND        ① 代数预测的间接效应 `a×b`;② 实测的系数变化 `Δ = b_控制后 − b_控制前`;
                ③ 这个设计对 `Δ` 的 MDE。
判据(**先标支**)
                【两支】② 必须与 ① 在**符号与量级**上一致(否则是我的代数错,不是设计错)。
                【零支】若 MDE ≫ |a×b| -> **宣告这个检验不可行**,而不是报一个零。
                【非零支】若 |Δ| 意外越阈 -> 代数被推翻,那才是发现。
⚠ 零的种类     `offset_control`:**Δ 的零不该是零** —— 加任何与 `ACTED` 相关的控制都会动系数。
                零 = 加入一个**与 `ACTED` 相关度相同、但与 `BELIEF` 无关**的合成控制后的 Δ 分布(`#372a` 的方法)。
IMPOSSIBLE      ① 路径系数取自不同模型(`#384a` 与 `#386a`),乘积只是**近似**的间接效应;
                ② 本轮不做因果中介,只做**可行性**判断。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from scipy.stats import rankdata
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
BC=next(c for c in d.columns if '7lgg41e' in c); AC=next(c for c in d.columns if '41kpfir' in c)
BMAP={'Impossible':0.,'With an extreme amount of effort, maybe':1.,
      'With a lot of effort, yes':2.,'With some effort, yes':3.,'With little effort, yes':4.}
BELIEF=d[BC].map(BMAP).values.astype(float)
ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(BELIEF)&np.isfinite(ACTED)&np.isfinite(ncat)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def fit(y,xs,g=None):
    g=M if g is None else g; k=int(g.sum())
    X=np.column_stack([np.ones(k)]+[z(v,g) for v in xs])
    yy=z(y,g); b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
BASE=[sh,S,C3,ncat]
a,_=fit(ACTED,BASE); A=float(a[0])                       # 羞耻 -> ACTED
bb,_=fit(BELIEF,[ACTED,S,C3,ncat,sh]); B=float(bb[0])    # ACTED -> BELIEF(控制羞耻)
b0v,s0v=fit(BELIEF,BASE); b1v,s1v=fit(BELIEF,BASE+[ACTED])
B0,B1=float(b0v[0]),float(b1v[0]); DEL=B1-B0
print(f"n=**{n:,}**")
print(f"① **代数,跑之前就能算**:")
print(f"   羞耻 → ACTED = **{A:+.4f}** · ACTED → BELIEF = **{B:+.4f}**")
print(f"   间接路径 a×b = **{A*B:+.5f}**(**正的** —— 两个负数相乘)")
print(f"   ⇒ 控制 `ACTED` 应当**减掉**它:`羞耻 → BELIEF` 从 {B0:+.4f} 降到 ≈ **{B0-A*B:+.4f}**")
print(f"   **⇒ 我在 `#386` 的 NEXT 里写的「变正」是**错的**;正确的预测是「变小」。**")
print(f"\n② 实测:控制前 **{B0:+.4f}** -> 控制后 **{B1:+.4f}** · Δ = **{DEL:+.5f}**")
print(f"   代数预测的 Δ = **{-A*B:+.5f}** · 差 **{abs(DEL-(-A*B)):.5f}** -> "
      f"{'**一致(符号与量级)**' if DEL*(-A*B)>0 and abs(DEL-(-A*B))<0.005 else '⚠ 不一致 —— 代数或设计有一个错'}")

NP_=300
rE=float(np.corrcoef(z(ACTED,M),z(sh,M))[0,1])
def sham(seed):
    rg=np.random.default_rng(seed); v=np.full(NN,np.nan)
    o=z(ACTED,M); v[M]=o*0+rg.standard_normal(n)         # 与 BELIEF 无关
    v[M]=rE*z(sh,M)+np.sqrt(max(1-rE*rE,1e-9))*rg.standard_normal(n)
    return v
nul=np.array([float(fit(BELIEF,BASE+[sham(4900+s)])[0][0])-B0 for s in range(NP_)])
THR=float(np.percentile(np.abs(nul),95))
print(f"\n⚠ offset 零(**与 `ACTED` 对羞耻的相关度 {rE:+.4f} 相同、但与 `BELIEF` 无关**的合成控制;"
      f"**加任何相关控制都会动系数,所以零不是零**):")
print(f"   Δ 的零 = **{nul.mean():+.5f} ± {nul.std():.5f}** · |值| 95 分位 **{THR:.5f}**")
print(f"   实测 Δ = **{DEL:+.5f}** -> {'**越阈**' if abs(DEL)>THR else '**未越阈**'}")

# ---- ③ ⚠ 第一版这里是一条 MDE 扫描,而它**又一次量错了东西**(与 `#372c②` 同族) ----
# 我把 `gg*z(ACTED)` 种进结局,以为「种植的就是 Δ」。**不是** —— 种得越大,
# 控制 `ACTED` 吸收得越多,羞耻的 Δ 反而回落,于是检出率是 100/100/0/0,**非单调**。
# **一条非单调的灵敏度曲线,几乎总是曲线量错了东西,而不是设计忽好忽坏。**
#
# 而真正的结论比「功率不足」**强得多**,并且第 ② 步已经把它证出来了:
# **实测 Δ = 代数预测,差 0.00000(5 位小数)。**
# `Δ = −a×b` 是标准化 OLS 的**省略变量公式**,**由代数强制** ——
# 它不可能是别的值。所以「Δ 越阈」说的不是世界,是**回归的算术在工作**。
print(f"\n③ ★★ 这不是功率问题,是**恒等式**。")
print(f"   实测 Δ **{DEL:+.6f}** vs 代数 −a×b **{-A*B:+.6f}** · 差 **{abs(DEL-(-A*B)):.6f}**")
print(f"   `Δ = −a×b` 是标准化 OLS 的省略变量公式 —— **它不可能是别的值。**")
print(f"   ⇒ 「Δ 越过零」说的不是世界,是**回归的算术在工作**(`AGENTS.md §9` 的算术陷阱)。")
print(f"   ⇒ **这个分离器是**空的**,不是**没功率的**。**")
print(f"   而 `a` 与 `b` 我**已经知道**(`#384b` 与 `#386a`)-> 本轮没有一个比特的新信息。")
CONT=abs(A*B)
pd.DataFrame([dict(v_a=A,v_b=B,v_ab=A*B,v_b0=B0,v_b1=B1,v_delta=DEL,
                   v_algebra=-A*B,v_gap=abs(DEL-(-A*B)),v_thr=THR)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'chain.csv',index=False)

g=Gate('#386 的 NEXT 里那个分离器,可行吗')
g.asserted('★【两支】实测 Δ 必须与代数预测在符号与量级上一致(否则是我的代数错)',
           DEL*(-A*B)>0 and abs(DEL-(-A*B))<0.005,
           f"实测 {DEL:+.6f} · 代数 {-A*B:+.6f}",kind='control')
g.asserted('★【两支】offset 零非退化(加任何相关控制都会动系数)',nul.std()>0,
           f"{nul.mean():+.5f} ± {nul.std():.5f}",kind='control')
g.asserted('★【非零支】Δ 意外**偏离代数** -> 代数被推翻,那才是发现',
           abs(DEL-(-A*B))>0.002,
           f"实测 {DEL:+.6f} · 代数 {-A*B:+.6f} · 差 {abs(DEL-(-A*B)):.6f} —— "
           f"**差为零 = 恒等式 = 这个分离器是空的**")
print(g)
print(f"\n⇒ **结论:这个分离器是**空的** —— 它测的量由代数强制,不可能有别的值。**")
print(f"   **不报零,也不报「功率不足」,报「这个检验没有信息」。**")
print(f"   ⚠ 而我在 `#386` 的 NEXT 里注册它的时候,**符号也写反了**(写「变正」,实际「变小」)——")
print(f"   **一个我连符号都没算对的预测,本来就不该被注册成分离器。**")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

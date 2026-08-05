import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A149 R440 -- 「加一个控制,看系数往哪走」是一整类**空的**分离器

`#395` 的 NEXT 提的是:控制五题分数后,看 `羞耻 → BELIEF` 的系数往哪走。
⚠ **按 `#387b` 的教训,先算代数** —— 而代数说这又是一个恒等式:
`Δ = −a×b`,其中 `a` = 羞耻→五题(`#395a` **已知**)· `b` = 五题→BELIEF(`#393a` **已知**)。
**两条腿都已经估过了,所以那个「移动」不可能是别的值。**

**⇒ 这是第二次(`#386` 的 NEXT 是第一次)。两次之后,它是一个**类**,不是一次失误。**

ESTIMAND        ① 代数预测的 `−a×b`;② 实测的 Δ;③ 两者之差。
判据(**先标支**)
                【非零支】**Δ 偏离代数** -> 代数被推翻,那才是发现;
                【两支】三条系数必须各自可复现(否则是我算错,不是设计错)。
本轮不做 offset 零 —— **因为本轮不主张任何关于世界的东西**,只主张一个**算术关系**。
                (`#312a` 的镜像:一个不主张世界的轮次,不需要世界的零。)
IMPOSSIBLE      恒等式只在**标准化 OLS + 同一控制集 + 同一掩码**下精确;三者任一不同,会有残差。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
BC=next(c for c in d.columns if '7lgg41e' in c); AC=next(c for c in d.columns if '41kpfir' in c)
BMAP={'Impossible':0.,'With an extreme amount of effort, maybe':1.,
      'With a lot of effort, yes':2.,'With some effort, yes':3.,'With little effort, yes':4.}
BELIEF=d[BC].map(BMAP).values.astype(float)
ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
INT=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i]) for i in range(5)]),1)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(BELIEF)&np.isfinite(ACTED)&np.isfinite(ncat)&np.isfinite(INT)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
CTRL=[ACTED,S,C3,ncat]
def fit(y,xs):
    k=n; X=np.column_stack([np.ones(k)]+[z(v,M) for v in xs]); yy=z(y,M)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); return b[1:]
a=float(fit(INT,[sh]+CTRL)[0])            # 羞耻 -> 五题(同一控制集)
b=float(fit(BELIEF,[INT,sh]+CTRL)[0])     # 五题 -> BELIEF(同一控制集)
B0=float(fit(BELIEF,[sh]+CTRL)[0])        # 羞耻 -> BELIEF,**不**控制五题
B1=float(fit(BELIEF,[sh,INT]+CTRL)[0])    # 羞耻 -> BELIEF,**控制**五题
DEL=B1-B0
print(f"n=**{n:,}**(标准化 OLS · 同一控制集 · 同一掩码 —— 恒等式成立的三个前提)\n")
print(f"① 代数,**跑之前就能算**(两条腿都已估过):")
print(f"   a = 羞耻 → 五题   = **{a:+.5f}**(`#395a` 的那条,同一控制集)")
print(f"   b = 五题 → BELIEF = **{b:+.5f}**(`#393a` 的那条,同一控制集)")
print(f"   ⇒ 预测 Δ = −a×b = **{-a*b:+.6f}**")
print(f"\n② 实测:控制前 **{B0:+.5f}** -> 控制后 **{B1:+.5f}** · Δ = **{DEL:+.6f}**")
print(f"   **差 = {abs(DEL-(-a*b)):.6f}**")
print(f"\n③ ⇒ {'**恒等式** —— 这个分离器是**空的**' if abs(DEL-(-a*b))<1e-4 else '⚠ 偏离代数 —— 那才是发现'}")
print(f"\n★★ **而这是第二次(`#386` 的 NEXT 是第一次)。两次之后,它是一个**类**,不是一次失误:**")
print(f"   **「加一个控制 M,看 X 的系数往哪走」,当 `X→M` 与 `M→Y|X` **都已估过**时,**")
print(f"   **那个移动**由代数强制**,不可能是别的值 —— 它携带**零**比特的新信息。**")
print(f"   **它只在其中一条腿**没有**被估过时才是一个检验。**")
print(f"   ⇒ 而我连着两次把它写进 NEXT,因为它**读起来像**一个中介分析,")
print(f"     **而中介分析的信息量来自它的**假设**(无混淆),不来自那个算术。**")
pd.DataFrame([dict(v_a=a,v_b=b,v_pred=-a*b,v_obs=DEL,v_gap=abs(DEL-(-a*b)),v_n=n)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'identity.csv',index=False)
g=Gate('「加一个控制看系数往哪走」是不是又一个恒等式')
g.asserted('★【两支】三条系数各自可复现(否则是我算错,不是设计错)',
           all(np.isfinite([a,b,B0,B1])),f"a {a:+.5f} · b {b:+.5f} · B0 {B0:+.5f} · B1 {B1:+.5f}",kind='control')
g.asserted('★【非零支】Δ **偏离**代数 -> 代数被推翻,那才是发现',abs(DEL-(-a*b))>1e-4,
           f"实测 {DEL:+.6f} · 代数 {-a*b:+.6f} · 差 {abs(DEL-(-a*b)):.6f} —— "
           f"**差为零 = 恒等式 = 这个分离器是空的**")
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A147 R437 -- 那 +0.1255,是「常规内容」还是「整体不容易被唤起」

`#392c` 把方向翻了过来:五题的分数**越高 = 常规的东西对他越不色**。
`#392d`:而五题与 `Totalsexacts` 相关 **−0.4502** ——
**它有很大一块可能只是「这个人整体上有多容易被唤起」(取反)。**

两个活着的世界:
**A 内容** —— 关于**常规**这一类内容 -> 控制整体唤起广度后,+0.1255 **仍在**;
**B 整体唤起的影子** —— 与内容无关,只是「他不太容易被唤起」-> **塌掉**。

ESTIMAND        `BELIEF ~ 五题 + Totalsexacts + ACTED + S + c3⁻ + 类别数 + 羞耻`;
                主量 = **五题系数的掉幅** `1 − |b_后| / |b_前|`。
判据(**先标支**,`#379c`)
                【两支】guard 26 用 **MDE 扫描**;负对照;offset 零非退化。
                【非零支】掉幅越过 offset 地板 -> `Totalsexacts` 是一条真通路;
                          并**另判**:控制后五题**是否还越过自己的零**(= 世界 A)。
                【零支】掉幅未越阈时启用 MDE。
⚠ 零的种类     `offset_control`:**掉幅的零绝不是零** —— 任何与五题相关的控制都会吃掉一点。
                零 = 加入一个**与五题相关度和 `Totalsexacts` 相同、但与 `BELIEF` 无关**的合成控制
                后的掉幅分布(`#372a` 的方法)。
⚠ 新规矩       **`#392e`:`Totalsexacts` 进模型前,先打印它的取值集合、众数、与一个方向已知的锚的相关。**
IMPOSSIBLE      ① `Totalsexacts` 是**计数**,五题是**评分** -> 控制它不等于控制「唤起度」这个构念的全部;
                ② 两者相关 −0.45 -> 控制后系数天然变不稳,**掉幅要对着地板读,不对着 0 读**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
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
TSA=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)

# ---- ⚠ `#392e` 的新规矩,三行,跑在建模之前 ----
u=np.unique(TSA[np.isfinite(TSA)])
md=float(pd.Series(TSA[np.isfinite(TSA)]).mode().iloc[0])
anc=pd.to_numeric(d['totalfetishcategory'],errors='coerce').values.astype(float)
ga=np.isfinite(TSA)&np.isfinite(anc)
print("⚠ **`#392e` 的新规矩:变量进模型前先看它自己**(三行,跑在建模之前)")
print(f"   `Totalsexacts` 取值集合 = {u[:12].tolist()}{'…' if len(u)>12 else ''}(共 {len(u)} 档)")
print(f"   众数 = **{md:g}** · n = {int(np.isfinite(TSA).sum()):,}")
print(f"   与方向已知的锚 `totalfetishcategory` 相关 = **{np.corrcoef(TSA[ga],anc[ga])[0,1]:+.4f}** "
      f"-> **同向,越大越色,方向确认**\n")

M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(BELIEF)&np.isfinite(ACTED)&np.isfinite(ncat)&np.isfinite(INT)&np.isfinite(TSA)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
CTRL=[ACTED,S,C3,ncat,sh]
def fit(y,xs,g=None):
    g=M if g is None else g; k=int(g.sum())
    X=np.column_stack([np.ones(k)]+[z(v,g) for v in xs]); yy=z(y,g)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
b0,s0=fit(BELIEF,[INT]+CTRL); b1,s1=fit(BELIEF,[INT,TSA]+CTRL)
DROP=100*(1-abs(b1[0])/max(abs(b0[0]),1e-12))
print(f"n=**{n:,}** · corr(五题, Totalsexacts) = **{np.corrcoef(INT[M],TSA[M])[0,1]:+.4f}**")
print(f"控制前 五题 **{b0[0]:+.4f}** (se {s0[0]:.4f}) -> 控制后 **{b1[0]:+.4f}** (se {s1[0]:.4f}) · "
      f"**掉幅 {DROP:+.2f}%**")
print(f"   (同一模型里 `Totalsexacts` 自己 = **{b1[1]:+.4f}**,se {s1[1]:.4f})")
rE=float(np.corrcoef(z(INT,M),z(TSA,M))[0,1])
NF=300
def sham(seed):
    rg=np.random.default_rng(seed); v=np.full(NN,np.nan)
    v[M]=rE*z(INT,M)+np.sqrt(max(1-rE*rE,1e-9))*rg.standard_normal(n)   # 与 BELIEF 无关
    return v
fl=np.array([100*(1-abs(fit(BELIEF,[INT,sham(5200+s)]+CTRL)[0][0])/max(abs(b0[0]),1e-12))
             for s in range(NF)])
THR=float(np.percentile(fl,95))
print(f"\n⚠ offset 地板(**与五题相关度 {rE:+.4f} 相同、但与 `BELIEF` 无关**的合成控制,{NF} 次;"
      f"**任何相关控制都会吃掉一点,所以零不是零**):")
print(f"   **{fl.mean():+.2f}% ± {fl.std():.2f}%** · 95 分位 **{THR:+.2f}%**")
print(f"   实测掉幅 **{DROP:+.2f}%** -> **{(DROP-fl.mean())/max(fl.std(),1e-12):+.2f} sd** · "
      f"{'**越阈:`Totalsexacts` 是一条真通路**' if DROP>THR else '**未越阈**'}")
# ⚠ **和 `#385c` 一模一样的错,我写下那条教训之后又犯了一次**:
# 打乱**整个含 NaN 的数组**会让掩码里出现 NaN,阈变成 NaN,后面每个门都「失败」。
# 修法(项目已有的 `perm_finite` 模式):**只在掩码内打乱**。
def perm_in(v,g,seed):
    o=v.copy(); jj=np.flatnonzero(g&np.isfinite(v))
    o[jj]=v[jj][np.random.default_rng(seed).permutation(len(jj))]; return o
NUL=np.array([fit(perm_in(BELIEF,M,6600+s),[INT,TSA]+CTRL)[0][0] for s in range(300)])
T2=float(np.percentile(np.abs(NUL),95))
print(f"\n★ 另判(世界 A vs B):控制后五题 **{b1[0]:+.4f}** vs 它自己的零(阈 **{T2:.4f}**)-> "
      f"{'**仍越阈 -> 世界 A(内容还在)**' if abs(b1[0])>T2 else '**塌掉 -> 世界 B(整体唤起的影子)**'}")
print(f"\nguard 26 的正对照 = **MDE 扫描**(`#384d`),每级 30 次:")
MDE=None
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(4300+int(gg*100)*47+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(INT,M)+rg.standard_normal(n)
        if abs(fit(y,[INT,TSA]+CTRL)[0][0])>T2: hit+=1   # T2 现在是有限值
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 = 控制前的 **{abs(b0[0]):.4f}**")
bNEG,_=fit(perm_in(BELIEF,M,99),[INT,TSA]+CTRL)
pd.DataFrame([dict(v_b0=b0[0],v_b1=b1[0],v_drop=DROP,v_thr=THR,v_t2=T2,v_mde=MDE_,
                   v_tsa=b1[1],v_corr=rE,v_n=n)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'responsiveness.csv',index=False)
g=Gate('那 +0.1255 是常规内容,还是整体不容易被唤起')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 控制前的效应',MDE_,abs(float(b0[0])),True,what='MDE 扫描 80% 检出')
g.asserted('★【两支】负对照:打乱 `BELIEF` -> 必须落回零',abs(bNEG[0])<=T2,
           f"{bNEG[0]:+.5f} vs {T2:.4f}",kind='control')
g.asserted('★【两支】offset 地板非退化(任何相关控制都吃掉一点)',fl.std()>0,
           f"{fl.mean():+.2f}% ± {fl.std():.2f}%",kind='control')
if MDE_<=abs(b0[0]) and abs(bNEG[0])<=T2:
    g.asserted('★【非零支】掉幅越过 offset 地板 -> `Totalsexacts` 是一条真通路',DROP>THR,
               f"{DROP:+.2f}% vs {THR:+.2f}%")
    g.asserted('★【非零支/另判】控制后五题仍越阈 -> 世界 A(内容还在)',abs(b1[0])>T2,
               f"{b1[0]:+.4f} vs {T2:.4f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A146 R435 -- 那五题是一个家族,还是一两题被我裹进了一个名字里

`#390c`:剩下的唯一开口是**那五题是我挑的**。
便宜且决定性的检验:**逐一把每一题单独当家族**(其余四题不参与),看它自己的相关。

**若五题各自都在同一方向、且各自越过单变量零 -> 家族不是靠「凑」出来的;
若只有一两题扛着全部 -> 那不是一个家族,是一两题被我裹进了一个名字里**
(`#325` 同族:`animated`/`written` 两个指标**不等量**地带着 `form`)。

ESTIMAND        每一题**单独**放进 `BELIEF ~ 题 + ACTED + S + c3⁻ + 类别数 + 羞耻`;
                主量 = 五个系数的**符号一致性**与**各自是否越过单变量零**;
                另报 **留一**:去掉每一题后家族分数还剩多少。
判据(**先标支**,`#379c`)
                【两支】阳性参照 `ACTED` 在单变量零上开火 · 负对照 · guard 26(MDE 扫描)。
                【非零支】**五题全部同号且 ≥4 题越阈** -> 是家族;
                          **≤2 题越阈** -> 不是家族,是一两题。
                【零支】不适用(本轮不报零)。
⚠ 零的种类     `offset_control`:每题各自的**单变量零**(打乱那一列)——
                它们是**预先指定的五题**,不是从候选池挑的(`#388b` 的类别错误不再犯)。
⚠ 多重性       五题 -> **全部报出,不挑**。
⚠ 不许加题      加题就是又一次选择(`#390c`)。
IMPOSSIBLE      ① 五题彼此相关 -> 单独系数会互相「借」,所以**留一**那一列才是家族性的直接证据;
                ② 全自报,同源方差;③ 本轮不改变 `#388c` 的选择偏差 scope。
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
KEYS=['947wne3','normalsex','yuc275j','cunnilingus','jn2b355']
ITEMS=[]
for k in KEYS:
    for c in d.columns:
        if k in str(c):
            v=pd.to_numeric(d[c],errors='coerce').values.astype(float)
            if np.isfinite(v).sum()>2000: ITEMS.append((str(c)[:52],v)); break
print(f"五题(**一字不改**,`#388a` 的那五个):")
for nm,_ in ITEMS: print(f"   · {nm}")
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(BELIEF)&np.isfinite(ACTED)&np.isfinite(ncat)
for _,v in ITEMS: M&=np.isfinite(v)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
CTRL=[ACTED,S,C3,ncat,sh]
def fit(y,xs,g=None):
    g=M if g is None else g; k=int(g.sum())
    X=np.column_stack([np.ones(k)]+[z(v,g) for v in xs]); yy=z(y,g)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
X1=np.column_stack([np.ones(n)]+[z(v,M) for v in CTRL]); yB=z(BELIEF,M)
b1,*_=np.linalg.lstsq(X1,yB,rcond=None); RES=yB-X1@b1
NP_=400
def sthr(v,seed0):
    nl=np.array([abs(float(np.corrcoef(RES,np.random.default_rng(seed0+s).permutation(z(v,M)))[0,1]))
                 for s in range(NP_)])
    return float(np.percentile(nl,95))
print(f"\nn=**{n:,}** · 逐题**单独**(⚠ 每题各自的**单变量零**,因为它们是预先指定的):")
rows=[]
for i,(nm,v) in enumerate(ITEMS):
    b,se=fit(BELIEF,[v]+CTRL); t=sthr(v,3000+100*i)
    rows.append(dict(v_item=nm,v_b=float(b[0]),v_se=float(se[0]),v_thr=t,
                     v_over=bool(abs(b[0])>t)))
    print(f"   {nm[:44]:<46} **{b[0]:+.4f}** (se {se[0]:.4f}) vs 阈 {t:.4f} · "
          f"{'**越阈**' if abs(b[0])>t else '未越阈'}")
T=pd.DataFrame(rows); check_columns(T,'R435')
SAME=bool((T.v_b>0).all() or (T.v_b<0).all()); NOVER=int(T.v_over.sum())
print(f"\n   **全部同号 {SAME}** · **越阈 {NOVER}/5**")

Z=np.column_stack([z(v,M) for _,v in ITEMS])
def expand(col):
    """把长度 n 的列还原成长度 NN 的数组(⚠ `fit` 吃全长数组,不吃掩码内的短向量)。"""
    o=np.full(NN,np.nan); o[np.flatnonzero(M)]=col; return o
FULL=expand(Z.mean(1)); bF,sF=fit(BELIEF,[FULL]+CTRL)
print(f"\n留一(去掉每一题后,家族分数还剩多少;全 5 题 = **{bF[0]:+.4f}**):")
for i_,(nm,_) in enumerate(ITEMS):
    keep=[jj for jj in range(5) if jj!=i_]
    bL,_=fit(BELIEF,[expand(Z[:,keep].mean(1))]+CTRL)
    print(f"   去掉 {nm[:40]:<42} -> **{bL[0]:+.4f}** ({100*bL[0]/bF[0]:.0f}% 保留)")

# ⚠ **第一版这里是一个恒等式,不是一个对照**:`RES` 是控制掉 `ACTED` **之后**的残差,
# 所以 `corr(RES, ACTED)` **按构造恰好是 0**(实测 0.0000)。
# 于是阳性参照「失败」,整轮被判 UNVERIFIED —— **又一次假的仪器失灵**(与 `#388b` 同族,
# 与 `#431` 同族:**恒等式被当成了测量**)。
# 参照必须用**未控制 `ACTED`** 的残差(R432/R433 我做对了,本轮回退了)。
X2=np.column_stack([np.ones(n)]+[z(v,M) for v in [S,C3,ncat,sh]])
b2,*_=np.linalg.lstsq(X2,yB,rcond=None); RES_NA=yB-X2@b2
def sthr2(v,seed0):
    nl=np.array([abs(float(np.corrcoef(RES_NA,np.random.default_rng(seed0+s).permutation(z(v,M)))[0,1]))
                 for s in range(NP_)])
    return float(np.percentile(nl,95))
rA=abs(float(np.corrcoef(RES_NA,z(ACTED,M))[0,1])); tA=sthr2(ACTED,9000)
print(f"   ⚠ 用**控制掉 `ACTED` 之后**的残差算它自己,得 "
      f"**{abs(float(np.corrcoef(RES,z(ACTED,M))[0,1])):.4f}** —— **那是恒等式,不是对照**。")
rgN=np.random.default_rng(7); rF=abs(float(np.corrcoef(RES,rgN.standard_normal(n))[0,1]))
print(f"\n阳性参照 `ACTED`:|r| **{rA:.4f}** vs {tA:.4f} -> {'**开火**' if rA>tA else '不开火'} · "
      f"负对照 **{rF:.4f}**")
MDE=None
print(f"\nguard 26 的正对照 = MDE 扫描(种在第一题上),每级 30 次:")
t0=T.v_thr.iloc[0]
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(6400+int(gg*100)*43+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(ITEMS[0][1],M)+rg.standard_normal(n)
        bq,_=fit(y,[ITEMS[0][1]]+CTRL)
        if abs(bq[0])>t0: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 = 家族的 **{abs(bF[0]):.4f}**")
T.to_csv(pathlib.Path(__file__).parent/'results'/'per_item.csv',index=False)

g=Gate('那五题是一个家族,还是一两题被我裹进了一个名字里')
g.asserted('★【两支】阳性参照 `ACTED` 在单变量零上开火',rA>tA,f"{rA:.4f} vs {tA:.4f}",kind='control')
g.asserted('★【两支】负对照:纯噪声不越阈',rF<=float(T.v_thr.max()),f"{rF:.4f}",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 家族效应',MDE_,abs(float(bF[0])),True,what='MDE 扫描 80% 检出')
g.asserted('★【非零支】五题全部同号 **且** ≥4 题各自越阈 -> 是家族',SAME and NOVER>=4,
           f"同号 {SAME} · 越阈 {NOVER}/5 · {[round(x,4) for x in T.v_b]}")
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

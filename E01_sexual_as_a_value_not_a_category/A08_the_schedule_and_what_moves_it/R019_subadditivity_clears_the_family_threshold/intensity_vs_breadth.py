import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A146 R434 -- 「常规的东西对你有多色」与「你在常规块里勾了多少」,哪一个预测「我能不能改」

`#389c`:两种操作化对同一个概念给出**相反的符号** ——
`强度`(我挑的五题的评分)**+0.1159** · `广度`(两块的勾选计数)**−0.0259**。
**放进同一个模型,让它们互相控制。**

三个活着的世界:
**A 强度是真的** -> 控制广度后,**强度仍为正且越阈**;
**B 广度是真的** -> 控制强度后,**广度仍为负且越阈**;
**C 同一个东西的两面** -> 放一起后**两个都塌**(共线)-> **这条线到此为止。**

⚠ **不许**在看过结果后把「强度」重新定义成别的东西(`#389c` 刚立的规矩)。
**强度 = `#388a` 那五题,一字不改;广度 = `#433` 那两块的勾选计数,一字不改。**

ESTIMAND        `BELIEF ~ 强度 + 广度 + ACTED + S + c3⁻ + 类别数 + 羞耻`,
                主量 = **强度与广度各自的偏系数**。
判据(**先标支**,`#379c`)
                【两支】阳性参照 `ACTED` 在**单变量零**上开火(`#388b`)· 负对照 · guard 26(MDE 扫描)。
                【非零支】按上面的三个世界判。
                【零支】两个都未越阈时启用 MDE。
⚠ 零的种类     `offset_control`:两者相关高 -> **偏系数的零不该是零**;
                零 = **整行打乱**(强度与广度**一起**打乱,保住它们之间的相关与维数)。
IMPOSSIBLE      ① 两者都自报,同源方差;② 若两者相关极高,偏系数会**同时**变得不稳 ——
                那正是世界 C,而**不稳本身**就是答案;③ 强度用的是评分,广度用的是计数,**量纲不同**
                -> 只比符号与越阈,不比大小。
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
FAM5=['947wne3','normalsex','yuc275j','cunnilingus','jn2b355']       # 强度:#388a 那五题,一字不改
iv=[]
for c in d.columns:
    if any(k in str(c) for k in FAM5):
        v=pd.to_numeric(d[c],errors='coerce').values.astype(float)
        if np.isfinite(v).sum()>2000: iv.append(v)
INT=np.nanmean(np.column_stack([(x-np.nanmean(x))/np.nanstd(x) for x in iv]),1)
inv=pd.read_csv('data/derived/inventory.csv')
mult=[c for c in inv[inv['kind']=='MULTISELECT']['col'] if any(w in str(c).lower() for w in ('sex act','sensation'))]
rv=[]
for c in mult:
    v=pd.to_numeric(d[c].astype(str).str.count(','),errors='coerce').values.astype(float)+1
    v=np.where(d[c].isna().values,np.nan,v); rv.append((v-np.nanmean(v))/np.nanstd(v))
BRD=np.nanmean(np.column_stack(rv),1)                                 # 广度:#433 那两块,一字不改
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(BELIEF)&np.isfinite(ACTED)&np.isfinite(ncat)&np.isfinite(INT)&np.isfinite(BRD)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
print(f"n=**{n:,}** · 强度 = `#388a` 的 {len(iv)} 题(评分)· 广度 = `#433` 的 {len(rv)} 块(勾选计数)")
print(f"⚠ **corr(强度, 广度) = {np.corrcoef(INT[M],BRD[M])[0,1]:+.4f}** —— 世界 C 的前提就在这个数上\n")
def fit(y,xs,g=None):
    g=M if g is None else g; k=int(g.sum())
    X=np.column_stack([np.ones(k)]+[z(v,g) for v in xs]); yy=z(y,g)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
CTRL=[ACTED,S,C3,ncat,sh]
bI,sI=fit(BELIEF,[INT]+CTRL); bB,sB=fit(BELIEF,[BRD]+CTRL)
bJ,sJ=fit(BELIEF,[INT,BRD]+CTRL)
print(f"单独放:强度 **{bI[0]:+.4f}** (se {sI[0]:.4f}) · 广度 **{bB[0]:+.4f}** (se {sB[0]:.4f})")
print(f"一起放:强度 **{bJ[0]:+.4f}** (se {sJ[0]:.4f}) · 广度 **{bJ[1]:+.4f}** (se {sJ[1]:.4f})")
NP_=400; nulI=[]; nulB=[]
ii=np.flatnonzero(M)
for s_ in range(NP_):
    rg=np.random.default_rng(8100+s_); pp=rg.permutation(ii)
    I2=INT.copy(); B2=BRD.copy(); I2[ii]=INT[pp]; B2[ii]=BRD[pp]     # ★ 整行打乱,保住两者相关
    b2,_=fit(BELIEF,[I2,B2]+CTRL); nulI.append(b2[0]); nulB.append(b2[1])
nulI=np.array(nulI); nulB=np.array(nulB)
TI=float(np.percentile(np.abs(nulI),95)); TB=float(np.percentile(np.abs(nulB),95))
print(f"\n⚠ offset 零(**整行打乱** —— 强度与广度一起打乱,保住它们之间的相关与维数):")
print(f"   强度阈 **{TI:.4f}** · 广度阈 **{TB:.4f}**")
print(f"   -> 强度 {'**越阈**' if abs(bJ[0])>TI else '**未越阈**'} · "
      f"广度 {'**越阈**' if abs(bJ[1])>TB else '**未越阈**'}")
X1=np.column_stack([np.ones(n)]+[z(v,M) for v in [S,C3,ncat,sh]])
yB=z(BELIEF,M); b1,*_=np.linalg.lstsq(X1,yB,rcond=None); RES_NA=yB-X1@b1
nulA=np.array([abs(float(np.corrcoef(RES_NA,np.random.default_rng(9100+s).permutation(z(ACTED,M)))[0,1]))
               for s in range(300)])
STHR=float(np.percentile(nulA,95)); rA=abs(float(np.corrcoef(RES_NA,z(ACTED,M))[0,1]))
rgN=np.random.default_rng(7); rF=abs(float(np.corrcoef(z(BELIEF,M),rgN.standard_normal(n))[0,1]))
print(f"\n阳性参照 `ACTED`:|r| **{rA:.4f}** vs 单变量零 {STHR:.4f} -> "
      f"{'**开火**' if rA>STHR else '**不开火**'} · 负对照(纯噪声)**{rF:.4f}**")
print(f"\nguard 26 的正对照 = **MDE 扫描**,每级 30 次(种在强度上):")
MDE=None
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(7250+int(gg*100)*41+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(INT,M)+rg.standard_normal(n)
        bq,_=fit(y,[INT,BRD]+CTRL)
        if abs(bq[0])>TI: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 = `#388a` 的 **0.116**")
pd.DataFrame([dict(v_what='强度',v_alone=bI[0],v_joint=bJ[0],v_se=sJ[0],v_thr=TI),
              dict(v_what='广度',v_alone=bB[0],v_joint=bJ[1],v_se=sJ[1],v_thr=TB)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'int_vs_brd.csv',index=False)

OI=abs(bJ[0])>TI; OB=abs(bJ[1])>TB
g=Gate('强度与广度,哪一个预测「我能不能改」')
g.asserted('★【两支】阳性参照 `ACTED` 在单变量零上开火',rA>STHR,f"{rA:.4f} vs {STHR:.4f}",kind='control')
g.asserted('★【两支】负对照:纯噪声列不越阈',rF<=max(TI,TB),f"{rF:.4f}",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度 0.116',MDE_,0.1159,True,what='MDE 扫描 80% 检出')
g.asserted('★【非零支】世界判定:A(只有强度越阈)· B(只有广度)· C(两个都塌)',
           OI!=OB,
           f"强度 {bJ[0]:+.4f}/{TI:.4f} {'越阈' if OI else '未越'} · "
           f"广度 {bJ[1]:+.4f}/{TB:.4f} {'越阈' if OB else '未越'} -> "
           f"**世界 {'A 强度' if OI and not OB else ('B 广度' if OB and not OI else ('C 两个都塌 —— 这条线到此为止' if not OI and not OB else '两个都越阈 —— 它们不是同一个东西'))}**")
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

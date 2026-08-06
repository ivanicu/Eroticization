import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A152 R444 -- 「灵活性」是 P1/P2 的共同因子,还是第三件事

`#399c` 的假设:P1(吸收对方的兴奋)与 P2(不需要对方的兴奋)都在描述
**「在不是自己默认的情形下也能被唤起」** —— 而 `corr(P1,P2) = 0.031` 说明**没有共同因子把它们绑在一起**。
**所以这个假设必须被测,而不是被读出来。**

灵活性的一个**可操作**版本:同一个人在**支配**与**臣服**两端**都**高。
`FLEX_min = min(z(支配), z(臣服))` —— **两端都要高才算高**(这是「灵活」的定义,不是「强度」)。

⚠ **两个判据,不是一个**(`#399b` 的教训:**名字要有它自己的判据**):
**① `FLEX` 自己越过单变量零** -> 它与「治疗性」有关;
**② 且它吃掉 P1/P2 的份额越过 offset 地板** -> 它是**它们的共同因子**。
**只有 ① 而无 ② -> 「灵活性」是第三件事,不是共同因子。**

⚠ `min()` 是**非线性**的 -> **同轮报 `min` 与 `均值` 两种构造**;
**若结论只在其中一种成立,那是构造在说话,不是数据。**

ESTIMAND        ① `治疗性 ~ FLEX + 羞耻 + S + c3⁻ + 类别数 + ACTED`;
                ② 加 `FLEX` 前后 P1/P2 系数的掉幅。
⚠ 零的种类     ① 用**单变量零**(`FLEX` 是预先指定的,不是从候选里挑的,`#388b`);
                ② 用 `offset_control` —— **掉幅的零绝不是零**(任何相关控制都吃一点),
                零 = `lib.nulls.sham_control`(与 `FLEX` 相关度相同、与结局无关)后的掉幅分布。
IMPOSSIBLE      ① 支配/臣服两题自报且互斥语义 -> `min` 会被**两端都低**的人拉低,而那不是「不灵活」是「都不感兴趣」;
                ② 同一份数据;③ 「共同因子」只在**线性**意义上被否定。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, sham_control, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
TC=next(c for c in d.columns if 'vmq8jqw' in c); THER=pd.to_numeric(d[TC],errors='coerce').values.astype(float)
AC=next(c for c in d.columns if '41kpfir' in c); ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
P1=pd.to_numeric(d[next(c for c in d.columns if 'If my partner is aroused by something' in str(c))],errors='coerce').values.astype(float)
P2=pd.to_numeric(d[next(c for c in d.columns if 'even if I know they' in str(c))],errors='coerce').values.astype(float)
DOM=pd.to_numeric(d[next(c for c in d.columns if '6w3xquw' in str(c))],errors='coerce').values.astype(float)
SUB=pd.to_numeric(d[next(c for c in d.columns if 'xem7hbu' in str(c))],errors='coerce').values.astype(float)
anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
print("⚠ **`#392e`:两道题进模型前先各自看清楚**")
for nm,v in (('支配',DOM),('臣服',SUB)):
    g=np.isfinite(v)&np.isfinite(anc)
    print(f"   {nm}:取值 {np.unique(v[np.isfinite(v)]).tolist()} · "
          f"众数 **{float(pd.Series(v[np.isfinite(v)]).mode().iloc[0]):g}** · 与锚相关 **{np.corrcoef(v[g],anc[g])[0,1]:+.4f}**")
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(THER)&np.isfinite(ACTED)&np.isfinite(P1)&np.isfinite(P2)&np.isfinite(DOM)&np.isfinite(SUB)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
zd=np.full(NN,np.nan); zs=np.full(NN,np.nan); zd[M]=z(DOM,M); zs[M]=z(SUB,M)
FMIN=np.full(NN,np.nan); FMIN[M]=np.minimum(zd[M],zs[M])
FAVG=np.full(NN,np.nan); FAVG[M]=(zd[M]+zs[M])/2
print(f"\nn=**{n:,}** · corr(支配, 臣服) = **{np.corrcoef(DOM[M],SUB[M])[0,1]:+.4f}** · "
      f"corr(FLEX_min, FLEX_avg) = **{np.corrcoef(FMIN[M],FAVG[M])[0,1]:+.4f}**")
CTRL=[sh,S,C3,ncat,ACTED]
def fit(y,xs):
    X=np.column_stack([np.ones(n)]+[z(v,M) for v in xs]); yy=z(y,M)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(n-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
rows=[]
for nm,F in (('FLEX_min',FMIN),('FLEX_avg',FAVG)):
    bF,sF=fit(THER,[F]+CTRL)
    nul=np.array([fit(perm_in(THER,M,5300+s),[F]+CTRL)[0][0] for s in range(400)])
    TF=float(np.percentile(np.abs(nul),95))
    bP,_=fit(THER,[P1,P2]+CTRL); bPF,_=fit(THER,[P1,P2,F]+CTRL)
    d1=100*(1-abs(bPF[0])/max(abs(bP[0]),1e-12)); d2=100*(1-abs(bPF[1])/max(abs(bP[1]),1e-12))
    rE=float(np.corrcoef(z(F,M),z(P1,M))[0,1])
    fl=np.array([100*(1-abs(fit(THER,[P1,P2,sham_control(P1,rE,M,6100+s)]+CTRL)[0][0])/max(abs(bP[0]),1e-12))
                 for s in range(200)])
    TD=float(np.percentile(fl,95))
    rows.append(dict(v_ctor=nm,v_b=float(bF[0]),v_se=float(sF[0]),v_thr=TF,
                     v_over=bool(abs(bF[0])>TF),v_dropP1=d1,v_dropP2=d2,v_dropthr=TD,
                     v_eats=bool(d1>TD or d2>TD)))
    print(f"\n**{nm}**:")
    print(f"   ① 自己:**{bF[0]:+.4f}** (se {sF[0]:.4f}) vs 单变量零 **{TF:.4f}** · "
          f"{'**越阈**' if abs(bF[0])>TF else '**未越阈**'}")
    print(f"   ② 吃掉 P1 **{d1:+.2f}%** · P2 **{d2:+.2f}%** vs offset 地板 95 分位 **{TD:+.2f}%** · "
          f"{'**越阈(是共同因子)**' if (d1>TD or d2>TD) else '**未越阈(不是共同因子)**'}")
T=pd.DataFrame(rows); check_columns(T,'R444')
T.to_csv(pathlib.Path(__file__).parent/'results'/'flex.csv',index=False)
negs=np.array([fit(perm_in(THER,M,80000+s),[FMIN]+CTRL)[0][0] for s in range(200)])
TF0=float(T.v_thr.iloc[0]); rate=float((np.abs(negs)>TF0).mean())
print(f"\n负对照(**越阈率**):**{100*rate:.1f}%**")
MDE=None
print(f"guard 26 = **MDE 扫描**,每级 30 次(种在 FLEX_min 上):")
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(900+int(gg*100)*71+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(FMIN,M)+rg.standard_normal(n)
        if abs(fit(y,[FMIN]+CTRL)[0][0])>TF0: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
CONS=bool(T.v_over.iloc[0]==T.v_over.iloc[1] and T.v_eats.iloc[0]==T.v_eats.iloc[1])
g=Gate('「灵活性」是 P1/P2 的共同因子,还是第三件事')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,max(abs(T.v_b)),True,what='MDE 扫描 80% 检出')
g.asserted('★【两支】两种构造(`min` / 均值)结论一致 —— 否则是构造在说话',CONS,
           f"越阈 {T.v_over.tolist()} · 吃掉 {T.v_eats.tolist()}",kind='control')
if 0.01<=rate<=0.12 and CONS:
    g.asserted('★【非零支】判据①:`FLEX` 自己越过单变量零',bool(T.v_over.all()),
               f"{T.v_b.tolist()} vs {T.v_thr.tolist()}")
    g.asserted('★【非零支】判据②:且它**吃掉** P1/P2 -> 才叫共同因子',bool(T.v_eats.all()),
               f"掉幅 P1 {T.v_dropP1.tolist()} · P2 {T.v_dropP2.tolist()} vs 地板 {T.v_dropthr.tolist()}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

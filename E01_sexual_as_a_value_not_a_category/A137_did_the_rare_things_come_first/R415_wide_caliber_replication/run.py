import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A137 R415 -- `#370c` 的跨人群复制:宽口径,同轮带窄口径参照臂

`#370c`:`ord ↔ 羞耻` = **−0.032**,越阈,但**低于自己的 MDE(0.040)** —— 强度不可报。
欠一次复制,而**正确形式不是重跑**(同一份数据同一个设计,那是同一个数)。

ESTIMAND        `#370c` 的同一个系数,在**放宽两个口径**后重估:
                块覆盖 **≥8 -> ≥4**(`CALIBER.md` 旋钮 ⑩)· `ord_i` 的类别数下限 **≥8 -> ≥6**。
                **同一次运行里带一条窄口径参照臂**(⑩ 的「必须同时报」条款)。
⚠ 这不是加功率  `#346b` 已证明宽口径与窄口径是**不同的人群**,不是同一群人加数据。
                **报告必须写成跨人群复制。**
KILL(条件式)  仅当对照过**且宽臂 MDE < |−0.032|** -> 判:宽臂系数是否同号且越阈。
                同号越阈 -> `#370c` 跨人群复制成功,可以上页面;
                不同号或未越阈 -> `#370c` 留在账上,不上页面。
POSITIVE / NEGATIVE  沿用 `#370d`(种植 0.12 / 打乱人),两臂各跑一次。
⚠ 顺序          **先算两臂各自的 MDE 再看系数**(`#369a`)。
                若宽臂 MDE 仍 > 0.032 -> **说这个设计还是看不见它**,不要因为符号一致就升级。
IMPOSSIBLE      放宽 `ord_i` 的下限会让秩相关更不稳(类别少的人噪声大),
                所以宽臂的**每个 `ord_i` 更差**,而**人更多** —— 两者反向,净效应事前不知道。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from scipy.stats import rankdata
_R414=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R414_within_person_order/run.py').read_text()
_BODY=_R414.split('"""',2)[2]
exec(_BODY.split('MINC=8')[0])          # 到 O · RAR · ncat · sh · ok · fit_apply 为止
COVB=cov.copy()                          # 块覆盖数(R333 定义 ok=cov>=8)
# ⚠ **guard 25 抓到的第一版缺陷**:第一版用 `S=fit_apply(...)`,而 `fit_apply` 内部硬编码
# `np.where(cv>=8, …)` —— 于是掩码里的 `isfinite(S)` **把 ≥8 又悄悄加了回来**,
# 「宽臂」只从 6,473 涨到 6,543(+70),而 `#350` 的同一个旋钮该给 ~8,000。
# **一个没有变宽的口径臂不是另一个人群的复制,而是同一群人换了一个更宽的标签** —— 讨好方向的 scope 主张。
# 修法:`S` 按**每条臂自己的覆盖阈**重建。
def make_S(K):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nb=M.sum(1)
        v=np.where(nb>0,(M@rr)/np.maximum(nb,1),np.nan)
        gq=np.isfinite(v); cv[ppl[gq]]+=1; ps[ppl[gq]]+=v[gq]
    return np.where(cv>=K,ps/np.maximum(cv,1),np.nan)
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)

def ordvec(Om,minc):
    out=np.full(NN,np.nan); F=np.isfinite(Om)
    for i in np.flatnonzero(F.sum(1)>=minc):
        j=np.flatnonzero(F[i]); a=Om[i,j]; b=RAR[j]
        if np.std(a)<1e-9 or np.std(b)<1e-9: continue
        out[i]=np.corrcoef(rankdata(a),rankdata(b))[0,1]
    return out

ARMS=[('窄(参照臂:覆盖≥8 · 类别≥8)',8,8),('宽(覆盖≥4 · 类别≥6)',4,6)]
RES={}
for label,mincov,minc in ARMS:
    okA=COVB>=mincov
    S=make_S(mincov)                     # ★ 每条臂自己的尺子口径
    OV=ordvec(O,minc)
    # ⚠ S 是在窄口径上估的尺子,两臂共用同一把尺(#293:换尺子不解决问题,但换尺子会换掉被测的东西)
    mA=okA&np.isfinite(OV)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(S)
    n=int(mA.sum())
    def coef(y,o=OV,g=mA):
        gg=g&np.isfinite(o)&np.isfinite(y); k=int(gg.sum())
        X=np.column_stack([np.ones(k),z(o,gg),z(ncat,gg),z(S,gg)])
        b,*_=np.linalg.lstsq(X,z(y,gg),rcond=None); r=z(y,gg)-X@b
        s2=float(r@r)/(k-4); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
        return float(b[1]),float(se[1]),k
    NUL=np.array([coef(sh[np.random.default_rng(3000+s).permutation(NN)])[0] for s in range(300)])
    THR=float(np.percentile(np.abs(NUL),95))
    MDE=None
    for g_ in (0.020,0.025,0.030,0.035,0.040,0.060):
        hit=0
        for s_ in range(40):
            rg=np.random.default_rng(9000+int(g_*1000)*7+s_)
            y=np.full(NN,np.nan); y[mA]=g_*z(OV,mA)+rg.standard_normal(n)
            if abs(coef(y)[0])>THR: hit+=1
        if MDE is None and hit>=32: MDE=g_
    MDE_=MDE if MDE else 0.10
    b,se,k=coef(sh)
    rg=np.random.default_rng(55)
    yp=np.full(NN,np.nan); yp[mA]=0.12*z(OV,mA)+rg.standard_normal(n)
    bp=coef(yp)[0]; bn=coef(sh[np.random.default_rng(99).permutation(NN)])[0]
    RES[label]=dict(n=k,b=b,se=se,thr=THR,mde=MDE_,bp=bp,bn=bn,
                    ordmean=float(np.nanmean(OV[mA])),nulsd=float(NUL.std()),nulmean=float(NUL.mean()))
    print(f"{label}")
    print(f"   n=**{k:,}** · mean(ord)= {np.nanmean(OV[mA]):+.4f} · "
          f"零 {NUL.mean():+.5f} ± {NUL.std():.5f} · 阈 {THR:.5f} · **MDE {MDE_:.3f}**")
    print(f"   系数 **{b:+.5f}** · se {se:.5f} · 95% CI [{b-1.96*se:+.5f}, {b+1.96*se:+.5f}] · "
          f"距零 **{(b-NUL.mean())/max(NUL.std(),1e-12):+.2f} sd** · "
          f"{'**越阈**' if abs(b)>THR else '**未越阈**'}")
    print(f"   正对照 {bp:+.5f} · 负对照 {bn:+.5f}\n")

T=pd.DataFrame([dict(v_arm=k,**{('v_'+kk):vv for kk,vv in v.items()}) for k,v in RES.items()])
check_columns(T,'R415'); T.to_csv(pathlib.Path(__file__).parent/'results'/'arms.csv',index=False)
NA,WA=ARMS[0][0],ARMS[1][0]
nn,ww=RES[NA],RES[WA]
print(f"⚠ **跨人群,不是加功率**(`#346b`):窄臂 n={nn['n']:,} -> 宽臂 n={ww['n']:,};"
      f"两臂的 mean(ord) 是 {nn['ordmean']:+.4f} vs {ww['ordmean']:+.4f} —— "
      f"{'差别可见,人群确实不同' if abs(nn['ordmean']-ww['ordmean'])>0.01 else '两臂的 ord 水平接近'}")
SAME=np.sign(nn['b'])==np.sign(ww['b'])
print(f"⇒ 同号 **{SAME}** · 宽臂 {'越阈' if abs(ww['b'])>ww['thr'] else '**未越阈**'} · "
      f"宽臂 MDE {ww['mde']:.3f} vs |窄臂系数| {abs(nn['b']):.3f} -> "
      f"{'**MDE 够小,可以判**' if ww['mde']<abs(nn['b']) else '**MDE 仍不够小 —— 这个设计还是看不见它**'}")

g=Gate('#370c 的跨人群复制')
CP=all(abs(v['bp'])>v['thr'] for v in RES.values())
CN=all(abs(v['bn'])<=v['thr'] for v in RES.values())
g.asserted('★ 正对照:两臂各种植 0.12 -> 都必须越阈',CP,
           ' · '.join(f"{k[:2]} {v['bp']:+.4f}/{v['thr']:.4f}" for k,v in RES.items()),kind='control')
g.asserted('★ 负对照:两臂各打乱人 -> 都必须是零',CN,
           ' · '.join(f"{k[:2]} {v['bn']:+.4f}" for k,v in RES.items()),kind='control')
g.asserted('★ ⑩ 的强制条款:同一次运行里带了窄口径参照臂',NA in RES,
           f"窄臂 n={nn['n']:,} · 宽臂 n={ww['n']:,}",kind='control')
g.relaxation_reached_the_population('★ guard 25:放宽的口径真的到达了人群',
                                    nn['n'],ww['n'],what='覆盖 ≥8 -> ≥4 · 类别 ≥8 -> ≥6')
if CP and CN:
    if ww['mde']<abs(nn['b']):
        g.asserted('★ 注册的 kill:宽臂同号且越阈 -> 复制成功',
                   SAME and abs(ww['b'])>ww['thr'],
                   f"宽臂 {ww['b']:+.5f} vs 阈 {ww['thr']:.5f} · 同号 {SAME}")
    else:
        g.asserted('★ 注册的 kill(宽臂 MDE ≥ |窄臂系数| -> 不判,不因符号一致而升级)',False,
                   f"宽臂 MDE {ww['mde']:.3f} ≥ {abs(nn['b']):.3f}")
else:
    g.asserted('★ 注册的 kill(对照未过 -> 不判)',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")

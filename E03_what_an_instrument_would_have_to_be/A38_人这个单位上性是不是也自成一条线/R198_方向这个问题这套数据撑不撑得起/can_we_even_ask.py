"""#757 · E03·A38·R198 —— 「方向」这个问题,这套数据撑不撑得起?

三具仪器都说:管性与要求服从是同一件事的两侧(`#753`·`#756`)。**下一个问题必然是方向** ——
是要求服从的人才管别人的性,还是管性的人才要求服从?

⚠ **而本轮的任务不是回答它,是先量这个设计有没有资格问它。**
`#746`① 定下的规矩:**先算 MDE,把它用在选设计上,而不是跑完再解释。**

⚠ **本轮结构上只此一具仪器,而这句话是量出来的**:本地 `data/external/` 下八个数据源
(brfss · dataverse/MFQ · dplace/SCCS · gss · ngram · nsfg · openpsych · yrbs)**全部是横截面,
没有任何面板** —— 而方向需要面板。**换不了仪器**,不是没去找。

**为什么人层不可识别**:GSS 是重复横截面,不是面板 —— 同一个人只被问一次,
**没有任何个体的时间先后**。⇒ 人层的方向**结构上不可得**,不是「计划中」(realstat §2)。

**唯一可识别的是年代层的领先-滞后**:`obey` 与 `premarsx` 都有 1986–2024 的年度均值。
⚠ 但两条序列都在强趋势上,**两条趋势序列的互相关由趋势主导,是伪的** ⇒ 必须先去趋势。
⚠ 而去趋势之后剩下的是噪声 + 抽样误差,**所以真正的问题是:23 个时点的残差,
   能不能分辨 lag −1 与 lag +1 的差别?** 这就是本轮要算的 MDE。

**预注册:**
  若在 80% 功效下能分辨的最小领先-滞后不对称 **小于** 一个合理的效应(比如 0.3 的互相关差),
  ⇒ 设计够用,**下一轮跑它**;
  若最小可分辨不对称 **大于** 1.0(即互相关差要超过理论上限才检得出),
  ⇒ **设计结构上撑不起,把「方向」写进页上的「这套数据做不到什么」,并且不许跑那个弱检验**;
  之间 ⇒ 报出这个数,让它自己说话。
"""
import pandas as pd, numpy as np, json, pathlib, sys
RNG=np.random.default_rng(198)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
d=pd.read_stata(ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta",
                columns=["year","obey","premarsx","homosex"],convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(lambda x:x>0) for c in ("obey","premarsx","homosex")})
M["homosex"]=M["homosex"].where(M["homosex"]<=4); M["year"]=d.year
Y=sorted(set(M.loc[M.obey.notna(),"year"]) & set(M.loc[M.premarsx.notna(),"year"]))
print(f"=== 硬规则①:两条序列共同的年份 {len(Y)} 个:{[int(y) for y in Y]} ===")
g=M[M.year.isin(Y)].groupby("year")
S=pd.DataFrame({"obey":g.obey.mean(),"premarsx":g.premarsx.mean(),"homosex":g.homosex.mean(),
                "n_obey":g.obey.count(),"n_pre":g.premarsx.count()})
S["se_obey"]=g.obey.std()/np.sqrt(S.n_obey); S["se_pre"]=g.premarsx.std()/np.sqrt(S.n_pre)
print(f"\n{'年':>6s} {'obey均值':>9s} {'se':>7s} {'premarsx均值':>12s} {'se':>7s}")
for y,r in S.iterrows(): print(f"{int(y):6d} {r.obey:9.4f} {r.se_obey:7.4f} {r.premarsx:12.4f} {r.se_pre:7.4f}")

t=np.arange(len(S))
def detrend(v):
    A=np.c_[np.ones(len(t)),t]; return v-A@np.linalg.lstsq(A,v,rcond=None)[0]
ro,rp=detrend(S.obey.to_numpy()),detrend(S.premarsx.to_numpy())
print(f"\n=== 去趋势 ===")
print(f"  obey     趋势斜率 {np.polyfit(t,S.obey,1)[0]:+.5f}/波 · 残差 sd {ro.std(ddof=1):.4f} · 平均抽样 se {S.se_obey.mean():.4f}")
print(f"  premarsx 趋势斜率 {np.polyfit(t,S.premarsx,1)[0]:+.5f}/波 · 残差 sd {rp.std(ddof=1):.4f} · 平均抽样 se {S.se_pre.mean():.4f}")
print(f"  ⚠ 残差 sd 与抽样 se 的比:obey {ro.std(ddof=1)/S.se_obey.mean():.2f}× · premarsx {rp.std(ddof=1)/S.se_pre.mean():.2f}×")
print(f"     -> 若接近 1,**去趋势后剩下的几乎全是抽样噪声,没有可供领先-滞后分析的信号**")

def xcorr(a,b,lag):
    if lag>0: a2,b2=a[lag:],b[:-lag]
    elif lag<0: a2,b2=a[:lag],b[-lag:]
    else: a2,b2=a,b
    return float(np.corrcoef(a2,b2)[0,1])
print(f"\n=== 观测的互相关(去趋势后)===")
obs={L:xcorr(rp,ro,L) for L in (-2,-1,0,1,2)}
for L,v in obs.items(): print(f"  lag {L:+d}: {v:+.4f}   ({'obey 领先' if L>0 else 'premarsx 领先' if L<0 else '同期'})")
asym=obs[1]-obs[-1]
print(f"  不对称 (lag+1 − lag−1) = {asym:+.4f}")

# --- MDE:在「无领先-滞后」的零下模拟,看这个设计能分辨多大的不对称 ---
NS=4000
def sim(rho_lead=0.0):
    """两条序列 = 共同 AR(1) 成分 + 各自噪声;rho_lead 注入一期领先。"""
    n=len(S); phi=0.5
    e=RNG.normal(size=n+2); c=np.zeros(n+2)
    for i in range(1,n+2): c[i]=phi*c[i-1]+e[i]
    a=c[2:]+RNG.normal(0,1,n)          # premarsx 残差
    b=(1-abs(rho_lead))*RNG.normal(0,1,n)+rho_lead*np.r_[c[1:-1]]  # obey 残差,可提前一期
    return detrend(a),detrend(b)
null=[ (lambda A,B: xcorr(A,B,1)-xcorr(A,B,-1))(*sim(0.0)) for _ in range(NS) ]
q=float(np.quantile(np.abs(null),.95))
print(f"\n=== MDE:无领先-滞后的零下,不对称的 95% 分位 = {q:.4f}(n={len(S)} 个时点,{NS} 次)===")
grid=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
print(f"  {'注入的领先强度':>12s} {'检出率(|不对称|>零95%)':>22s}")
mde=None
for g_ in grid:
    p=np.mean([abs((lambda A,B: xcorr(A,B,1)-xcorr(A,B,-1))(*sim(g_)))>q for _ in range(800)])
    print(f"  {g_:12.2f} {p:22.3f}")
    if mde is None and p>=0.80: mde=g_
print(f"\n  80% 功效需要的最小领先强度 = {mde if mde else '>0.8 —— 网格内达不到'}")
print("\n"+"="*66)
if mde is None or mde>0.7:
    v=(f"**结构上撑不起:23 个时点、去趋势后残差 sd/抽样 se = {rp.std(ddof=1)/S.se_pre.mean():.2f}×,"
       f"要 80% 功效需要领先强度 {'>0.8' if mde is None else mde} —— **「方向」写进「这套数据做不到什么」,不跑那个弱检验**")
elif mde<=0.3: v=f"**设计够用:80% 功效只需领先强度 {mde} ⇒ 下一轮跑它**"
else: v=f"**边缘:80% 功效需要领先强度 {mde} —— 报出这个数,不假装它是零或有**"
print(v)
json.dump(dict(years=[int(y) for y in Y],obs_xcorr=obs,asym=asym,null_q95=q,mde=mde,
    resid_sd_over_se=dict(obey=float(ro.std(ddof=1)/S.se_obey.mean()),premarsx=float(rp.std(ddof=1)/S.se_pre.mean())),
    verdict=v),open(OUT/"can_we_even_ask.json","w"),ensure_ascii=False,indent=1)

"""E03·A19·R101 —— 一件事变得更少见,它就更被谴责吗

**类型:FRONTIER**。**这一页开场那句话从来没有被测过。**
E02 的立论是:**「E01 的 corr(rarity, shame)=0.758 已从一个机制降级为『两者都是谴责的下游』」**。
**那次降级本身,从未在任何单位上被直接检验。** 而 GSS 同时有**谴责**与**行为发生率**,同一具仪器。

⚠ **BASIN**:连着五轮都在测我自己的仪器 —— 这是盆地。
  **本轮下注 W2/W3**,即**下注这一页开场当成已定的那句话是未经检验就采用的**。

W1 **降级对了** —— 两条序列的一阶差分不相关 ⇒ 在年代这个单位上,稀有度不驱动谴责。
W2 **降级下早了** —— 差分正相关 ⇒ 一件事变少的那些年,它也更被谴责。
W3 **反向** —— 谴责的变化**领先**发生率的变化 ⇒ 谴责在上游,这才是那句话的**强形式**,
   而它会是一个**发现**,不只是一个零。
**Meta-separator**:若在年代这个单位上两条序列是**同一条序列加噪声**,
  **那么「稀有度」与「谴责」根本不是两个变量**,我的分解就是错的 —— 由「同期相关 ≈ 1」抓。

## 硬规则①(已跑,写在这里)
`xmarsex` n=46,266 · **30 年 1973–2024** · `evstray` n=34,737 · **17 年 1991–2022**。
**共同年份 17 个,联合 n=20,408,每年都 >=200。** ⇒ **16 个一阶差分。**
`evstray` 码:`1=yes` `2=no` **`3=never married`** —— **码 3 必须剔除**,他们从不在风险集里。

## G1 ESTIMAND(先于方法)
逐年两条序列:**谴责** = `xmarsex` 答「always wrong」的份额;**发生率** = 已婚过者里 `evstray=1` 的份额。
主量 = **两条序列一阶差分的 Pearson 相关**(差分,不是水平 —— 水平会被共同趋势制造出来,
本页第七件已量过「看起来像共同气候的东西大半是一条共同趋势」)。
## G2 CONTROLS
**正对照**:`premarsx` 与 `homosex` 的差分必须正相关(两道题已知同向变软)。
**安慰剂**:`xmarsex` 差分 × `cappun`(死刑)差分 —— 与性无关的道德序列,必须 ≈0。
**⚠ 而零该不该是零?** 安慰剂**该**趋零(两件无关的事),用 `negative_control`。
**MDE 必须报** —— n=16 的零是弱的,不报 MDE 的零只是沉默(本页已记过这条)。
## G3/G4
网格:{全体已婚过 / 仅当前已婚} × {谴责=份额 / 谴责=均值} × {同期 / 谴责领先 1 期 / 发生率领先 1 期}。
## KILL(条件式,预注册)
if 正对照>0 and |安慰剂| < 0.5*|主量|:
  同期差分相关的 bootstrap 区间**含零 且 |r| < MDE** -> W1(而必须同时说「只排除了大于 MDE 的效应」)
  区间不含零且为正 -> W2 · 领先项显著而同期不显著 -> W3
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
16 个差分 ⇒ **功率极低,这是结构性的** · 无干预 · 自报 · `evstray` 的分母是「已婚过」,
**随离婚率变动**(网格里控制) · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import pearsonr, spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
P="data/external/gss/GSS_stata/gss7224_r3a.dta"
df,_=pyreadstat.read_dta(P, usecols=["year","evstray","xmarsex","premarsx","homosex","cappun","marital"], encoding="latin1")

def series(mask=None):
    d=df if mask is None else df[mask]
    g=d[d.evstray.isin([1,2])].groupby("year")
    prev=g.apply(lambda x:(x.evstray==1).mean()); n1=g.size()
    g2=d[d.xmarsex.isin([1,2,3,4])].groupby("year")
    cond=g2.apply(lambda x:(x.xmarsex==1).mean()); n2=g2.size()
    m=pd.DataFrame(dict(prev=prev,cond=cond,n1=n1,n2=n2)).dropna()
    return m[(m.n1>=200)&(m.n2>=200)]

S=series()
print("=== 逐年两条序列(硬规则①:先打印)===")
print(f"  {'年':>6s}{'发生率':>9s}{'谴责':>9s}{'n(行为)':>9s}{'n(态度)':>9s}")
for y,r in S.iterrows(): print(f"  {int(y):>6d}{r.prev:>9.4f}{r.cond:>9.4f}{int(r.n1):>9d}{int(r.n2):>9d}")
print(f"  年数 = **{len(S)}** ⇒ 差分 **{len(S)-1}** 个")

dp=np.diff(S.prev.values); dc=np.diff(S.cond.values)
r_same=float(pearsonr(dp,dc).statistic); rs=float(spearmanr(dp,dc).statistic)
lvl=float(pearsonr(S.prev.values,S.cond.values).statistic)
print(f"\n=== 主量 ===\n  同期差分相关 r = **{r_same:+.4f}** (Spearman {rs:+.4f}) · n = {len(dp)}")
print(f"  ⚠ 对照用的**水平**相关 = {lvl:+.4f} —— 差分与水平的差距就是共同趋势那一份")
print(f"  meta-separator:水平相关若 ≈1 则两者不是两个变量 -> |{lvl:.3f}| {'<' if abs(lvl)<0.9 else '>='} 0.9")

def boot(x,y,n=2000):
    rng=np.random.default_rng(SEEDS[0]); out=[]
    for _ in range(n):
        i=rng.integers(0,len(x),len(x))
        if len(set(i))<3: continue
        v=pearsonr(x[i],y[i]).statistic
        if np.isfinite(v): out.append(v)
    return np.array(out)
bs=boot(dp,dc); lo,hi=np.quantile(bs,[.025,.975])
print(f"  95% CI = [{lo:+.4f}, {hi:+.4f}]  -> {'**含零**' if lo<0<hi else '**不含零**'}")
# MDE:n=16 时 80% 功率能测到多大的 |r|
from scipy.stats import norm
n=len(dp); z=(norm.ppf(0.975)+norm.ppf(0.80))/np.sqrt(n-3)
mde=float(np.tanh(z))
print(f"  **MDE(n={n}, 80% 功率, 双侧 0.05)= |r| = {mde:.4f}** —— 小于它的效应本设计看不见")

print("\n=== G4 网格 ===")
grid={}
for tag,mask in [("全体已婚过",None),("仅当前已婚",df.marital==1)]:
    T=series(mask)
    if len(T)<6: grid[tag]=None; print(f"  {tag}: 年数 {len(T)} 不足"); continue
    a=np.diff(T.prev.values); b=np.diff(T.cond.values)
    row={"同期":float(pearsonr(a,b).statistic)}
    row["谴责领先1期"]=float(pearsonr(a[1:],b[:-1]).statistic)
    row["发生率领先1期"]=float(pearsonr(a[:-1],b[1:]).statistic)
    grid[tag]=row
    print(f"  {tag:10s} n={len(a):2d} · " + " · ".join(f"{k} {v:+.4f}" for k,v in row.items()))

print("\n=== 控制 ===")
def dseries(col, val=1):
    g=df[df[col].isin([1,2,3,4])].groupby("year")
    s=g.apply(lambda x:(x[col]==val).mean()); nn=g.size()
    s=s[nn>=200]
    return s
pm,hm=dseries("premarsx"),dseries("homosex")
cy=sorted(set(pm.index)&set(hm.index))
pc=float(pearsonr(np.diff(pm[cy].values),np.diff(hm[cy].values)).statistic)
cp=dseries("cappun"); cy2=sorted(set(S.index)&set(cp.index))
pl=float(pearsonr(np.diff(S.cond[cy2].values),np.diff(cp[cy2].values)).statistic) if len(cy2)>=6 else np.nan
print(f"  正对照 premarsx×homosex 差分 r = **{pc:+.4f}** (n={len(cy)-1}) · 必须 > 0")
print(f"  安慰剂 xmarsex×cappun 差分 r = **{pl:+.4f}** (n={len(cy2)-1}) · 应 ≈ 0")

G=Gate("一件事变得更少见,它就更被谴责吗")
p1=G.positive_control("premarsx×homosex 的差分必须正相关",planted=pc,floor=0.0,spread=0.05)
p2=G.negative_control("安慰剂:谴责婚外性 × 支持死刑,两件无关的道德序列",null=abs(pl),effect=abs(r_same),
                      null_spread=0.05,null_kind="与性无关的道德态度序列,同样做一阶差分")
if p1 and p2:
    if lo<0<hi and abs(r_same)<mde:
        verdict=(f"**W1 —— 差分不相关(r={r_same:+.4f},CI 含零,|r| < MDE {mde:.3f})。"
                 f"⚠ 而这只排除了大于 {mde:.3f} 的效应,没有排除小效应。**")
    elif not (lo<0<hi) and r_same>0: verdict=f"**W2 —— 降级下早了:差分正相关 {r_same:+.4f}**"
    elif not (lo<0<hi): verdict=f"**方向相反:{r_same:+.4f}**"
    else: verdict=f"**判不了 —— CI 含零但 |r|={abs(r_same):.3f} >= MDE {mde:.3f},功率不足以称零**"
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(years=[int(x) for x in S.index],prev=S.prev.tolist(),cond=S.cond.tolist(),
               r_diff=r_same,r_spearman=rs,r_level=lvl,ci=[float(lo),float(hi)],mde=mde,n_diff=int(n),
               grid=grid,positive=pc,placebo=pl,verdict=verdict,unchallenged=True),
          open(OUT/"condemnation_vs_prevalence.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'condemnation_vs_prevalence.json'}")

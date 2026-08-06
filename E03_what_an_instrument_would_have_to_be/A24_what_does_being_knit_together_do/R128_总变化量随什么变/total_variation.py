"""E03·A24·R128 —— 那个更基本的量:一个人的答案总共变化多少,随什么变

**类型:FRONTIER。** `#685` 的 NEXT。上一轮杀掉「分辨 vs 一刀切」这个二分,
**顺手立起了一个更基本的量:该人答案的总变化量** —— 而它在这条线上从没被当主角看过一眼。

## 硬规则①(已跑),两条停止条件都没触发
**① `corr(总变化量, 答题数) = +0.2089`** ⇒ **不是完成度指标**(预注册的停止条件是 >0.5 当场停)。
**② 题集随年份剧变**(每年可用题数 8–33)⇒ **人内标准差跨年份不可比** ——
**题集不变的最长连续区间 = 1988–2018,18 波,n = 42,829**,本轮**全部计算限制在这个窗口内**。

## G1 ESTIMAND
`totvar` = 该人在 33 题上标准化答案的人内标准差(仅 1988–2018)。**主量 = `ρ(educ, totvar)`,n = 人数。**
按 `#683` 的做法在 **cohort 五分位内**分层重算。
## ⑧ 判据(`#685` 在跑之前写死,不得改)
**cohort 内五格同号 且 中位 ≥ 0.5 × 整体** ⇒ **是教育**;
**任一格变号,或中位 < 0.5 × 整体** ⇒ **记「不是教育」。**
## G2 CONTROLS
**正对照(`#685` 指定)**:该量必须复现 `#682` 的题层结果 ——
**高教育端容忍题的熵更低、性题的熵更高**,否则量错了当场停。
**零**:打乱 `educ`。**这个零该不该是零?** 该 ⇒ `negative_control`。
**混淆控制**:偏相关控掉**答题数**(`nans`),因为答得多的人 sd 估得更稳。
## KILL(条件式)
if 正对照复现 and 观测超过打乱 educ 的零: evaluate(判据⑧) else UNVERIFIED
## IMPOSSIBLE(不写 planned)
1988–2018 之外的年份**结构性无法纳入**(题集不同);
`totvar` 混合了「立场极端」与「跨题不一致」两种来源,**本设计分不开**,只报合量;
跨仪器:MFQ 题集固定但无 cohort;因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
BAT={"性四题":["premarsx","xmarsex","homosex","teensex"],
     "堕胎":["abrape","abhlth","abdefect","abpoor","abnomore","absingle","abany"],
     "自杀":["suicide1","suicide2","suicide3","suicide4"],
     "容忍·言论":["spkath","spkrac","spkcom","spkmil","spkhomo"],
     "容忍·任教":["colath","colrac","colcom","colmil","colhomo"],
     "容忍·藏书":["libath","librac","libcom","libmil","libhomo"],
     "性别角色":["fefam","fepol","fepresch"]}
ALL=sorted({c for v in BAT.values() for c in v})
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ","cohort"]+ALL, apply_value_formats=False, encoding="latin1")
d=df[(df.year>=1988)&(df.year<=2018)&df.educ.notna()&df.cohort.notna()].copy()
Z=pd.DataFrame({c:(d[c]-d[c].mean())/d[c].std() for c in ALL})
d["nans"]=d[ALL].notna().sum(1); d["totvar"]=Z.std(1)
d=d[d.nans>=5].dropna(subset=["totvar"])
print(f"窗口 1988–2018 · 至少答 5 题 · 有 educ 与 cohort ⇒ **n = {len(d):,}**")
def H(s):
    p=s.value_counts(normalize=True).to_numpy(); p=p[p>0]; k=len(p)
    return float(-(p*np.log(p)).sum()/np.log(k)) if k>1 else 0.0
lo=d[d.educ<=12]; hi=d[d.educ>=16]
tol=[c for c in ALL if c[:3] in ("spk","col","lib")]; sex=BAT["性四题"]
dtol=float(np.median([H(hi[c].dropna())-H(lo[c].dropna()) for c in tol]))
dsex=float(np.median([H(hi[c].dropna())-H(lo[c].dropna()) for c in sex]))
print(f"正对照(复现 #682 的题层结果):容忍题 Δ熵中位 **{dtol:+.4f}**(应<0) · "
      f"性题 Δ熵中位 **{dsex:+.4f}**(应>0)  {'✅ 复现' if dtol<0<dsex else '⛔ 不复现 —— 量错了,当场停'}")
def rc(a,b):
    m=np.isfinite(a)&np.isfinite(b)
    return float(np.corrcoef(pd.Series(a[m]).rank(),pd.Series(b[m]).rank())[0,1])
def part(a,b,*ctrl):
    A=pd.Series(a).rank().to_numpy(); B=pd.Series(b).rank().to_numpy()
    C=np.column_stack([pd.Series(c).rank().to_numpy() for c in ctrl]+[np.ones(len(A))])
    ra=A-C@np.linalg.lstsq(C,A,rcond=None)[0]; rb=B-C@np.linalg.lstsq(C,B,rcond=None)[0]
    return float(np.corrcoef(ra,rb)[0,1])
E=d.educ.to_numpy(float); T=d.totvar.to_numpy(float); N=d.nans.to_numpy(float)
r=rc(E,T); pn=part(E,T,N)
rng=np.random.default_rng(20260806)
nul=np.array([abs(rc(rng.permutation(E),T)) for _ in range(300)]); q=float(np.quantile(nul,.95))
print(f"\n整体 ρ(educ, totvar) = **{r:+.4f}** · 控答题数 **{pn:+.4f}** · 打乱 educ 零 95% **{q:.4f}** "
      f"{'✅ 超零' if abs(r)>q else '⛔ 在零里'}")
print(f"\ncohort 五分位内:")
q5=pd.qcut(d.cohort,5,duplicates='drop'); ds=[]
for lv,g in d.groupby(q5,observed=True):
    rr=rc(g.educ.to_numpy(float),g.totvar.to_numpy(float)); ds.append(rr)
    print(f"   {str(lv):22s} n={len(g):>6,}  ρ = **{rr:+.4f}**")
ds=np.array(ds); mid=float(np.median(ds)); same=len(set(np.sign(ds)))==1
print(f"\n   五格中位 **{mid:+.4f}** · 同号 {'✅ 5/5' if same else f'⛔ 只有 {int(max((np.sign(ds)==1).sum(),(np.sign(ds)==-1).sum()))}/5'} · "
      f"占整体 **{mid/r*100:.0f}%** {'✅ ≥50%' if abs(mid)>=0.5*abs(r) else '⛔ <50%'}")
G=Gate("一个人的答案总共变化多少,随什么变")
p1=G.positive_control("复现 #682 题层结果:容忍题熵降、性题熵升",
                      planted=1.0 if (dtol<0<dsex) else 0.0,floor=0.5,spread=0.01)
p2=G.negative_control("打乱 educ 后关系应消失",null=q,effect=abs(r),null_spread=0.005,
                      null_kind="人层打乱 educ —— 若教育与总变化量无关,打乱后应无差别")
if p1 and p2:
    v=(f"**是教育:整体 ρ = {r:+.4f}(控答题数 {pn:+.4f}),cohort 内五格同号,中位 {mid:+.4f} = 整体的 {mid/r*100:.0f}%**"
       if (same and abs(mid)>=0.5*abs(r)) else
       f"**不是教育:五格{'变号' if not same else ''}{'、' if not same and abs(mid)<0.5*abs(r) else ''}"
       f"{'中位仅为整体的 %.0f%%'%(mid/r*100) if abs(mid)<0.5*abs(r) else ''} ⇒ 判据⑧ 命中否定支**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(d)),rho=r,part_nans=pn,null_q95=q,strata=[float(x) for x in ds],
               median_strata=mid,same_sign=bool(same),d_tol=dtol,d_sex=dsex,verdict=v,unchallenged=True),
          open(OUT/"total_variation.json","w"),indent=1,ensure_ascii=False)

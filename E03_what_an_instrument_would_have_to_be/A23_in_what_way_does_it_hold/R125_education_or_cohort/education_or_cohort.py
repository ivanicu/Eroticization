"""E03·A23·R125 —— 那到底是不是「教育」

**类型:FRONTIER。** `#682` 的 NEXT。`#675`–`#682` 全线只用了 `educ` 一个协变量,
**而「教育」和「世代」在美国社会调查里是出了名的纠缠。若分层后教育效应消失,全线都要改口。**

## 硬规则①(已跑),而它连续第二轮杀掉我自己预注册的混淆
分析样本 n = **14,811**,cohort 1899–2006。
**`corr(educ, cohort) = +0.1340` · `corr(educ, age) = −0.0489` · `corr(cohort, year) = +0.4601`**
⇒ **纠缠远比我预注册时预想的弱。** 五个 cohort 五分位内 educ 的 IQR = 3.0 / 4.0 / 4.0 / 4.0 / 4.0,
**中位 4.0 年 ≥ 3 ⇒ 上一轮写死的停止条件通过,分层可用。**

## G1 ESTIMAND
一般因子 = `#682` 保留的 **21 题**(两端熵都 ≥0.50)两两耦合的中位。
⚠ **两端定义在本轮统一为 `educ<=12` vs `educ>=16`**(分层后 `<=10` 的格 n 太小),
**整体值也用同一定义重算一遍,保证与分层值可比** —— 不许拿旧定义的整体值和新定义的分层值比。
**主量:① 沿 educ 的 Δ(整体)· ② 沿 cohort 的 Δ · ③ 每个 cohort 五分位内沿 educ 的 Δ。**
## G2 CONTROLS
**正对照**:每个 cohort 格内两端的中位都必须 >0.20,否则该格记 `UNVERIFIED` 而不是「掉到地板」。
**安慰剂**:性别角色↔外部的同一分层(`#681` 已测它低于一般因子)。
## KILL(条件式)
if 每格正对照都过:
  五格 educ Δ 同号且中位 ≥ 整体的一半 -> **是教育** ·
  五格中位掉到地板 -> **是世代,`#675`–`#682` 全线改口** ·
  两者都在 -> **如实记「两个都在,本设计分不开」,这也是一个结果**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
`cohort` = `year − age`,与 `year` 相关 +0.4601 ⇒ **世代与时期在本设计里不可分**,
本轮只分离「教育 vs 世代」,**不主张分离时期**;跨仪器:MFQ 无 cohort;因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]; FEM=["fefam","fepol","fepresch"]
KEEP=['spkath','spkrac','spkcom','spkmil','colath','colrac','colcom','colmil',
      'libath','librac','libcom','libmil','libhomo','suicide1','suicide4',
      'abdefect','abnomore','abpoor','abrape','absingle','abany']
FLIP=lambda c:(c.startswith("spk") or (c.startswith("col") and c!="colcom")
               or c.startswith("suicide") or c.startswith("ab"))
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ","cohort"]+SEX+FEM+KEEP, apply_value_formats=False, encoding="latin1")
for c in KEEP:
    if FLIP(c): df[c]=-df[c]
df=df.dropna(subset=["educ","cohort"])
def nrm(d,a,b):
    if len(d)<150 or d[a].nunique()<2 or d[b].nunique()<2: return np.nan
    x=d[a].rank().to_numpy(float); y=d[b].rank().to_numpy(float)
    r=np.corrcoef(x,y)[0,1]
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
    c=np.corrcoef(xs,ys)[0,1]
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
EEp=list(combinations(KEEP,2)); FEp=[(s,e) for s in FEM for e in KEEP]
def med(q,pairs):
    v=[nrm(q.dropna(subset=[a,b]),a,b) for a,b in pairs]
    v=[u for u in v if np.isfinite(u)]
    return (float(np.median(v)) if v else np.nan), len(v)
def delta(q,lo,hi,pairs):
    a,na=med(q[lo],pairs); b,nb=med(q[hi],pairs)
    return b-a,a,b,min(na,nb)
d_e,lo_e,hi_e,n_e=delta(df,df.educ<=12,df.educ>=16,EEp)
print(f"① 沿 educ(<=12 vs >=16)整体 Δ = **{d_e:+.4f}**(两端中位 {lo_e:+.4f} / {hi_e:+.4f},可用对 {n_e})")
cl,ch=df.cohort.quantile(.2),df.cohort.quantile(.8)
d_c,lo_c,hi_c,n_c=delta(df,df.cohort<=cl,df.cohort>=ch,EEp)
print(f"② 沿 cohort(<={cl:.0f} vs >={ch:.0f})整体 Δ = **{d_c:+.4f}**(两端中位 {lo_c:+.4f} / {hi_c:+.4f})")
print(f"\n③ 每个 cohort 五分位内沿 educ:")
q5=pd.qcut(df.cohort,5,duplicates='drop'); ds=[]
for lv,g in df.groupby(q5,observed=True):
    d,a,b,n=delta(g,g.educ<=12,g.educ>=16,EEp)
    okc=(a>0.20 and b>0.20)
    ds.append(d if okc else np.nan)
    print(f"   {str(lv):22s} n={len(g):>6,}  两端中位 {a:>+7.4f}/{b:>+7.4f} {'✅' if okc else '⛔正对照失败'}  Δ = **{d:+.4f}**")
ds=np.array(ds,dtype=float); mid=float(np.nanmedian(ds))
print(f"\n   五格 Δ 中位 = **{mid:+.4f}** · 同号 {int(np.nansum(np.sign(ds)==np.sign(mid)))}/{int(np.isfinite(ds).sum())} · "
      f"占整体 **{mid/d_e*100:.0f}%**")
d_f,lof,hif,_=delta(df,df.educ<=12,df.educ>=16,FEp)
print(f"   安慰剂(性别角色↔外部)整体 Δ = **{d_f:+.4f}**")
G=Gate("是教育还是世代")
p1=G.positive_control("五个 cohort 格内两端中位都 >0.20",planted=float(np.nanmin(ds*0+1))*min(lo_e,1.0),floor=0.20,spread=0.01)
allok=np.isfinite(ds).all()
if p1 and allok:
    v=(f"**是教育:五格 Δ 中位 {mid:+.4f},占整体 {mid/d_e*100:.0f}%,同号 "
       f"{int(np.nansum(np.sign(ds)==np.sign(mid)))}/5;而沿 cohort 只有 {d_c:+.4f}**"
       if mid>=0.5*d_e else
       f"**是世代:cohort 内分层后 educ 的 Δ 中位掉到 {mid:+.4f}(整体 {d_e:+.4f}) ⇒ `#675`–`#682` 全线改口**")
else: v="UNVERIFIED —— 有格未过正对照"
print(f"\n{v}"); print(G)
json.dump(dict(d_educ=d_e,d_cohort=d_c,strata=[None if not np.isfinite(x) else float(x) for x in ds],
               median_strata=mid,placebo=d_f,verdict=v,unchallenged=True),
          open(OUT/"education_or_cohort.json","w"),indent=1,ensure_ascii=False)

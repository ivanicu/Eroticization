"""E03·A23·R123 —— 「性与堕胎绑紧」里有多少只是一个一般因子

**类型:FRONTIER。** `#680` 的 NEXT,而它是这条线上第一个可能把整条线降级的检验。
`#680` 报「教育让性、堕胎、安乐死一起绑紧」。**但若教育让任意两道道德题都绑得更紧,
那「性与它们绑紧」就不是关于性的。**

## 硬规则①(已跑)
对数:**外部↔外部 325 对 · 性↔外部 104 对**(全部可用)。对内 n 中位 5,355(低端)/ 9,049(高端)。
低端:外部↔外部 **+0.2753** · 性↔外部 **+0.2553**;高端:**+0.4302** · **+0.4660**。
⇒ **外部↔外部 Δ = +0.1549,性↔外部 Δ = +0.2107,差只有 +0.0558。**

## G1 ESTIMAND
**主量 = (性↔外部的 Δ 中位) − (外部↔外部的 Δ 中位)**,一个**差**,不是比(`#679`:比值不单调)。
## G2 CONTROLS
**正对照**:两端的两个中位都必须为正且 >0.20(对齐正确的证据,`#680` 已建)。
**⚠ 自助的单位是「题」不是「对」** —— 325 对之间共享题,不独立;重抽 26 个外部题(有放回),
每次重算两个中位与差。**这是预注册⑤要求的那件事。**
**安慰剂**:性别角色三题↔外部的同一个差(`#676` 已测它内部一致性不动)。
## KILL(条件式)
if 正对照过:
  差的 95% 自助区间**不含零** -> **性与外部的联结特别** ·
  含零 -> **是一个一般因子;`#680` 那句话缩成「教育让所有道德判断彼此绑紧,性不例外」**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
26 个外部题**不是道德题的随机样本**(GSS 问什么就有什么)⇒ 「一般因子」只在这 26 题的范围内成立;
`nat*` 无宽容轴已排除;跨仪器:MFQ 无容忍题组;因果:横断面无干预。`[unchallenged]`
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
EXT=["spkath","spkrac","spkcom","spkmil","spkhomo","colath","colrac","colcom","colmil","colhomo",
     "libath","librac","libcom","libmil","libhomo","suicide1","suicide2","suicide3","suicide4",
     "abdefect","abnomore","abhlth","abpoor","abrape","absingle","abany"]
FLIP=lambda c:(c.startswith("spk") or (c.startswith("col") and c!="colcom")
               or c.startswith("suicide") or c.startswith("ab"))
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["educ"]+SEX+FEM+EXT, apply_value_formats=False, encoding="latin1")
for c in EXT:
    if FLIP(c): df[c]=-df[c]
def nrm(d,a,b):
    if len(d)<200 or d[a].nunique()<2 or d[b].nunique()<2: return np.nan
    x=d[a].rank().to_numpy(float); y=d[b].rank().to_numpy(float)
    r=np.corrcoef(x,y)[0,1]
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
    c=np.corrcoef(xs,ys)[0,1]
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
LO=df[df.educ<=10]; HI=df[df.educ>=16]
def matrix(q,A,Bs):
    M={}
    for a in A:
        for b in Bs:
            if a==b: continue
            k=tuple(sorted((a,b)))
            if k in M: continue
            M[k]=nrm(q.dropna(subset=[a,b]),a,b)
    return M
EEl,EEh=matrix(LO,EXT,EXT),matrix(HI,EXT,EXT)
SEl,SEh=matrix(LO,SEX,EXT),matrix(HI,SEX,EXT)
FEl,FEh=matrix(LO,FEM,EXT),matrix(HI,FEM,EXT)
def med(M,keys): 
    v=[M[k] for k in keys if k in M and np.isfinite(M[k])]
    return float(np.median(v)) if v else np.nan
def stats(items):
    ee=[tuple(sorted(p)) for p in combinations(items,2)]
    se=[tuple(sorted((s,e))) for s in SEX for e in items]
    fe=[tuple(sorted((s,e))) for s in FEM for e in items]
    dee=med(EEh,ee)-med(EEl,ee); dse=med(SEh,se)-med(SEl,se); dfe=med(FEh,fe)-med(FEl,fe)
    return dee,dse,dfe
dee,dse,dfe=stats(EXT)
print(f"外部↔外部 Δ = **{dee:+.4f}** · 性↔外部 Δ = **{dse:+.4f}** · 性别角色↔外部 Δ = **{dfe:+.4f}**")
print(f"主量 差 = **{dse-dee:+.4f}** · 安慰剂差 = **{dfe-dee:+.4f}**")
rng=np.random.default_rng(20260806); bs=[];bp=[]
for _ in range(3000):
    it=list(rng.choice(EXT,size=len(EXT),replace=True))
    if len(set(it))<4: continue
    a,b,c=stats(sorted(set(it)))
    if np.isfinite(a) and np.isfinite(b): bs.append(b-a)
    if np.isfinite(a) and np.isfinite(c): bp.append(c-a)
bs=np.array(bs); bp=np.array(bp)
ci=(float(np.quantile(bs,.025)),float(np.quantile(bs,.975)))
cp=(float(np.quantile(bp,.025)),float(np.quantile(bp,.975)))
print(f"\n按题自助(B={len(bs)}):主量差 95% 区间 **[{ci[0]:+.4f}, {ci[1]:+.4f}]**  "
      f"{'✅ 不含零' if ci[0]*ci[1]>0 else '⛔ 含零'}")
print(f"                      安慰剂差 95% 区间 **[{cp[0]:+.4f}, {cp[1]:+.4f}]**  "
      f"{'不含零' if cp[0]*cp[1]>0 else '含零'}")
pos=min(med(EEl,[tuple(sorted(p)) for p in combinations(EXT,2)]),
        med(SEl,[tuple(sorted((s,e))) for s in SEX for e in EXT]))
G=Gate("性与外部的联结是特别的,还是一个一般因子")
p1=G.positive_control("两端两个中位都为正且 >0.20(对齐正确的证据)",planted=float(pos),floor=0.20,spread=0.01)
if p1:
    v=(f"**性与外部的联结特别:差 {dse-dee:+.4f},95% 区间 [{ci[0]:+.4f}, {ci[1]:+.4f}] 不含零 —— "
       f"但一般因子占 {dee/dse*100:.0f}%,必须同报**" if ci[0]*ci[1]>0 else
       f"**是一个一般因子:差 {dse-dee:+.4f} 的 95% 区间 [{ci[0]:+.4f}, {ci[1]:+.4f}] 含零 ⇒ "
       f"`#680` 缩成「教育让所有道德判断彼此绑紧,性不例外」**")
else: v="UNVERIFIED —— 正对照失败"
print(f"\n{v}"); print(G)
json.dump(dict(d_ext=dee,d_sex=dse,d_fem=dfe,contrast=dse-dee,ci=ci,plac_ci=cp,
               general_share=dee/dse if dse else None,verdict=v,unchallenged=True),
          open(OUT/"general_factor.json","w"),indent=1,ensure_ascii=False)

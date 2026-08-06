"""E03·A25·R130 —— 这条线量到的是立场,还是结构

**类型:FRONTIER,而且是这条线上唯一一个能把自己十二轮整体降级的检验。**
`#687` 暴露的窟窿:**`ρ(educ, stance) = +0.3555` 比这条线量过的任何压缩效应大一倍以上,
而 `#675`–`#686` 十二轮从没把它放进过模型** ⇒ 所有「教育 → 耦合/压缩」都可能只是「教育 → 立场」的影子。

## ⑤ 最强混淆先算(`#687` 预注册):均值与标准差在有界刻度上必然相关
`stance` = 33 题标准化答案的**均值**,`totvar` = **同一批**的标准差。
**把每一题在人之间独立打乱**(边际完全保留,跨题结构全毁),重算 `corr(mean, sd)` ——
**这就是「同边际下必然的那部分」,它必须作为偏相关的基线,不许算进 `stance` 的解释力。**

## G1 ESTIMAND
**① `ρ(educ, totvar)` 控掉 `stance` 后的偏相关**(`#686` 的量)。
**② `#681` 的一般因子:26 题两两耦合的两端中位差,在「stance 五分位内」重算后取中位**
   —— 一般因子是**题对层**的量,控 `stance` 只能靠分层,不能靠人层偏相关。
## ⑧ 判据(`#687` 在跑之前写死,不得改)
**控掉 `stance` 后仍保号 且 ≥ 0.5 × 原值 ⇒ 教育有独立于立场的作用;
收缩到 0.5 以下或变号 ⇒ 记「这条线量到的是立场,不是结构」,`#675`–`#686` 全部降级。**
## G2 CONTROLS
**正对照**:容忍题熵降 / 性题熵升(沿用 `#686`)。
**安慰剂(④ 预注册)**:用 **`cohort` 代替 `stance`** 做同一控制 ——
`#683` 已证 cohort 只吃掉约 7%,**若控 cohort 也吃掉同样多,说明偏相关做法本身有问题,当场停**。
**这个零该不该是零?** 不该(cohort 已知吃掉约 7%)⇒ `offset_control`,
**零的种类 = 用一个已知只吃 7% 的协变量做同一控制时的收缩量。**
## KILL(条件式)
if 正对照复现 and 控 cohort 的收缩量确实小: evaluate(判据⑧) else UNVERIFIED
## IMPOSSIBLE(不写 planned)
`stance` 与 `totvar` 同源 ⇒ **「控掉 stance」必然连着结构一起控掉一部分**,合成基线只能定量下界,
**不能声称已经分开**;一般因子是题对层的量 ⇒ 只能分层,**分层格内 n 更小、精度更低**;
跨仪器:MFQ 题集不同且无 cohort;因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
BAT={"性四题":["premarsx","xmarsex","homosex","teensex"],
     "堕胎":["abrape","abhlth","abdefect","abpoor","abnomore","absingle","abany"],
     "自杀":["suicide1","suicide2","suicide3","suicide4"],
     "容忍·言论":["spkath","spkrac","spkcom","spkmil","spkhomo"],
     "容忍·任教":["colath","colrac","colcom","colmil","colhomo"],
     "容忍·藏书":["libath","librac","libcom","libmil","libhomo"],
     "性别角色":["fefam","fepol","fepresch"]}
FLIP=lambda c:(c.startswith("spk") or (c.startswith("col") and c!="colcom")
               or c.startswith("suicide") or c.startswith("ab"))
ALL=sorted({c for v in BAT.values() for c in v})
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ","cohort"]+ALL, apply_value_formats=False, encoding="latin1")
for c in ALL:
    if FLIP(c): df[c]=-df[c]
d=df[(df.year>=1988)&(df.year<=2018)&df.educ.notna()&df.cohort.notna()].copy()
Z=pd.DataFrame({c:(d[c]-d[c].mean())/d[c].std() for c in ALL},index=d.index)
d["nans"]=d[ALL].notna().sum(1); d["totvar"]=Z.std(1); d["stance"]=Z.mean(1)
keep=(d.nans>=5)&d.totvar.notna()&d.stance.notna()
d=d[keep].reset_index(drop=True); Z=Z[keep.values].reset_index(drop=True)
print(f"n = **{len(d):,}**")
rc=lambda a,b:float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
obs_ms=rc(d.stance.to_numpy(float),d.totvar.to_numpy(float))
rng=np.random.default_rng(20260806); sims=[]
for _ in range(50):
    Zp=pd.DataFrame({c:rng.permutation(Z[c].to_numpy()) for c in ALL})
    sims.append(rc(Zp.mean(1).to_numpy(),Zp.std(1).to_numpy()))
forced=float(np.median(sims))
print(f"\n⑤ 合成基线:每题独立打乱(边际保留、跨题结构全毁)⇒ **corr(mean, sd) = {forced:+.4f}**")
print(f"   实测 corr(stance, totvar) = **{obs_ms:+.4f}** ⇒ "
      f"**必然的那部分占 {abs(forced)/abs(obs_ms)*100:.0f}%**,其余 {100-abs(forced)/abs(obs_ms)*100:.0f}% 才是真结构")
def H(s):
    p=s.value_counts(normalize=True).to_numpy(); p=p[p>0]; k=len(p)
    return float(-(p*np.log(p)).sum()/np.log(k)) if k>1 else 0.0
lo,hi=d[d.educ<=12],d[d.educ>=16]
tol=[c for c in ALL if c[:3] in ("spk","col","lib")]
dtol=float(np.median([H(hi[c].dropna())-H(lo[c].dropna()) for c in tol]))
dsex=float(np.median([H(hi[c].dropna())-H(lo[c].dropna()) for c in BAT["性四题"]]))
print(f"\n正对照:容忍题 Δ熵 **{dtol:+.4f}** · 性题 **{dsex:+.4f}** {'✅' if dtol<0<dsex else '⛔ 当场停'}")
def part(a,b,*ctrl):
    A=pd.Series(a).rank().to_numpy(); B=pd.Series(b).rank().to_numpy()
    C=np.column_stack([pd.Series(c).rank().to_numpy() for c in ctrl]+[np.ones(len(A))])
    ra=A-C@np.linalg.lstsq(C,A,rcond=None)[0]; rb=B-C@np.linalg.lstsq(C,B,rcond=None)[0]
    return float(np.corrcoef(ra,rb)[0,1])
E=d.educ.to_numpy(float); T=d.totvar.to_numpy(float); S=d.stance.to_numpy(float); CO=d.cohort.to_numpy(float)
r0=rc(E,T); rs=part(E,T,S); rcoh=part(E,T,CO)
print(f"\n=== ① `#686` 的量:ρ(educ, totvar) ===")
print(f"  原值 **{r0:+.4f}** · **控 stance {rs:+.4f}(保留 {rs/r0*100:.0f}%)** · 控 cohort(安慰剂) {rcoh:+.4f}(保留 {rcoh/r0*100:.0f}%)")
EXT=['spkath','spkrac','spkcom','spkmil','colath','colrac','colcom','colmil',
     'libath','librac','libcom','libmil','libhomo','suicide1','suicide4',
     'abdefect','abnomore','abpoor','abrape','absingle','abany']
def nrm(q,a,b):
    m=q[[a,b]].dropna()
    if len(m)<150 or m[a].nunique()<2 or m[b].nunique()<2: return np.nan
    x=m[a].rank().to_numpy(float); y=m[b].rank().to_numpy(float)
    r=np.corrcoef(x,y)[0,1]
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
    c=np.corrcoef(xs,ys)[0,1]
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
PAIR=list(combinations(EXT,2))
def genfac(q):
    a=[nrm(q[q.educ<=12],x,y) for x,y in PAIR]; b=[nrm(q[q.educ>=16],x,y) for x,y in PAIR]
    a=[v for v in a if np.isfinite(v)]; b=[v for v in b if np.isfinite(v)]
    return (float(np.median(b))-float(np.median(a))) if a and b else np.nan
g0=genfac(d)
q5=pd.qcut(d.stance,5,duplicates='drop'); gs=[genfac(g) for _,g in d.groupby(q5,observed=True)]
gs=[v for v in gs if np.isfinite(v)]; gmid=float(np.median(gs))
q5c=pd.qcut(d.cohort,5,duplicates='drop'); gc=[genfac(g) for _,g in d.groupby(q5c,observed=True)]
gc=[v for v in gc if np.isfinite(v)]; gcmid=float(np.median(gc))
print(f"\n=== ② `#681` 的一般因子(题对层,只能分层)===")
print(f"  原值 **{g0:+.4f}** · **stance 五分位内中位 {gmid:+.4f}(保留 {gmid/g0*100:.0f}%)** "
      f"· cohort 五分位内中位(安慰剂) {gcmid:+.4f}(保留 {gcmid/g0*100:.0f}%)")
print(f"  五格分别 stance: {[f'{v:+.4f}' for v in gs]}")
G=Gate("量到的是立场还是结构")
p1=G.positive_control("复现 #686:容忍题熵降、性题熵升",planted=1.0 if dtol<0<dsex else 0.0,floor=0.5,spread=0.01)
p2=G.offset_control("控 stance 的收缩必须超出「控一个已知只吃 7% 的协变量」的收缩",
                    effect=abs(1-rs/r0),offset=abs(1-rcoh/r0),spread=0.02,
                    null_kind="用 cohort(#683 已证只吃约 7%)做同一控制时的收缩量 —— 若它也吃掉同样多,是偏相关做法本身的问题")
keep1=(np.sign(rs)==np.sign(r0)) and abs(rs)>=0.5*abs(r0)
keep2=(np.sign(gmid)==np.sign(g0)) and abs(gmid)>=0.5*abs(g0)
if p1:
    v=(f"**教育有独立于立场的作用:`#686` 的量保留 {rs/r0*100:.0f}%,一般因子保留 {gmid/g0*100:.0f}%,两者都 ≥50% 且保号**"
       if (keep1 and keep2) else
       f"**这条线量到的是立场,不是结构:`#686` 保留 {rs/r0*100:.0f}%、一般因子保留 {gmid/g0*100:.0f}% ⇒ "
       f"`#675`–`#686` 全部降级**" if (not keep1 and not keep2) else
       f"**一半一半,照登:`#686` 保留 {rs/r0*100:.0f}%({'过' if keep1 else '不过'})、"
       f"一般因子保留 {gmid/g0*100:.0f}%({'过' if keep2 else '不过'})—— 两个量对立场的依赖不同,不许合并成一句话**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(d)),forced_corr=forced,obs_corr=obs_ms,r0=r0,r_stance=rs,r_cohort=rcoh,
               g0=g0,g_stance=gmid,g_cohort=gcmid,g_strata=gs,verdict=v,unchallenged=True),
          open(OUT/"stance_or_structure.json","w"),indent=1,ensure_ascii=False)

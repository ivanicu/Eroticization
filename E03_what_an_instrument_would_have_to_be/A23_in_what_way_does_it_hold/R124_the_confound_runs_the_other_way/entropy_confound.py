"""E03·A23·R124 —— 我预注册的混淆量出来是反的

**类型:FRONTIER。** `#681` 的 NEXT。预注册的怀疑是「高教育端答案更分散 ⇒ 更少天花板压制 ⇒ 耦合虚高」。
**硬规则①在跑之前就把它量死了,并且量出了它的反面。**

## 硬规则①(已跑)
30 题在 educ 两端的标准化熵:**Δ熵中位 −0.1156,只有 12/30 上升**
⇒ **高教育端并没有更分散。预注册的那条混淆当场被量掉,不是被论证掉的。**
⚠ **而分布是系统性的,方向相反:**
容忍题**全部大跌**(`spkhomo` −0.6373 · `colhomo` −0.5094 · `libath` −0.4635)⇒ 高教育端**近乎一致**;
性四题**上升**(`homosex` +0.2729 · `xmarsex` +0.2557 · `teensex` +0.2184)⇒ 高教育端**更分裂**;
安乐死 2/3/4 **+0.3177…+0.3713**。
**⇒ 真正的混淆是反的:近乎一致的题天花板会塌,而「除天花板」= 除一个很小的数 ⇒ 归一化会放大噪声。**

## G1 ESTIMAND(预注册那条 + 反向那条,两条都跑)
**① 预注册**:把每对耦合对「该对两题在该端的熵之积」回归,取**残差**的 Δ 中位。
**② 反向(本轮新增,因为混淆方向变了)**:只保留**两端熵都 ≥0.50** 的题(天花板不塌的那批),重算一般因子。
**两条规格都登,不挑一条。**
## G2 CONTROLS
**正对照**:限制后剩余题数必须 ≥8 且两端中位仍 >0.20,否则**记「限制之后没有仪器了」并停**。
**安慰剂**:性别角色↔外部的同一处理(`#681` 已测它低于一般因子)。
## G3/G4:线性与秩两种回归规格都报(⑤:熵与耦合的关系不是线性的)。
## KILL(条件式)
if 正对照过:
  残差 Δ 与限制后的一般因子**都**仍显著为正 -> **一般因子是实的** ·
  任一掉到地板 -> **一般因子里有归一化伪影,`#681` 那句话要再缩一层**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
熵与天花板不是同一个量(熵是单题的,天花板是成对的)⇒ **回归掉熵不等于回归掉天花板**,只是它的代理;
跨仪器:MFQ 无容忍题组;因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
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
LO=df[df.educ<=10]; HI=df[df.educ>=16]
def H(s):
    p=s.value_counts(normalize=True).to_numpy(); p=p[p>0]; k=len(p)
    return float(-(p*np.log(p)).sum()/np.log(k)) if k>1 else 0.0
Hlo={c:H(LO[c].dropna()) for c in SEX+FEM+EXT}; Hhi={c:H(HI[c].dropna()) for c in SEX+FEM+EXT}
def nrm(d,a,b):
    if len(d)<200 or d[a].nunique()<2 or d[b].nunique()<2: return np.nan
    x=d[a].rank().to_numpy(float); y=d[b].rank().to_numpy(float)
    r=np.corrcoef(x,y)[0,1]
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
    c=np.corrcoef(xs,ys)[0,1]
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
def build(q,Hm,pairs):
    v=[];h=[]
    for a,b in pairs:
        u=nrm(q.dropna(subset=[a,b]),a,b)
        if np.isfinite(u): v.append(u); h.append(Hm[a]*Hm[b])
    return np.array(v),np.array(h)
EEp=list(combinations(EXT,2)); SEp=[(s,e) for s in SEX for e in EXT]; FEp=[(s,e) for s in FEM for e in EXT]
def resid_delta(pairs,mode):
    out={}
    for tag,q,Hm in [("lo",LO,Hlo),("hi",HI,Hhi)]:
        v,h=build(q,Hm,pairs)
        if mode=="rank": x=pd.Series(h).rank().to_numpy(float); y=pd.Series(v).rank().to_numpy(float)
        else: x,y=h,v
        A=np.vstack([x,np.ones_like(x)]).T
        beta=np.linalg.lstsq(A,y,rcond=None)[0]
        out[tag]=float(np.median(y-A@beta))
        out[tag+"_raw"]=float(np.median(v))
    return out["hi"]-out["lo"], out["hi_raw"]-out["lo_raw"]
print("=== ① 预注册:回归掉熵之积后的残差 Δ 中位 ===")
res={}
for mode in ("linear","rank"):
    de,rawe=resid_delta(EEp,mode); ds,raws=resid_delta(SEp,mode); dfm,rawf=resid_delta(FEp,mode)
    res[mode]=dict(ext=de,sex=ds,fem=dfm,raw_ext=rawe)
    print(f"  {mode:7s} 外部↔外部残差 Δ **{de:+.4f}**(生 {rawe:+.4f}) · 性↔外部 **{ds:+.4f}** · 性别角色 **{dfm:+.4f}**")
KEEP=[c for c in EXT if Hlo[c]>=0.50 and Hhi[c]>=0.50]
print(f"\n=== ② 反向:只留两端熵都 ≥0.50 的题 ⇒ **{len(KEEP)}/{len(EXT)} 题** {KEEP} ===")
ok=len(KEEP)>=8
if ok:
    kp=list(combinations(KEEP,2)); ksp=[(s,e) for s in SEX for e in KEEP]
    vlo,_=build(LO,Hlo,kp); vhi,_=build(HI,Hhi,kp)
    slo,_=build(LO,Hlo,ksp); shi,_=build(HI,Hhi,ksp)
    gen=float(np.median(vhi)-np.median(vlo)); sx=float(np.median(shi)-np.median(slo))
    print(f"  限制后 两端中位 {np.median(vlo):+.4f} / {np.median(vhi):+.4f}  一般因子 Δ = **{gen:+.4f}**")
    print(f"  限制后 性↔外部 Δ = **{sx:+.4f}**  差 = **{sx-gen:+.4f}**")
    # ⚠ 自助里绝不能重算相关:预算好每对的两端值,自助只重抽题并查表。
    #   第一版每次自助都从原始数据重算 325 对 -> 2 分钟超时。这是成本表(door ⑦),不是算法问题。
    PV={}
    for a,b in combinations(KEEP,2):
        u=nrm(LO.dropna(subset=[a,b]),a,b); w=nrm(HI.dropna(subset=[a,b]),a,b)
        if np.isfinite(u) and np.isfinite(w): PV[(a,b)]=(u,w)
    rng=np.random.default_rng(20260806); bs=[]
    for _ in range(3000):
        it=sorted(set(rng.choice(KEEP,size=len(KEEP),replace=True)))
        if len(it)<4: continue
        vv=[PV[k] for k in combinations(it,2) if k in PV]
        if len(vv)<3: continue
        bs.append(float(np.median([z[1] for z in vv])-np.median([z[0] for z in vv])))
    bs=np.array(bs); ci=(float(np.quantile(bs,.025)),float(np.quantile(bs,.975)))
    print(f"  一般因子 95% 自助区间(按题) **[{ci[0]:+.4f}, {ci[1]:+.4f}]** {'✅ 不含零' if ci[0]*ci[1]>0 else '⛔ 含零'}")
else:
    gen=sx=np.nan; ci=(np.nan,np.nan); print("  ⛔ 限制之后剩余题数 <8 —— 记「限制之后没有仪器了」并停")
G=Gate("我预注册的混淆量出来是反的")
p1=G.positive_control("限制后剩余题数 ≥8 且两端中位仍 >0.20",
                      planted=float(np.median(vlo)) if ok else 0.0,floor=0.20,spread=0.01)
if p1:
    v=(f"**一般因子是实的:回归掉熵后残差 Δ 线性 {res['linear']['ext']:+.4f} / 秩 {res['rank']['ext']:+.4f};"
       f"只留不塌的 {len(KEEP)} 题后仍 {gen:+.4f},区间 [{ci[0]:+.4f}, {ci[1]:+.4f}] 不含零**"
       if (res['linear']['ext']>0 and res['rank']['ext']>0 and ci[0]*ci[1]>0) else
       f"**一般因子里有归一化伪影:残差 Δ {res['linear']['ext']:+.4f}/{res['rank']['ext']:+.4f},"
       f"限制后 {gen:+.4f} 区间 [{ci[0]:+.4f}, {ci[1]:+.4f}] ⇒ `#681` 再缩一层**")
else: v="UNVERIFIED —— 限制之后没有仪器了"
print(f"\n{v}"); print(G)
json.dump(dict(resid=res,keep=KEEP,gen_restricted=gen,sex_restricted=sx,ci=ci,
               entropy_lo=Hlo,entropy_hi=Hhi,verdict=v,unchallenged=True),
          open(OUT/"entropy_confound.json","w"),indent=1,ensure_ascii=False)

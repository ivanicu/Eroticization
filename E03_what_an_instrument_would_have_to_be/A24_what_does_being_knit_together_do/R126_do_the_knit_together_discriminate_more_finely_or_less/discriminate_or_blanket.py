"""E03·A24·R126 —— 连成一片的人,分辨得更细,还是一刀切

**类型:FRONTIER。A23 关弧,A24 开弧**(新决定:连成一片对这个人本身做了什么 —— 「算子」那一层)。

## ⚠ 两处在跑之前就改掉的东西,理由与时点一并入账

**① `#683` 的 NEXT ② 在人层不可识别。** 它写「主量 = `ρ(该人的块内耦合, 他的情境敏感度)`」——
**而耦合是跨人的相关,一个人没有相关。** 这与 `#675` 第一版是同一类错(`#674` 已点名)。
G1 要求先识别再谈功率 ⇒ **人层「连成一片」重定义为:同一人在 7 个题组上的标准化位置的跨组标准差**
(低 = 连成一片),`n` = 人数。**这是本轮唯一可识别的版本,不是一个更方便的版本。**

**② 硬规则①逼停了预注册的「情境敏感度 = 秩相关」。**
实测:堕胎七题皆答 n = 36,840,而**其中 17,075 人(46.3%)七题答案全同 ⇒ 秩相关无定义**。
预注册写的处置是「排除」,**而排除他们等于把「一刀切」这个被研究的现象本身排除掉** ——
**在看到任何结果之前**改为对所有人有定义的两个量:
  **`discrim`** = 该人七题答案的标准差(全同 = 0)· **`guttman`** = 1 − 该人相对总体情境序的 Guttman 误差比。

## 硬规则①(已跑)
堕胎七题同意率:`abhlth` **89.8%** · `abrape` 81.7% · `abdefect` 79.7% · `abpoor` 48.2% ·
`abnomore` 45.8% · `absingle` 45.2% · `abany` **43.3%**(各题 n 41,398–48,488)。
**正对照 ✅**:前三名 = {健康, 强奸, 畸形},与文献情境序一致 ⇒ 仪器没错。
七个题组全部答全 + `educ` 的 **n = 8,229**。

## G1 ESTIMAND
**主量 = `ρ(knit, discrim)` 与 `ρ(knit, guttman)`**,`knit = −sd_across_batteries`,**n = 人数**。
## G2 CONTROLS
**正对照**:情境序复现(已过,见上)。
**最强混淆 = 极端性**:跨组标准差与「答案极端」机械相关 ⇒ **同时报偏相关,控制 |跨组平均标准化位置|**。
**零**:把情境序**打乱**后重算 `guttman` —— **这个零该不该是零?** 该:若人们不按情境序作答,
打乱后应无差别 ⇒ `negative_control`,零的种类 = **随机情境序下的同一个 Guttman 一致度**。
**安慰剂**:同一流水线用在**自杀四题**(也有情境序)与**容忍·藏书五题**(是对象不是情境)。
## G3/G4:两个主量 × 三个题组 × {原始, 偏相关} = 12 格,全部照登;BH 覆盖整格。
## KILL(条件式)
if 情境序正对照过 and 打乱序的零确实为零:
  `ρ(knit, discrim)` 与 `ρ(knit, guttman)` **同为正** -> **W1 连成一片的人分辨得更细** ·
  **同为负** -> **W2 一刀切(不受欢迎的那个结果,而设计允许它)** ·
  两者反号 -> **W3 或二分本身错了,记「这个二分不成立」**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
七题的情境序是**总体的**,不是每人自己的 ⇒ Guttman 一致度衡量的是「与多数人同序」,不是「内在一致」;
`knit` 与 `discrim` 都来自同一批答案 ⇒ **共同方法方差无法用本设计消掉**,只能用偏相关缩小;
跨仪器:MFQ 无情境梯度题组;因果:横断面无干预。`[unchallenged]`
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
FLIP=lambda c:(c.startswith("spk") or (c.startswith("col") and c!="colcom")
               or c.startswith("suicide") or c.startswith("ab"))
ALL=sorted({c for v in BAT.values() for c in v})
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["educ"]+ALL, apply_value_formats=False, encoding="latin1")
for c in ALL:
    if FLIP(c): df[c]=-df[c]
RAW,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["educ"]+ALL, apply_value_formats=False, encoding="latin1")   # 未翻极性的副本,用于算同意率
J=df.dropna(subset=ALL+["educ"]).copy(); JR=RAW.loc[J.index]
print(f"分析样本 n = **{len(J):,}**")
Z={k:((J[v]-J[v].mean())/J[v].std()).mean(1) for k,v in BAT.items()}
ZM=pd.DataFrame(Z)
knit=(-ZM.std(1)).to_numpy(float)          # 低跨组 sd = 连成一片
extremity=ZM.mean(1).abs().to_numpy(float) # 最强混淆
def person_metrics(items, order):
    V=J[items].to_numpy(float)
    d=V.std(1)
    o=[items.index(c) for c in order]
    g=[]
    for row in V:
        seq=row[o]                      # 按总体情境序排好
        errs=sum(1 for i in range(len(seq)-1) for j in range(i+1,len(seq)) if seq[i]<seq[j])
        mx=len(seq)*(len(seq)-1)/2
        g.append(1-errs/mx)
    return d, np.array(g)
def rc(a,b):
    m=np.isfinite(a)&np.isfinite(b)
    if m.sum()<200 or np.std(a[m])<1e-9 or np.std(b[m])<1e-9: return np.nan
    return float(np.corrcoef(pd.Series(a[m]).rank(),pd.Series(b[m]).rank())[0,1])
def partial(a,b,c):
    rab,rac,rbc=rc(a,b),rc(a,c),rc(b,c)
    den=np.sqrt((1-rac**2)*(1-rbc**2))
    return (rab-rac*rbc)/den if den>1e-9 else np.nan
rows=[]; rng=np.random.default_rng(20260806)
for bn in ["堕胎","自杀","容忍·藏书"]:
    items=BAT[bn]
    # ⚠ 第一版写的是 `J[c].rank(pct=True).mean()` —— 对任何列都 ≈0.5,**退化统计量**,
    #   给出的序是 ['abrape','abnomore','abhlth',...],而真序是 ['abhlth','abrape','abdefect',...]。
    #   同意率必须在**未翻极性**的原始码上算,而且是「取最宽容那一档的比例」。
    rates={c: float((JR[c]==JR[c].min()).mean()) if bn!="容忍·藏书" else float((JR[c]==JR[c].max()).mean())
           for c in items}
    order=sorted(items,key=lambda x:-rates[x])
    d,g=person_metrics(items,order)
    r1,r2=rc(knit,d),rc(knit,g)
    p1_,p2_=partial(knit,d,extremity),partial(knit,g,extremity)
    nul=[]
    for _ in range(200):
        o2=list(rng.permutation(items)); _,g2=person_metrics(items,o2); nul.append(abs(rc(knit,g2)))
    nul=np.array(nul); q=float(np.nanquantile(nul,.95))
    rows.append(dict(batt=bn,rho_discrim=r1,rho_guttman=r2,part_discrim=p1_,part_guttman=p2_,
                     null_q95=q,order=order))
    print(f"\n=== {bn} · 情境序 {order} ===")
    print(f"  ρ(knit, discrim)  = **{r1:+.4f}**   偏相关(控极端性) **{p1_:+.4f}**")
    print(f"  ρ(knit, guttman)  = **{r2:+.4f}**   偏相关(控极端性) **{p2_:+.4f}**")
    print(f"  打乱情境序的零 95% 分位 **{q:.4f}**  {'✅ 观测超零' if abs(r2)>q else '⛔ 观测在零里'}')" if False else
          f"  打乱情境序的零 95% 分位 **{q:.4f}**  {'✅ 观测超零' if abs(r2)>q else '⛔ 观测在零里'}")
ab=[r for r in rows if r["batt"]=="堕胎"][0]
G=Gate("连成一片的人分辨得更细还是一刀切")
# ⚠ 第一版这里传的是字面量 1.0 —— **一道结构上不可能失败的正对照**,而它本该抓住上面那个退化统计量。
#   改成:用**脚本自己算出的那个序**去比对预注册的期望,命中数即种入值。
_top3=set(ab["order"][:3]); _hit=len(_top3 & {"abhlth","abrape","abdefect"})/3
print(f"\n正对照(用脚本自己算出的序):前三 = {ab['order'][:3]} · 命中 {_hit*3:.0f}/3")
pc=G.positive_control("堕胎七题的情境序必须复现文献(健康/强奸/畸形居前)",planted=_hit,floor=0.6,spread=0.01)
nc=G.negative_control("打乱情境序后 Guttman 一致度与 knit 的相关应消失",
                      null=ab["null_q95"],effect=abs(ab["rho_guttman"]),null_spread=0.01,
                      null_kind="随机情境序下的同一个 Guttman 一致度 —— 若人们不按情境序作答,打乱后应无差别")
if pc and nc:
    s1,s2=np.sign(ab["rho_discrim"]),np.sign(ab["rho_guttman"])
    v=("**W1 连成一片的人分辨得更细**" if s1>0 and s2>0 else
       "**W2 连成一片的人更一刀切 —— 这是我预注册的「不受欢迎的正结果」**" if s1<0 and s2<0 else
       f"**二分不成立:两个主量反号({ab['rho_discrim']:+.4f} / {ab['rho_guttman']:+.4f})—— 记「这个二分本身错了」**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(J)),rows=rows,verdict=v,unchallenged=True),
          open(OUT/"discriminate_or_blanket.json","w"),indent=1,ensure_ascii=False)

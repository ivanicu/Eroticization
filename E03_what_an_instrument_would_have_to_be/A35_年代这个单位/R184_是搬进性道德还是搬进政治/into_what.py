"""E03·A35·R184 —— 同性恋是搬进了性道德,还是搬进了政治

**类型:FRONTIER。这是 `#740`① —— 而它是这三轮里最难切、也最可能推翻结论的那一刀。**

**心理学的那一句(本轮要判的):三十六年里,对同性恋的判断和对婚前性的判断融合了。
但如果它同时也在和**警察该不该打人**、**该不该允许堕胎**、**该往哪里花钱**一起靠拢,
那它搬进的就不是「性道德」,是「政治」——而页上那句新话就是过强的。**

## ⚠ 我把预注册的检验设计得更狠了,理由写在这里
`#740`① 写的是「拿 `homosex × 警察四` 与 `homosex × 堕胎四` 跑同一条世代内斜率」。
**三个手挑的靶子回答不了元分离器**:也许 GSS 里**每一对**态度都在这几十年里靠拢(大众排队)。
⇒ **本轮跑 23 个题的全部 253 对**,三个靶子作为其中的**具名格**,
而**253 条斜率的分布**同时充当「一对典型的态度会怎么动」的参照。**这不是偏离预注册,是把它包住。**

## 硬规则①(已跑)
23 题全部满足「≥15 个调查年、档数 ≤6」:性 4(4–5 档)· 堕胎 4(**二值**)· 警察 4(**二值**)·
自杀 4(**二值**)· 性别角色 3(2–4 档)· 支出 4(3 档)。⚠ **档数从 2 到 5**,
而 `#731`/`#737` 已证天花板对比值影响极大 ⇒ **生相关与天花板归一各报一次(G4)。**

## G1 ESTIMAND
对每一对,**世代内合并斜率**(与 `#740` 逐字同一条路径:四个过地板的世代,按 n 加权)。
然后三个量:
- **S_sex** = `homosex` 与其余**三个性题**的斜率中位;
- **S_non** = `homosex` 与 **19 个非性题**的斜率中位;
- **S_all** = 全部 **253** 对的斜率中位(**典型的一对**)。

## W1–W4(三分 + 元分离器,双边)
| 世界 | 判据 | 读法 |
|---|---|---|
| **W1 专属于性** | `S_non ≤ 0.25 × S_sex` | 页上那句站得住 |
| **W2 搬进政治** | `S_non ≥ 0.60 × S_sex` | **「搬进性道德」改成「搬进政治」** |
| **W3 之间** | 0.25–0.60 | 报份额,不报判决 |
| **W4 元分离器** | `S_all ≥ 0.40 × S_sex` | **每一对都在靠拢 ⇒「搬进 X」这个提法本身是错的**,只能报名次 |

⚠ **W2/W3/W4 的正结果全都不利于页上最新那一行 —— 这正是本轮设计成能出它们的理由。**

## G2 CONTROLS
**④ 正对照**:`premarsx × homosex` 必须复现 `#740` 的 **+0.00672/年**(容差 0.0005)。
**零** = `negative_control`,**零的种类 = 在每个世代内部打乱受访者的年份标签 ——
保住世代构成、每年 n 与作答分布,只毁掉「同一世代里谁属于哪一年」**(与 `#740` 同一具零)。
⚠ **零只对关键那一对跑(B=200)** —— **253 对各跑一次置换在算力上做不到,如实说明**;
**其余的参照是那 253 条斜率自己的分布,而它是一个更贴题的参照(「一对典型的态度会怎么动」)。**
## G3:253 对全算,按族汇总全报,含不支持结论的。G4:生 / 天花板归一 各一次。
## ⑤ 停止条件(**双边**,跑之前写死)
- **`premarsx × homosex` 复现不到 0.0005 ⇒ UNVERIFIED 并停。**
- 依次判 **W4 → W2 → W1 → W3**(元分离器优先:若全体都在靠拢,后面三个都不成立)。
## IMPOSSIBLE(不写 planned)
23 题是我挑的(六个族),**不是 GSS 全部态度题的随机样本** ⇒ **`S_all` 是这 23 题的典型,不是 GSS 的典型**;
仍是**重复横断面**;**换不了仪器**。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
FAM={"性":["premarsx","xmarsex","homosex","teensex"],
     "堕胎":["abdefect","abnomore","abpoor","abrape"],
     "警察":["polabuse","polmurdr","polescap","polattak"],
     "自杀":["suicide1","suicide2","suicide3","suicide4"],
     "性别角色":["fefam","fepresch","fepol"],
     "支出":["natspac","natenvir","natheal","natcrime"]}
ALL=[c for v in FAM.values() for c in v]; FAMOF={c:f for f,v in FAM.items() for c in v}
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year","cohort"]+ALL,encoding="latin1")
g=g.dropna(subset=["cohort"]).copy()
g["gen"]=pd.cut(g.cohort,[1880,1928,1946,1965,1981,1997,2010],
                labels=["前1929","1929–45","婴儿潮46–64","X 65–80","千禧81–96","Z 97+"])
GENS=["1929–45","婴儿潮46–64","X 65–80","千禧81–96"]      # #740 的四个过地板世代
FN,FY=150,5
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def pair_slope(a,b,frame,norm=False):
    """世代内合并斜率(按 n 加权),与 `#740` 同一条路径。"""
    num=[];wt=[]
    for gn in GENS:
        fr=frame[frame.gen==gn][["year",a,b]].dropna()
        pts=[]
        for y,sub in fr.groupby("year"):
            if len(sub)<FN or sub[a].nunique()<2 or sub[b].nunique()<2: continue
            x=sub[a].to_numpy(float); yv=sub[b].to_numpy(float); r=sp(x,yv)
            if norm:
                xs=np.sort(x); ys=np.sort(yv); ys=ys if r>0 else ys[::-1]
                c=abs(sp(xs,ys))
                if c<1e-9: continue
                r=r/c
            pts.append((float(y),r,len(sub)))
        if len(pts)>=FY:
            num.append(float(np.polyfit([p[0] for p in pts],[p[1] for p in pts],1)[0])); wt.append(sum(p[2] for p in pts))
    if not num: return np.nan
    return float(np.average(num,weights=wt))
print("=== ④ 正对照:`premarsx × homosex` 必须复现 `#740` 的 +0.00672/年(容差 0.0005)===")
key=pair_slope("premarsx","homosex",g)
print(f"  实测 **{key:+.5f}** · 账本 +0.00672 · 差 {abs(key-0.00672):.5f} {'✅' if abs(key-0.00672)<=0.0005 else '⛔ ⑤ 触发'}")
if abs(key-0.00672)>0.0005:
    print("⛔ 停"); json.dump(dict(stop="旧值不可复现",key=key),open(OUT/"into.json","w"),indent=1,ensure_ascii=False); sys.exit(0)
print(f"\n=== G3:253 对全算(生相关)===")
S={}
for a,b in itertools.combinations(ALL,2): S[(a,b)]=pair_slope(a,b,g)
vals=np.array([v for v in S.values() if np.isfinite(v)])
print(f"  可算 {vals.size}/253 · 斜率中位 **{np.median(vals):+.5f}** · "
      f"5–95% [{np.quantile(vals,.05):+.5f}, {np.quantile(vals,.95):+.5f}] · 最大 {vals.max():+.5f}")
def hs(x): return S.get(("homosex",x),S.get((x,"homosex"),np.nan))
sexs=[hs(c) for c in FAM["性"] if c!="homosex"]
nons=[hs(c) for c in ALL if FAMOF[c]!="性"]
S_sex=float(np.nanmedian(sexs)); S_non=float(np.nanmedian(nons)); S_all=float(np.median(vals))
print(f"\n=== 三个量 ===")
print(f"  S_sex(homosex × 其余三个性题)中位 **{S_sex:+.5f}**  逐个:"+" ".join(f"{v:+.5f}" for v in sexs))
print(f"  S_non(homosex × 19 个非性题)中位 **{S_non:+.5f}** ⇒ 占 S_sex 的 **{S_non/S_sex:.0%}**")
print(f"  S_all(全部 {vals.size} 对)中位 **{S_all:+.5f}** ⇒ 占 S_sex 的 **{S_all/S_sex:.0%}**")
print(f"\n=== `homosex` × 各族的中位(全报,含不支持结论的)===")
for f,cs in FAM.items():
    v=[hs(c) for c in cs if c!="homosex"]
    v=[x for x in v if np.isfinite(x)]
    if v: print(f"  homosex × {f:8s} 中位 **{np.median(v):+.5f}**  ({len(v)} 对)")
rank=int((vals>S.get(("premarsx","homosex"),np.nan)).sum())+1
print(f"\n  `premarsx × homosex` = {key:+.5f} ⇒ **在 {vals.size} 对里排第 {rank}**")
rng=np.random.default_rng(20260806); nul=[]
for _ in range(200):
    P=g.copy(); P["year"]=P.groupby("gen",observed=True)["year"].transform(lambda s: rng.permutation(s.to_numpy()))
    nul.append(pair_slope("premarsx","homosex",P))
nul=np.array([x for x in nul if np.isfinite(x)]); q=np.quantile(nul,[0.025,0.975])
print(f"\n=== 零(世代内打乱年份,B={nul.size},只对关键那一对 —— 253 对各跑一次算力上做不到)===")
print(f"  95% 区间 [{q[0]:+.5f}, {q[1]:+.5f}] · 实测 {key:+.5f} ⇒ {'✅ 在零之外' if not (q[0]<=key<=q[1]) else '⚠ 落在零里'}")
print(f"\n=== G4:天花板归一版(档数 2–5 差异大,`#731`/`#737` 已证天花板影响极大)===")
Sn={}
for a,b in itertools.combinations(ALL,2): Sn[(a,b)]=pair_slope(a,b,g,norm=True)
vn=np.array([v for v in Sn.values() if np.isfinite(v)])
def hsn(x): return Sn.get(("homosex",x),Sn.get((x,"homosex"),np.nan))
Sn_sex=float(np.nanmedian([hsn(c) for c in FAM["性"] if c!="homosex"]))
Sn_non=float(np.nanmedian([hsn(c) for c in ALL if FAMOF[c]!="性"]))
print(f"  归一:S_sex {Sn_sex:+.5f} · S_non {Sn_non:+.5f}(占 {Sn_non/Sn_sex:.0%})· 全体中位 {np.median(vn):+.5f}"
      f"(占 {np.median(vn)/Sn_sex:.0%})")
G=Gate("同性恋是搬进了性道德还是搬进了政治")
p1=G.positive_control("`premarsx × homosex` 必须复现 #740(容差 0.0005)",planted=float(0.0005-abs(key-0.00672)),floor=0.0,spread=0.00002)
p2=G.negative_control("世代内打乱年份后关键那一对的斜率应回到零",null=float(max(abs(q[0]),abs(q[1]))),
    effect=abs(key),null_spread=0.00002,
    null_kind="在每个世代内部打乱受访者的年份标签 —— 保住世代构成、每年 n 与作答分布,只毁掉「同一世代里谁属于哪一年」")
r_non,r_all=S_non/S_sex,S_all/S_sex
if not p1: v="**UNVERIFIED:旧值不可复现**"
elif r_all>=0.40: v=f"**W4:全部 {vals.size} 对的中位就占 S_sex 的 {r_all:.0%} ⇒ 每一对都在靠拢,「搬进 X」这个提法本身是错的**"
elif r_non>=0.60: v=f"**W2:S_non 占 S_sex 的 {r_non:.0%} ⇒ 「搬进性道德」要改成「搬进政治」**"
elif r_non<=0.25: v=f"**W1:S_non 只占 S_sex 的 {r_non:.0%} ⇒ 专属于性,页上那句站得住**"
else: v=f"**W3:S_non 占 {r_non:.0%} ⇒ 报份额,不报判决**"
print(f"\n{v}"); print(G)
json.dump(dict(key=key,S_sex=S_sex,S_non=S_non,S_all=S_all,rank=rank,n_pairs=int(vals.size),
  null_ci=[float(q[0]),float(q[1])],norm=dict(S_sex=Sn_sex,S_non=Sn_non,S_all=float(np.median(vn))),
  by_family={f:float(np.nanmedian([hs(c) for c in cs if c!="homosex"])) for f,cs in FAM.items() if any(c!="homosex" for c in cs)},
  verdict=v,unchallenged=True),open(OUT/"into.json","w"),indent=1,ensure_ascii=False)

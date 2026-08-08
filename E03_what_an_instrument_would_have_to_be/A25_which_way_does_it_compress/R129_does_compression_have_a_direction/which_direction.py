"""E03·A25·R129 —— 压缩到底收到哪一头去

**类型:FRONTIER。A24 关弧,A25 开弧。**
**⚠ BASIN(诚实标注):连续五轮都在确认「教育压缩」。本轮不受欢迎的正结果 = W2**——
若压缩完全是单侧的,「压缩」就只是「受教育者更自由派」的换词,**这条线降级为对已知事实的重新发现**
(P4 / L21:来得容易的想法多半已在数据里)。**设计必须允许 W2 出现,而它在描述层已经先露头了。**

## 硬规则①(已跑)
窗口 1988–2018,至少答 5 题,**n = 38,284**。全部题已按 `#680` 的极性表对齐到「高 = 更宽容」。
**`corr(平均立场, totvar) = −0.2189`** ⇒ |r| ≤ 0.3,**⑤ 预注册的改道条件没触发,分侧可用**。
`totvar` **最低的一成**(n = 3,907)里:**宽容侧 79.5% · 严厉侧 20.5%**(全样本 50/50)
⇒ 两侧都 ≥10%,**可双侧报**。
⚠ **但这 79.5% 里有一部分是构造出来的**(立场与 totvar 相关 −0.2189)⇒ **安慰剂必须能吸收这部分。**

## G1 ESTIMAND
`stance` = 该人 33 题标准化答案的均值(高 = 更宽容);`totvar` = 同一批答案的人内标准差。
**主量:`ρ(educ, totvar)` 在「宽容侧」与「严厉侧」各算一遍**(按全样本 `stance` 中位分侧),n = 人数。
附:`ρ(educ, stance)` 整体值。
## ⑧ 判据(`#686` 在跑之前写死,不得改)
**两侧的 `ρ(educ, totvar)` 同号 且 量级比在 [0.5, 2] 内 ⇒ 压缩与方向无关(教育压缩,不选边);
一侧显著更强 ⇒ 压缩是有方向的,那是完全不同的一句话。**
## G2 CONTROLS
**正对照(沿用 `#686`)**:高教育端**容忍题熵降、性题熵升**,否则量错了当场停。
**安慰剂 A(纯随机)**:把「立场」换成**随机分配的侧**,两侧应无差别。
**安慰剂 B(相关度匹配的 sham,本轮新增)**:构造一个**与 `totvar` 的相关等于 −0.2189、
但不携带任何心理含义**的变量来分侧 —— **它复制了分侧带来的选择结构,却没有立场的内容**。
**这个零该不该是零?** **不该**:B 已知会制造一点两侧差异 ⇒ `offset_control`,
**零的种类 = 与 totvar 相关度匹配的随机分侧下的同一个两侧比值。**
## KILL(条件式)
if 正对照复现 and 安慰剂 A 两侧无差别: evaluate(判据⑧,并把观测比值对安慰剂 B 的比值作差) else UNVERIFIED
## IMPOSSIBLE(不写 planned)
`stance` 与 `totvar` 出自同一批答案 ⇒ **分侧永远不可能与 totvar 完全正交**,只能用 sham B 定量吸收;
1988–2018 之外题集不同 ⇒ 无法纳入;跨仪器:MFQ 无 cohort 且题集不同 ⇒ 拿不到第二具仪器的同一分侧;
因果:横断面无干预、教育非随机。`[unchallenged]`
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
      usecols=["year","educ","cohort"]+ALL, apply_value_formats=False, encoding="latin1")
for c in ALL:
    if FLIP(c): df[c]=-df[c]
d=df[(df.year>=1988)&(df.year<=2018)&df.educ.notna()&df.cohort.notna()].copy()
Z=pd.DataFrame({c:(d[c]-d[c].mean())/d[c].std() for c in ALL})
d["nans"]=d[ALL].notna().sum(1); d["totvar"]=Z.std(1); d["stance"]=Z.mean(1)
d=d[d.nans>=5].dropna(subset=["totvar","stance"]).reset_index(drop=True)
def H(s):
    p=s.value_counts(normalize=True).to_numpy(); p=p[p>0]; k=len(p)
    return float(-(p*np.log(p)).sum()/np.log(k)) if k>1 else 0.0
lo,hi=d[d.educ<=12],d[d.educ>=16]
tol=[c for c in ALL if c[:3] in ("spk","col","lib")]
dtol=float(np.median([H(hi[c].dropna())-H(lo[c].dropna()) for c in tol]))
dsex=float(np.median([H(hi[c].dropna())-H(lo[c].dropna()) for c in BAT["性四题"]]))
print(f"n = {len(d):,} · 正对照:容忍题 Δ熵 **{dtol:+.4f}** · 性题 **{dsex:+.4f}** "
      f"{'✅ 复现' if dtol<0<dsex else '⛔ 当场停'}")
E=d.educ.to_numpy(float); T=d.totvar.to_numpy(float); S=d.stance.to_numpy(float)
rc=lambda a,b:float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
print(f"\nρ(educ, stance) 整体 = **{rc(E,S):+.4f}**(正 = 受教育者更宽容)")
def two_sided(side):
    a=rc(E[side],T[side]); b=rc(E[~side],T[~side]); return a,b,int(side.sum()),int((~side).sum())
m=np.median(S); obs=two_sided(S>m)
print(f"\n=== 主量:按立场分侧 ===")
print(f"  宽容侧 n={obs[2]:,}  ρ(educ, totvar) = **{obs[0]:+.4f}**")
print(f"  严厉侧 n={obs[3]:,}  ρ(educ, totvar) = **{obs[1]:+.4f}**")
ratio=abs(obs[0])/abs(obs[1]) if obs[1] else np.inf
print(f"  同号 {'✅' if obs[0]*obs[1]>0 else '⛔'} · 量级比 **{ratio:.3f}** "
      f"{'✅ 在 [0.5,2] 内' if 0.5<=ratio<=2 else '⛔ 在 [0.5,2] 外'}")
rng=np.random.default_rng(20260806)
ra=[]
for _ in range(300):
    s=rng.permutation(S>m); a,b,_,_=two_sided(s)
    ra.append(abs(a)/abs(b) if b else np.nan)
ra=np.array(ra,dtype=float)
print(f"\n安慰剂 A(纯随机分侧):比值中位 **{np.nanmedian(ra):.3f}** · 95% 区间 "
      f"[{np.nanquantile(ra,.025):.3f}, {np.nanquantile(ra,.975):.3f}]")
target=rc(S,T); rb=[]
for _ in range(300):
    g=rng.standard_normal(len(T))
    tr=pd.Series(T).rank().to_numpy(); tr=(tr-tr.mean())/tr.std()
    sham=target*tr+np.sqrt(max(0.0,1-target**2))*g     # 与 totvar 相关度匹配、无心理含义
    a,b,_,_=two_sided(sham>np.median(sham))
    rb.append(abs(a)/abs(b) if b else np.nan)
rb=np.array(rb,dtype=float)
print(f"安慰剂 B(与 totvar 相关度匹配 {target:+.4f} 的 sham 分侧):比值中位 **{np.nanmedian(rb):.3f}** · "
      f"95% 区间 [{np.nanquantile(rb,.025):.3f}, {np.nanquantile(rb,.975):.3f}]")
G=Gate("压缩是有方向的吗")
p1=G.positive_control("复现 #686:容忍题熵降、性题熵升",planted=1.0 if dtol<0<dsex else 0.0,floor=0.5,spread=0.01)
p2=G.negative_control("纯随机分侧下两侧应无差别(比值≈1)",
                      null=abs(float(np.nanmedian(ra))-1.0),effect=abs(ratio-1.0),null_spread=0.02,
                      null_kind="随机分侧下的同一个两侧比值 —— 若分侧本身不制造差异,它应当等于 1")
p3=G.offset_control("观测比值须超出「与 totvar 相关度匹配的 sham 分侧」制造的比值",
                    effect=abs(ratio-1.0),offset=abs(float(np.nanmedian(rb))-1.0),spread=0.02,
                    null_kind="与 totvar 相关度匹配的随机分侧下的同一个两侧比值 —— 它复制选择结构而不携带立场内容")
if p1 and p2:
    v=(f"**W1 压缩与方向无关:两侧 {obs[0]:+.4f} / {obs[1]:+.4f},同号,量级比 {ratio:.3f} 在 [0.5,2] 内**"
       if (obs[0]*obs[1]>0 and 0.5<=ratio<=2 and not p3) else
       f"**W2 压缩是有方向的:两侧 {obs[0]:+.4f} / {obs[1]:+.4f},量级比 {ratio:.3f}"
       f"{'(超出 sham 制造的 %.3f)'%float(np.nanmedian(rb)) if p3 else ''} ⇒ 这是完全不同的一句话**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(d)),rho_educ_stance=rc(E,S),lib=obs[0],harsh=obs[1],n_lib=obs[2],n_harsh=obs[3],
               ratio=float(ratio),placeboA=float(np.nanmedian(ra)),placeboB=float(np.nanmedian(rb)),
               corr_stance_totvar=target,d_tol=dtol,d_sex=dsex,verdict=v,unchallenged=True),
          open(OUT/"which_direction.json","w"),indent=1,ensure_ascii=False)

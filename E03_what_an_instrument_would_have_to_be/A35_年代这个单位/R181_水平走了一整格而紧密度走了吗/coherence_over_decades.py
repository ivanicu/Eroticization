"""E03·A35·R181 —— 水平走了一整格,而紧密度走了吗

**类型:FRONTIER。新弧 `A35` —— 年代,是 loop 点名而这一整段一次没碰过的那个单位。**

**心理学的那一句(本轮要判的):1988 到 2024,美国人对性的判断整体松了一整格。
那么在同一段时间里,这四个判断彼此**绑得更紧**了、**更散**了,还是**没动**?
「变宽容」与「变成一个立场」是不是同一件事,在年代这个单位上第一次被问。**

## ⚠ 换方向的理由先写下
`#735`/`#736`/`#737` 连着三轮在确认「我的估计量偏多少」——**那是盆地,而且是舒服的那一种。**
BASIN 规则要求换方向;而 `#737`① 是一个**caveat 的量级**问题,本轮是一句**关于人**的问题。
**`#737`① 不作废,记在 NEXT 里,并在本轮顺手交回它要的一半(逐年天花板)。**

## 硬规则①(已跑)
四题同时非缺失 **n=15,056 · 1988–2024 · 21 个调查年**,逐年 n **428–1,302,中位 731**。
水平(1=always wrong … 4=not wrong at all)**1.758 → 2.370,移动 +0.612**;
逐题 `premarsx +0.486 · xmarsex +0.200 · homosex +1.209 · teensex +0.554`。

## G1 ESTIMAND(先命名,后选统计量)
**逐年的「紧密度」= 六对天花板归一相关的均值**(不是最弱一环:min 在 n≈700 上噪声太大,
**这个选择在看到任何趋势之前写下**;最弱一环进 G4)。
**趋势统计量 = 紧密度对年份的 OLS 斜率(每年)。**

## W1 / W2 / W3 / W4
| 世界 | 斜率 | 读法 |
|---|---|---|
| **W1 越松越成一块** | **> 0** | 「变宽容」与「变成一个立场」同时发生 —— `#676` 的教育发现在**年代**上复制(跨单位) |
| **W2 越松越散** | **< 0** | 松开的同时判断也彼此脱钩,宽容**不是**一个立场 |
| **W3 没动** | 落在零里 | 年代这个单位上只有**水平**在动,结构不动 |
| **W4 判不了** | 零的半宽 > 36 年的总变化 | **这具仪器在年代单位上分辨不了紧密度** —— 而这本身是结论 |

⚠ **W1 是我想要的那个(它给页上添一条跨单位复制)** ⇒ **所以本轮必须让 W2/W3/W4 真的赢得了。**
判据因此**写双边**,且 W4 的地板**用年内劈半量出来,不假设**。

## G2 CONTROLS
- **④ 正对照**:**水平**的斜率必须显著为正(已知 +0.612/36 年)。
  ⚠ **且必须在 g=0 时失败**:打乱年份标签后,水平的斜率必须回到零。
  **测不出已知的水平趋势,就测不了紧密度趋势。**
- **零** = `negative_control`,**零的种类 = 打乱受访者的年份标签 ——
  保住每年的 n 与全体的作答分布,只毁掉「谁属于哪一年」。**
- **噪声地板(测出来,不假设)**:每年随机劈半,两半各算紧密度,取 |差| 的分布。
- **PLACEBO**:警察四题在同一批年份上的紧密度斜率 —— 一个没有理由随年代变的题组。
## G3:21 年全报。G4:均值/最弱一环 × 归一/生 × 全样本/劈半。
## ⑤ 停止条件(**双边**,跑之前写死 —— 兑现 `#736`②)
- **水平的斜率不显著为正,或打乱年份后仍显著 ⇒ UNVERIFIED 并停。**
- **零的 95% 半宽 × 36 年 > 观测的 36 年总变化 ⇒ 判 W4「这具仪器分辨不了」并停。**
- **斜率 > +0.0020/年 且在零之外 ⇒ W1;斜率 < −0.0020/年 且在零之外 ⇒ W2;
  在零之内 ⇒ W3。**(±0.0020/年 × 36 年 = ±0.072,相对紧密度 ~0.5–0.7 的量级)
- ⚠ **两侧都写了阈,两侧都能裁决 —— 不许出现「只有一侧有阈」的哑口。**
## IMPOSSIBLE(不写 planned)
GSS 是**重复横断面**,不是面板 ⇒ **不能说「同一个人变了」,只能说「这一年的美国人」**;
**换不了仪器**:没有第二份跨 1988–2024 问同一批性道德题的公开数据(`#700`/`#732` 已枚举);
1988 之前四题不齐(`teensex` 始于 1986)⇒ **1972–1987 结构上进不来。**`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]; POL=["polabuse","polmurdr","polescap","polattak"]
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+SEX+POL,encoding="latin1")
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def coh(fr,items,stat="mean"):
    """六对天花板归一相关的均值(或最弱一环)。返回 (值, 天花板中位)。"""
    v=[];c=[]
    for a,b in itertools.combinations(items,2):
        m=fr[[a,b]].dropna()
        if len(m)<60 or m[a].nunique()<2 or m[b].nunique()<2: return np.nan,np.nan
        x=m[a].to_numpy(float); y=m[b].to_numpy(float); r=sp(x,y)
        xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        cc=abs(sp(xs,ys))
        if cc<1e-9: return np.nan,np.nan
        v.append(r/cc); c.append(cc)
    return (float(np.mean(v)) if stat=="mean" else float(np.min(v))), float(np.median(c))
def slope(ys,xs): return float(np.polyfit(np.asarray(xs,float),np.asarray(ys,float),1)[0])
J=g.dropna(subset=SEX); yrs=sorted(J.year.unique())
print(f"=== G3 逐年全报({len(yrs)} 年)===")
print(f"{'年':>6s}{'n':>6s}{'水平':>8s}{'紧密度(均值)':>13s}{'最弱一环':>10s}{'天花板中位':>11s}")
rows={}
for y in yrs:
    fr=J[J.year==y]; cm,ce=coh(fr,SEX,"mean"); cw,_=coh(fr,SEX,"min")
    lv=float(fr[SEX].mean().mean()); rows[y]=dict(n=int(len(fr)),level=lv,coh=cm,weak=cw,ceil=ce)
    print(f"{int(y):>6d}{len(fr):>6d}{lv:>8.3f}{cm:>13.4f}{cw:>10.4f}{ce:>11.4f}")
Y=np.array(yrs,float); CO=np.array([rows[y]["coh"] for y in yrs]); LV=np.array([rows[y]["level"] for y in yrs])
WK=np.array([rows[y]["weak"] for y in yrs]); CE=np.array([rows[y]["ceil"] for y in yrs])
s_lv,s_co,s_wk=slope(LV,Y),slope(CO,Y),slope(WK,Y)
print(f"\n斜率(每年):水平 **{s_lv:+.5f}** · 紧密度 **{s_co:+.5f}** · 最弱一环 {s_wk:+.5f} · 天花板 {slope(CE,Y):+.5f}")
print(f"36 年总变化:水平 {s_lv*36:+.3f} · 紧密度 **{s_co*36:+.4f}** · 最弱一环 {s_wk*36:+.4f}")
rng=np.random.default_rng(20260806)
print("\n=== ④ 正对照:水平的斜率必须显著为正,且打乱年份后回到零 ===")
nl=[]
for _ in range(1000):
    p=J.assign(year=rng.permutation(J.year.to_numpy()))
    nl.append(slope([float(p[p.year==y][SEX].mean().mean()) for y in yrs],Y))
q=np.quantile(nl,[0.025,0.975])
pc=(s_lv>0) and (s_lv>q[1])
print(f"  水平斜率 {s_lv:+.5f} · 打乱年份后的零 95% 区间 [{q[0]:+.5f}, {q[1]:+.5f}] ⇒ {'✅' if pc else '⛔ ⑤ 触发'}")
print("\n=== 零:打乱受访者的年份标签(保住每年 n 与全体作答分布,只毁掉谁属于哪一年)===")
nc=[]
for _ in range(1000):
    p=J.assign(year=rng.permutation(J.year.to_numpy()))
    nc.append(slope([coh(p[p.year==y],SEX,"mean")[0] for y in yrs],Y))
qc=np.quantile(nc,[0.025,0.975]); half=(qc[1]-qc[0])/2
print(f"  紧密度斜率的零:95% 区间 [{qc[0]:+.5f}, {qc[1]:+.5f}] · 半宽 {half:.5f} ⇒ 36 年 **±{half*36:.4f}**")
print("\n=== 噪声地板(年内随机劈半,5 个种子,测出来不假设)===")
fl=[]
for s in range(5):
    r=np.random.default_rng(700+s)
    for y in yrs:
        fr=J[J.year==y]; idx=r.permutation(len(fr)); h=len(fr)//2
        a,_=coh(fr.iloc[idx[:h]],SEX,"mean"); b,_=coh(fr.iloc[idx[h:2*h]],SEX,"mean")
        if np.isfinite(a) and np.isfinite(b): fl.append(abs(a-b))
print(f"  年内劈半 |差| 中位 **{np.median(fl):.4f}** · 95% 分位 {np.quantile(fl,0.95):.4f}(n≈365/半)")
print("\n=== PLACEBO:警察四题在同一批年份上的紧密度斜率 ===")
P=g.dropna(subset=POL); yp=[y for y in yrs if len(P[P.year==y])>=200]
cp=[coh(P[P.year==y],POL,"mean")[0] for y in yp]
ok=[(y,c) for y,c in zip(yp,cp) if np.isfinite(c)]
s_pl=slope([c for _,c in ok],[y for y,_ in ok]) if len(ok)>5 else np.nan
print(f"  可算年 {len(ok)}/{len(yp)} · 斜率 **{s_pl:+.5f}**(36 年 {s_pl*36:+.4f})")
G=Gate("水平走了一整格,而紧密度走了吗")
p1=G.positive_control("水平的斜率必须显著为正且打乱年份后回到零",planted=1.0 if pc else 0.0,floor=0.0,spread=0.1)
p2=G.negative_control("打乱年份标签后紧密度的斜率应回到零",null=float(max(abs(qc[0]),abs(qc[1]))),
    effect=abs(s_co),null_spread=0.00002,
    null_kind="打乱受访者的年份标签 —— 保住每年的 n 与全体的作答分布,只毁掉「谁属于哪一年」")
tot=abs(s_co*36)
if not p1: v="**UNVERIFIED:正对照没过**"
elif half*36>tot: v=(f"**W4 判不了:零的 36 年半宽 ±{half*36:.4f} 大过观测的总变化 {s_co*36:+.4f} "
                     f"⇒ 这具仪器在年代单位上分辨不了紧密度**")
elif s_co>0.0020 and s_co>qc[1]: v=f"**W1:紧密度斜率 {s_co:+.5f}/年(36 年 {s_co*36:+.4f})—— 越松越成一块**"
elif s_co<-0.0020 and s_co<qc[0]: v=f"**W2:紧密度斜率 {s_co:+.5f}/年(36 年 {s_co*36:+.4f})—— 越松越散**"
elif qc[0]<=s_co<=qc[1]: v=f"**W3:紧密度斜率 {s_co:+.5f} 落在零的 95% 区间内 —— 年代上只有水平在动,结构不动**"
else: v=f"**斜率 {s_co:+.5f} 在零之外但未过 ±0.0020 的实质阈 ⇒ 方向可报,量级不可报**"
print(f"\n{v}"); print(G)
json.dump(dict(rows={int(y):rows[y] for y in yrs},slope_level=s_lv,slope_coh=s_co,slope_weak=s_wk,
  slope_ceiling=slope(CE,Y),null_ci=[float(qc[0]),float(qc[1])],noise_floor_median=float(np.median(fl)),
  placebo_slope=float(s_pl),verdict=v,unchallenged=True),open(OUT/"decades.json","w"),indent=1,ensure_ascii=False)

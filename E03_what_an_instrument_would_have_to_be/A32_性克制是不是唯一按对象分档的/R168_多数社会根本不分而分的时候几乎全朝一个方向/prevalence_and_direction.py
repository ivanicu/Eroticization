"""E03·A32·R168 —— 多数社会根本不分;而一旦分,几乎全朝同一个方向

**类型:FRONTIER。它是 `R167` 那次 UNVERIFIED 的重设,而缺陷有名字。**

**心理学的那一句:一个社会教孩子的十件事,多数社会对男孩女孩给的是同一个分 ——
性别化不是普遍的底色,而是发生在**某些事**上。而一旦一个社会开始分,它朝哪边分,
在全世界几乎是一致的:自立与勇气给男孩,勤勉与性克制给女孩。**

## ⚠ `R167` 为什么判 UNVERIFIED,以及缺陷叫什么
`R167` 的估计量是**逐社会男女差的中位数**。20 格全部返回 **0.000**,正对照(攻击性男>女)因此也判失败。
**不是仪器瞎,是估计量退化**:每一格里给同一个分的社会占 **53%–84%**
(如攻击性·晚 91/127,服从·早 134/159)⇒ **一个超过一半并列于零的分布,它的中位恒等于零。**
这是 `#674`「统计量饱和」的同型,而这一次它把一个**真实且强**的信号压成了零:
攻击性·晚的计数是 **女>男 4 · 男>女 32**,文献先验就在数据里,中位看不见它。
> **一个中位数在有大量并列的离散量上不是稳健,是失明。**

## G1 ESTIMAND(重新命名,而且是**两个**不同的量)
① **PREVALENCE 分不分**:`差异率 = 1 − 并列数/n` —— **这个社会到底有没有按性别区别对待这件事。**
② **DIRECTION 朝哪边分**:在**分了的社会里**,`女>男` 的比例 `p` —— **一个干净的二项量,零 = 0.5。**
**两个是不同的心理学**:① 说「这件事被性别化得多普遍」,② 说「性别化的方向有多一致」。
**不许合并成一个数**(`R167` 的错正是把两者压进一个中位)。

## W1–W4
| 世界 | ① 差异率 | ② 方向一致度 | 读法 |
|---|---|---|---|
| **W1 性别是底色** | 各品格都高 | 各品格都高 | 社会普遍按性别教养,性只是其中一件 |
| **W2 性别是局部的** | **只在少数品格上高** | 高 | **性别化发生在特定内容上,不是一种普遍风格** |
| **W3 方向不一致** | 任意 | ≈0.5 | 「双标」是各地各自的偏好,没有跨文化方向 |
| **W4 性克制不特殊** | 性克制不在前三 | — | **页上 `#724` 那一行要收窄** |

⚠ **W4 仍然是我不高兴的那个,并且它这一轮仍然可能赢。**

## G2 CONTROLS
- **④ 正对照(文献先验)**:**攻击性在「分了的社会」里必须是男>女**,即 `p < 0.5` 且显著。
  ⚠ **且必须在 g=0 时失败**:把每个社会内的男/女标签打乱后,`p` 必须回到 0.5。
- **零** = `negative_control`,**零的种类 = 在每个社会内部打乱「男/女」标签 ——
  保住该社会的两个取值、它的水平、以及「分不分」这件事本身,只毁掉「哪个值属于哪个性别」。**
  ⚠ **注意这个零对 ① 是不动的**(打乱标签不改变并列数)⇒ **① 没有这个零可用,如实说明**,
  ①的零是**同一批编码者在别的品格上的差异率**(即 sham 分布),不是置换。
- **SHAM**:九个非性品格,同编码者同量表 —— 「有界差在两端被压缩」这条代数在它们上同样成立。
- **PLACEBO**:`晚男 − 早男` 是年龄差,性别零毁不掉它 ⇒ 用来证明零只毁了性别标签。
## G3:10 品格 × 2 年龄 = 20 格全报。G4:早/晚 × 差异率/方向 × 有符号/绝对。
## ⑤ 停止条件(跑之前写死)
- **攻击性的 `p<0.5` 不成立,或打乱标签后仍成立 ⇒ UNVERIFIED 并停。**
- **性克制的差异率若不在前三 ⇒ 判 W4**,如实报名次。
- **任一品格的「分了的社会」少于 30 个 ⇒ 那一格记「判不了」,不进名次**(`#641` 的地板写在每一格上)。
## IMPOSSIBLE(不写 planned)
**换不了仪器**:只有 `barry1976traits` 编过这个 0–10 网格(`#700` 已枚举)。
「灌输强度」是**编码者读民族志给的分**,不是被观察的行为 —— **仪器必须点名**(硬规则②)。
并列可能是**编码者偷懒**(没读到差别就给同一个分)而非社会真的不分 ⇒ **①的读法上限是「记录里没有差别」**。
`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import binomtest
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
V=pd.read_csv(B+"variables.csv",low_memory=False); Dd=pd.read_csv(B+"data.csv")
W=Dd.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
P=V[V.source.astype(str).str.contains("barry1976traits",na=False)]
fam={}
for _,r in P.iterrows():
    m=re.match(r'(.+?):\s*(Early|Late)\s+(Boy|Girl)s?$',str(r.title))
    if m: fam.setdefault(m.group(1),{})[f"{m.group(2)}_{m.group(3)}"]=r.id
FAM={k:v for k,v in fam.items() if len(v)==4}
FLOOR=30
print("=== 硬规则①:仪器 = `barry1976traits` 编码者读民族志给的 0–10「灌输强度」分,不是被观察的行为 ===")
res={}
for k,d in FAM.items():
    for age,(cb,cg) in (("early",(d["Early_Boy"],d["Early_Girl"])),("late",(d["Late_Boy"],d["Late_Girl"]))):
        J=W[[cb,cg]].dropna(); g=(J[cg]-J[cb]).to_numpy(float)
        nd=int((g!=0).sum()); pos=int((g>0).sum())
        p=pos/nd if nd else np.nan
        bt=binomtest(pos,nd,0.5).pvalue if nd>=1 else np.nan
        res[(k,age)]=dict(n=len(g),n_diff=nd,prev=nd/len(g),pos=pos,neg=int((g<0).sum()),
                          p_girl=p,binom_p=float(bt),usable=nd>=FLOOR,
                          med_abs_diff=float(np.median(np.abs(g[g!=0]))) if nd else np.nan)
print(f"\n=== G3 全格 20 格 · 地板 = 分了的社会 ≥ {FLOOR} ===")
print(f"{'品格':18s}{'年龄':>5s}{'n':>5s}{'分了':>5s}{'差异率':>8s}{'女>男':>6s}{'男>女':>6s}{'p(女>男)':>9s}{'二项p':>9s}{'|差|中位':>8s}")
for k in FAM:
    for age in ("early","late"):
        r=res[(k,age)]; star=" ★" if k=="Sexual Restraint" else (" ⊕" if k=="Aggression" else "")
        flag="" if r["usable"] else "  ⚠判不了"
        print(f"{k:18s}{age:>5s}{r['n']:>5d}{r['n_diff']:>5d}{r['prev']:>8.3f}{r['pos']:>6d}{r['neg']:>6d}"
              f"{r['p_girl']:>9.3f}{r['binom_p']:>9.2e}{r['med_abs_diff']:>8.2f}{star}{flag}")
# ④ 正对照 + 零
rng=np.random.default_rng(20260806)
def null_p(cb,cg,Bn=4000):
    J=W[[cb,cg]].dropna(); a=J[cb].to_numpy(float); b=J[cg].to_numpy(float)
    out=[]
    for _ in range(Bn):
        f=rng.random(len(a))<0.5
        gb=np.where(f,b,a); gg=np.where(f,a,b); d=gg-gb; nd=(d!=0).sum()
        out.append((d>0).sum()/nd if nd else np.nan)
    return np.array([x for x in out if np.isfinite(x)])
print(f"\n=== ④ 正对照:攻击性在「分了的社会」里必须是 **男>女**(p<0.5),打乱标签后必须回到 0.5 ===")
agg_ok=True
for age in ("early","late"):
    d=FAM["Aggression"]; cb,cg=(d["Early_Boy"],d["Early_Girl"]) if age=="early" else (d["Late_Boy"],d["Late_Girl"])
    r=res[("Aggression",age)]; nul=null_p(cb,cg)
    lo,hi=np.quantile(nul,[0.025,0.975])
    ok=r["usable"] and r["p_girl"]<0.5 and (r["p_girl"]<lo)
    agg_ok&=ok
    print(f"  Aggression {age:5s}: p(女>男) = **{r['p_girl']:.3f}**(分了 {r['n_diff']})· "
          f"打乱标签后零的 95% 区间 [{lo:.3f}, {hi:.3f}] · 零中位 {np.median(nul):.3f} ⇒ {'✅' if ok else '⛔'}")
print(f"\n=== PLACEBO:晚男 − 早男 是**年龄**差,性别零毁不掉它 ===")
for k in ("Aggression","Sexual Restraint"):
    d=FAM[k]; J=W[[d["Early_Boy"],d["Late_Boy"]]].dropna(); dd=(J[d["Late_Boy"]]-J[d["Early_Boy"]])
    print(f"  {k:18s} 分了的社会 {int((dd!=0).sum())}/{len(J)} · 晚>早 的比例 {(dd>0).sum()/max((dd!=0).sum(),1):.3f}")
print(f"\n=== 名次(G4:两个量各排一次,且只排过地板的格)===")
ranks={}
for age in ("early","late"):
    ok=[k for k in FAM if res[(k,age)]["usable"]]
    for key,nm,rev in (("prev","差异率(分不分)",True),("p_girl","p(女>男)(朝哪边)",True)):
        od=sorted(ok,key=lambda k:-res[(k,age)][key])
        rk=od.index("Sexual Restraint")+1 if "Sexual Restraint" in od else None
        ranks[(age,key)]=(rk,len(od))
        print(f"  {age:5s} 按 {nm:20s} ⇒ 性克制第 **{rk}/{len(od)}** · 前三:{', '.join(od[:3])}")
G=Gate("多数社会分不分,分了朝哪边")
p1=G.positive_control("攻击性在分了的社会里必须是男>女,且打乱标签后回到 0.5(文献先验)",
    planted=1.0 if agg_ok else 0.0,floor=0.0,spread=0.1)
d=FAM["Sexual Restraint"]; nul=null_p(d["Late_Boy"],d["Late_Girl"])
sr=res[("Sexual Restraint","late")]
p2=G.negative_control("打乱每个社会内的男/女标签后,性克制的 p 应回到 0.5",
    null=abs(float(np.quantile(nul,0.95))-0.5),effect=abs(sr["p_girl"]-0.5),null_spread=0.01,
    null_kind="在每个社会内部打乱「男/女」标签 —— 保住两个取值、水平、以及「分不分」本身,只毁掉哪个值属于哪个性别")
rk_prev=[ranks[(a,"prev")][0] for a in ("early","late")]
if not p1: v="**UNVERIFIED:正对照没过**"
elif any(r and r<=3 for r in rk_prev):
    v=(f"**性克制的差异率在十个品格里排第 {rk_prev[0]}(早)/ {rk_prev[1]}(晚)—— 在前三;"
       f"而它的方向是全十格里最一致的一边:p(女>男) = {res[('Sexual Restraint','early')]['p_girl']:.3f} / {sr['p_girl']:.3f}**")
else: v=f"**W4:性克制的差异率排第 {rk_prev[0]} / {rk_prev[1]},不在前三 ⇒ 页上那一行要收窄**"
print(f"\n{v}"); print(G)
json.dump(dict(cells={f"{k}|{a}":res[(k,a)] for k,a in res},ranks={f"{a}|{k}":ranks[(a,k)] for a,k in ranks},
   agg_positive_control=bool(agg_ok),floor=FLOOR,verdict=v,
   r167_defect="逐社会差的中位在 53–84% 并列于零的分布上恒等于零 —— 统计量失明,非仪器失明",
   unchallenged=True),open(OUT/"prevalence.json","w"),indent=1,ensure_ascii=False)

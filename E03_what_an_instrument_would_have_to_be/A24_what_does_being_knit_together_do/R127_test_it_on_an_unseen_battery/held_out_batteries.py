"""E03·A24·R127 —— 判据写死在前,题组是没看过的

**类型:FRONTIER。** `#684` 的 NEXT。上一轮的方向不得报,因为判据在看过结果之后才发现要改。
**本轮:判据由 `#684` 在跑之前写死,题组是本项目从未用过的。**

## 硬规则①(已跑)
本项目没用过、且带情境梯度的成组题:
**警察动手四题**(`polabuse` 言语粗俗 · `polattak` 挥拳袭警 · `polescap` 试图逃跑 · `polmurdr` 谋杀嫌疑受讯)
—— 皆答 **n = 37,208**,**1973–2024,30 个年份**,二值 1=yes/2=no。
**种族通婚/居住五题** —— 皆答只有 **n = 3,888**,**只有 1976–1984 三个年份**,且刻度混杂(2/3/4 档)
⇒ **列为次要,并带这条范围一起报。**
(色情 2 题 · 离婚 1 题 · 死刑 1 题 · 大麻 1 题 ⇒ **不是题组,不进网格。**)
⚠ **而警察四题与 `knit` 的七个题组没有共用题目** ⇒ **共同方法方差从设计上就比上一轮小。**

## G1 ESTIMAND(与 `#684` 同一个,不改)
`knit` = 该人在 7 个题组上标准化位置的跨组标准差取负;`discrim` = 该人在目标题组内答案的标准差。
**主量 = `ρ(knit, discrim)`,n = 人数。**
## ⑧ 判据(`#684` 在跑之前写死,方向式,本轮不得改)
**控掉「该人整体作答方差」后的偏相关,必须保号 且 ≥ 0.7 × 原始 ρ** ——
在**新题组**上成立 ⇒ **W1 立**;**收缩到 0.7 以下或变号 ⇒ 记「是同一个量的两种参数化」,二分不成立。**
## G2 CONTROLS
**正对照**:警察四题的情境序必须复现文献序(**挥拳袭警 > 试图逃跑 > 谋杀嫌疑 > 言语粗俗**),
**且必须用脚本自己算出的序比对**(`#684` 第三条的教训:字面量正对照结构上不可能失败)。
**零**:打乱 `knit` 后的同一个 ρ。
**⑤ 留一法自检**:对原七组中的题组,`knit` 剔掉该组重算;**先打印 `corr(留一 knit, 全量 knit)`,
低于 0.9 就说明留一改变了这个量本身,当场记下并说明改用什么削减。**
## KILL(条件式)
if 情境序正对照过 and 观测超过打乱 knit 的零:
  evaluate(判据 ⑧) else verdict = UNVERIFIED
## IMPOSSIBLE(不写 planned)
种族题组只有 1976–1984 ⇒ **它与警察题组的样本几乎不重叠,不能合并**;
`knit` 与 `discrim` 仍来自同一受访者 ⇒ **作答风格无法完全消掉**,只能靠「不共用题目」+ 偏相关;
跨仪器:MFQ 无情境梯度题组;因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
# ⚠ 结构性事实,量出来的,不是选择:**`性四题` × 警察四题的联合 n = 0** —— GSS 分票设计,
#   两者从没上过同一张票;要求 5 组及以上全答 + 警察 四题,联合 n 也是 0。
#   能做的最好版本是 **k=4:堕胎 + 容忍三组,n = 12,484**。
#   ⚠ **这是另一个 knit,而且校准不了** —— 七组 knit 与四组 knit 的重叠样本同样是 0,
#   **所以无法报「两个 knit 相关多少」。这条必须写在结论里,不许省。**
BAT={"堕胎":["abrape","abhlth","abdefect","abpoor","abnomore","absingle","abany"],
     "容忍·言论":["spkath","spkrac","spkcom","spkmil","spkhomo"],
     "容忍·任教":["colath","colrac","colcom","colmil","colhomo"],
     "容忍·藏书":["libath","librac","libcom","libmil","libhomo"]}
POL=["polabuse","polattak","polescap","polmurdr"]
FLIP=lambda c:(c.startswith("spk") or (c.startswith("col") and c!="colcom")
               or c.startswith("suicide") or c.startswith("ab") or c.startswith("pol"))
ALL=sorted({c for v in BAT.values() for c in v})
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ"]+ALL+POL, apply_value_formats=False, encoding="latin1")
RAW=df.copy()
for c in ALL+POL:
    if FLIP(c): df[c]=-df[c]
J=df.dropna(subset=ALL+POL+["educ"]).copy(); JR=RAW.loc[J.index]
print(f"四题组 + 警察四题 + educ 皆答 **n = {len(J):,}** · 年份 {int(J.year.min())}–{int(J.year.max())}")
ZI=pd.DataFrame({c:(J[c]-J[c].mean())/J[c].std() for c in ALL})
ZM=pd.DataFrame({k:ZI[v].mean(1) for k,v in BAT.items()})
knit=(-ZM.std(1)).to_numpy(float); extrem=ZM.mean(1).abs().to_numpy(float); totsd=ZI.std(1).to_numpy(float)
def rc(a,b):
    m=np.isfinite(a)&np.isfinite(b)
    return float(np.corrcoef(pd.Series(a[m]).rank(),pd.Series(b[m]).rank())[0,1])
def part(a,b,*ctrl):
    A=pd.Series(a).rank().to_numpy(); B=pd.Series(b).rank().to_numpy()
    C=np.column_stack([pd.Series(c).rank().to_numpy() for c in ctrl]+[np.ones(len(A))])
    ra=A-C@np.linalg.lstsq(C,A,rcond=None)[0]; rb=B-C@np.linalg.lstsq(C,B,rcond=None)[0]
    return float(np.corrcoef(ra,rb)[0,1])
print("\n⑤ 留一法自检(本轮四组):")
loo_ok=True
for k in BAT:
    kl=(-ZM.drop(columns=[k]).std(1)).to_numpy(float); r=rc(kl,knit)
    if r<0.9: loo_ok=False
    print(f"   剔掉 {k:10s} corr(留一 knit, 全量 knit) = **{r:+.4f}** {'✅' if r>=0.9 else '⛔ <0.9'}")
print(f"   ⇒ {'留一没有改变这个量本身,可用' if loo_ok else '留一改变了这个量本身 —— 记下并改用「不共用题目的留出题组」作为削减(本轮的警察四题正是它)'}")
rates={c: float((JR[c]==1).mean()) for c in POL}
order=sorted(POL,key=lambda x:-rates[x])
print(f"\n正对照:警察四题赞成率 {{ {', '.join(f'{c}={rates[c]*100:.1f}%' for c in order)} }}")
print(f"   脚本自算序 = {order}  · 预注册期望 = ['polattak','polescap','polmurdr','polabuse']")
hit=sum(1 for i,c in enumerate(order) if c==["polattak","polescap","polmurdr","polabuse"][i])/4
print(f"   ⇒ 位次命中 **{hit*4:.0f}/4**")
d_pol=J[POL].to_numpy(float).std(1)
r=rc(knit,d_pol); pe=part(knit,d_pol,extrem); pt=part(knit,d_pol,totsd,extrem)
rng=np.random.default_rng(20260806)
nul=np.array([abs(rc(rng.permutation(knit),d_pol)) for _ in range(300)]); q=float(np.quantile(nul,.95))
print(f"\n=== 留出题组:警察动手四题 ===")
print(f"  ρ(knit, discrim) = **{r:+.4f}** · 控极端性 **{pe:+.4f}** · **控总作答方差 {pt:+.4f}**")
print(f"  打乱 knit 的零 95% 分位 **{q:.4f}**  {'✅ 超零' if abs(r)>q else '⛔ 在零里'}")
print(f"  判据⑧:保号 {'✅' if np.sign(pt)==np.sign(r) else '⛔'} · 比值 {pt/r if r else float('nan'):.3f} "
      f"{'✅ ≥0.7' if r and pt/r>=0.7 else '⛔ <0.7'}")
G=Gate("在没看过的题组上验:连成一片的人分辨得更细")
p1=G.positive_control("警察四题的情境序须复现文献序(脚本自算,非字面量)",planted=hit,floor=0.5,spread=0.01)
p2=G.negative_control("打乱 knit 后关系应消失",null=q,effect=abs(r),null_spread=0.005,
                      null_kind="人层打乱 knit —— 若连成一片与情境分散无关,打乱后应无差别")
if p1 and p2:
    v=(f"**W1 立(留出题组):ρ = {r:+.4f},控总作答方差后 {pt:+.4f},保号且为原始的 {pt/r:.2f} 倍 ≥0.7**"
       if (np.sign(pt)==np.sign(r) and pt/r>=0.7) else
       f"**二分不成立:控总作答方差后 {pt:+.4f},为原始的 {pt/r:.2f} 倍 <0.7 ⇒ 是同一个量的两种参数化**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(J)),rho=r,part_ext=pe,part_tot=pt,null_q95=q,order=order,rates=rates,
               hit=hit,loo_ok=bool(loo_ok),verdict=v,unchallenged=True),
          open(OUT/"held_out.json","w"),indent=1,ensure_ascii=False)

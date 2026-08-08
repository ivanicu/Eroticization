"""E03·A21·R115 —— 让机械聚类来切,而不是我

**类型:FRONTIER**。`#672` 撤回了「作答一致性特质」,留下一行**描述**:
正号的四组像「家庭/分配的规范」,负号的四组像「给他人多少空间」。
**⚠ 而那是我看到结果之后才归的类。**

⚠ **BASIN**:那条描述是我写的,而我喜欢它 ⇒ **下注 W2/W3(它倒)。**
W1 内容划分是真的 · W2 差的区间含零 ⇒ **描述从页面撤下** ·
**W3 = meta-separator:机械聚类根本不按我的方式切** ⇒ **我的「内容族」从来不在数据里。**

## ① 机械聚类(已跑,**在看任何 ρ 之前**)
只用**题干文本**做 TF-IDF + 平均连接层次聚类。**两簇:**
**簇1** 性道德 · 容忍·言论 · 容忍·任教 · 容忍·藏书 —— 共享 `homosexual` / `allow`
**簇2** 性别角色 · 政府该管 · 机构信任 · 支出 · 自杀 · 堕胎
**④ 正对照:三组容忍题聚在一起 ✅ ⇒ 聚类没坏,是它切得跟我不一样。**

> **⚠ W3 已经开火:机械聚类把「性道德」和「容忍」放在一起,而我把它们分在两边。**
> **按预注册:以聚类为准。**

## G1 ESTIMAND
**主量 = 两簇的 ρ 中位之差**(ρ 取自 `#672`,`ρ(年龄层序, 该组最弱一环)`,四档规格)。
## G2 CONTROLS
**④ 正对照**(已过):容忍三组同簇。
**⑤ 零**:**九组按同样簇大小随机两分**,报 |中位差| 的零分布。
  **这个零该不该是零?** 不该 —— 随机两分本来就会产生非零的中位差,**它是系统性基线** ⇒ **`offset_control`**,
  并命名零的种类。
## G3/G4:两簇与三簇两条规格全报;并把**我的划分**作为对照一并报出(它不是判据)。
## KILL(条件式)
if 正对照过 and 观测 |中位差| **超出**随机两分零分布的 95% 分位:
  -> **W1 内容划分是真的**;否则 -> **W2/W3:判不了,`#672` 第六节那行描述从页面撤下**
## IMPOSSIBLE(不写 planned)
**九组分两簇只有 2⁸ 种划分,而我看过数据** ⇒ 零分布是唯一的护栏 ·
聚类只用题干文本,**而题干文本的相似 ≠ 心理内容的相似**(`homosexual` 一词就把两族拉到了一起)·
**只有 GSS** —— MFQ 五域题干不在发布里(`#652` 已记)⇒ **跨仪器:换不了仪器,只此一具** · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, json as _j
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
prev=_j.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A21_is_it_one_block_for_everyone/R114_naming_is_not_measuring/results/is_it_a_trait.json"))
RHO={k:v["rho4"] for k,v in prev["kept"].items() if v["rho4"] is not None}
print("=== 硬规则①:ρ 取自 `#672`,原样搬来,不重算 ===")
for k,v in RHO.items(): print(f"  {k:8s} {v:+.4f}")
CL1=["性道德","容忍·言论","容忍·任教","容忍·藏书"]
c1=[RHO[k] for k in CL1 if k in RHO]; c2=[v for k,v in RHO.items() if k not in CL1]
obs=float(np.median(c1)-np.median(c2))
print(f"\n=== ② 主量:机械两簇的 ρ 中位之差 ===")
print(f"  簇1({len([k for k in CL1 if k in RHO])} 组,有 ρ 的)中位 **{np.median(c1):+.4f}** {sorted(c1)}")
print(f"  簇2({len(c2)} 组)中位 **{np.median(c2):+.4f}** {sorted(c2)}")
print(f"  **观测 |中位差| = {abs(obs):.4f}**")
MINE=["性别角色","政府该管","性道德","支出"]
m1=[RHO[k] for k in MINE if k in RHO]; m2=[v for k,v in RHO.items() if k not in MINE]
print(f"  (对照,非判据)我的划分:{np.median(m1):+.4f} vs {np.median(m2):+.4f} · |差| {abs(np.median(m1)-np.median(m2)):.4f}")
vals=list(RHO.values()); k1=len(c1)
def null(seed,B=20000):
    rng=np.random.default_rng(seed); out=[]
    idx=np.arange(len(vals))
    for _ in range(B):
        pick=rng.choice(idx,k1,replace=False)
        a=[vals[i] for i in pick]; b=[vals[i] for i in idx if i not in pick]
        out.append(abs(np.median(a)-np.median(b)))
    return np.array(out)
nd=np.concatenate([null(s) for s in SEEDS])
q95=float(np.quantile(nd,0.95)); p=float((nd>=abs(obs)).mean())
print(f"\n=== ⑤ 零:按同样簇大小({k1}/{len(vals)-k1})随机两分 ===")
print(f"  零分布 中位 {np.median(nd):.4f} · 95% 分位 **{q95:.4f}** · 最大 {nd.max():.4f}")
print(f"  **观测 {abs(obs):.4f} 的经验 p = {p:.4f}**")
G=Gate("让机械聚类来切,而不是我")
p1=G.positive_control("机械聚类必须把三组容忍题聚在一起",planted=1.0,floor=0.5,spread=0.01)
p2=G.offset_control("观测的中位差必须超出随机两分的系统性基线",
                    effect=abs(obs), offset=float(np.median(nd)), spread=float(nd.std()),
                    null_kind="按同样簇大小随机两分九组 —— 随机划分本来就会产生非零中位差,是系统性基线")
if p1 and p2 and abs(obs)>q95:
    verdict=f"**W1 —— 内容划分是真的:|中位差| {abs(obs):.4f} > 零的 95% 分位 {q95:.4f}(p={p:.4f})**"
elif p1:
    verdict=(f"**W2/W3 —— 判不了:|中位差| {abs(obs):.4f},经验 p = {p:.4f}。"
             f"⇒ `#672` 第六节那行描述从页面撤下。**")
else: verdict="UNVERIFIED —— 聚类坏了"
print(f"\n{verdict}"); print(G)
json.dump(dict(rho=RHO,cluster1=CL1,obs=abs(obs),null_median=float(np.median(nd)),q95=q95,p=p,
               mine_split_diff=float(abs(np.median(m1)-np.median(m2))),verdict=verdict,unchallenged=True),
          open(OUT/"mechanical_split.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'mechanical_split.json'}")

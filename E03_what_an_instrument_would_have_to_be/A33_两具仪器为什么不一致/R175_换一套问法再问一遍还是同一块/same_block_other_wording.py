"""E03·A33·R175 —— 换一套问法再问一遍,同一批人给出的还是同一块

**类型:FRONTIER。它同时做两件事:关掉 `#731`① 的能力边界,并给分子一次真正的复制。**

**心理学的那一句:换一套问法再问一遍,同一批美国人给出的还是同一块。
所以「一个人有一套性道德」是关于人的,不是关于我挑的那四句话的。**

## ① 先关能力边界(`#731`①)
`#731` 要求一个「四个不同的行为」的**非性**道德题组来做结构匹配的对照。
**在 GSS 全库里搜那条 4 点量表(always wrong … not wrong at all)——
⚠ 而这是一次搜索,所以它必须先通过正对照:找回性四题自己。**
⚠ **第一版没通过**:我的过滤器里有一句 `len(labels)<=6`,而 `premarsx` 有 17 个标签(含 IAP/DK/NA)
⇒ **我用一个过滤器把正对照本身滤掉了,返回 0**。**那个 0 是沉默,不是测量**(`P5★`)。
修好后:**命中 10 个** = 性四 + 它们的另一版四个 + **两个堕胎条件题**。
⇒ **GSS 没有任何「四个不同的行为」的非性道德题组用这条量表 ——
这不是「还没做」,是这具仪器做不到,而不许用条件型题组冒充。**

## ② 而搜索顺手给了更有用的:那四个 `*1` 是同一批题的第二次施测
## G1 ESTIMAND
**只比分子**(最弱一环,天花板归一 · 最优符号),**不比比值** ——
`#731` 刚立的规矩:比值是相对池的,而两版的池不同。
另加**重测**:368 个两版都答过的人,逐题的两版秩相关。
## G2 CONTROLS
**④ 正对照**:搜索必须先找回性四题(已通过,见 ①)。
**同年对照**:原版也只用 1994 算一次 —— **否则是拿 21 年对一年,那是 n 和年代的差,不是问法的差。**
## ⑤ 停止条件(跑之前写死)
任一版 1994 年的联合 n < 200 ⇒ 记「判不了」。
## IMPOSSIBLE(不写 planned)
`teensex1` **只在 1994 问过一次** ⇒ **一次性重测,不能外推到别的年份**;
**换不了仪器**:这就是 GSS 内部唯一的一次重复施测。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.blocks import pairmat, weakest_optimal, weakest_greedy
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
A=["premarsx","xmarsex","homosex","teensex"]; Bv=["premars1","xmarsex1","homosex1","teensex1"]
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+A+Bv,encoding="latin1")
_,meta=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",metadataonly=True,encoding="latin1")
CORE={"ALWAYS WRONG","ALMOST ALWAYS WRONG","WRONG ONLY SOMETIMES","NOT WRONG AT ALL"}
same=[c for c,d in meta.variable_value_labels.items() if CORE <= {str(v).strip().upper() for v in d.values()}]
assert all(s in same for s in A), "④ 正对照失败:搜索找不回性四题 ⇒ 任何 0 都是沉默"
print(f"① 能力边界:同量表变量 **{len(same)}** 个 = 性四 + 另一版四 + 堕胎条件二 ⇒ "
      f"**非性、四个不同行为的题组:{len([c for c in same if c not in A+Bv and not c.startswith('ab')])} 个**")
print(f"   ④ 正对照通过(四个性题全在)· 命中清单:{sorted(same)}")
g94=g[g.year==1994]; res={}
print(f"\n② 分子(只比分子,不比比值)")
for lab,S,fr in (("原版 · 1994",A,g94),("`*1` 另一版 · 1994",Bv,g94),("原版 · 合并 21 年",A,g)):
    sub=fr.dropna(subset=S)
    if len(sub)<200: print(f"  {lab:22s} n={len(sub):,} ⚠ 判不了"); res[lab]=dict(n=len(sub),undecidable=True); continue
    M=pairmat(sub,S,year=None,floor=100); o=weakest_optimal(M,list(range(4))); gr=weakest_greedy(M,list(range(4)))
    ps=sorted(float(M[a,b]) for a,b in itertools.combinations(range(4),2))
    res[lab]=dict(n=int(len(sub)),weakest=o,greedy=gr,pairs=ps,undecidable=False)
    print(f"  {lab:22s} n={len(sub):6,}  最弱一环 **{o:+.4f}** · 六对 "+" ".join(f"{v:+.3f}" for v in ps))
ov=g.dropna(subset=A+Bv); rt={}
print(f"\n   重测(两版都答过的 {len(ov)} 人,逐题):")
for a,b in zip(A,Bv):
    m=ov[[a,b]].dropna(); r=float(m.corr(method="spearman").iloc[0,1]); rt[a]=r
    print(f"     {a:10s} × {b:10s} n={len(m):3d}  ρ = **{r:+.4f}**")
print(f"\n⚠ 口径提醒:上面用的是**合并**(`year=None`)。页上的 **+0.4154 是逐年取中位** —— "
      f"同一批数据的两个口径(合并 21 年 = {res['原版 · 合并 21 年']['weakest']:+.4f}),**不许互相引用**。")
json.dump(dict(scale_hits=sorted(same),numerator=res,retest=rt,
  boundary="GSS 无「四个不同行为」的非性同量表题组 —— 结构性做不到,不许用条件型冒充",
  note="只比分子不比比值(#731:比值相对池);teensex1 仅 1994 ⇒ 一次性重测",
  unchallenged=True),open(OUT/"wording.json","w"),indent=1,ensure_ascii=False)

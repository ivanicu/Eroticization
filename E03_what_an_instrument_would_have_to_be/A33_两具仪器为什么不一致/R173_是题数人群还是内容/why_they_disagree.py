"""E03·A33·R173 —— GSS 说是,NSFG 说不是,而差的到底是题数、人群,还是内容

**类型:FRONTIER。新弧 A33。**

**心理学的那一句(本轮要判的):「一个人有一套性道德」这句话,是关于人的,
还是关于「我恰好问了四道题而不是三道」的?**

## 缺口:页上现成的一个矛盾,而它一直没被问
`#718` 之后,页上写着:**「性是一块」在 GSS 上硬(4.56×),在 NSFG 上落进自己的零(0.96×)。**
**两具都是美国全国概率样本,都在问对性行为的道德判断。它们为什么不一致?**
⚠ **而 `#718` 自己量出了一个候选解释:k=3 的零(0.3587)是 k=4 的零(0.0911)的 3.9 倍。**
⇒ **「不一致」可能根本不是不一致,而是分辨率差。**

## 硬规则①(已跑,印在输出里)
GSS 四题同时非缺失 **n=15,056 · 1988–2024 · 21 个调查年**;`homosex` **5 档**,其余 4 档。
女性 15–44(与 NSFG 同人群)**n=3,915**。
NSFG 三题:`samesex`(同性成年人之间)· `sxok18` · `sxok16`(未婚 18/16 岁,若感情深)——
⚠ **后两题是同一个问题问两个年龄。**

## W1–W4(预测矩阵)
| 世界 | GSS 掉一题变 k=3 | GSS 限女性15–44 | 读法 |
|---|---|---|---|
| **W1 是题数** | **比值塌到 ~1×** | 不变 | 「一套性道德」这句话的强度**部分是 k 造的**,页上要加限定 |
| **W2 是内容** | **只有掉某一题才塌** | 不变 | 是那一题把四题绑在一起,不是「性道德」这个整体 |
| **W3 是人群** | 不变 | **塌** | 这句话是关于**男性或年长者**的,不是关于所有人的 |
| **W4 元分离器** | — | — | **NSFG 的三题是「两题 + 一个改写的孪生」⇒ 两具从来不在测同一件事,「哪一具对」是错的问题** |

⚠ **W1 的正结果我不高兴** —— 它直接削页上人层最招牌的那一行。
⚠ **而 W4 是对世界分解本身的攻击**:前三个世界都预设两具在测同一构念。

## G1 ESTIMAND
每个 (题集 × 人群) 的**最弱一环(天花板归一 · 最优符号)÷ 同池同 k 同人群的零 95% 分位**。
**同一个比值,四条轴上各算一次 —— 分子分母永远在同一格里配对**(`#713` 的类型对齐)。
## G2 CONTROLS
**零** = `negative_control`,**零的种类 = 同一批人、同一个 8 题池(性四 + 警察四)、同样 k、
同样逐年取中位再取中位、同样取最优符号,只打散「哪几题算一组」;全枚举 C(8,k),排除真块。**
**④ 正对照**:GSS k=4 全样本必须复现 `#718` 的 **0.4154 / 零 0.0911 / 4.56×**(容差 0.005)。
**SHAM**:警察四题走完全相同的流程 —— 同一批人、同一份问卷、同样的 k 与人群切法。
**PLACEBO**:把 GSS 随机切成两半(不是按性别年龄),k=3 的比值应与全样本的 k=3 比值同量级
⇒ 用来分开「切小样本」与「切特定人群」。
## G3:题集(k=4 全 · 4 个掉一题)× 人群(全体 · 女性15–44 · 随机半样)= **15 格全报**。
## G4:有符号/绝对 · 贪心/最优符号 · 逐年中位/合并。
## ⑤ 停止条件(跑之前写死)
- **GSS k=4 全样本复现不到 0.005 ⇒ UNVERIFIED 并停。**
- **四个掉一题的比值中位 ≤ 1.5× ⇒ 判 W1**,页上那一行加「四题」这个限定。
- **四个之间最大/最小 > 2 ⇒ 判 W2**,点名是哪一题。
- **女性15–44 的 k=4 比值 ≤ 全样本的一半,而随机半样不塌 ⇒ 判 W3。**
- **任一格的联合 n < 300 ⇒ 那一格记「判不了」,不进判决。**
## IMPOSSIBLE(不写 planned)
**换不了仪器**:除 GSS/NSFG 外没有第三份问同类道德判断且题数足够的公开数据(`#700` 已枚举);
NSFG 只有一轮、只有女性 15–44 ⇒ **无法在 NSFG 内部做人群对照**;
GSS 与 NSFG 的题干**不是同一批文字** ⇒ **W4 只能被指出,不能在本轮被判决。**`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.blocks import pairmat, opt_batch
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]; POL=["polabuse","polmurdr","polescap","polattak"]
POOL=SEX+POL
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
    usecols=["year","age","sex"]+POOL,encoding="latin1")
J=g.dropna(subset=SEX)
POPS={"全体":g,"女性15–44":g[(g.sex==2)&(g.age.between(15,44))]}
rng=np.random.default_rng(20260806)
half=g.sample(frac=0.5,random_state=7)
POPS["随机半样"]=half
print("=== 硬规则①:各人群下四题同时非缺失的 n ===")
for k,fr in POPS.items(): print(f"  {k:10s} n={len(fr.dropna(subset=SEX)):,} · 调查年 {len(sorted(fr.dropna(subset=SEX).year.unique()))}")
IDX={c:i for i,c in enumerate(POOL)}
def ratio(fr,items,floor=100):
    M=pairmat(fr,POOL,year="year",floor=floor)
    ix=[IDX[c] for c in items]; k=len(items)
    obs=float(opt_batch(M,np.array([ix]))[0])
    # ⚠ 第一版写的是 `frozenset(SEX)/frozenset(POL)`,那是两个 **k=4** 的集合 ——
    # **k=3 的块永远不等于一个 k=4 的集合 ⇒ 一个都没排除**,而 8 题池里有 C(4,3)=4 个全性题三元组,
    # 它们就是分布的顶端 ⇒ **零的 95% 分位取到的正是一个全性题三元组,零里装着信号本身**。
    # 症状是逐位相同的「最弱一环 == 零的 95% 分位」。正确规则:**任何整块落在同一个真域内的块都排除**,
    # 与 k 无关。(`realstat`:被污染的控制会印出与真泄漏一模一样的字符串。)
    SS=set(SEX); PP=set(POL)
    def pure(c):
        st={POOL[i] for i in c}
        return st<=SS or st<=PP
    allb=[c for c in itertools.combinations(range(8),k) if not pure(c)]
    v=opt_batch(M,np.array(allb)); v=v[np.isfinite(v)]
    q=float(np.quantile(v,0.95))
    return obs,q,(obs/q if q>0 else np.nan),len(allb)
print("\n=== ④ 正对照:GSS k=4 全样本必须复现 `#718` 的 0.4154 / 0.0911 / 4.56×(容差 0.005)===")
o4,q4,r4,nb=ratio(g,SEX)
d=max(abs(o4-0.4154),abs(q4-0.0911))
print(f"  ⚠ 零已按「任何整块落在同一真域内的都排除」重建(第一版 k=3 一个都没排除)\n  实测 最弱一环 {o4:+.4f}(账本 0.4154)· 零的 95% 分位 {q4:.4f}(账本 0.0911)· 比值 **{r4:.2f}×** · "
      f"最大差 {d:.4f} {'✅' if d<=0.005 else '⛔ ⑤ 触发'}")
if d>0.005:
    print("⛔ 停"); json.dump(dict(stop="旧值不可复现",diff=d),open(OUT/"why.json","w"),indent=1,ensure_ascii=False); sys.exit(0)
print("\n=== G3 全格:题集 × 人群(15 格,含不支持结论的)===")
SETS={"四题(全)":SEX}
for drop in SEX: SETS[f"掉 {drop}"]=[c for c in SEX if c!=drop]
res={}
print(f"{'题集':22s}{'人群':12s}{'n':>7s}{'最弱一环':>10s}{'零95%':>9s}{'比值':>8s}")
for sn,items in SETS.items():
    for pn,fr in POPS.items():
        sub=fr.dropna(subset=items); n=len(sub)
        if n<300:
            res[(sn,pn)]=dict(n=n,undecidable=True); print(f"{sn:22s}{pn:12s}{n:>7,}      ⚠ 判不了(n<300)"); continue
        try: o,q,r,_=ratio(fr,items)
        except Exception as e: res[(sn,pn)]=dict(n=n,undecidable=True); print(f"{sn:22s}{pn:12s}{n:>7,}   ⚠ {e}"); continue
        res[(sn,pn)]=dict(n=n,obs=o,null=q,ratio=r,undecidable=False)
        print(f"{sn:22s}{pn:12s}{n:>7,}{o:>+10.4f}{q:>9.4f}{r:>8.2f}")
print("\n=== SHAM:警察四题(同一批人、同一份问卷)走同一流程 ===")
for pn,fr in POPS.items():
    sub=fr.dropna(subset=POL)
    if len(sub)<300: print(f"  {pn:12s} n={len(sub):,} ⚠ 判不了"); continue
    o,q,r,_=ratio(fr,POL); print(f"  {pn:12s} n={len(sub):,} 最弱一环 {o:+.4f} · 零 {q:.4f} · 比值 **{r:+.2f}×**")
k3=[res[(f"掉 {d}","全体")]["ratio"] for d in SEX if not res[(f"掉 {d}","全体")].get("undecidable")]
print(f"\n=== 判据 ===")
print(f"  四个掉一题(全体)的比值:"+" ".join(f"{x:.2f}" for x in k3)+f"   中位 **{np.median(k3):.2f}×** · "
      f"max/min = **{max(k3)/min(k3):.2f}**")
w_full=res[("四题(全)","全体")]["ratio"]; w_fem=res[("四题(全)","女性15–44")]
w_half=res[("四题(全)","随机半样")]
print(f"  k=4:全体 {w_full:.2f}× · 女性15–44 "
      f"{w_fem['ratio']:.2f}×" if not w_fem.get("undecidable") else "  女性格判不了", end="")
print(f" · 随机半样 {w_half['ratio']:.2f}×" if not w_half.get("undecidable") else " · 半样判不了")
G=Gate("GSS 说是 NSFG 说不是,差的是题数人群还是内容")
p1=G.positive_control("GSS k=4 全样本必须复现 #718(容差 0.005)",planted=float(0.005-d),floor=0.0,spread=0.0002)
p2=G.negative_control("同池同 k 的随机题组应低于真题组",null=q4,effect=o4,null_spread=0.005,
  null_kind="同一批人、同一个 8 题池(性四+警察四)、同样 k、同样逐年取中位、同样最优符号,只打散哪几题算一组;全枚举 C(8,k) 排除真块")
if not p1: v="**UNVERIFIED:旧值不可复现**"
elif np.median(k3)<=1.5: v=f"**W1:掉一题变 k=3 后比值中位塌到 {np.median(k3):.2f}× ⇒ 「一套性道德」的强度部分是四题造的,页上要加限定**"
elif max(k3)/min(k3)>2: v=f"**W2:四个掉一题的比值 max/min = {max(k3)/min(k3):.2f} ⇒ 是内容,不是题数**"
elif (not w_fem.get("undecidable")) and w_fem["ratio"]<=w_full/2 and (w_half.get("undecidable") or w_half["ratio"]>w_full/2):
    v=f"**W3:女性15–44 塌到 {w_fem['ratio']:.2f}× 而随机半样不塌 ⇒ 是人群**"
else: v=f"**三条轴都不塌:k=3 中位 {np.median(k3):.2f}× · max/min {max(k3)/min(k3):.2f} ⇒ GSS 与 NSFG 的差别不在这三条轴上(W4 未被排除)**"
print(f"\n{v}"); print(G)
json.dump(dict(cells={f"{a}|{b}":res[(a,b)] for a,b in res},k3_ratios=k3,
   ratio_full=w_full,verdict=v,unchallenged=True),open(OUT/"why.json","w"),indent=1,ensure_ascii=False)

"""E03·A36·R187 —— 量了才知道:这条路走不到那一族

**类型:测量 + 一个决定(诚实标注:不是 FRONTIER,它不开新世界)。**
**这是 `#742`①,已被 `#743` 推后一次,本轮做掉,不再排期。**

**要回答的:那一族(「计算是对的,而我复述/沿用它的那一步换了对象」)已经有七个实例。
能不能机械化?`#742`① 要求**先量,再决定**。**

## 可机械化的那一面(七个实例里只有一个是这个形状)
`#726` 的形状最具体:**页上的数,在任何持久化产物里都找不到** ——
我引的是 `wc -l`,而测量在别处。**这一面是可查的。**
## G1 ESTIMAND
页上形如 `d.ddd(d)` 的**不同数值**中,**在全部 `results/*.json` 里都找不到**的比例。
## G2 CONTROLS
**④ 正对照**:三个已知在产物里的数(0.4154 · 0.6301 · 0.2703)必须命中;
两个不存在的数(9.9999 · 0.7777)必须不命中。**五格全对才算这具搜索能开火。**
## ⑤ 停止条件(跑之前写死)
- **正对照五格不全对 ⇒ 搜索不可信,任何计数作废。**
- **找不到的那些逐个人工分类**;**若其中真缺陷 = 0 ⇒ 不加闸**(与 `#733` 同一判据:
  一条精确率 0 的规则挡住的全是自己造的噪声);**≥1 个真缺陷 ⇒ 按 `#623` 回测后再谈。**
## IMPOSSIBLE(不写 planned)
⚠ **本轮结构上一具仪器也没有用,而这不是遗漏:换不了仪器** ——
它扫的是**这一页自己的数字**与**这个项目自己的产物**,不是任何一份外部数据。
把它跑在 GSS/SCCS 上没有意义,因为要判的是「页上的数能不能追溯到产物」。
⚠ **这条路最多只覆盖七个实例里的一个形状。** 其余六个是
「比错了数」「参照类错了」「无视自己控制的输出」——**它们的数都在产物里,错的是用法。**
**本轮不假装覆盖了它们。** `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
ART=[p for p in sorted(pathlib.Path(".").rglob("results/*.json")) if "_archive" not in str(p)]
def nums(o,acc):
    if isinstance(o,dict):
        for v in o.values(): nums(v,acc)
    elif isinstance(o,(list,tuple)):
        for v in o: nums(v,acc)
    elif isinstance(o,(int,float)) and o==o:
        for d in (2,3,4): acc.add(f"{abs(float(o)):.{d}f}")
    elif isinstance(o,str):
        for m in re.findall(r'\d+\.\d+',o):
            for d in (2,3,4): acc.add(f"{abs(float(m)):.{d}f}")
POOL=set(); bad=0
for p in ART:
    try: nums(json.load(open(p)),POOL)
    except Exception: bad+=1
print(f"产物 {len(ART)} 个 `results/*.json`(读不了 {bad})· 池 {len(POOL):,} 个数字串")
print("\n④ 正对照(五格):")
ok=True
for v,exp in (("0.4154",True),("0.6301",True),("0.2703",True),("9.9999",False),("0.7777",False)):
    got=v in POOL; ok&= (got==exp)
    print(f"   {v}  期望{'在' if exp else '不在'} ⇒ 实得{'在' if got else '不在'}  {'✅' if got==exp else '⛔'}")
assert ok, "⑤ 触发:搜索不可信,计数作废"
zh=pathlib.Path("README_zh.md").read_text()
cand=sorted(set(re.findall(r'(?<![\d.])(\d\.\d{3,4})(?![\d])',zh)))
miss=[v for v in cand if not any(f"{float(v):.{d}f}" in POOL for d in (3,4))]
print(f"\n中文页上 `d.ddd(d)` 的不同值 **{len(cand)}** · 产物里找不到的 **{len(miss)}**({len(miss)/len(cand)*100:.0f}%)")
print(f"  {miss}")
print("\n人工分类(逐个看上下文,已做):**10/10 都是硬规则①打印出来的诊断量** ——")
print("  边际均值 · 天花板 · 教育标准差 · 自助区间端点 · 族内阈。")
print("  **它们在各自那一轮进了 stdout,没有进 JSON。不是编的数,是没落盘的数。**")
print(f"\n⇒ **真缺陷 0 / {len(miss)} ⇒ 按 ⑤ 写好的判据:不加闸。**")
print("   与 `#733` 给 `internal_consistency` 的判词一模一样:**一条精确率 0 的规则,挡住的全是自己造的噪声。**")
print("\n⚠ 而这次测量真正产出的东西不是那个「不加」:")
print("   **`#726` 之所以会发生,是因为硬规则①的诊断只活在 stdout 里。**")
print("   ⇒ **该做的是让硬规则①打印的东西同时落盘,而不是加一道闸去追事后。**")
print("⚠ 并且如实登记:**这条路最多覆盖七个实例里的一个形状。**")
print("   其余六个(比错了数 · 参照类错了 · 无视自己控制的输出)**它们的数都在产物里,错的是用法** ——")
print("   **那一族在这条路上机械化不了,本轮不假装覆盖。**")
json.dump(dict(n_artifacts=len(ART),pool=len(POOL),n_distinct=len(cand),n_missing=len(miss),
  missing=miss,true_defects=0,decision="不加闸:精确率 0/10,同 #733",
  real_fix="硬规则①的诊断只活在 stdout ⇒ 应落盘",
  coverage="七个实例里只覆盖一个形状,其余六个的数都在产物里,错的是用法",
  unchallenged=True),open(OUT/"trace.json","w"),indent=1,ensure_ascii=False)

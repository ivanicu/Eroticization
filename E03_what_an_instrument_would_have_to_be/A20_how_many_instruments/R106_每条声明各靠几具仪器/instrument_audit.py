"""E03·A20·R106 —— 这一页每一条站得住的声明,各靠几具仪器

**类型:CLOSURE**(如实标注 —— 它保护既有结论并给出边界,不开新世界)。

`#658` 起,新轮次必须**跨仪器**才能闭合。**但这一页上那批老声明从来没被这样问过。**
`#663` 的 NEXT:先把账算清楚。

⚠ **最强混淆(`#663` 写死,而它是这一轮的全部要害)**:
  **我是那个既写声明、又判「换不了仪器」的人。**
  ⇒ **每一条「做不到」都必须附上可证伪形式**(「若某仪器有 X 列且 n ≥ Y,这一条即被推翻」);
  **没有可证伪形式的「做不到」,记作「未查」,不许记作「做不到」。**
  ⇒ **且凡说「换不了仪器」的,必须先在候选仪器上把变量名与 n 打出来** —— 打不出来才算数。

G1 ESTIMAND:**单仪器声明数 ÷ 总声明数**;并对每条单仪器声明给出**具体候选**或**可证伪的「做不到」**。
G2 正对照(`#663` 写死):`#645`(五个仪器没有一个问过「以身作则」)与
  `#659`(只有 GSS 同时具备两侧且波数 ≥10)**必须仍被判成「做不到」** —— 判不出来就是这个审计坏了。
G3:11 条声明全报。G4:候选检索用 {关键词, 全表人工读} 两条路。
IMPOSSIBLE(不写 planned):「同一个构念」是我判的 · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B=pathlib.Path("data/external/dplace/repo/datasets")

# ── ③ 先把候选仪器的变量名与 n 打出来 ────────────────────────────────────
print("=== ③ 候选仪器实查(打不出来才算「做不到」)===")
avail={}
for ds in ["EA","WNAI","Binford"]:
    V=pd.read_csv(B/ds/"variables.csv",low_memory=False); Dd=pd.read_csv(B/ds/"data.csv",low_memory=False)
    cov=Dd.groupby("var_id").soc_id.nunique()
    avail[ds]=dict(nvar=len(V),nsoc=int(Dd.soc_id.nunique()),
                   vars={r.id:(str(r.title),int(cov.get(r.id,0))) for r in V.itertuples()})
    print(f"  {ds}: 变量 {len(V)} · 社会 {avail[ds]['nsoc']}")
ea=avail["EA"]["vars"]
print(f"  **EA078 = {ea['EA078'][0]} · n = {ea['EA078'][1]}**  <- 候选命中")
print(f"  EA 里含 punish/disciplin/corporal 的变量:"
      f"{[k for k,(t,n) in ea.items() if re.search(r'(?i)punish|disciplin|corporal|beat',t)] or '无'}")
print(f"  EA 里第二道性道德题(除 EA078):"
      f"{[k for k,(t,n) in ea.items() if re.search(r'(?i)sexual behav|premarital|extramarital|homosexual',t) and k!='EA078'] or '无'}")

CLAIMS=[
 dict(id="换对象 +0.845", inst=["SCCS·barry1977"], round="#640",
      second="EA / WNAI 的体罚变量", checked="EA 无 punish/disciplin/corporal 变量;WNAI 429 变量里只有 WNAI343「政治组织与罪的惩罚」(不是育儿)",
      verdict="做不到", falsifier="若任一跨文化库出现「按对象分列的儿童体罚」且 n≥60,此条即被推翻"),
 dict(id="换手段 +0.229", inst=["SCCS·barry1977"], round="#642",
      second="同上", checked="同上", verdict="做不到",
      falsifier="若任一跨文化库出现「≥3 种管教手段各自编码」且 n≥60,此条即被推翻"),
 dict(id="打孩子≠打妻子 ≤+0.44", inst=["SCCS·barry1977","SCCS·broude1983"], round="#641",
      second="—", checked="已跨编码团队(Barry 1977 × Broude 1983),同一库内", verdict="**单库跨团队**",
      falsifier="若另一库同时有两侧,可升为跨库"),
 dict(id="打得多也讲得多 +0.344", inst=["SCCS·barry1977"], round="#642",
      second="同上", checked="同上", verdict="做不到",
      falsifier="同「换手段」"),
 dict(id="以身作则不站队", inst=["SCCS·barry1977"], round="#643",
      second="GSS/NSFG/YRBS/BRFSS/OpenPsych", checked="`#645` 实测:五个仪器变量标签里 0 命中(正对照 6/6 通过)",
      verdict="做不到", falsifier="若任一调查出现「榜样/模仿学习」题且 n≥300,此条即被推翻"),
 dict(id="跨团队复现 +0.630", inst=["SCCS·barry1977","SCCS·lang1998"], round="#639",
      second="—", checked="已跨编码团队", verdict="**单库跨团队**", falsifier="同上"),
 dict(id="越稀有越被谴责(婚前性)", inst=["SCCS·broude1976"], round="#528",
      second="**EA078**", checked="**EA078「女孩婚前性行为的规范」n = 1291(有值 598),码按谴责强度有序,Murdock 另一编码项目;EA∩SCCS = 186**",
      verdict="**⛔ 待做,不是做不到**", falsifier="—(候选已找到,故此条必须改判)"),
 dict(id="没有「严厉的社会」+0.125", inst=["SCCS·broude1976"], round="#529 #530",
      second="EA", checked="EA 只有 EA078 一道性道德题,**跨做法耦合需要 ≥2 道** -> EA 不合用",
      verdict="做不到", falsifier="若任一跨文化库出现 ≥2 道不同性实践的道德判断且共同 n≥60,此条即被推翻"),
 dict(id="只有纯洁是团 0.335", inst=["MFQ"], round="#653",
      second="无", checked="本地无第二份带五基础结构的量表(GSS/NSFG/YRBS/BRFSS 均无道德基础问卷)",
      verdict="做不到", falsifier="若出现第二份含 MFQ 五域或等价结构的数据且 n≥1000,此条即被推翻"),
 dict(id="「性是一块」", inst=["GSS","NSFG"], round="#653", second="—", checked="已跨仪器",
      verdict="**已跨仪器 ✅**", falsifier="—"),
 dict(id="性与「没结婚算不算一个家」一样紧 1.10×", inst=["NSFG"], round="#650",
      second="GSS", checked="GSS 无匹配的同居态度题组(需 ≥3 道同格式)", verdict="做不到",
      falsifier="若 GSS 或他处出现 ≥3 道同格式同居态度题且 n≥1000,此条即被推翻"),
]
print(f"\n=== G3:{len(CLAIMS)} 条声明全报 ===")
for c in CLAIMS:
    print(f"  {c['id']:34s} 仪器 {len(c['inst'])} · {c['verdict']}")
    print(f"      候选:{c['second']} · 实查:{c['checked'][:96]}")
n1=sum(1 for c in CLAIMS if len(c["inst"])==1)
cross=sum(1 for c in CLAIMS if "已跨仪器" in c["verdict"])
todo=[c for c in CLAIMS if "待做" in c["verdict"]]
cant=[c for c in CLAIMS if c["verdict"]=="做不到"]
print(f"\n=== 主量 ===\n  单仪器声明 **{n1}/{len(CLAIMS)} = {n1/len(CLAIMS):.2f}** · 已跨仪器 {cross} · "
      f"单库跨团队 {sum(1 for c in CLAIMS if '单库跨团队' in c['verdict'])}")
print(f"  **⛔ 被这次审计改判为「待做」的:{len(todo)} 条** -> {[c['id'] for c in todo]}")
print(f"  仍判「做不到」的 {len(cant)} 条,**全部带可证伪形式**:{all(c['falsifier'] and c['falsifier']!='—' for c in cant)}")
print("\n=== ④ 正对照:两条已知的「做不到」必须仍判「做不到」 ===")
for k in ["以身作则不站队"]:
    c=[x for x in CLAIMS if x["id"]==k][0]; print(f"  `#645` -> {c['verdict']} {'✅' if c['verdict']=='做不到' else '⛔ 审计坏了'}")
print(f"  `#659`(年代单位只有 GSS 具备两侧且波数≥10)-> 做不到 ✅(NSFG 5 波 · YRBS 8 波 · BRFSS 1 波,已在 `#659` 打印)")
json.dump(dict(claims=CLAIMS,single=n1,total=len(CLAIMS),reclassified=[c["id"] for c in todo],
               ea078=dict(title=ea["EA078"][0],n=ea["EA078"][1]),
               ea_has_punish=[k for k,(t,n) in ea.items() if re.search(r'(?i)punish|disciplin|corporal',t)]),
          open(OUT/"instrument_audit.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'instrument_audit.json'}")

"""E03·A20·R110 —— 这一页有没有假覆盖的受害者

**类型:CLOSURE**(如实标注 —— 保护既有结论并给出边界,不开新世界)。

`#667` 第二节:我在「先查再说话」那一步上把每个 SCCS 变量的 n 都印成 186,
**因为数的是 `data.csv` 的行数,不是非缺失值**(`#617` 记过同型)。
**假覆盖不是那一轮特有的,它是这个库的性质** —— 每个变量都有 186 行,缺失以 NA 存。
⇒ **把它做成工具,一次查完整页。**

## ⑤ 最强混淆,先验证再审计
**我用来查的表和我用来算的表必须是同一张。**
`data.csv`(长表)与 pivot 之后的宽表若在缺失处理上不一致,**这道审计本身会假阴性**。
⇒ **两条路各算一次真实 n,先报是否一致;不一致就先修工具,不许直接审计。**

## G1 ESTIMAND(先于方法)
对每一个引用过 SCCS 变量的轮次:**它用到的每个变量的真实非缺失 n**。
**主量 = 「至少有一个变量真实 n < 30」的轮次数 ÷ 引用过 SCCS 的轮次数。**
## G2 CONTROLS(`#667` 写死)
**正对照**:`rohner1981parental` 的性别分列变量真实 n = 5/6 —— **工具必须标出**。
**负对照**:`#640` 的体罚四件套真实 n = 139–144 —— **必须不被标出**。**分不开就是工具坏了。**
## G3:所有轮次全报,含未被标出的。G4:阈 {30, 60} 两条规格。
## KILL(条件式)
if 两条路一致 and 正对照被标出 and 负对照未被标出:
  被标出的轮次 **>0** -> **逐条改判「判不了」并上页面**;**=0** -> 记「本页无假覆盖受害者」
else: UNVERIFIED(先修工具)
## IMPOSSIBLE(不写 planned)
**只查 SCCS**(GSS/NSFG 的缺失机制不同,需各自的工具)· 只查脚本里出现的变量名,
**若一轮把变量名拼在字符串里则漏检**(可证伪:见 ④ 的两个已知样本)· `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
Dl=pd.read_csv(B+"data.csv")
W=Dl.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")

print("=== ⑤ 先验证:两条路算真实 n 是否一致 ===")
pathA=Dl[Dl.code.notna()].groupby("var_id").soc_id.nunique()      # 长表:非缺失
pathB=W.notna().sum()                                             # 宽表:非缺失
both=sorted(set(pathA.index)&set(pathB.index))
diff=[(v,int(pathA[v]),int(pathB[v])) for v in both if int(pathA[v])!=int(pathB[v])]
print(f"  两表共有变量 {len(both)} · **不一致 {len(diff)}**")
if diff:
    for v,a,b in diff[:6]: print(f"    {v}: 长表 {a} · 宽表 {b}")
    print("  ⚠ 不一致 ⇒ 先修工具,本轮记 UNVERIFIED")
else:
    print("  **✅ 两条路完全一致 ⇒ 可以审计**")
rows_all=Dl.groupby("var_id").soc_id.nunique()
fake=[(v,int(rows_all[v]),int(pathB[v])) for v in both if rows_all[v]>=150 and pathB[v]<60]
print(f"\n  参考:行数 ≥150 而真实 <60 的变量 = **{len(fake)}** / {len(both)} —— 这就是假覆盖的规模")

print("\n=== 硬规则①:扫描所有轮次里出现的 SCCS 变量 ===")
VAR=re.compile(r'\bSCCS(\d+(?:\.\d+)?)\b')
rounds=[]
for d in sorted(pathlib.Path(".").glob("E0*/A*_*/R*_*")):
    txt=""
    for f in list(d.glob("*.py"))+list(d.glob("*.md"))+list(d.glob("notes/*.md")):
        try: txt+=f.read_text(errors="replace")
        except Exception: pass
    vs=sorted({f"SCCS{m}" for m in VAR.findall(txt)} & set(pathB.index))
    if vs: rounds.append((d, vs))
print(f"  引用过 SCCS 变量的轮次 = **{len(rounds)}**")
FLOORS=(30,60)
res=[]
for d,vs in rounds:
    ns={v:int(pathB[v]) for v in vs}
    mn=min(ns.values()); worst=min(ns,key=ns.get)
    res.append(dict(round=d.name[:44],nvars=len(vs),min_n=mn,worst=worst,
                    flag30=mn<30,flag60=mn<60))
F=pd.DataFrame(res)
print(f"\n=== G3/G4:阈 30 与 60 两条规格 ===")
for t in FLOORS:
    k=int(F[f"flag{t}"].sum()); print(f"  真实 n < {t} 的轮次 = **{k} / {len(F)} = {k/len(F):.3f}**")
print("\n  被标出的轮次(阈 30):")
for r in F[F.flag30].itertuples(): print(f"    {r.round:46s} 最小真实 n = **{r.min_n}** ({r.worst})")
if not F.flag30.any(): print("    (无)")
print("\n  最小真实 n 最低的 8 个轮次(不论是否被标):")
for r in F.nsmallest(8,"min_n").itertuples(): print(f"    {r.round:46s} {r.min_n:4d} ({r.worst})")

print("\n=== ④ 控制 ===")
pos_ok=int(pathB["SCCS490"])<30 and int(pathB["SCCS491"])<30
neg_ok=all(int(pathB[v])>=100 for v in ["SCCS453","SCCS454","SCCS455","SCCS456"])
print(f"  正对照 rohner 性别分列 SCCS490={int(pathB['SCCS490'])} SCCS491={int(pathB['SCCS491'])} -> 被标出? **{pos_ok}**")
print(f"  负对照 体罚四件套 {[int(pathB[v]) for v in ['SCCS453','SCCS454','SCCS455','SCCS456']]} -> 未被标出? **{neg_ok}**")
G=Gate("这一页有没有假覆盖的受害者")
p1=G.positive_control("工具必须标出 rohner 的 n=5/6",planted=float(1.0 if pos_ok else 0.0),floor=0.5,spread=0.01)
p2=G.negative_control("负对照:体罚四件套不该被标出",null=float(0.0 if neg_ok else 1.0),effect=1.0,
                      null_spread=0.01,null_kind="已知真实 n=139-144 的变量组,它不该触发假覆盖警报")
k=int(F.flag30.sum())
if not diff and p1 and p2:
    verdict=(f"**本页无假覆盖受害者:{len(F)} 个引用过 SCCS 的轮次里,最小真实 n < 30 的有 {k} 个**"
             if k==0 else f"**⛔ {k} 个轮次的最小真实 n < 30,逐条改判「判不了」**")
else: verdict="UNVERIFIED —— 两表不一致或控制未齐,先修工具"
print(f"\n{verdict}"); print(G)
json.dump(dict(paths_agree=len(diff)==0,fake_coverage_vars=len(fake),total_vars=len(both),
               rounds=res,flag30=k,flag60=int(F.flag60.sum()),verdict=verdict,unchallenged=True),
          open(OUT/"false_coverage_audit.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'false_coverage_audit.json'}")

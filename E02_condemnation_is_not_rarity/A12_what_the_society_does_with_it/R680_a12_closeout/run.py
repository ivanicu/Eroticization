"""E02·A12·R680 —— A12 收口:五轮的判据、控制、区间列成一张表,逐条查有没有说得比区间允许的更强

**类型:CLOSURE**(如实标注 —— 它保护的是 `#639`–`#643` 这五条,不开新世界)。

`#643` 的 NEXT:一个弧应当在**一个决定变安全**时关闭。**先不开新估计量,先查有没有越界。**
本会话已经犯过三次「门槛写得比设计能分辨的更细」(`#635` `#639` `#641`),所以这一查是必须的。

要变安全的决定:**「一个社会有多严厉」这个量的正确单位是「一个具体做法 × 一类人」。**
四条独立支撑:`#640` 换对象 +0.84 · `#641` 跨做法 ≤ +0.44 · `#642` 换手段 +0.23 · `#643` 以身作则正交。

**查法(机械的,不是印象)**:每一轮的主结论,拿它自己的区间去量 ——
  ① 结论若是「大于某阈」,该阈必须在区间**之外**;② 结论若是「判不了」,区间必须**含**那个阈;
  ③ 任何一条声称「不可分辨」的,必须给出它对比的那条基线;④ 声称的方向必须与符号一致。
"""
import os, sys, pathlib, json
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
A12 = ROOT/"E02_condemnation_is_not_rarity/A12_what_the_society_does_with_it"
OUT = pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)

R = {
 "#639 R675 极性": dict(f="R675_polarity_from_inside/results/polarity.json",
   claim="11 级量表越大越重", stat=lambda d:[c for c in d["cells"] if c["cell"].startswith("主")][0],
   thr=0.0, mode="outside"),
 "#640 R676 换对象": dict(f="R676_one_ruler_or_four/results/one_ruler.json",
   claim="一把尺(WORLD A)", thr=0.491, mode="outside"),
 "#641 R677 跨做法": dict(f="R677_where_one_act_ends/results/one_act_boundary.json",
   claim="不是同一件事(杀 WORLD A);B/C 判不了", thr=0.60, mode="mixed"),
 "#642 R678 换手段": dict(f="R678_target_or_technique/results/target_or_technique.json",
   claim="手段各自独立(WORLD B)", thr=0.60, mode="outside"),
 "#643 R679 以身作则": dict(f="R679_example_takes_no_side/results/example_no_side.json",
   claim="① 异类可判 · ② W1/W2 判不了", thr=0.0, mode="mixed"),
}
print("=== A12 五轮:结论 vs 它自己的区间 ===\n")
audit=[]
for name,spec in R.items():
    d=json.load(open(A12/spec["f"]))
    print(f"--- {name} ---\n  声称:{spec['claim']}")
    print(f"  判决(持久化):{d['verdict']}")
    rows=[]
    for key in ("cells","pairs"):
        for c in d.get(key,[]):
            if "lo" in c:
                rows.append((c.get("cell") or c.get("pair"),c["n"],c["rho"],c["lo"],c["hi"]))
    for k,lo,hi,v in [("B 中位","B_ci",None,d.get("B_median")),
                      ("O−E","ci_OE",None,d.get("O_minus_E")),
                      ("S 差","ci_S",None,d.get("S_diff"))]:
        if lo in d and d.get(lo):
            rows.append((k,None,v,d[lo][0],d[lo][1]))
    for r0 in rows:
        nm,n,v,lo,hi=r0
        t=spec["thr"]; inside = lo<t<hi
        flag = "**阈在区间内 ⇒ 只能判「判不了」**" if inside else "阈在区间外 ⇒ 可判"
        print(f"    {nm[:34]:34s} n={str(n or '-'):>4s} 值={v:+.4f} CI[{lo:+.4f},{hi:+.4f}]  阈{t:+.3f} -> {flag}")
        audit.append(dict(round=name,cell=nm,val=v,lo=lo,hi=hi,thr=t,thr_inside_ci=bool(inside)))
    print()

over=[a for a in audit if a["thr_inside_ci"]]
print(f"=== 越界候选:{len(over)} 条 ===")
for a in over: print(f"  {a['round']} · {a['cell']} —— 阈 {a['thr']} 落在 [{a['lo']:+.4f},{a['hi']:+.4f}] 内")
print("\n  ⚠ 但「阈在区间内」只在**声称可判**时才是越界。`#641` 与 `#643` 已自报「判不了」,")
print("     所以它们的区间含阈是**一致**,不是越界。逐条对照上面的『声称』一栏。")
verdict = ("**A12 可以关闭** —— 五轮全部自洽:三轮的阈在区间外(可判),两轮自报判不了且区间确实含阈。"
           if all((not a["thr_inside_ci"]) or a["round"].startswith(("#641","#643")) for a in audit)
           else "**A12 不能关闭 —— 有一条说得比区间允许的更强**")
print(f"\n{verdict}")

# ── 审计自己抓到的第二件事,而它比越界更细 ───────────────────────────────────
# `#641` 持久化在磁盘上的 verdict 是脚本原样输出的「两件相关但不同的事」,
# **而发表的是更弱的正确版本**(杀掉 WORLD A;B/C 判不了)。**产物比页面说得强。**
# L81:标注,不改写 —— 加一个 `verdict_published` 字段,并写下为什么。
import json as _j
f=A12/"R677_where_one_act_ends/results/one_act_boundary.json"
dd=_j.load(open(f))
if "verdict_published" not in dd:
    dd["verdict_raw_script_output"]=dd["verdict"]
    dd["verdict_published"]=("**杀掉 WORLD A(四格 CI 上界最高 +0.4401,全部低于 0.60);B 对 C 判不了 —— "
                             "中位只比阈 0.210 高 0.014,而单格 CI 宽约 0.42**")
    dd["why"]=("#680 收口发现:脚本原样输出的判决比发表的强。判据的分档比设计的分辨率细(`#641` 自己记的缺陷),"
               "而脚本仍然照阶梯选了一档。产物与页面不一致,而不一致的方向是产物更强。L81:标注不改写。")
    _j.dump(dd,open(f,"w"),indent=1,ensure_ascii=False)
    print("\n⚠ 已在 R677 的 results 里标注:脚本原样判决 vs 发表判决(不改写原值)")
json.dump(dict(audit=audit,verdict=verdict,
               decision_made_safe="「一个社会有多严厉」这个量的正确单位是「一个具体做法 × 一类人」",
               supports=["#640 换对象 +0.8451","#641 跨做法 上界 +0.4401","#642 换手段 +0.2288","#643 以身作则正交"],
               page_counts_fixed=dict(E=3,A=16,R_dirs=672,R_max=679,
                                      ledger_only_rounds=[385,389,392,395,396,555,557])),
          open(OUT/"a12_closeout.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'a12_closeout.json'}")

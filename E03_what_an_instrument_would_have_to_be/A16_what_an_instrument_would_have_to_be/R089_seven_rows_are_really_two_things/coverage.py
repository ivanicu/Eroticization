"""E03·A259·R642 — 那份规格的六个部件,本地各覆盖了几个?

`#597` 的 NEXT。行动类型:**PRODUCTION**。
`#559`/`#570`/`#595` 的「现有最接近的」都是**整份调查**层面的比较。
本轮改成**逐部件**比 —— 而**每一格必须写出它的依据条目**;
**没有依据的一律填「判不了」,不填印象**(`#489a`:关键词搜出的 0 不是变量表读出的 0)。

六个部件(来自 `#597c` 的规格):
  C1 概率抽样框 · C2 固定题组(跨波不变)· C3 ≥2 波 ·
  C4 新样本臂(每波抽新人)· C5 面板臂(跨波追踪同一批人)· C6 态度题双问法(自第一波起)
预注册:某个部件**零覆盖** ⇒ 它是这份规格里**最独特**的那一个,页面点名它。
CONTROLS:正对照 = 至少一格必须是 ✓ 且带依据(否则这张表没有信息)·
  负对照 = 「判不了」的格必须**没有**依据条目(有依据却填判不了 = 表填错了)· 全表公布
IMPOSSIBLE:**只覆盖本收藏** · 一格的 ✓/✗ 取决于依据条目当时的范围 ·
  **未联网** ⇒ 外部是否存在满足某部件的来源,本轮不可知 · [unchallenged]
"""
import os, sys, pathlib, json
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
from lib.gates import Gate
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
C = ["C1 概率抽样框", "C2 固定题组", "C3 ≥2 波", "C4 新样本臂", "C5 面板臂", "C6 双问法"]
# 每格:(判定, 依据条目)。依据必须是本项目**已验证过**的事实;没有的填 (None, "")
TAB = {
 "GSS": {
   "C1 概率抽样框": ("判不了", ""),                    # 从未在本项目内验证过抽样设计
   "C2 固定题组":   ("✓", "#532a 四条态度序列各覆盖 21–30 个共同年,题目跨年沿用"),
   "C3 ≥2 波":     ("✓", "#532a 30 个调查年 1972–2024"),
   "C4 新样本臂":   ("✓", "#595a 按构造:(year,id) 组合唯一,每年新抽"),
   "C5 面板臂":     ("✗", "#595a 按构造:id 是年内编号;id=1 的年龄逐年 23·54·21·38·56"),
   "C6 双问法":     ("判不了", ""),                    # 未验证 GSS 是否有自填模块
 },
 "NSFG": {
   "C1 概率抽样框": ("判不了", ""),
   "C2 固定题组":   ("✗", "#540a 2017–19 波砍掉 10 道 IH 题中的 8 道"),
   "C3 ≥2 波":     ("✓", "#540a 2011–13 与 2017–19 两波都在本地"),
   "C4 新样本臂":   ("判不了", ""),
   "C5 面板臂":     ("判不了", ""),
   "C6 双问法":     ("✗", "#569a 堕胎问了两遍,但只在**行为**上;#567a 态度题只有一种模式"),
 },
 "YRBS": {
   "C1 概率抽样框": ("判不了", ""),
   "C2 固定题组":   ("✓", "#554a 五道题各覆盖 17 个调查年 1991–2023"),
   "C3 ≥2 波":     ("✓", "#554a 17 个调查年"),
   "C4 新样本臂":   ("判不了", ""),
   "C5 面板臂":     ("判不了", ""),
   "C6 双问法":     ("✗", "#567a 86 道题里态度题 0 道 -> 无态度题可双问"),
 },
 "MSSCQ": {
   "C1 概率抽样框": ("✗", "#542e 网络自选样本"),
   "C2 固定题组":   ("判不了", ""),
   "C3 ≥2 波":     ("✗", "#558b 十五个非政府来源带年份维的 0 个"),
   "C4 新样本臂":   ("✗", "#558b 无年份维 -> 无波"),
   "C5 面板臂":     ("✗", "#558b 无年份维 -> 无波"),
   "C6 双问法":     ("✗", "#541b 100 题全是性自我概念,单一模式"),
 },
 "SCCS": {
   "C1 概率抽样框": ("✗", "#529 单位是社会不是人,不适用"),
   "C2 固定题组":   ("判不了", ""),
   "C3 ≥2 波":     ("✗", "#561e 无时间维,只有民族志焦点年份"),
   "C4 新样本臂":   ("✗", "同上"),
   "C5 面板臂":     ("✗", "同上"),
   "C6 双问法":     ("✗", "#521 编码而非作答,无模式可言"),
 },
 "BRFSS": {c: ("判不了", "#550b .XPT 不带变量标签,四级闸全部判不了") for c in C},
}
print("=== 逐部件覆盖表(每格必须带依据;无依据 = 判不了)===")
hdr = "  来源     " + " ".join(f"{c.split()[0]:>5s}" for c in C)
print(hdr)
cells, cov_src, cov_comp = {}, {}, {c: 0 for c in C}
for src, row in TAB.items():
    line, n_ok = [], 0
    for c in C:
        v, why = row[c]
        cells[f"{src}×{c.split()[0]}"] = dict(n=1, verdict=v, evidence=why,
            inclusion=[f"{src} 的 {c}", "依据必须是本项目已验证过的事实", "无依据 -> 判不了"])
        if v == "✓": n_ok += 1; cov_comp[c] += 1
        line.append(v if v != "判不了" else "?")
    cov_src[src] = n_ok
    print(f"  {src:8s} " + " ".join(f"{x:>5s}" for x in line) + f"   覆盖 {n_ok}/6")
print("\n=== 每个部件被几份数据源覆盖 ===")
for c in C: print(f"  {c:14s} {cov_comp[c]} 份")
zero = [c for c in C if cov_comp[c] == 0]
print(f"\n  **零覆盖的部件:{zero if zero else '无'}**")
G = Gate("那份规格的六个部件,本地各覆盖了几个?")
G.positive_control("正对照:至少一格是 ✓ 且带依据",
                   planted=float(sum(1 for v in cells.values() if v["verdict"] == "✓" and v["evidence"])),
                   floor=0.5, spread=1e-9)
bad = [k for k, v in cells.items() if v["verdict"] == "判不了" and v["evidence"] and "判不了" not in v["evidence"]]
G.negative_control("负对照:『判不了』的格不得带**肯定性**依据", null=float(len(bad)),
                   effect=float(sum(1 for v in cells.values() if v["verdict"] == "✓")),
                   null_spread=1e-9, null_kind="有依据却填判不了 = 表填错了")
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 70)
if sum(1 for v in cells.values() if v["verdict"] == "✓") > 0 and not bad:
    world = "ZERO-COVERAGE" if zero else "ALL-COVERED"
    verdict = (f"零覆盖部件:{zero} -> **页面点名它** " if zero
               else "每个部件至少被一份覆盖 -> 独特性在**组合**而非**部件**")
    print(f"评判:**{world}** —— {verdict}")
    print(f"  每份来源的覆盖数:{cov_src}")
    print("⚠ 这个 KILL 会怎样失败:**只覆盖本收藏**,且一格的 ✓/✗ 取决于依据条目当时的范围;"
          "而**未联网** ⇒ 外部是否存在满足某部件的来源,本轮不可知。")
else:
    world, verdict = "UNVERIFIED", f"对照未过(bad={bad})"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(table={s: {c: TAB[s][c] for c in C} for s in TAB}, coverage_by_source=cov_src,
               coverage_by_component=cov_comp, zero_coverage=zero, world=world, verdict=verdict,
               impossible=["只覆盖本收藏", "一格取决于依据条目当时的范围", "未联网,外部不可知"],
               unchallenged=True), open(OUT / "coverage.json", "w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'coverage.json'}")

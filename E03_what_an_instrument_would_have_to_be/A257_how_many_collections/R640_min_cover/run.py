"""E03·A257·R640 — 这七个缺口,最少需要几次采集?

`#595` 的 NEXT。行动类型:**PRODUCTION**(结构分析,**不提议采数据**)。
七条规格的「差在哪」现在散在一张表里,读起来像**七个独立的洞**。
**而「它们是不是七件事」是一个有可失败答案的问题。**

G1 ESTIMAND(先于方法):**最小采集设计数** = 把七条规格分成若干组,
   使**每一组内部两两兼容**,所需的**最少组数**(最小团覆盖)。
**兼容的判据(先于分组写死,机械应用):**
   两条规格**不兼容**,当且仅当它们要求
   ① **不同的观察单位**(人 vs 社会),或 ② **不同的原材料**(新访谈 vs 既有民族志)。
   ⚠ **时间跨度不算不兼容** —— 一条要等二十年才能满足的规格,
     与一条立刻能满足的规格,**仍在同一次采集里**;这一点会低估**时间成本**,必须单独报。
预注册:
   最小组数 **≤ 2** ⇒ 这七条其实是**两件事**,页面应当这样说;
   **≥ 4** ⇒ 它们是彼此独立的七个缺口,页面维持现状;
   **= 3** ⇒ 报三组各是什么。
CONTROLS:正对照 = 一对**明显不兼容**的规格(R2 民族志编码 × R7 个人面板)必须判不兼容 ·
   负对照 = 一条规格与**它自己**必须判兼容 · 全矩阵公布
IMPOSSIBLE:兼容判据由我定(两条,写在前面)· **不估成本、不估可行性**(`#559c` 的纪律)·
   **时间跨度不进判据**,所以「一次采集」低估了它的时间成本 · [unchallenged]
"""
import os, sys, pathlib, json, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
from lib.gates import Gate
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
# 七条规格,逐条标注它的**观察单位**与**原材料** —— 这两栏是判据的全部输入
SPECS = {
 "R1 抽样框 + 第二波": dict(unit="人", source="新访谈",
   gap="A sampling frame, and a second wave", horizon="立即起,第二波需数年"),
 "R2 第二组编码者": dict(unit="社会", source="既有民族志",
   gap="A second set of coders", horizon="立即"),
 "R3 下一波留住十道题": dict(unit="人", source="新访谈",
   gap="Keeping all ten in the next wave", horizon="下一波"),
 "R4 跨波条目稳定": dict(unit="人", source="新访谈",
   gap="Item stability across waves", horizon="跨两波"),
 "R5 第二条跨年序列": dict(unit="人", source="新访谈",
   gap="A second series like it", horizon="**二十年以上**"),
 "R6 态度题的双问法": dict(unit="人", source="新访谈",
   gap="The same double-asking, applied to the attitude items", horizon="立即"),
 "R7 面板维度": dict(unit="人", source="新访谈",
   gap="The entire panel dimension", horizon="跨两波"),
}
def compatible(a, b):
    A, B = SPECS[a], SPECS[b]
    if A["unit"] != B["unit"]: return False, f"观察单位不同({A['unit']} vs {B['unit']})"
    if A["source"] != B["source"]: return False, f"原材料不同({A['source']} vs {B['source']})"
    return True, "单位与原材料都相同"
K = list(SPECS)
print("=== 七条规格的两栏输入(判据的全部依据)===")
for k in K: print(f"  {k:16s} 单位={SPECS[k]['unit']:2s} 原材料={SPECS[k]['source']:6s} 时间跨度={SPECS[k]['horizon']}")
print("\n=== 兼容矩阵(全矩阵公布)===")
M, cells = {}, {}
print("      " + " ".join(f"{k.split()[0]:>4s}" for k in K))
for a in K:
    row = []
    for b in K:
        ok, why = compatible(a, b)
        M[(a, b)] = ok
        cells[f"{a.split()[0]}×{b.split()[0]}"] = dict(n=1, compatible=bool(ok), why=why,
            inclusion=["判据:单位相同且原材料相同", "机械应用,无人工例外"])
        row.append("✓" if ok else "✗")
    print(f"  {a.split()[0]:>4s}  " + " ".join(f"{x:>4s}" for x in row))
# 最小团覆盖:七个点,兼容关系是等价关系(单位×原材料),所以直接按 (unit,source) 分组
groups = {}
for k in K: groups.setdefault((SPECS[k]["unit"], SPECS[k]["source"]), []).append(k)
print(f"\n=== 分组(兼容关系由两个属性决定,故为等价关系 -> 组 = 等价类)===")
for (u, s), items in groups.items():
    print(f"  ({u} · {s}):{len(items)} 条 —— {[x.split()[0] for x in items]}")
n_groups = len(groups)
print(f"\n  **最小采集设计数 = {n_groups}**(预注册:≤2 -> 两件事 · =3 -> 报三组 · ≥4 -> 七个独立缺口)")
G = Gate("这七个缺口,最少需要几次采集?")
ok_pos, why_pos = compatible("R2 第二组编码者", "R7 面板维度")
G.positive_control("正对照:R2(民族志编码)× R7(个人面板)必须判不兼容",
                   planted=float(not ok_pos), floor=0.5, spread=1e-9)
selfok = all(compatible(k, k)[0] for k in K)
G.negative_control("负对照:每条规格与它自己必须兼容", null=float(not selfok), effect=1.0,
                   null_spread=1e-9, null_kind="同一条规格与自身,判据上必然兼容")
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 70)
if (not ok_pos) and selfok:
    if n_groups <= 2:
        world = "TWO-THINGS"; verdict = f"最小组数 {n_groups} ≤ 2 -> **这七条其实是两件事**"
    elif n_groups == 3:
        world = "THREE"; verdict = "三组,逐组报"
    else:
        world = "SEVEN"; verdict = f"最小组数 {n_groups} ≥ 4 -> **彼此独立的缺口**"
    print(f"评判:**{world}** —— {verdict}")
    longest = max(SPECS.values(), key=lambda v: ("二十年" in v["horizon"], v["horizon"]))
    print(f"⚠ 而「一次采集」低估了时间成本:同一组里最长的一条是 **{longest['horizon']}** "
          f"({[k for k in K if SPECS[k]['horizon']==longest['horizon']][0]}) —— "
          f"**判据里没有时间,所以这个 {n_groups} 是「几种采集」,不是「多久」。**")
else:
    world, verdict = "UNVERIFIED", "对照未过"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(specs=SPECS, matrix={f"{a}|{b}": M[(a, b)] for a, b in itertools.product(K, K)},
               groups={f"{u}·{s}": v for (u, s), v in groups.items()}, n_groups=n_groups,
               world=world, verdict=verdict,
               impossible=["兼容判据由我定(两条,写在前面)", "不估成本不估可行性",
                           "时间跨度不进判据 -> 低估时间成本"],
               unchallenged=True), open(OUT / "min_cover.json", "w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'min_cover.json'}")

"""#854 · E03·A88·R293 —— 我一轮前写下的那堵墙:把盘上每一具仪器都点一遍名

**Frontier。** `#853` 在「结构性做不到」里写下:
**「九十年代复现不了 —— NSFG 的字典只存在于 2011–2013 与 2017–2019。」**
⇒ 而 `#844` 刚教过:**三堵墙里倒的那一堵,恰好是我最有把握的那一堵。**
`#843` 更狠 —— 我登记的「结构性改不了」,**下一轮用一条 `git show` 就推翻了。**
**这堵墙才一轮大,而且它挡着本项目的头条。⇒ 现在查它,不是以后。**

**⚠⚠ 而更要紧的是:`#853` 给的理由太窄。** 它只说了 NSFG 的字典 ——
**而「九十年代能不能跨仪器复现」这个问题的总体,是盘上的每一具仪器,不是我恰好看过的那一具。**
`#849`① 点了 NSFG 的名,于是我只查了 NSFG。**这正是 `feedback_check_only_as_good_as_its_population`:
一个由我随手选定的总体,让一次客观检查变回自我报告。**

`G1` **估计量(先于方法命名)**:
   **盘上每一具仪器,在三根「承载该主张所必需」的轴上各自的通过情况** ——
   ① **态度轴**:有没有一道问「对同性性关系持什么态度」的题(**不是**问受访者自己的身份或行为);
   ② **宗教轴**:有没有任何宗教度测量;
   ③ **时间轴**:数据是否覆盖 1990–1999。
   **主张要能被第二具仪器检验,需要某一具仪器三轴全过。**

**⚠ 这是一次「搜索即仪器」的测量,所以它必须有正控**
(`realstat` 点名的那条:**一个从没返回过非零的搜索给的零是沉默,不是无罪**):
   **正控**:搜索必须在**已知三轴俱全**的地方开火 —— GSS(`homosex` + `attend`/`reliten`/`fund` + 1972–2024)
   与 **NSFG-2017**(`samesex` + `attndnow`/`reldlife`,`#853` 刚亲手验过)。
   **若搜索在这两处看不见,它在别处给的零一律作废。**
   **负控**:一根**不存在的轴**(一个胡编的词根)必须在每一具仪器上都返回 0 ——
   ⚠ **「这个零该不该是零?」该**:盘上没有任何数据集会包含一个我刚编出来的词,**按定义**。

三个世界:
   A **有某具仪器三轴全过** ⇒ **墙倒了** ⇒ 九十年代可以跨仪器复现,`#853` 那句要撤。
   B **每具仪器都至少缺一轴** ⇒ **墙立着,而理由比 `#853` 写的宽也更精确** ——
     不是「NSFG 没字典」,是**盘上没有任何仪器同时具备这三轴**。
   C **搜索在正控上就看不见** ⇒ **UNVERIFIED**,本轮什么都没测到(而这不是墙立着)。

预测矩阵:
   | 世界 | 现在 | 有仪器三轴全过 | 每具都缺轴 | 正控失败 |
   | A 墙倒     | 0.25 | **0.85** | 0.05 | 0.10 |
   | B 墙立(更宽) | 0.60 | 0.05 | **0.85** | 0.10 |
   | C 搜索不合用 | 0.15 | 0.10 | 0.10 | **0.80** |

预注册判词(条件式):
  if 正控在 GSS 与 NSFG-2017 上都开火 and 负控在每具仪器上都为 0:
      任一仪器三轴全过 -> A(撤 `#853` 那句)
      每具都缺 ≥1 轴   -> B(墙立着,并把理由改写成「所有仪器」)
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**词表决定结论。**
  「态度轴」的词表若写宽(含 `gay`/`lesbian`),会把 **YRBS 的 `sexid`(受访者自己的性身份)**
  当成态度题命中 —— **那是身份,不是态度,构念完全不同。**
  ⇒ 控制:**态度轴与身份/行为轴分开数、分开报**,并把命中的**变量名与标签原样印出来**供复核;
  **判据只认态度轴,而身份轴的命中数一起发表**,让读者看见我没有把它算进去。
⚠ 本轮换不了仪器 —— **对象就是「盘上有哪些仪器」**,这是唯一的总体。
⚠ 标注 **Frontier**:它有一个能倒的墙,而墙倒会撤掉一句已发表的话。
"""
import json, pathlib, re, sys, zipfile
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
EXT = ROOT / "data/external"

# ── 三根轴的词表:态度 / 身份行为 / 宗教;外加一根不存在的轴做负控 ────────────
ATT = re.compile(r"(all right|wrong|approve|disapprove|acceptab|attitude|opinion|should be)"
                 r"[^\n]{0,80}(same.?sex|homosex|gay|lesbian)"
                 r"|(same.?sex|homosex|gay|lesbian)[^\n]{0,80}"
                 r"(all right|is wrong|approve|disapprove|acceptab|attitude|opinion|should)", re.I)
IDB = re.compile(r"(sexual identity|sexid|ever had sexual contact|partners.{0,20}same.?sex|"
                 r"orientation)", re.I)
REL = re.compile(r"(religio|church|attend[^\n]{0,20}service|how often[^\n]{0,20}pray|"
                 r"fundamentalis|denominat)", re.I)
FAKE = re.compile(r"zqxwjv[a-z]*", re.I)          # 负控:一根不存在的轴

def text_of(paths, cap=4_000_000):
    """只读**字典/码本/文本**,不读 .dat 主体(成本表:.dat 最大 141 MB)。"""
    buf = []
    for p in paths:
        try:
            if p.suffix.lower() == ".zip":
                with zipfile.ZipFile(p) as z:
                    for nm in z.namelist()[:40]:
                        if nm.lower().endswith((".sas", ".dct", ".do", ".txt", ".csv", ".sps")):
                            buf.append(z.read(nm)[:cap].decode("utf-8", "replace"))
            elif p.suffix.lower() in (".sas", ".dct", ".do", ".txt", ".csv", ".sps", ".json", ".md"):
                buf.append(p.read_text(errors="replace")[:cap])
        except Exception:
            pass
        if sum(len(b) for b in buf) > cap: break
    return "\n".join(buf)

# ⚠⚠ **第一版的年份轴是从字典散文里刮 4 位数字 —— 而那不是覆盖,是标签里出现过的年份。**
#    实测:NSFG 因此被判为「九十年代 有[1990,1991,1992]」,而那是**出生年份区间的标签**,
#    NSFG 的字典只到 2011。**一个把「提到过 1990」当成「覆盖 1990」的代理,方向是系统性多报。**
#    ⇒ 改成:**覆盖 = 可读的覆盖** —— 由「存在数据文件」**且**「存在能解它的字典」共同决定,
#      逐仪器显式列出,不从散文刮。
READABLE = {
 "GSS":       [(1972, 2024)],                 # 单个 .dta,标签在二进制里(见下)
 "NSFG":      [(2011, 2013), (2017, 2019)],   # ⚠ 只有这两轮有 .dct;1988/1995 的 .dat 无字典
 "YRBS":      [(1991, 2023)],                 # SADC 合并档 + SAS input program
 "BRFSS":     [(2023, 2023)],                 # 盘上只有 LLCP2023
 "SCCS":      [],                             # 无时间维:一次性民族志编码
 "OpenPsych": [],                             # 无时间序列
 "Dataverse": [],                             # 单次研究复制包
}
def has_nineties(name):
    return any(a <= 1999 and b >= 1990 for a, b in READABLE[name])

# 每具仪器:名字 -> (要读的字典/码本文件, 已知年份补充)
INSTR = {
 "GSS":   (sorted((EXT/"gss").rglob("*.dct")) + sorted((EXT/"gss").rglob("*.do"))[:3], range(1972, 2025)),
 "NSFG":  (sorted((EXT/"nsfg/setup").glob("*")), (2012, 2018)),
 "YRBS":  (sorted((EXT/"yrbs").glob("*.sas")), ()),
 "BRFSS": (sorted((EXT/"brfss").glob("*")), ()),
 "SCCS":  (sorted((EXT/"dplace").rglob("variables.csv"))[:1], ()),
 "OpenPsych": (sorted((EXT/"openpsych").rglob("*.txt"))[:6] + sorted((EXT/"openpsych").rglob("*.csv"))[:4], ()),
 "Dataverse": (sorted((EXT/"dataverse").rglob("*.do"))[:6] + sorted((EXT/"dataverse").rglob("*.R"))[:4], ()),
}
print("=== ⓪ 成本表 + 每具仪器读了哪些字典(只读字典,不读 .dat 主体)===")
ROWS = {}
for name, (paths, extra) in INSTR.items():
    txt = text_of(paths)
    # ⚠ GSS 的字典在 `.dta` 二进制里,没有独立的 .dct/.do ⇒ 第一版读到 0 字符,**正控当场失败**。
    #   那正是这条控制存在的理由:我差点在「搜索看不见 GSS」的情况下报出「只有 NSFG 三轴全过」,
    #   而 GSS 正是本项目每一个数字的来源。⇒ 从 `.dta` 里读变量标签。
    if name == "GSS":
        import pandas as _pd
        _r = _pd.io.stata.StataReader(EXT/"gss/GSS_stata/gss7224_r3a.dta")
        txt = "\n".join(f"{k} {v}" for k, v in _r.variable_labels().items())
    nineties = [a for a, b in READABLE[name] if a <= 1999 and b >= 1990]
    att = ATT.findall(txt); idb = IDB.findall(txt); rel = REL.findall(txt); fk = FAKE.findall(txt)
    ROWS[name] = dict(n_files=len(paths), chars=len(txt),
                      attitude=len(att), identity_behaviour=len(idb), religion=len(rel),
                      fake=len(fk), readable_spans=READABLE[name], has_nineties=has_nineties(name))
    print(f"  {name:10s} 文件 {len(paths):>3} · 字典字符 {len(txt):>9,} · "
          f"**态度 {len(att):>4}** · 身份/行为 {len(idb):>4} · **宗教 {len(rel):>5}** · "
          f"**可读年段** {READABLE[name] or '无时间维'} ⇒ 九十年代 "
          f"{'**有**' if has_nineties(name) else '**无**'}")
print("  ⚠ **态度轴与身份/行为轴分开数** —— YRBS 的 `sexid` 是**受访者自己的性身份**,"
      "**不是对同性性关系的态度**,构念完全不同,不许算进态度轴")

print("\n=== ① 三轴同时满足的仪器 ===")
# ⚠⚠⚠ **第一版的 kill 写成「任一仪器三轴全过」,于是 GSS 自己过了,判词印出「墙倒了」——
#    而 GSS 正是那条主张已经住在里面的仪器。问的是**第二具**仪器能不能承载它,
#    而我把在位者算进了 kill 的总体。**判据的尺是对的,判据的总体是错的。**
#    `#836`① 让我命名尺,而这一次错的不是尺 —— 是我从没问过「这条判据该在谁身上评」。
INCUMBENT = "GSS"
def passes(r): return r["attitude"] > 0 and r["religion"] > 0 and r["has_nineties"]
ok_all = [n for n, r in ROWS.items() if passes(r)]
ok = [n for n in ok_all if n != INCUMBENT]          # ⇐ 排除在位者,这才是「第二具仪器」
print(f"  ⚠ 在位者 `{INCUMBENT}` 三轴全过 = {INCUMBENT in ok_all} —— **这是同义反复,不计入 kill**;"
      f"kill 的总体是**除它以外**的 {len(ROWS)-1} 具")
for n, r in ROWS.items():
    miss = [a for a, v in (("态度", r["attitude"] > 0), ("宗教", r["religion"] > 0),
                           ("九十年代", r["has_nineties"])) if not v]
    print(f"  {n:10s} {'**三轴全过**' if not miss else '缺:' + '·'.join(miss)}")
print(f"  ⇒ **三轴全过的仪器:{ok or '无'}**")

print("\n=== ② 控制 ===")
pc_gss = ROWS["GSS"]["attitude"] > 0 and ROWS["GSS"]["religion"] > 0
pc_nsfg = ROWS["NSFG"]["attitude"] > 0 and ROWS["NSFG"]["religion"] > 0
print(f"  正控:搜索必须在**已知三轴俱全**处开火 ⇒ GSS 态度 {ROWS['GSS']['attitude']}/宗教 "
      f"{ROWS['GSS']['religion']} ⇒ **{pc_gss}** · NSFG 态度 {ROWS['NSFG']['attitude']}/宗教 "
      f"{ROWS['NSFG']['religion']} ⇒ **{pc_nsfg}**")
print(f"     ⚠ **若搜索在这两处看不见,它在别处给的零一律作废**(`realstat`:一个从没返回过非零的"
      f"搜索给的零是沉默,不是无罪)")
nc = all(r["fake"] == 0 for r in ROWS.values())
print(f"  负控:一根**不存在的轴**(胡编词根)必须在每具仪器上都返回 0 ⇒ "
      f"{ {n: r['fake'] for n, r in ROWS.items()} } ⇒ **{nc}**")
print(f"     ⚠ **「这个零该不该是零?」该** —— 盘上没有任何数据集会包含一个我刚编出来的词,**按定义**")

G = Gate("#854 · 把盘上每一具仪器都点一遍名")
G.asserted("① 前提(跑前写下的最强混淆):**词表决定结论** —— 「态度」词表若写宽会把 YRBS 的 "
           "`sexid`(**受访者自己的性身份**)当成态度题命中,**那是身份不是态度,构念完全不同** ⇒ "
           "**态度轴与身份/行为轴分开数、分开报**,判据只认态度轴",
           bool(all("identity_behaviour" in r for r in ROWS.values())),
           " · ".join(f"{n}:态度{r['attitude']}/身份{r['identity_behaviour']}" for n, r in ROWS.items()),
           kind="control")
G.asserted("② 正控:搜索必须在**已知三轴俱全**的 GSS 与 NSFG-2017 上开火 —— "
           "**否则它在别处给的零一律作废**",
           bool(pc_gss and pc_nsfg), f"GSS {pc_gss} · NSFG {pc_nsfg}", kind="control")
G.asserted("③ 负控:一根不存在的轴必须在每具仪器上返回 0"
           "(⚠ **这个零该是零**:盘上不会有我刚编出来的词,按定义)",
           nc, str({n: r["fake"] for n, r in ROWS.items()}), kind="control")
G.asserted("④ 前提(总体):**总体是盘上的每一具仪器,不是我恰好看过的那一具** —— "
           "`#853` 只查了 NSFG,因为 `#849`① 点了它的名"
           "(`feedback_check_only_as_good_as_its_population`)",
           bool(len(ROWS) >= 6), f"点名的仪器 {len(ROWS)} 具:{sorted(ROWS)}", kind="control")
G.asserted("⑤ kill(预注册):「九十年代那条缝可以跨仪器复现」要成立,需**至少一具「非 GSS」的仪器三轴全过** "
           "—— ⚠ **在位者 GSS 自己过关是同义反复,不计入**(第一版忘了排除它,判词因此印出「墙倒了」)",
           bool(ok), f"非 GSS 三轴全过 {ok or '无'}(在位者 GSS 全过 {INCUMBENT in ok_all},不计)",
           kind="kill",
           yardstick="每具仪器在三根必需轴上的通过情况(态度 / 宗教 / 覆盖 1990–1999)",
           yardstick_noise=0.0)
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif ok:
    VERD = (f"**A 墙倒了 —— 非 GSS 的 {ok} 三轴全过 ⇒ `#853` 那句「九十年代复现不了」要撤。**")
else:
    misses = {n: [a for a, v in (("态度", r["attitude"] > 0), ("宗教", r["religion"] > 0),
                                 ("九十年代", r["has_nineties"])) if not v]
              for n, r in ROWS.items() if n != INCUMBENT}
    VERD = (f"**B 墙立着,而理由比 `#853` 写的宽,也更精确。**\n"
            f"  `#853` 说的是「NSFG 的字典只到 2011」——**那只是一具仪器的一个理由**;\n"
            f"  点完名之后,真正的理由是:**盘上没有任何一具仪器同时具备这三轴。**\n"
            + "\n".join(f"  · **{n}** 缺:{'·'.join(m)}" for n, m in misses.items() if m) + "\n"
            f"  ⇒ **一句关于人的话,而这一轮它是关于「我们目前只能从一个地方知道这件事」的:\n"
            f"  九十年代那条缝到底是不是真的张开了,现在只有一份问卷说过。\n"
            f"  不是因为没人去查别的问卷 —— 是因为要查它,一份问卷得同时问「你怎么看同性关系」、\n"
            f"  「你信不信教」,而且得在九十年代就开始问。盘上这几份,没有一份三样都占。**")
print(VERD)
print(f"\n⚠ **本轮量的是「字典里有没有这样的题」,不是「题问得好不好」** —— "
      f"构念效度需要外部金标准,本站没有(`realstat` §2)。")
print(f"⚠ 且 **YRBS 的身份/行为命中 {ROWS['YRBS']['identity_behaviour']} 处被如实印出、"
      f"但不计入态度轴** —— 让读者看见我没有把它算进去,而不是只看见我的结论。")
json.dump(dict(instruments=ROWS, three_axis_pass_excluding_incumbent=ok,
               three_axis_pass_all=ok_all, incumbent=INCUMBENT,
               pos_control=dict(gss=pc_gss, nsfg=pc_nsfg), neg_control=nc,
               scope="dictionaries/codebooks only; construct validity needs an external gold standard",
               admissible=adm, verdict=VERD, gate_ok=G.verdict()),
          open(OUT / "enumerate_instruments.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'enumerate_instruments.json'}")

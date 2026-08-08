"""#843 · E03·A81·R282 —— `#798` 的 `naked==0` 重算:而「当时错没错」和「今天对不对」是两个问题

**还 `#842`①。** `#842` 登记过一条「改不了的」:
**那三轮的判词要重算就得整个重跑,而它们依赖当时的页面状态,今天的页面已多了四十多轮
⇒ 重跑给的是「今天的答案」,不是「当时是否算错了」。**

**⚠⚠ 而那句「改不了」是错的,本轮当场推翻它 —— 因为 `git` 里存着当时的页面。**
`git log -S"Entry 798 ·" -- RETRACTIONS.md` 定位到引入 `#798` 的那次提交,
`git show <c>:README.md` 就是**当时那一版页面**。
⇒ **两个问题可以分开回答,而且必须分开:**
   **(a) 当时错没错** —— 用**当时的页面**,比较旧切法(它当时用的)与新切法给出的判词;
   **(b) 今天对不对** —— 用**今天的页面**,新切法下判词是什么。
⚠ **`#842` 把这两个问题合并成「做不到」,而做不到的只是把它们混在一起做。**
**一条「结构性做不到」的登记,本身就是一个该被攻击的对象**
(`realstat` 的 *a wall never checked*:三条「永久限制」里有一条只是没跑过的查询)。

G1 估计量:`#798` 的判词 `naked == 0` —— **使用侧(owner < 797)且提到三堆结构、
却没有撤回指针的单位数** —— 在 **{当时页面, 今天页面} × {旧切法, 新切法}** 四格上各是多少。

四个世界:
   A **当时对、今天也对** ⇒ 切法缺陷没碰到这条判词 ⇒ `#798` 完好。
   B **当时对、今天错** ⇒ 判词没错,但**页面在这四十多轮里退化了** ⇒ 要修的是页面不是账。
   C **当时就错** ⇒ **`#798` 的「撤回落地」是切法造出来的** ⇒ 必须撤 `#798` 的判词。
   D **旧新切法在当时页面上给同一个答案** ⇒ 这条判词对切法不敏感 ⇒ `#833` 的影响面比想的小。

预注册判词(条件式):
  if 正控开火(**在当时页面里植入一条裸着的使用侧单位,两种切法都必须数到它**)
     and 负控开火(**给那条植入的单位补上撤回指针,两种切法都必须不再数它**):
      当时(旧) != 当时(新)              -> C(**撤 `#798` 的判词**)
      当时(旧) == 当时(新) 且 今天(新)>0 -> B(**修页面**)
      四格全 0                            -> A
      当时两格相同而与今天不同            -> D 并入 A/B 判定
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**`#798` 的谓词里 `PTR` 含 `#797`,而「当时的页面」正是刚写完 `#797`
  的那一版** —— 若我从今天的文件读谓词却套到当时的页面上,**谓词与页面不同代**。
  ⇒ 控制:**谓词也从当时那次提交的脚本里读**(`git show <c>:<script>`),
  **让谓词与页面同代**;并把两代谓词是否相同印出来。
⚠ 本轮换不了仪器(对象是页面自己)。⚠ 标注 **Production + 推翻 `#842` 的一条「改不了」**。
"""
import ast, json, pathlib, re, subprocess, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.page_units import units
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
SCRIPT = ("E03_what_an_instrument_would_have_to_be/A50_撤回要走到它引发的物件上/"
          "R237_how_many_page_places_still_use_a_dead_structure/verify_the_retraction_landed.py")

def git(*a):
    return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout

COMMIT = git("log", "--format=%H", "-S", "Entry 798 ·", "--", "RETRACTIONS.md").strip().split("\n")[-1]
print(f"=== ⓪ `#798` 是哪一次提交引入的(从对象读,不靠记忆)===")
print(f"  commit **{COMMIT[:12]}** · {git('log','-1','--format=%ad %s','--date=short',COMMIT).strip()[:90]}")

def patterns(src):
    out = {}
    for n in ast.parse(src).body:
        if (isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)
                and isinstance(n.value, ast.Call) and getattr(n.value.func, "attr", None) == "compile"):
            try: out[n.targets[0].id] = eval(compile(ast.Expression(n.value), "<p>", "eval"), {"re": re})
            except Exception: pass
    return out

P_now = patterns((ROOT/SCRIPT).read_text(encoding="utf-8"))
P_then = patterns(git("show", f"{COMMIT}:{SCRIPT}"))
same = {k: (P_now[k].pattern == P_then.get(k, re.compile("")).pattern) for k in ("CLUMP", "PTR")}
print(f"\n=== ① 跑前混淆的控制:谓词与页面必须同代 ===")
for k in ("CLUMP", "PTR"):
    print(f"  {k}: 今天版与当时版{'**相同**' if same[k] else '**不同** ⚠'}")
print(f"  ⇒ 下面**一律用当时那一版谓词配当时的页面,今天的谓词配今天的页面**。")

PAGES = {"当时": {f: git("show", f"{COMMIT}:{f}") for f in ("README_zh.md", "README.md")},
         "今天": {f: (ROOT/f).read_text(encoding="utf-8") for f in ("README_zh.md", "README.md")}}
PATS = {"当时": P_then, "今天": P_now}

def naked(page_txt, pats, fixed):
    """`#798` 的判词:使用侧(owner<797)提到三堆结构、却没有撤回指针的单位数。"""
    n, rows = 0, []
    for k, owner, body in units(page_txt, fixed=fixed):
        if not pats["CLUMP"].search(body): continue
        if owner >= 797: continue                       # 撤回侧
        if pats["PTR"].search(body): continue           # 带撤回指针
        n += 1; rows.append((k, owner))
    return n, rows

print(f"\n=== ② `#798` 的判词 `naked == 0` 在四格上各是多少 ===")
G4 = {}
for era in ("当时", "今天"):
    for cut, fx in (("旧切法", False), ("新切法", True)):
        tot, det = 0, []
        for f in ("README_zh.md", "README.md"):
            n, rows = naked(PAGES[era][f], PATS[era], fx)
            tot += n; det += [(f, *r) for r in rows]
        G4[(era, cut)] = dict(naked=tot, rows=det[:8])
        print(f"  {era} × {cut}:**naked = {tot}**" + (f"  ⇒ {det[:4]}" if tot else ""))
then_old, then_new = G4[("当时", "旧切法")]["naked"], G4[("当时", "新切法")]["naked"]
now_new = G4[("今天", "新切法")]["naked"]

print(f"\n=== ③ 控制 ===")
INJ = ("\n本段是植入的对照:三堆结构在这里被当作事实断言。 `[#700「植入」]`\n")
inj_page = PAGES["当时"]["README_zh.md"] + INJ
pc = {cut: naked(inj_page, PATS["当时"], fx)[0] - naked(PAGES["当时"]["README_zh.md"], PATS["当时"], fx)[0]
      for cut, fx in (("旧", False), ("新", True))}
print(f"  正控:在当时页面末尾植入**一条裸着的使用侧单位**(owner=700<797,提三堆,无指针)"
      f"⇒ 旧切法 +{pc['旧']} · 新切法 +{pc['新']} —— 两种切法都必须 +1")
INJ2 = ("\n本段是植入的对照:三堆结构在这里被当作事实断言,而它已在 `#797` 撤回。 `[#701「植入带指针」]`\n")
nc = {cut: naked(PAGES["当时"]["README_zh.md"]+INJ2, PATS["当时"], fx)[0]
             - naked(PAGES["当时"]["README_zh.md"], PATS["当时"], fx)[0]
      for cut, fx in (("旧", False), ("新", True))}
print(f"  负控:同一条单位**补上撤回指针** ⇒ 旧切法 +{nc['旧']} · 新切法 +{nc['新']} —— "
      f"两种切法都必须 **+0**(⚠ **「这个零该不该是零?」该**:带了指针就不算裸着,按定义)")

G = Gate("#843 · #798 的 naked==0 重算:当时错没错 vs 今天对不对")
G.asserted("① 推翻 `#842` 登记的那条「改不了」:**当时的页面在 `git` 里** —— "
           "`git log -S\"Entry 798 ·\"` 定位提交,`git show <c>:README.md` 取回当时那一版 ⇒ "
           "**两个问题可以分开回答,而做不到的只是把它们混在一起做**",
           bool(len(COMMIT) == 40 and all(len(PAGES["当时"][f]) > 1000 for f in PAGES["当时"])),
           f"commit {COMMIT[:12]} · 当时中文页 {len(PAGES['当时']['README_zh.md']):,} 字符",
           kind="control")
G.asserted("② 前提(跑前写下的最强混淆):`PTR` 含 `#797`,而当时的页面正是刚写完 `#797` 那一版 ⇒ "
           "**谓词必须与页面同代** ⇒ 谓词也从那次提交的脚本里读",
           bool(len(P_then) >= 2), f"当时版谓词 {sorted(P_then)} · 与今天版是否相同 {same}",
           kind="control")
G.asserted("③ 正控:在当时页面植入**一条裸着的使用侧单位**,两种切法都必须数到它(各 +1)",
           bool(pc["旧"] == 1 and pc["新"] == 1), f"旧 +{pc['旧']} · 新 +{pc['新']}", kind="control")
G.asserted("④ 负控:同一条单位补上撤回指针后,两种切法都必须**不再**数它(各 +0)"
           "(⚠ **这个零该是零**:带指针就不算裸着,按定义)",
           bool(nc["旧"] == 0 and nc["新"] == 0), f"旧 +{nc['旧']} · 新 +{nc['新']}", kind="control")
G.asserted("⑤ kill(预注册):「`#798` 当时没算错」要成立,需**当时页面上旧切法与新切法给出同一个 naked**",
           bool(then_old == then_new), f"当时:旧 {then_old} vs 新 {then_new}", kind="kill",
           yardstick="同一页面同一谓词下,两种切法各自的 naked 计数", yardstick_noise=0.0)
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif then_old != then_new:
    V = (f"**C `#798` 当时就算错了 ⇒ 它的判词要撤。** 当时页面上,旧切法 naked={then_old}、"
         f"新切法 naked={then_new} —— **「撤回落地」这个结论是切法造出来的。**")
elif now_new > 0:
    V = (f"**B 当时对,今天错 ⇒ 要修的是页面,不是账。** 当时两种切法都给 naked={then_old};"
         f"而**今天的页面在新切法下 naked={now_new}** ⇒ 这四十多轮里,页面上又长出了"
         f"{now_new} 处在断言已撤结构而不带指针的地方。\n"
         f"  ⇒ **一句关于方法的话:一条「已落地」的更正不是一次性的事件,是一个会退化的状态 ——\n"
         f"  它在写下的那天为真,而没有任何机制让它明天还为真。**")
else:
    V = (f"**A 当时对、今天也对。** 四格里 naked 全为 {then_old}/{then_new}/{now_new} ——\n"
         f"  **切法缺陷没有碰到这条判词**,而这一点现在是量出来的,不是猜的。\n"
         f"  ⇒ **一句关于方法的话:`#842` 说这两个问题「结构性分不开」,而它们只是没被分开做过。\n"
         f"  一条写在账本里的「改不了」,和一条写在论文里的「不可能」,是同一种东西 ——\n"
         f"  它让停下来显得有理由,所以从来没人回去查它。**")
print(V)
json.dump(dict(commit=COMMIT, grid={f"{e}|{c}": v for (e, c), v in G4.items()},
               patterns_same_generation=same, pos_control=pc, neg_control=nc,
               admissible=adm, verdict=V, gate_ok=G.verdict(),
               action="Production + 推翻 #842 登记的一条「改不了」"),
          open(OUT/"then_versus_now.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'then_versus_now.json'}")

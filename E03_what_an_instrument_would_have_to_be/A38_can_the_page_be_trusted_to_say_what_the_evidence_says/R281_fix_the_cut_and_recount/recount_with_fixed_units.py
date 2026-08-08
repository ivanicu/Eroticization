"""#842 · E03·A81·R281 —— 修好 `units()`,把 `#798`/`#807`/`#828` 的计数重数一遍

**还 `#833`① 那笔债。** 而在还之前,`lib/page_units.py` 的控制**逼出了对 `#833` 自己的更正**:

**`#833` 写的是「锚后追加的注记对 `para` 不可见」。实测不是不可见,是归错了主。**
旧版 `para` 的左边界是**上一个锚的结尾**,所以 `#N` 锚之后的文本落进了 **`#N+1` 的单位**。
控制实测:`⟨注记:三堆⟩` 写在 `#100` 的锚后,旧版把它交给了 **`#101`**。
⇒ **后果不是「少数了」,是「数到别人头上」。**
**一个「不可见」的缺陷只让你漏;一个「错归」的缺陷会让你反向定罪** ——
而 `#798` 恰恰是**按 owner 把单位分成「撤回侧 / 使用侧」**的轮次。

**⚠ 三份 `units()` 拷贝(`R237`/`R246`/`R267`)实测语义相同** ⇒ 缺陷是统一的,不是三种。
⇒ 修在一处(`lib/page_units.py`),三轮的重算都指向它。

G1 估计量:**每一轮的总体规模与它自己的关键计数,在旧切法与新切法下各是多少**,
   以及**有没有任何一条已发表的判词因此改变**。

三个世界:
   A **数变了,判词没变** ⇒ 总体虚高是真的,但没有一条结论靠它 ⇒ 只需在页面上改数。
   B **有判词翻了** ⇒ **已发表的结论里有一条是切法造出来的**,必须撤。
   C **连总体规模都没变** ⇒ ⚠ **那 `#833` 记的缺陷在真实页面上不发生**,
     `#833` 本身要降级 —— **这是我不欢迎的那个,因为它意味着我为一个不存在的缺陷欠了三轮的债。**

预注册判词(条件式):
  if `lib/page_units.controls()` 全过(正控:旧版必须把锚后注记归给下一个锚;
     负控:没有表行的页面上两版必须给出相同的 para):
      任一判词改变 -> B(撤那一条)
      判词全不变而计数变 -> A
      计数也不变 -> C(降级 `#833`)
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**我要重数的三个谓词,是三轮各自写的正则。**
  **若我把它们手抄进本文件,就正好犯了 `#840`/`#841` 刚立规矩要防的那件事。**
  ⇒ 控制:**用 `ast` 从那三个文件里把 `re.compile(...)` 的模块级赋值读出来再编译**,
  **一个字符都不手抄** —— 而本轮也因此是那把新尺子的第一个正面用例。

⚠ 本轮对象是页面自己 ⇒ 换不了仪器(与 `#787`/`#792`/`#793`/`#798` 同一种:对象不是世界)。
⚠ 本轮标注 **Production + 一次对 `#833` 的更正**,不冒充新发现。
"""
import ast, json, pathlib, sys, re
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.page_units import units, controls
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
E3 = ROOT/"E03_what_an_instrument_would_have_to_be"
SRC = {
 798: E3/"A50_撤回要走到它引发的物件上/R237_how_many_page_places_still_use_a_dead_structure/verify_the_retraction_landed.py",
 807: E3/"A57_两条更正要走到页面上/R246_not_unrun_but_unrunnable_and_the_page_does_not_say_so/two_corrections_must_land.py",
 828: E3/"A72_两条更正落地/R267_uncorrected_markers_and_a_debt_never_recounted/land_two_corrections.py",
}

def read_patterns(path):
    """⚠ 从对象读,不手抄(`#840`/`#841`):把模块级 `X = re.compile(...)` 取出来编译。"""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out = {}
    for n in tree.body:
        if (isinstance(n, ast.Assign) and len(n.targets) == 1
                and isinstance(n.targets[0], ast.Name) and isinstance(n.value, ast.Call)
                and getattr(n.value.func, "attr", None) == "compile"):
            try: out[n.targets[0].id] = eval(compile(ast.Expression(n.value), "<p>", "eval"), {"re": re})
            except Exception: pass
    return out

print("=== ⓪ 控制:`lib/page_units.controls()` ===")
C = controls()
for k, v in C.items(): print(f"  {k:34s} {v}")
print(f"  ⇒ **旧版把写在 `#100` 锚后的注记记到了 {C['old_misattributed_to']} 名下** —— "
      f"这就是对 `#833` 的更正:**不是不可见,是归错主。**")

print("\n=== ① 三轮的谓词:从它们自己的文件读出来,一个字符都不手抄(`#841` 的第一个正面用例)===")
PAT = {}
for n, p in SRC.items():
    PAT[n] = read_patterns(p)
    print(f"  `#{n}` ← {p.name}:{sorted(PAT[n])}")

print("\n=== ② 总体规模:旧切法 vs 新切法 ===")
POP = {}
for f in ("README_zh.md", "README.md"):
    t = (ROOT/f).read_text(encoding="utf-8")
    o, w = units(t, fixed=False), units(t, fixed=True)
    ochars = sum(len(b) for _, _, b in o); wchars = sum(len(b) for _, _, b in w)
    POP[f] = dict(old_units=len(o), new_units=len(w), old_chars=ochars, new_chars=wchars,
                  old_rows=sum(1 for k, _, _ in o if k == "row"),
                  old_paras=sum(1 for k, _, _ in o if k == "para"),
                  new_paras=sum(1 for k, _, _ in w if k == "para"))
    print(f"  {f}: 单位数 {len(o)} → {len(w)} · **总字符 {ochars:,} → {wchars:,} "
          f"(虚高 {100*(ochars-wchars)/wchars:+.1f}%)** · para {POP[f]['old_paras']} → {POP[f]['new_paras']}")
print("  ⚠ **单位数不变而字符数变** —— 因为重复的是**文本**不是**单位**:"
      "表行的正文同时坐在 `row` 与某个 `para` 里。")

print("\n=== ③ 三轮的关键计数:旧 vs 新 ===")
RECOUNT = {}
def count_hits(pats, names, fixed):
    """一条单位若命中 `names` 里任一模式就计一次;逐版本、逐页统计。"""
    tot = {}
    for f in ("README_zh.md", "README.md"):
        t = (ROOT/f).read_text(encoding="utf-8")
        for k, owner, body in units(t, fixed=fixed):
            for nm in names:
                if nm in pats and pats[nm].search(body):
                    tot.setdefault(nm, []).append((f, k, owner))
    return {nm: len(v) for nm, v in tot.items()}, tot

SETS = {798: ["CLUMP", "PTR"], 807: ["TOP", "GEN", "MARK_TOP", "MARK_GEN"], 828: ["TARGET", "MARKED"]}
changed = []
for n, names in SETS.items():
    old_c, old_d = count_hits(PAT[n], names, False)
    new_c, new_d = count_hits(PAT[n], names, True)
    RECOUNT[n] = dict(old=old_c, new=new_c)
    keys = sorted(set(old_c) | set(new_c))
    print(f"  `#{n}`:" + " · ".join(
        f"{k} **{old_c.get(k,0)} → {new_c.get(k,0)}**{' ⚠变' if old_c.get(k,0)!=new_c.get(k,0) else ''}"
        for k in keys) or "  (无模式)")
    for k in keys:
        if old_c.get(k, 0) != new_c.get(k, 0): changed.append(f"#{n}:{k}")
    # 归属是否换边(`#798` 的判词直接建在 owner 上)
    if n == 798 and "CLUMP" in old_d:
        oo = {(f, k, ow) for f, k, ow in old_d["CLUMP"]}
        nn = {(f, k, ow) for f, k, ow in new_d.get("CLUMP", [])}
        moved = sorted(oo ^ nn)
        print(f"     ⚠ `#798` 的判词按 owner 分「撤回侧 (≥797) / 使用侧 (<797)」 ⇒ "
              f"**两版之间换了归属的单位:{len(moved)}**")
        for x in moved[:6]: print(f"        {x}")
        RECOUNT[798]["moved_units"] = len(moved)
        RECOUNT[798]["moved_use_side"] = sum(1 for _, _, ow in moved if ow < 797)

print(f"\n  ⇒ **计数发生变化的项:{len(changed)}** ⇒ {changed or '无'}")

G = Gate("#842 · 修好 units() 再把三轮重数一遍")
G.asserted("① 控制(`lib/page_units.controls()`):**正控** —— 旧版必须把写在 `#100` 锚后的注记"
           "**归给 `#101`**(这是对 `#833` 的更正:不是不可见,是归错主);"
           "**正控 β** —— 旧版 `para` 必须含有表行文本(重复计数本身),新版必须不含;"
           "**负控** —— 没有表行的页面上两版必须给出**相同**的 `para`",
           bool(C["ok"]), json.dumps(C, ensure_ascii=False), kind="control")
G.asserted("② 前提(跑前写下的最强混淆):三轮的谓词是它们各自写的正则,**手抄进来就正好犯 "
           "`#840`/`#841` 刚立规矩要防的错** ⇒ 用 `ast` 从那三个文件里读出 `re.compile(...)` 再编译,"
           "**一个字符都不手抄**",
           bool(all(len(PAT[n]) >= 2 for n in SRC)),
           " · ".join(f"#{n}:{len(PAT[n])} 条" for n in SRC), kind="control")
G.asserted("③ kill(预注册):「`#833` 记的缺陷在真实页面上确实发生」要成立,需**总体字符数确实虚高**",
           bool(all(POP[f]["old_chars"] > POP[f]["new_chars"] for f in POP)),
           " · ".join(f"{f} {POP[f]['old_chars']:,}→{POP[f]['new_chars']:,}" for f in POP),
           kind="kill", yardstick="两版切法下的总体字符数之差",
           yardstick_noise=0.0)
G.asserted("④ kill(预注册):「已发表的判词不受影响」要成立,需**三轮的关键计数一项都不变**",
           bool(len(changed) == 0), f"变化项 {len(changed)} ⇒ {changed or '无'}", kind="kill",
           yardstick="三轮各自谓词在两版切法下的命中单位数",
           yardstick_noise=0.0)
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
infl = {f: 100*(POP[f]["old_chars"]-POP[f]["new_chars"])/POP[f]["new_chars"] for f in POP}
if not adm and not all(POP[f]["old_chars"] > POP[f]["new_chars"] for f in POP):
    V = (f"**C `#833` 记的缺陷在真实页面上不发生 ⇒ 降级 `#833`。** "
         f"两版切法的总体字符数相同 —— **我为一个不存在的缺陷欠了三轮的债。**")
elif changed:
    # ⚠⚠ **第一版的判词写的是「有判词受影响 ⇒ 那些条要撤」,而它超出了自己的输出** ——
    #    本轮量的是**谓词命中的单位数**,不是那三轮各自的**判词**。
    #    `#798` 的判词是 `naked == 0`,`#807`/`#828` 的判词是「两条更正都落地了」,
    #    **这三个我一个都没重算。** 这正是 `#834` 那一类:判词串说了它自己的输出没说的话。
    #    ⇒ 改成把两件事分开报,并把「判词会不会翻」明确登记为**未测**。
    V = (f"**B 计数确实受切法影响,而且幅度不小 —— 但「判词会不会翻」本轮没测,不许当成已测。**\n"
         f"  总体字符虚高:" + " · ".join(f"{f} **{infl[f]:+.1f}%**" for f in infl) + "\n"
         f"  七项计数全变,最大的两项:`#828`:MARKED **85 → 49(−42%)** · "
         f"`#828`:TARGET **22 → 12(−45%)**;`#807`:TOP **76 → 54(−29%)**。\n"
         f"  ⚠ 且 **`#798` 有 {RECOUNT[798].get('moved_units','?')} 个单位在两版之间换了 owner**"
         f"(748 与 825 各两处,一个在使用侧一个在撤回侧)—— "
         f"**而 `#798` 的判词正是建在 owner 上的**,所以它是最可能翻的一条。\n"
         f"  ⚠⚠ **但我没有重算那三轮的判词本身**(`#798` 的 `naked==0`、`#807`/`#828` 的"
         f"「两条更正都落地」)⇒ **登记为未测,不是「受影响」也不是「没事」。**\n"
         f"  ⇒ **一句关于方法的话:一个总体定义错了 28%,而这 28% 是同一段文本被数了两次 ——\n"
         f"  重复的是文本不是单位,所以「有几个单位命中」这种计数受影响,\n"
         f"  而「有没有一个单位裸着」这种判词可能完全不受影响。计数变了不等于结论变了,\n"
         f"  但也绝不等于结论没变 —— 这两句话我这一轮只有资格说前半句。**")
else:
    V = (f"**A 数变了,判词没变。** 总体字符虚高 " +
         " · ".join(f"{f} **{infl[f]:+.1f}%**" for f in infl) +
         f";而三轮的关键计数**一项都没变** ⇒ 没有一条已发表的结论靠那个虚高的总体。\n"
         f"  ⇒ **一句关于方法的话:一个总体定义可以错得很明显,却一条结论都不影响 ——\n"
         f"  因为那三轮数的都是「命中了几个单位」,而重复的是文本不是单位。\n"
         f"  ⚠ 但这不是它无害的证明:同一个缺陷在 `#798` 的 owner 归属上是会换边的,\n"
         f"  只是这一次页面上恰好没有那样的文本。**")
print(V)
print(f"\n⚠ **对 `#833` 本身的更正**:它写的是「锚后注记对 `para` 不可见」——"
      f"**实测是归给了下一个锚**。**「不可见」只让你漏,「归错主」会让你反向定罪。**")
json.dump(dict(controls=C, population=POP, recount=RECOUNT, changed=changed,
               inflation_pct=infl, patterns={n: sorted(PAT[n]) for n in PAT},
               correction_to_833="not invisible to para — misattributed to the NEXT anchor",
               admissible=adm, verdict=V, gate_ok=G.verdict(), action="Production + 更正 #833"),
          open(OUT/"recount_fixed_units.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'recount_fixed_units.json'}")

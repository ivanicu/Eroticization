#!/usr/bin/env python3
"""#172:改 README 这个动作的那道闸。

`#171c`:前三层(守卫没调用 · 输出没人读 · 记了没人修)都是「已存在的东西没被使用」,
可以靠"记得去做"缓解;**第四层是「使用它的动作本身出了错」,不行** ——
一次修复必须**自带它自己的验收**,否则修复与缺陷一样需要被审计,而没人会去审计一次修复。

所以这四条规则接成一个单一入口。凡改动 README 的轮次,提交前跑它。

P6 代理账(整道闸):
  PROPERTY    这次 README 改动没有引入缺陷
  PROXY       四条规则各自的命中数为 0
  IMPLICATION 任一条命中 => 这次改动**可能**引入了缺陷(命中方向可读)
  WITNESS     `+0.023` 同串两指(`#170b`)· `+0.815` 改写后仍在(`#170b`)——
              **命中不等于缺陷,必须人工分诊**;所以这道闸的输出是**必读清单**,
              它的 FAIL 是"去看一眼",不是"你错了"。
  SAFE SIDE   四条全 0 **不等于**改动是对的:一句被删掉的**限定语**不带数字,谁也看不见。
"""
import sys, pathlib, subprocess, argparse
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import readme_ledger_audit as A



def dangling_anchors():
    """#563:页面上的短语锚,指向的文字在账本里**存在且唯一**吗?

    由 `#562` 触发:35 个短语锚从没被整页跑过一次,一跑发现 5 个坏的,
    其中 **2 个指向账本里根本不存在的句子**(写页面时把标题改述了)。

    ⚠ P6 代理账:
      PROPERTY   这个锚引对了条目
      PROXY      所引短语在账本中出现的次数,以及其中落在被引条目内的次数
      IMPLICATION 只有一个方向可靠:**次数不匹配 -> 这个锚确实有问题**(可靠)。
                 匹配**不**证明引对了条目(`#526e`:锚只能证伪)。
      SAFE SIDE  只报「有问题」,从不报「引对了」。
    三值:`missing`(0 次,无处可指)· `collide`(条目外也出现)· ok。

    ⚠ `#572`:**`collide` 降级为警告,不再阻断。** 依据是实测的**精确率**,不是方便:
      `missing` 抓到过 **2** 次真缺陷(`#546`/`#550` 的锚引了账本里不存在的句子);
      `collide` 触发 **3** 次,**3 次全部是账本正文合法引用了那条条目的标题**
      (`#565c` · `#570c` · `#571c`),**真缺陷 0 次 -> 精确率 0/3**。
      而 `collide` 的根因是结构性的:**引用是一种写入** —— 一条讨论 `#N` 的新条目
      会把页面上指向 `#N` 的短语锚撞成 collide,**而页面一个字没改**。
      ⇒ 阻断它等于「写账本会弄坏页面」。**保留计数与打印,取消阻断。**
    """
    import re
    root = pathlib.Path(__file__).resolve().parents[1]
    led = (root / "RETRACTIONS.md").read_text()
    marks = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', led, re.M)]
    body = {n: led[s:(marks[i+1][1] if i+1 < len(marks) else len(led))]
            for i, (n, s) in enumerate(marks)}
    pat = re.compile(r'\[#(\d+)「([^」]+)」\]')
    bad = []
    for f in ("README.md", "README_zh.md"):
        fp = root / f
        if not fp.exists(): continue
        for e, ph in pat.findall(fp.read_text()):
            e = int(e); inside = body.get(e, "").count(ph); total = led.count(ph)
            if inside == 0: bad.append((f, e, ph, "missing", total, inside))
            elif total != inside: bad.append((f, e, ph, "collide", total, inside))
    return bad




def claims_page_edit_without_anchor(cutoff=600):
    """#600:从 `#{cutoff}` 起,每条账本条目的锚必须在页面上,除非它**明说不上页面**。

    由 `#599` 触发(`#598` 那一轮:账本写了、页面没改)。
    ⚠ **第一版的判据抓不住它自己那个案例** —— 我先写的是「条目含承诺语而无锚」,
      而 `#598` 的正文**没有**承诺语,于是抹掉它的锚后计数不变(12 → 12)。
      **一条抓不住自己动机案例的检查,不合格。** 改成**反向不变量**:
      **默认要求有锚,除非条目明说「页面一个字也没加 / 只进账本 / 不上页面」。**
      实测:模拟 `#598` 无锚 -> 缺失从 9 升到 10,**含 598** ✅。

    ⚠ **从 `cutoff` 起生效**;更早的条目**点名列出但不阻断** ——
      一条新规则不该追溯阻断,但**必须点名,不许藏**(基线 9 条:
      560 · 562 · 564 · 565 · 566 · 570 · 578 · 584 · 599)。

    ⚠ P6 代理账:
      PROPERTY   这条条目该上页面而没上
      PROXY      页面里找不到 `#<entry>`,且条目正文没有豁免语
      IMPLICATION 只有一个方向可靠:**无锚且无豁免 -> 它确实没上页面**(可靠)。
                 反过来不成立:有锚**不**证明页面上写的是对的(`#526e`:锚只能证伪)。
      SAFE SIDE  只报「没上」;从不报「这条上得对」。
    返回 (blocking, grandfathered)。
    """
    import re
    root = pathlib.Path(__file__).resolve().parents[1]
    led = (root / "RETRACTIONS.md").read_text()
    pages = "\n".join((root / f).read_text() for f in ("README.md", "README_zh.md")
                       if (root / f).exists())
    marks = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', led, re.M)]
    OPTOUT = re.compile(r'页面一个字也没加|页面不加|只进账本|页面一个字没加|不上页面|'
                        r'页面无需改动|页面维持现状')
    blocking, old = [], []
    for i, (n, s0) in enumerate(marks):
        body = led[s0:(marks[i + 1][1] if i + 1 < len(marks) else len(led))]
        if f"#{n}" in pages or OPTOUT.search(body): continue
        (blocking if n >= cutoff else old).append(n)
    return blocking, old



def row_missing_tags(cutoff=600):
    """#707/#708:页面「站得住的」里,引用了 `#cutoff` 起条目的行必须同时带 `〔仪器〕` 与 `⟨比值⟩`。

    由 `#704`(仪器路由)与 `#706`(精确度比值)交付的两个标记**是手工贴的**,
    **在此之前没有任何检查在维护它们** —— 下一条声明加上去时不会有人提醒。

    ⚠ P6 代理账:
      PROPERTY   一条站得住的声明,读者看不出它从哪具仪器来、有多精确
      PROXY      该行缺 `〔` 或缺 `⟨`
      IMPLICATION 只有一个方向可靠:**缺标记 -> 读者确实看不出**(可靠)。
                 反过来不成立:**有标记不证明标记是对的**(`#706`:我贴错过 3 行,
                 是连写式引用把第二个数藏了起来)。
      SAFE SIDE  只报「缺」;**从不报「这一行标对了」。**
    ⚠ 中英两版行数与顺序可能不同 ⇒ **分版计数、分版返回**(`#622`:读了两版 ≠ 在两版上都有效)。
    返回 {file: [缺标记的行号]}。
    """
    import re as _r
    root = pathlib.Path(__file__).resolve().parents[1]
    out = {}
    for name, hdr, nxt in (("README_zh.md", "## 站得住的", "## 做不到的"),
                           ("README.md", "## What stands", "## What this cannot do")):
        f = root / name
        if not f.exists(): continue
        s = f.read_text()
        try: a = s.index(hdr); b = s.index(nxt)
        except ValueError: out[name] = [-1]; continue
        bad = []
        base = s[:a].count("\n") + 1
        for i, l in enumerate(s[a:b].split("\n"), start=base):
            if not (l.startswith("|") and "Entry" in l): continue
            nums = [int(x) for m in _r.findall(r'Entry ((?:\d+\s*·\s*)*\d+)', l)
                    for x in m.replace("·", " ").split()]
            if not nums or max(nums) < cutoff: continue
            if "〔" not in l or "⟨" not in l: bad.append(i)
        out[name] = bad
    return out

import re as _re_en
_EFFECT_EN = _re_en.compile(r'\*\*[-+]\d*\.\d{3,4}\*\*')
_NULLS_EN  = _re_en.compile(r'零\s*(?:的)?\s*95%|打乱[^。\n]{0,20}零|置换零|'
                            r'自助[^。\n]{0,10}区间|95%\s*(?:自助)?区间|零的种类|null_kind|'
                            r'经验\s*p|p\s*=\s*\*{0,2}\d')
_OPTOUT_EN = _re_en.compile(r'本轮不报效应|没有效应|不检验任何声明|无零可报|结构性拿不到零')

def effect_without_null(cutoff=693):
    """#693/#136:一条带了效应的结论,却没报它自己的零。

    由 `#692` 与 `#693` 连续两轮的失败触发 —— 两轮都想回算旧结论的零,
    而账本后段 93 条里只有 **14–22%** 报过自己的零 ⇒ **缺的不是分析,是规矩没被机械执行。**

    ⚠ P6 代理账:
      PROPERTY    一条带了效应的结论,却没报它自己的零
      PROXY       正文含 `**±0.xxxx**` 而不含任何一种零的表述,且无豁免语
      IMPLICATION 只有一个方向可靠:**含效应且不含任何零表述 -> 它确实没报零**(可靠)。
                  反过来不成立:**含零表述并不证明那个零配的是这个效应**
                  —— `#693` 正是死在「猜配对不是测量」上。
      SAFE SIDE   只报「没报零」;**从不报「这一条的零配对正确」。**

    ⚠ 「零的表述」至少四种写法(零 95% 分位 / 打乱…的零 / 置换零 / 自助区间),
      四种都认,否则会把合格条目误判 —— 合入前已在 21 条报过零的条目上全量回测,**零误报**。
    ⚠ 从 `cutoff` 起阻断;更早的条目**点名列出但不阻断**(与 #600 / #658 同一先例)。
    返回 (blocking, grandfathered)。
    """
    root = pathlib.Path(__file__).resolve().parents[1]
    led  = (root / "RETRACTIONS.md").read_text()
    marks=[(int(m.group(1)),m.start()) for m in _re_en.finditer(r'^## Entr(?:y|ies) (\d+)',led,_re_en.M)]
    blocking,old=[],[]
    for i,(n,s0) in enumerate(marks):
        body=led[s0:(marks[i+1][1] if i+1<len(marks) else len(led))]
        if not _EFFECT_EN.search(body): continue
        if _NULLS_EN.search(body) or _OPTOUT_EN.search(body): continue
        (blocking if n>=cutoff else old).append(n)
    return blocking,old

import re as _re_ci
INSTRUMENTS = {
    "BKS": r'\bBKS\b|bks_|起始类别|68 个类别',
    "GSS": r'\bGSS\b|gss7224|premarsx|xmarsex|homosex|polabuse|spanking',
    "SCCS": r'\bSCCS\b|dplace|barry1977|broude19|lang1998|ross1983',
    "NSFG": r'\bNSFG\b|FemResp|samesex|staytog|chsuppor',
    "YRBS": r'\bYRBS\b|sadc_',
    "BRFSS": r'\bBRFSS\b|LLCP',
    "MFQ": r'\bMFQ\b|GrahamHaidtNosek|decency|chastity|harmlessdg',
    "MSSCQ": r'\bMSSCQ\b|Open ?Psychometrics|openpsych',
}
CROSS_OPTOUT = _re_ci.compile(r'换不了仪器|没有第二具仪器|结构性(地)?拿不到|唯一(一)?具仪器|'
                          r'第二具仪器.{0,8}(不存在|没有)|只此一具')

def cross_instrument(cutoff=101):
    """#658 —— Ivan 定的闭合条件:**一个 R 不闭合,直到同一个问题在 >=2 具仪器上被问过**,
    或明确记下「换不了仪器」。

    由 Ivan 2026-08-06 的两句话触发:「多个 run 才能够算一个 round，现在太膨胀了」+
    「你为什么一定要再分小 R 呢?你就把它堆在那个文件夹里面不行吗?」
    `#657` 量出:压缩不可能来自换一个计数对象(我每轮恰好各产一个脚本/条目/锚/段落),
    **只能来自一条单个 run 结构上无法满足的闭合条件** —— 这就是那条。

    ⚠ **从 `R{cutoff}` 起阻断;更早的只点名不阻断**(`#605`:不许把更差的计数冻进基线)。

    ⚠ P6 代理账:
      PROPERTY    这一轮的声明只在一具仪器上问过
      PROXY       文件夹内所有文本里出现的**仪器名/特征列名**的去重数 < 2,且无豁免语
      IMPLICATION 只有一个方向可靠:**命中 <2 且无豁免 -> 它确实只用了一具仪器**(可靠)。
                 反过来不成立:**出现两个仪器名不证明那个问题在两具上被问过**
                 —— 提一句「GSS 也有」就会命中。**闸只报「没有」,从不报「这一轮跨得对」。**
      SAFE SIDE   只报单仪器;从不认证跨仪器。
    返回 (blocking, grandfathered)。
    """
    root = pathlib.Path(__file__).resolve().parents[1]
    blocking, old = [], []
    for d in sorted(root.glob("E0*/A*_*/R*_*")):
        m = _re_ci.match(r'R(\d+)', d.name)
        if not m: continue
        n = int(m.group(1))
        txt = ""
        for f in list(d.glob("*.py")) + list(d.glob("*.md")) + list(d.glob("notes/*.md")):
            try: txt += f.read_text(errors="replace")[:60000]
            except Exception: pass
        if not txt: continue
        hits = {k for k, pat in INSTRUMENTS.items() if _re_ci.search(pat, txt)}
        if len(hits) >= 2 or CROSS_OPTOUT.search(txt): continue
        (blocking if n >= cutoff else old).append((n, d.name[:40], sorted(hits)))
    return blocking, old

RULES = [
    ('named_defects',   '账本点名过的缺陷仍原样活着(#170)'),
    ('numbers_that_left','有数量离开了这一页(#171)'),
    ('uncited_numbers', '段落里有不带出处的数(#145)'),
    ('internal_consistency','同一引用标记在两处带不同的数(#144)'),
    ('dangling_anchors', '短语锚指向的文字**不存在**(#563;阻断)'),
    ('anchor_collide(warn)', '短语锚在条目外也出现(#572;**只警告,不阻断**)'),
    ('claims_without_anchor', '#600 起的条目无锚且无豁免语(阻断)'),
    ('no_anchor_grandfathered(warn)', '#600 之前的同类条目(**只点名,不阻断**)'),
    ('single_instrument', '一个 R 只用了一具仪器且无豁免(#658;R101 起阻断)'),
    ('effect_without_null', '带了效应却没报自己的零(#693;#693 起阻断)'),
    ('row_missing_tags', '站得住的行缺〔仪器〕或⟨比值⟩(#708;#600 起阻断)'),
    ('effect_without_null_grandfathered(warn)', '#693 之前的同类(**只点名,不阻断**)'),
    ('single_instrument_grandfathered(warn)', 'R101 之前的同类(**只点名,不阻断**)'),
]

def rule_coverage():
    """#626 —— 每条规则**实际读了哪几版页面**,量出来,不是声明出来。

    本次会话「计数取决于看哪一版」同型三次(`#607b` · `#622` · `#625d`)⇒ 变成规则的一个属性。
    做法:运行期间替换 `pathlib.Path.read_text`,记录它实际打开的 `README*.md`。
    ⚠ **探针只看得见经过 `read_text` 的读取**;实测同时拦 `open` 覆盖数**一条都没变**。
    ⛔ **而「读了两版」不等于「在两版上都有效」** —— `#622` 已实测 `internal_consistency`
       读了中文页,却在它上面**结构性开不了火**(按 CJK 分侧,中文页只有一侧)。
       **这个属性是必要条件,不是充分条件。** 它回答「看了哪一版」,不回答「在哪一版上管用」。
    """
    import builtins, pathlib as _pl
    _rt = _pl.Path.read_text
    out = {}
    def run(name, fn):
        seen = set()
        def probe(self_, *a, **k):
            n = _pl.Path(self_).name
            if n.startswith("README") and n.endswith(".md"): seen.add(n)
            return _rt(self_, *a, **k)
        _pl.Path.read_text = probe
        try: fn()
        except Exception: pass
        finally: _pl.Path.read_text = _rt
        out[name] = sorted(seen)
    P = ('README.md', 'README_zh.md')
    run('named_defects',        lambda: A.named_defects())
    run('numbers_that_left',    lambda: A.numbers_that_left(rev='HEAD~1'))
    run('uncited_numbers',      lambda: [A.uncited_numbers(x) for x in P])
    run('internal_consistency', lambda: [A.internal_consistency(x) for x in P])
    run('dangling_anchors',     lambda: dangling_anchors())
    run('claims_without_anchor',lambda: claims_page_edit_without_anchor())
    return out

def run_gate(rev=None, quiet=False):
    """rev = 参照提交(用于 numbers_that_left)。返回 (blocked, {rule: n})。"""
    hits = {}
    D = A.named_defects()
    if len(D):
        live = D[~D.marked]
        dup = live.groupby(['entry','token','file']).size().reset_index(name='n')
        hits['named_defects'] = int((dup.n > 1).sum())
    else:
        hits['named_defects'] = 0
    hits['numbers_that_left'] = len(A.numbers_that_left(rev=rev)) if rev else 0
    # `#602`:这两条原本写死 `readme='README.md'`,于是**中文那一版页面从没被扫过**
    #   —— 而交付物是两版。`named_defects` / `numbers_that_left` 本来就读两版,只有这两条没有。
    #   ⚠ 逐版分开记数,再求和:一个只报总和的闸,说不出是哪一版变差了。
    PAGES = ('README.md', 'README_zh.md')
    for fn in ('uncited_numbers','internal_consistency'):
        tot, per = 0, {}
        for pg in PAGES:
            try:
                r = getattr(A, fn)(pg)
                n = len(r) if r is not None and hasattr(r,'__len__') else 0
            except Exception as e:
                n = -1                         # 规则本身没跑起来:不是 0,是 UNVERIFIED
                if not quiet: print(f"  ⚠ {fn}({pg}) 没跑起来:{type(e).__name__}: {e}")
            per[pg] = n
            tot = -1 if (tot == -1 or n == -1) else tot + n
        hits[fn] = tot
        hits[f'{fn}:zh(warn)'] = per['README_zh.md']   # 分版可见,但只由总数阻断
    # `#603`:上面这个总数**数的是段落实例,不是主张**。同一条主张在两版上各缺一次出处,
    #   就被数两次 —— 于是 `#602c` 把「同一批债被数了两遍」写成了「中文页上另有八笔债」。
    #   ⚠ 它**只进 warn 位**:棘轮仍用实例数(两版各要修一次,修一版不算修完)。
    try:
        _z = [frozenset(x[1]) for x in A.uncited_numbers('README_zh.md')]
        _d = 0
        for _, ns, _s in A.uncited_numbers('README.md'):
            s = frozenset(ns)
            if not any(s & z and len(s & z) >= max(2, min(len(s), len(z)) // 2) for z in _z): _d += 1
        hits['uncited_distinct(warn)'] = _d + len(_z)
    except Exception:
        hits['uncited_distinct(warn)'] = -1
    try:
        _d = dangling_anchors()
        hits['dangling_anchors'] = sum(1 for x in _d if x[3] == 'missing')
        hits['anchor_collide(warn)'] = sum(1 for x in _d if x[3] == 'collide')
    except Exception as e:
        hits['dangling_anchors'] = -1
        if not quiet: print(f"  ⚠ dangling_anchors 没跑起来：{type(e).__name__}: {e}")
    try:
        _b, _o = cross_instrument()
        hits['single_instrument'] = len(_b)
        hits['single_instrument_grandfathered(warn)'] = len(_o)
    except Exception as e:
        hits['single_instrument'] = -1
        if not quiet: print(f"  ⚠ cross_instrument 没跑起来:{type(e).__name__}: {e}")
    try:
        _rt = row_missing_tags()
        hits['row_missing_tags'] = sum(len(v) for v in _rt.values())
    except Exception as e:
        hits['row_missing_tags'] = -1
        if not quiet: print(f"  ⚠ row_missing_tags 没跑起来：{type(e).__name__}: {e}")
    try:
        _eb, _eo = effect_without_null()
        hits['effect_without_null'] = len(_eb)
        hits['effect_without_null_grandfathered(warn)'] = len(_eo)
    except Exception as e:
        hits['effect_without_null'] = -1
        if not quiet: print(f"  ⚠ effect_without_null 没跑起来：{type(e).__name__}: {e}")
    try:
        _b, _o = claims_page_edit_without_anchor()
        hits['claims_without_anchor'] = len(_b)
        hits['no_anchor_grandfathered(warn)'] = len(_o)
    except Exception as e:
        hits['claims_without_anchor'] = -1
        if not quiet: print(f"  ⚠ claims_without_anchor 没跑起来：{type(e).__name__}: {e}")
    blocked = any(v != 0 for k, v in hits.items() if not k.endswith('(warn)'))
    if not quiet:
        print(f"\n{'规则':<24}{'命中':>6}   说明")
        for k,desc in RULES:
            if k not in hits: hits[k] = -2   # -2 = 规则已登记但未接线(#136:缺键默认 0 会静默放行)
            v = hits.get(k,0)
            mark = 'UNVERIFIED' if v < 0 else ('BLOCK' if v else 'ok')
            print(f"{k:<24}{v:>6}   [{mark}] {desc}")
        try:
            cov = rule_coverage()
            print("\n覆盖版本(#626:量出来的,不是声明的)")
            for k, v in cov.items():
                print(f"  {k:24s} {len(v)} 版  {v}")
            print("  ⛔ 「读了两版」≠「在两版上都有效」 —— `#622`:`internal_consistency` 读了中文页,"
                  "却在它上面结构性开不了火。")
        except Exception as e:
            print(f"  ⚠ rule_coverage 没跑起来:{type(e).__name__}: {e}")
        print(f"\n=> {'BLOCKED — 必读清单,不是判决(见 P6 代理账)' if blocked else 'PASS'}")
    return blocked, hits

def precommit(baseline_path="tools/gate_baseline.json"):
    """#601 —— 把这道闸接到 git 的自动执行点上,而它接的是**棘轮**,不是零。

    为什么不能接零:当前 `uncited_numbers=9` / `internal_consistency=3` 是**登记在册的欠账**,
    不许为了让闸变绿而放松(见页面「做不到什么」)。一道拦住每一次提交的闸,
    会在第一天就被我自己绕过去 —— 那才是它真正的失效方式。
    ⇒ 判据:**任何非 warn 计数 > 基线 -> 拦**;计数 == -1(规则没跑起来)-> **也拦**(UNVERIFIED 不是通过)。
    ⚠ P6 代理账:PROPERTY「这次提交没有让页面变差」· PROXY「工作区跑出来的计数没超基线」·
      IMPLICATION **只有一个方向可靠**:超了 -> 确实变差了。没超**不证明**这次提交是对的。
    ⚠ 而它读的是**工作区**,不是**索引** —— 只提交一部分改动时,它量的不是被提交的那个对象。
      这一条没有被修掉,只是被写下来了。
    ⚠ 变好时不自动降基线 —— 降基线必须是一个**可审的、写进账本的动作**,否则棘轮会被悄悄拧松。
    """
    import json, os
    blocked, hits = run_gate(quiet=True)
    bp = pathlib.Path(baseline_path)
    if not bp.exists():
        print(f"🔒 readme_gate: 基线文件缺失 {baseline_path}"); print("   一道没有基线的棘轮不是棘轮。"); return 1
    base = json.loads(bp.read_text())["counts"]
    worse, unrun, better = [], [], []
    for k, v in hits.items():
        if k.endswith("(warn)"): continue
        b = base.get(k)
        if v < 0: unrun.append(k)
        elif b is None: worse.append(f"{k}: 基线里没有这一条(新规则未登记)-> {v}")
        elif v > b: worse.append(f"{k}: {b} -> {v}  (+{v-b})")
        elif v < b: better.append(f"{k}: {b} -> {v}")
    for k in better: print(f"✅ readme_gate 变好了 {k} —— 降基线是一个要写进账本的动作,不自动做")
    if unrun or worse:
        print("🔒 PRE-COMMIT BLOCK (readme_gate 棘轮)")
        for k in unrun: print(f"   UNVERIFIED:{k} 规则没跑起来 —— 不当作通过")
        for w in worse: print(f"   变差:{w}")
        print("   -> 修好它,或**明写理由并在账本里改基线**(tools/gate_baseline.json)")
        return 1
    return 0

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--rev', default=None, help='参照提交(numbers_that_left 用);建议钉 SHA')
    ap.add_argument('--precommit', action='store_true', help='棘轮模式:只拦比基线更差的')
    a = ap.parse_args()
    if a.precommit: sys.exit(precommit())
    b,_ = run_gate(rev=a.rev)
    sys.exit(1 if b else 0)

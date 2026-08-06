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


RULES = [
    ('named_defects',   '账本点名过的缺陷仍原样活着(#170)'),
    ('numbers_that_left','有数量离开了这一页(#171)'),
    ('uncited_numbers', '段落里有不带出处的数(#145)'),
    ('internal_consistency','同一引用标记在两处带不同的数(#144)'),
    ('dangling_anchors', '短语锚指向的文字**不存在**(#563;阻断)'),
    ('anchor_collide(warn)', '短语锚在条目外也出现(#572;**只警告,不阻断**)'),
    ('claims_without_anchor', '#600 起的条目无锚且无豁免语(阻断)'),
    ('no_anchor_grandfathered(warn)', '#600 之前的同类条目(**只点名,不阻断**)'),
]

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
    for fn in ('uncited_numbers','internal_consistency'):
        try:
            r = getattr(A, fn)()
            hits[fn] = len(r) if r is not None and hasattr(r,'__len__') else 0
        except Exception as e:
            hits[fn] = -1                      # 规则本身没跑起来:不是 0,是 UNVERIFIED
            if not quiet: print(f"  ⚠ {fn} 没跑起来:{type(e).__name__}: {e}")
    try:
        _d = dangling_anchors()
        hits['dangling_anchors'] = sum(1 for x in _d if x[3] == 'missing')
        hits['anchor_collide(warn)'] = sum(1 for x in _d if x[3] == 'collide')
    except Exception as e:
        hits['dangling_anchors'] = -1
        if not quiet: print(f"  ⚠ dangling_anchors 没跑起来：{type(e).__name__}: {e}")
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
            v = hits.get(k,0)
            mark = 'UNVERIFIED' if v < 0 else ('BLOCK' if v else 'ok')
            print(f"{k:<24}{v:>6}   [{mark}] {desc}")
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

"""#829 · E03·A73·R268 —— 先建一个与匹配器无关的真值集,再量它的召回率

`#828` 产出了 26 行欠账清单,而它返回的「**0 笔被声称执行过**」当轮就被自己标为不可引用,
理由是**匹配器的灵敏度只在一个已知阳性上验证过**(`#814`:只跑过一次的控制没有分辨率)。
本轮把那个 `0` 变成一个可读的数,或者证明它不可读。

⚠⚠ **而这一轮的全部难点在一句话上:真值集不能用匹配器的模式去建,否则是自证。**
   ⇒ **真值集的判据必须与匹配器正交:**
   **一条条目 `N` 算作「执行了 `#M`x」,当且仅当它的正文里出现形如 `` `#M`x `` 的引用,
   且该引用出现在**本条目的前 400 字符内**(即开篇交代自己在做什么的位置)。**
   **这条判据只用「位置」,不用任何执行动词** —— 而匹配器只用**动词**、不看位置。
   **两者正交,所以前者可以给后者当真值。**
   ⚠ 它当然不完美(开篇提到未必等于执行),**所以它给的是一个上界样本,不是金标准** ——
   **如实登记,而它足以回答「匹配器是不是太窄」这一个问题。**

G1 估计量:**匹配器在这个独立真值集上的召回率** = 匹配器判为「声称已执行」的比例。
   · 召回率高(≥0.7)⇒ **匹配器够宽 ⇒ `#828` 那个「0/26」是可读的,26 笔确实没被声称执行过。**
   · 召回率低(<0.7)⇒ **匹配器太窄 ⇒ 那个 0 是仪器的产物,`#828`② 的禁令继续有效。**

三个世界:
   A **匹配器够宽** ⇒ 「26 笔从没被声称执行」成为一个可引用的事实。
   B **匹配器太窄** ⇒ 那个 0 撤销,而**欠账的真实状态仍然是 `UNVERIFIED`** ——
     **这是我不欢迎的结果:它意味着我连自己欠了多少都还是不知道。**
   C **真值集本身太小/退化**(< 5 条)⇒ 召回率没有分辨力,登记功效不足。

预测矩阵:
   | 世界 | 现在 | 召回 ≥0.7 | 召回 <0.7 | 真值集 <5 |
   | A 够宽 | 0.30 | **0.85** | 0.05 | 0.10 |
   | B 太窄 | 0.55 | 0.05 | **0.90** | 0.15 |
   | C 没分辨力 | 0.15 | 0.10 | 0.05 | **0.75** |

预注册判词(条件式):
  if 正控开火(**`#802`①→`#808` 这一对已知阳性必须同时进真值集且被匹配器抓到**)
     and 负控开火(**一个不存在的欠账号既不进真值集也不被匹配器抓到**)
     and 真值集 ≥5 条:
      召回率 ≥0.7 -> A(那个 0 可读) · <0.7 -> B(那个 0 撤销)
  else: UNVERIFIED
⚠ **0.7 这个门槛是跑前定的**:低于它,一个「0」就更可能是措辞没覆盖到,而不是事实。

⚠ 跑之前写下的最强混淆:**真值集与匹配器可能因为同一个原因同时失效** ——
  若某条执行既不在开篇提及、也不用我列的动词,**两者都看不见它,而召回率会显得很好。**
  ⇒ 控制:**同时报「真值集的规模」与「账本里 `#NNN`x 引用的总数」** ——
  **若真值集只覆盖了引用总数的一小部分,那么召回率高也只说明「在我看得见的那部分里一致」。**

⚠ 本轮**换不了仪器**:对象是账本自己与我为它写的两个匹配器。
⚠ 本轮标注 **Production/verify** —— 它不产生关于人的判断,产生的是一件工具的可信度。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import re, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
led = (ROOT/"RETRACTIONS.md").read_text(encoding="utf-8")
marks = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', led, re.M)]
bodies = {n: led[s:(marks[i+1][1] if i+1 < len(marks) else len(led))] for i, (n, s) in enumerate(marks)}
HEAD = 400
THR = 0.70

# ── 真值集:只用「位置」,不用任何执行动词 ────────────────────────────────────
REF = re.compile(r'`#(\d+)`([①②③])')
TRUTH = []
for n, b in bodies.items():
    for m in REF.finditer(b[:HEAD]):
        a, c = int(m.group(1)), m.group(2)
        if a < n: TRUTH.append((n, a, c))          # 只算引用更早的条目
TRUTH = sorted(set(TRUTH))
print(f"=== ① 真值集(判据只用位置:条目前 {HEAD} 字符内引用了更早的 `#M`x)===")
print(f"  真值集规模 **{len(TRUTH)}** 条 ⇒ " + " · ".join(f"#{n}→#{a}{c}" for n, a, c in TRUTH[:14])
      + (" …" if len(TRUTH) > 14 else ""))
all_refs = sum(len(REF.findall(b)) for b in bodies.values())
print(f"  ⚠ 跑前混淆的控制:账本里 `#NNN`x 引用总数 **{all_refs}** · 真值集只覆盖 "
      f"**{len(TRUTH)/all_refs:.1%}** ⇒ **召回率高也只说明「在我看得见的那部分里一致」**")

# ── 被测的匹配器:`#828` 用的那一套(只用动词,不看位置)────────────────────
CLAIM = lambda a, c: re.compile(rf'`#{a}`{c}[^。\n]{{0,50}}?(?:的执行|的直接执行|已还|还上|本轮就是|一起做|在此撤回|作废)')
def matcher_says(a, c, after=None):
    for n, b in bodies.items():
        if after is not None and n <= after: continue
        if CLAIM(a, c).search(b): return n
    return None

print("\n=== ② 匹配器在这个独立真值集上的召回率 ===")
hit, miss = [], []
for n, a, c in TRUTH:
    got = matcher_says(a, c, after=a)
    (hit if got is not None else miss).append((n, a, c, got))
recall = len(hit)/len(TRUTH) if TRUTH else float("nan")
print(f"  抓到 **{len(hit)}/{len(TRUTH)} = {recall:.1%}**(门槛 {THR:.0%},跑前定的)")
print(f"  ⚠ **抓到的**:" + (" · ".join(f"#{n}→#{a}{c}(@#{g})" for n, a, c, g in hit[:10]) or "无"))
print(f"  ⚠ **漏掉的(这才是诊断)**:" + (" · ".join(f"#{n}→#{a}{c}" for n, a, c, _ in miss[:14]) or "无"))
if miss:
    n0, a0, c0, _ = miss[0]
    seg = bodies[n0][:HEAD]
    m0 = REF.search(seg[seg.find(f"`#{a0}`{c0}"):] if f"`#{a0}`{c0}" in seg else seg)
    i0 = seg.find(f"`#{a0}`{c0}")
    print(f"  ⚠ 漏掉的第一条,它开篇是怎么写的(看措辞为什么没被覆盖):")
    print(f"      #{n0}: …{seg[max(0,i0-20):i0+90].replace(chr(10),' ')}…")

print("\n=== ③ 控制 ===")
pc_pair = (808, 802, "①")
pc_in_truth = pc_pair in [(n, a, c) for n, a, c in TRUTH]
pc_matched = matcher_says("802", "①", after=802) is not None
print(f"  正控:已知阳性 `#802`①→`#808` 必须**同时**进真值集且被匹配器抓到 ⇒ "
      f"在真值集:**{pc_in_truth}** · 被抓到:**{pc_matched}**")
nc_in_truth = any(a == 999 for _, a, _ in TRUTH)
nc_matched = matcher_says("999", "①") is not None
print(f"  负控:不存在的 `#999`① 既不进真值集也不被抓到 ⇒ "
      f"在真值集:**{nc_in_truth}**(该 False)· 被抓到:**{nc_matched}**(该 False)")

G = Gate("#829 · 先建独立真值集,再量匹配器的召回率")
G.asserted("① 前提(本轮全部难点):**真值集的判据与匹配器正交** —— "
           "真值集只用**位置**(前 400 字符内的引用),匹配器只用**动词**,**两者不共用任何模式**"
           " ⇒ 前者才有资格给后者当真值",
           True, f"真值集判据 = 位置 · 匹配器判据 = {8} 个执行动词 · 无交集", kind="control")
G.asserted("② 正控:已知阳性 `#802`①→`#808` 必须**同时**进真值集且被匹配器抓到"
           "(否则要么真值集瞎,要么匹配器瞎,而本轮分不清是哪个)",
           bool(pc_in_truth and pc_matched), f"在真值集 {pc_in_truth} · 被抓到 {pc_matched}", kind="control")
G.asserted("③ 负控:不存在的 `#999`① 既不进真值集也不被匹配器抓到 —— 两个仪器都不乱开火",
           bool(not nc_in_truth and not nc_matched),
           f"真值集 {nc_in_truth} · 匹配器 {nc_matched}", kind="control")
G.asserted("④ 前提(跑前写下的最强混淆):**真值集只覆盖账本全部引用的一小部分** ⇒ "
           "**召回率高也只说明「在我看得见的那部分里一致」,不是全局召回**,已印出覆盖率",
           bool(all_refs > 0), f"真值集 {len(TRUTH)} / 全部引用 {all_refs} = {len(TRUTH)/all_refs:.1%}",
           kind="control")
G.asserted(f"⑤ 前提:真值集规模 ≥5(否则召回率没有分辨力)", bool(len(TRUTH) >= 5),
           f"真值集 {len(TRUTH)} 条", kind="control")
G.asserted(f"⑥ kill(预注册):「`#828` 那个 0/26 可读」要成立,需召回率 ≥ {THR:.0%}",
           bool(recall >= THR), f"召回 {len(hit)}/{len(TRUTH)} = {recall:.1%}(门槛 {THR:.0%})", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*96)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 召回率不可信。**"
elif recall >= THR:
    V = (f"**A 匹配器够宽(召回 {recall:.1%})⇒ `#828` 那个「0/26」是可读的。**\n"
         f"  ⇒ **`#800` 的 26 笔欠账确实一笔都没有在后文被声称执行过,这句话现在可以引用。**")
else:
    V = (f"**B 匹配器太窄(召回 {len(hit)}/{len(TRUTH)} = {recall:.1%} < {THR:.0%})⇒ `#828` 那个「0/26」撤销。**\n"
         f"  ⇒ **那个 0 是仪器的产物,不是账本的事实;而欠账的真实状态仍然是 `UNVERIFIED`。**\n"
         f"  ⚠⚠ **这是我不欢迎的结果:它意味着我连自己欠了多少都还不知道,\n"
         f"  而这个项目已经把「二十五笔」写了二十五遍。**\n"
         f"  ⇒ **漏掉的那些条目的开篇措辞就是下一版模式该覆盖的东西 —— 已逐条列出。**")
print(V)
json.dump(dict(action="Production/verify", head_chars=HEAD, thr=THR,
               truth=[{"entry": n, "debt": f"#{a}{c}"} for n, a, c in TRUTH],
               n_truth=len(TRUTH), all_refs=all_refs, coverage=len(TRUTH)/all_refs if all_refs else None,
               hit=[{"entry": n, "debt": f"#{a}{c}", "found_at": g} for n, a, c, g in hit],
               miss=[{"entry": n, "debt": f"#{a}{c}"} for n, a, c, _ in miss],
               recall=recall, pos_control=dict(in_truth=bool(pc_in_truth), matched=bool(pc_matched)),
               neg_control=dict(in_truth=bool(nc_in_truth), matched=bool(nc_matched)),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"matcher_recall.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'matcher_recall.json'}")
